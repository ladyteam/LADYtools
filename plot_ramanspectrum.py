#!/usr/bin/python

def strType(var):
    try:
        if int(var) == float(var):
            return 'int'
    except:
        try:
            float(var)
            return 'float'
        except:
            return 'str'

def isNumeric(var):
    isnum=strType(var)
    if (isnum == 'int' or isnum == 'float'):
        return (True)
    return(False)

# Lorenz A -- integral, max I=y0+2A/(pi*w) FWHM=w 
def Lorenz(y0,x,xc,w,A):
    if(4*(x-xc)**2+w**2==0):
        print('Filed calculate lorenc at x=%f w=%f' % (x,w))
        return(0)
    return(y0+A*2/pi*w/( 4*(x-xc)**2+w**2 ))

def Gauss(y0,x,A,w,xc):
    return (y0+ A/(w*sqrt(pi/2))*exp(-2*(x-xc)**2/w**2) )


# $hkt=1.0545919/(1.3806221*$T)*10**(-11)
# $C=($nuex-$nu)**4/(  (1 - exp(-1*$hkt*$nu*2*$pi*$cs)) * 30*$nu );
def ramanpref(wl,wi,T):
    hkc=-1.4394
    cm1toau=4.55634*10**-6
    c=137.036
    try:
        return ( ((wl-wi)**4 /  ((1.0 - exp( hkc*wi/T ))* wi))/(2*c**4))*cm1toau**3
    except:
        print("error in ramanpref wi=%f T=%f" % (wi,T))
        return(0)


#def Voigtian()
#Gauss+Lorenz
import sys
from math import pi
from math import exp
from math import log
from math import sqrt
import numpy as np
import scipy as sc
from scipy import linalg
import argparse
import re

nm2cm1=10000000

parser = argparse.ArgumentParser(description='Script to plot scpectrum from activity data. !!! Please note be carefull with -l flag')

parser.add_argument("-i", "--input", action="store", type=str, dest="input_fn",  help="Input data filename")
parser.add_argument("-o", "--output", action="store", type=str, dest="output_fn",  help="Output data filename")
parser.add_argument("-c", "--comment", action="store", type=str, dest="comment",  help="Comments to header of output file")
parser.add_argument("-s", "--step", action="store", type=float, dest="step", default=0.1, help="Step (resolution) of output spectra")
parser.add_argument("-r", "--range", action="store", type=str, dest="range", default='0 1000',  help="Range from to. Two float numbers space separated")
parser.add_argument("-T", "--temp", action="store", type=int, dest="temp", default=0, help="Temperature of Raman spectra experiment")
parser.add_argument("-e", "--excline", action="store", type=float, dest="lexc", default=514.8, help="Excitation line wavenumber in nm")
parser.add_argument("-l", "--linear", action='store_true', dest="islin", default=False, help="Set the flag if the components of tensor are not in square degree")
parser.add_argument("-m", "--mult", action="store", type=float, dest="mult", default=1, help="Multiplyer")
parser.add_argument("-n", "--norm", action='store_true', dest="norm", default=False, help="Set the flag to make normilized spectra")

args = parser.parse_args()
nulaser=nm2cm1/args.lexc
print('nulaser=%f' % nulaser)
if (args.input_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

try:
    input_fh = open(args.input_fn, 'r')
except IOError:
    print("ERROR Couldn't open input file, exiting...\n")
    sys.exit(1)

if (args.output_fn == None):
    print('Error. No output filename was given.')
    sys.exit(1)

try:
    output_fh = open(args.output_fn, 'w')
except IOError:
    print("ERROR Couldn't open input file, exiting...\n")
    sys.exit(1)


if (isNumeric(args.range.split()[0]) and isNumeric(args.range.split()[1])):
    spefrom=float(args.range.split()[0])
    speto=float(args.range.split()[1])
else:
    print('Error. No spectrum range was provided')
    sys.exit(1)


xc=[]
A=[]
w=[]

if(args.islin):
    print("Warning. Linear components provided not intencity. Components of raman tensor will be squared")

for line in input_fh:
    if (re.search('^\s*$',line)):
        continue
    if (re.search('^\s*#',line)):
        continue
    if (len(line.split())<3):
        print('Warning the line %s contains less data than needed' % line)
        continue
    if (isNumeric(line.split()[0]) and isNumeric(line.split()[1]) and isNumeric(line.split()[2])):
        xc.append(float(line.split()[0]))
        if(args.islin):
            A.append(float(line.split()[1])**2)
        else:
            A.append(float(line.split()[1]))
        w.append(float(line.split()[2]))


np=int(abs(speto-spefrom)/args.step)
x=[spefrom+args.step*i for i in range(np)]
y=[0 for i in range(np)]
#Add each lorenzian sequentially
for i in range(len(xc)):
#Scan whole range of spectra
    for k in range(np):
        if(args.temp == 0):
            y[k]=y[k]+Lorenz(0,x[k],xc[i],w[i],A[i])*args.mult
        else:
            y[k]=y[k]+Lorenz(0,x[k],xc[i],w[i],A[i])*ramanpref(nulaser,xc[i],args.temp)
if (args.norm):
    max=0
    for i in range(np):
        if(max<abs(y[i])):
            max=abs(y[i])
    for i in range(np):
        y[i]=y[i]/max

if (args.comment):
    output_fh.write(args.comment+'\n')
for i in range(np):
    output_fh.write('% 10.7f\t% 14.10f\n' % (x[i],y[i]*args.mult))
