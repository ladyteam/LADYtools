#!/usr/bin/env python3.4
# -*- coding: utf-8 -*-
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Shift atom positions in cif file
#
# Author: Eugene Roginskii


def chunks(l, n):
    if n < 1:
        n = 1
    return [l[i:i + n] for i in range(0, len(l), n)]

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


import numpy as np
import re
from math import sqrt
from math import ceil
import sys
from math import pi
from shutil import move
import os
import datetime
import time
import argparse
import xml.etree.cElementTree as etree
import yaml


parser = argparse.ArgumentParser(description='Shift atom positions in cif file')

parser.add_argument("-i", "--input", action="store", type=str, dest="input_fn", default='findsym.cif', help="Inputfile in cif format")
parser.add_argument("-o", "--output", action="store", type=str, dest="output_fn", default='out.cif', help="Inputfile in cif format")
parser.add_argument("-s", "--shift", action="store", type=str, dest="shift_fn", default='shift.dat', help="The file name of shift atoms data")

args = parser.parse_args()

try:
    input_fh = open(args.input_fn, 'r')
except IOError:
    print ("ERROR Couldn't open input file, exiting...\n")
    sys.exit(1)

try:
    shift_fh = open(args.shift_fn, 'r')
except IOError:
    print ("ERROR Couldn't open file with shift of atoms data, exiting...\n")
    sys.exit(1)

shiftdat=[]
for line in shift_fh:
# skip empty lines and comments
    if (re.search('^\s*$',line)):
        continue
    if (re.search('^\s*#',line)):
        continue
    elif (len(line.split())<4):
        print('The data in shift file wrong format. Required: label x y z')
        sys.exit(1)

    elif (not isNumeric(line.split()[1]) or not isNumeric(line.split()[2]) or not isNumeric(line.split()[3])):
        print('The data in shift file wrong format. Required: label x y z')
        sys.exit(1)

    else:
        shiftdat.append([line.split()[0],float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])

#print(shiftdat)
#loop_   
#_atom_site_label
#_atom_site_type_symbol
#_atom_site_symmetry_multiplicity
#_atom_site_Wyckoff_label        
#_atom_site_fract_x             
#_atom_site_fract_y             
#_atom_site_fract_z             
#_atom_site_occupancy

lines=input_fh.readlines()
ii=0
for i in range(len(lines)):
    if 'loop_' in lines[i]:
        for k in range(1,len(lines)-i):
            if ( not re.search('^\s*_',lines[i+k]) ):
                break
            elif ('loop' in lines[i+k]):
                break
            elif ('_site_fract_x' in lines[i+k]):
                ii=i
                break
    if (ii>0):
        break



if (ii==0):
    print('Error, not found fractional coordinates in cif file. Exit')
    sys.exit(1)

colx=-1
coly=-1
colz=-1
collabel=-1
k=0
for i in range(ii,len(lines)):
    if ( re.search('^\s*_',lines[i]) ):
        k+=1
    if ('_site_fract_x' in lines[i]):
        colx=k-1
    elif ('_site_fract_y' in lines[i]):
        coly=k-1
    elif ('_site_fract_z' in lines[i]):
        colz=k-1
    elif ('_site_label' in lines[i]):
        collabel=k-1
    if (colx>=0 and coly>=0 and colz>=0 and collabel>=0):
#        print("%d %d %d" % (colx, coly, colz))
        break

if (colx<0 or coly<0 or colz<0 or collabel<0):
    print('Error. tag fract coordinates not found')
    sys.exit(1)

if(os.path.isfile(args.output_fn)):
    print('Error output file already exist')
    sys.exit(1)

try:
    out_fh = open(args.output_fn, 'w')
except IOError:
    print ("ERROR Couldn't open input file, exiting...\n")
    sys.exit(1)

for i in range(ii+1):
    out_fh.write(lines[i])

for i in range(ii+1,len(lines)):
    if ( re.search('^\s*_',lines[i]) ):
        out_fh.write(lines[i])
        continue
    if (len(lines[i].split())<max(colx,coly,colz)):
        out_fh.write(lines[i])
        continue
    for k in range(len(lines[i].split())):
        if(k==colx):
            if(isNumeric(lines[i].split()[colx])):
                for j in range(len(shiftdat)):
                    if(shiftdat[j][0].strip()==lines[i].split()[collabel].strip()):
                        break
                out_fh.write('%8.5f ' % (float(lines[i].split()[colx])+shiftdat[j][1]))
            else:
                print('Error. The coordinates is non numeric')
                sys.exit(1)
        elif(k==coly):
            if(isNumeric(lines[i].split()[coly])):
                for j in range(len(shiftdat)):
                    if(shiftdat[j][0].strip()==lines[i].split()[collabel].strip()):
                        break
                out_fh.write('%8.5f ' % (float(lines[i].split()[coly])+shiftdat[j][2]))
            else:
                print('Error. The coordinates is non numeric')
                sys.exit(1)
        elif(k==colz):
            if(isNumeric(lines[i].split()[colz])):
                for j in range(len(shiftdat)):
                    if(shiftdat[j][0].strip()==lines[i].split()[collabel].strip()):
                        break
                out_fh.write('%8.5f ' % (float(lines[i].split()[colz])+shiftdat[j][3]))
            else:
                print('Error. The coordinates is non numeric')
                sys.exit(1)
        else:
            out_fh.write('%s ' % lines[i].split()[k]) 
    out_fh.write('\n')


