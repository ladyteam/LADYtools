#!/usr/bin/python
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
# Fit Energy dependence on volume by Murnaghan equation
#
# Author: Eugene Roginskii

#from pylab import * #this includes numpy as np!

from scipy.optimize import leastsq
import numpy as np
import fileinput
import re
import pyx as px
import argparse


parser = argparse.ArgumentParser(description='''Fit Energy dependence on volume by Murnaghan equation.
                                                 Multiplyer calculated at the end.''')

parser.add_argument("-i", "--input", action="store", type=str, dest="dfn",  help='''Filename
        with data E(V)''')
parser.add_argument("-p", "--poscar", action="store", type=str, dest="poscar_fn", default='POSCAR',
        help="POSCAR fullpath")

args = parser.parse_args()

if (args.dfn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

try:
    dfh = open(args.dfn, 'r')
except IOError:
    print "ERROR Couldn't open init file, exiting...\n"
    sys.exit(1)

try:
    psfh = open(args.poscar_fn, 'r')
except IOError:
    print "ERROR Couldn't open POSCAR file, exiting...\n"
    sys.exit(1)

rprim=[]
psfh.readline()
psfh.readline()
for j in range(3):
    line=psfh.readline()
    rprim.append([float(line.split()[i]) for i in range(3)])

basis=np.array(rprim).reshape(3,3)
cvol=np.dot(basis[0],np.cross(basis[1],basis[2]))

# raw data from 2.2.3-al-analyze-eos.py
#v = np.array([13.72, 14.83, 16.0, 17.23, 18.52])
#e = np.array([-56.29, -56.41, -56.46, -56.46, -56.42])
#read data from stdin or file
v=np.array([])
e=np.array([])
for line in dfh:
    if (re.match('^\s*#',line)):
        continue
    if (re.search('\d+.*\s+[+|-]?\d+', line) != None):
        vals=line.split()
        v=np.append(v,float(vals[0]))
        e=np.append(e,float(vals[1]))


#make a vector to evaluate fits on with a lot of points so it looks smooth
vfit = np.linspace(min(v),max(v),100)

### fit a parabola to the data
# y = ax^2 + bx + c
a,b,c = np.polyfit(v,e,2) #this is from pylab

'''
the parabola does not fit the data very well, but we can use it to get
some analytical guesses for other parameters.

V0 = minimum energy volume, or where dE/dV=0
E = aV^2 + bV + c
dE/dV = 2aV + b = 0
V0 = -b/2a

E0 is the minimum energy, which is:
E0 = aV0^2 + bV0 + c

B is equal to V0*d^2E/dV^2, which is just 2a*V0

and from experience we know Bprime_0 is usually a small number like 4
'''

#now here are our initial guesses.
v0 = -b/(2*a)
e0 = a*v0**2 + b*v0 + c
b0 = 2*a*v0
bP = 4

#now we have to create the equation of state function
def Murnaghan(parameters,vol):
    '''
    given a vector of parameters and volumes, return a vector of energies.
    equation From PRB 28,5480 (1983)
    '''
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    
    E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)

    return E

# and we define an objective function that will be minimized
def objective(pars,y,x):
    #we will minimize this function
    err =  y - Murnaghan(pars,x)
    return err

x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function

murnpars, ier = leastsq(objective, x0, args=(e,v)) #this is from scipy

#now we make a figure summarizing the results
y1=np.array([])
y2=np.array([])
for x in vfit.tolist():
    y1=np.append(y1,a*x**2 + b*x + c)
    y2=np.append(y2,Murnaghan(murnpars,x))

px.text.set(mode="latex")
g=px.graph.graphxy(width=8,key=px.graph.key.key(pos="tr", dist=0.1),
                   x=px.graph.axis.linear(title=r"Volume, $\mbox{\normalfont\AA}^3$"),
                   y=px.graph.axis.linear(title="E, eV"))

data1=px.graph.data.values(title="parabolic fit", x=vfit.tolist(), y=y1.tolist())
data2=px.graph.data.values(title="Murnanghan fit", x=vfit.tolist(), y=y2.tolist())
data3=px.graph.data.values(title="Calculation", x=v.tolist(), y=e.tolist())
g.plot(data1,
    styles=[px.graph.style.line([px.color.rgb.blue, px.style.linestyle.dashed, px.style.linewidth.thick])])
g.plot(data2,
    styles=[px.graph.style.line([px.color.rgb.red, px.style.linestyle.solid, px.style.linewidth.thick])])
g.plot(data3,
    styles=[px.graph.style.symbol(px.graph.style.symbol.circle, symbolattrs=[px.deco.filled,px.color.rgb.green])])

g.text(g.width/2, g.height + 0.2, "Optimal Volume Murnaghan equation", 
       [px.text.halign.center, px.text.valign.bottom, px.text.size.Large])

g.text(0.2,0.2, r'\parbox{14em}{Min Volume = %1.5f $\mbox{\normalfont\AA}^3$ \\ Bulk modulus = %1.2f GPa }' % (murnpars[3], murnpars[1]*160.21773),
       [px.text.halign.left, px.text.valign.bottom, px.text.size.small])
g.writePDFfile("fit")

#plot(v,e,'ro')
#plot(vfit, a*vfit**2 + b*vfit + c,'--',label='parabolic fit')
#plot(vfit, Murnaghan(murnpars,vfit), label='Murnaghan fit')
#xlabel('Volume ($\AA^3$)')
#ylabel('Energy (eV)')
#legend(loc='best')

#add some text to the figure in figure coordinates
#ax = gca()
#text(0.4,0.5,'Min volume = %1.2f $\AA^3$' % murnpars[3],
#     transform = ax.transAxes)
#text(0.4,0.4,'Bulk modulus = %1.2f eV/$\AA^3$ = %1.2f GPa' % (murnpars[1],
#                                                             murnpars[1]*160.21773)
#    , transform = ax.transAxes)
#savefig('a-eos.png')
#show()


print('initial guesses                    : ',x0)
print('fitted parameters(V0 is last param): ', murnpars)
print('%f' % (float(murnpars.tolist()[3])/1.0))
print('Ucell volume: % 12.8f Multiplyer: % 12.8f' %  (cvol, (float(murnpars.tolist()[3])/cvol)**(1.0/3)))