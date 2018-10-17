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
# This program calculates reciprocal vectors and plots the distances in reciprocal
# space between points given in separate file. An finally the program suggest the Number
# of points in each sector to preserve the scale factor.
#
#
# Author: Eugene Roginskii


import sys
from math import pi
import numpy as np
import scipy as sc
from scipy import linalg
import argparse
import re

Angst2Bohr=1.889725989

parser = argparse.ArgumentParser(description='''This program calculates reciprocal vectors and plots
        the distances in reciprocal space''')

parser.add_argument("-i", "--input", action="store", type=str, dest="abinit_fn",  help="Abinit filename")
parser.add_argument("-p", "--points", action="store", type=str, dest="points_fn",  help="Filename with point in reciprocal space provided")
parser.add_argument("-n", "--minpoints", action="store", type=int, dest="minpt", default=6, help="Number of points in minimal section")

args = parser.parse_args()

if (args.abinit_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

try:
    abinit_fh = open(args.abinit_fn, 'r')
except IOError:
    print "ERROR Couldn't open abinit file, exiting...\n"
    sys.exit(1)

# abinit input file process
acell=[]
rprim=[]
scalecart=[1.0,1.0,1.0]
while True:
    line=abinit_fh.readline() 
    if not line: break
    if (re.match('^\s*#',line)):
        continue
    if 'acell' in line:
        acell=[float(line.split()[i+1]) for i in range(3)] 
        print('acell = %s' % acell)
    elif 'rprim' in line:
        try:
            for coor in line.split()[1:]:
                rprim.append(float(coor))
            while True:
                subline=abinit_fh.readline()
                if not subline: break
                if((re.match('^\s*-?\d*\.',subline) is None) and (re.match('^\s*-?\d',subline) is None)):
                    break
                else:
                    for coor in subline.split():
                        rprim.append(float(coor))
        except:
            while True:
                subline=abinit_fh.readline()
                if not subline: break
                if((re.match('^\s*-?\d*\.',subline) is None) and (re.match('^\s*-?\d',subline) is None)):
                    break
                else:
                    for coor in subline.split():
                        rprim.append(float(coor))
    elif 'scalecart' in line:
        scalecart=[float(line.split()[i+1]) for i in range(3)] 

if (len(rprim) == 0):
    rprim=[1,0,0,0,1,0,0,0,1]

rprimd=[]

for i in range(3):
    rprimd.append([scalecart[k]*rprim[i*3+k]*acell[i]/Angst2Bohr for k in range(3)])

basis=np.array(rprimd).reshape(3,3)

print('basis: ')
for i in range(3):
    print('\t'.join('% 6.4f' % basis[i][k] for k in range(3)))

print('scalecart = %s' % scalecart)


#A = sc.mat('[1 3 5; 2 5 1; 2 3 8]')
#old basis

cvol=np.dot(basis[0],np.cross(basis[1],basis[2]))

b1=np.cross(basis[1],basis[2])/cvol*pi
b2=np.cross(basis[2],basis[0])/cvol*pi
b3=np.cross(basis[0],basis[1])/cvol*pi

print('unit cell volume=%7.4f' % cvol)
print('pi units:')
print(b1)
print(b2)
print(b3)
print('reciprocal units:')
print(b1/pi)
print(b2/pi)
print(b3/pi)

if (args.points_fn == None):
    print('Error. No filename with points data was given.')
    sys.exit(1)
try:
    points_fh = open(args.points_fn, 'r')
except IOError:
    print "ERROR Couldn't open file with points data, exiting...\n"
    sys.exit(1)

points=[]
for line in points_fh.readlines():
    if (re.match('^\s*#',line)):
        continue
    if ('#' in line):
        line=re.sub('\s*#.*','',line)
    if (len(line.split()) < 4):
        print('Not enought data in line: %sSkip it' % line)
        continue
    try:
        points.append([[float(line.split()[i]) for i in range(3)],line.split()[3]])
    except:
        print('Error reading points data')

vec=[]
for i in range(len(points)-1):
    vec1=points[i][0][0]*b1+points[i][0][1]*b2+points[i][0][2]*b3
    vec2=points[i+1][0][0]*b1+points[i+1][0][1]*b2+points[i+1][0][2]*b3

    vec.append(vec1-vec2)

dist=[]
print('Distances:')
for i in range(len(points)-1):
    dist.append(np.linalg.norm(vec[i]))
    print('Vector between points %s-%s %s, norm: %5.3f' % (points[i][1],points[i+1][1]," ".join('% 5.3f' % v for v in vec[i]), dist[i])  )

mind=dist[0]
for d in dist:
    if(mind>d):
        mind=d

mult=mind/float(args.minpt)

print('Suggested number in each section is:')
ntot=0
for i in range(len(points)-1):
    npt=dist[i]/mult
    print('Section %s-%s, npt: %d' % (points[i][1],points[i+1][1], int(round(npt)))  )
    ntot+=int(round(npt))
print('Total number of points: %d' % ntot)
