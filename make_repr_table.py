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
# Exctract character table from anaddb utility of ABINIT package
#
# Author: Eugene Roginskii
#

import sys
import string
import re
import argparse
import numpy as np
import scipy as sc
from scipy import linalg

parser = argparse.ArgumentParser(description='Plot a table of characters')
parser.add_argument('-i', '--input', required=True, type=argparse.FileType('r'), help='input file name')
parser.add_argument('-o', '--output', required=True, type=argparse.FileType('w'), help='output file name')
args = parser.parse_args()


for line in args.input:
    if (" Phonon frequencies in cm-1" in line):
        break
freq=[]
for line in args.input:
    if('-' == line[0]):
        for i in range(1,len(line.split())):
            print(line.split()[i])
            freq.append(float(line.split()[i]))
    else:
        break

args.input.seek(0)
# Analysis of degeneracies and characters (maximum tolerance=1.00E-06 a.u.)                                                                                      
#  Symmetry characters of vibration mode #   1                                                                                                                    
#         degenerate with vibration modes #   2 to    3 
modedeg=0

#init character table
chrt=[]

for line in args.input:
    if ('Analysis of degeneracies and characters' in line):
        break

k=len(freq)
for spgline in args.input:
    if (spgline.find("=======================") != -1):
      break
    if(k<=0):
        break
    num=re.search('^\s*Symmetry characters of vibration mode #\s+(\d+)', spgline)
    if (num != None):
      mode=int(num.group(1))
      next
    numdeg=re.search('^\s*degenerate with vibration modes? #\s+(\d+)(\s+to\s+(\d+))?', spgline)
    if (numdeg != None):
      if(numdeg.group(3) != None):
          modedeg=int(numdeg.group(3))
      else:
          modedeg=int(numdeg.group(1))
      next
# read line of characters
    char=re.search('^\s*([+|-]?\d\.\d)+', spgline)
    if (char != None):
        if (modedeg !=0):
            chrt.append([freq[mode-1],(modedeg-mode+1),spgline.rstrip('\n')])
            print('%d-%d %s' % (mode, modedeg, spgline.rstrip('\n')))
            k-=(modedeg-mode+1)
            modedeg=0
        else:
            chrt.append([freq[mode-1],1,spgline.rstrip('\n')])
            k-=1
            print('%d %s' % (mode, spgline.rstrip('\n')))

#print(chrt)
args.output.write("# N freq,cm-1  deg    Characters\n")

if(len(chrt) <= 0):
    for i in range(len(freq)):
        args.output.write("%3d %8.3f\n" % (i+1, freq[i]))

else:
    for i in range(len(chrt)):
        args.output.write("%3d %8.3f %4d %s\n" % (i+1, chrt[i][0], chrt[i][1], chrt[i][2]))
