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
# Using character table in format of make_repr_table.py script generated from anaddb 
# utility of ABINIT package and crystallographic data in format (comment out first string,
# only representation names and characters allowed): 
#
# The head of file generated with make_repr_table.py Should be modifyed by adding 
# crystallographic data included inside block started by "REPRES_START" and ended with
# REPRES_STOP. Symmetry operations order in crystallographic block should be exactly the same 
# as anaddb output file character table. The symmetry operations used by abinit could be translated 
# from matrix form to international symbols with PHONOPY plot representation feature applied to
# FC_CONSTANTS file extracted from abinit output with ladytools utils.
#
# Example of crystallographic block:
# D2h(mmm)  #   1  -1   2y   my   2z  mz   2x   mx
#
#REPRES_START
# Ag    Γ1+     1   1    1   1    1   1    1    1 
# B1g   Γ3+     1   1   -1  -1    1   1   -1   -1 
# B2g   Γ2+     1   1    1   1   -1  -1   -1   -1 
# B3g   Γ4+     1   1   -1  -1   -1  -1    1    1 
# Au    Γ1-     1  -1    1  -1    1  -1    1   -1 
# B1u   Γ3-     1  -1   -1   1    1  -1   -1    1 
# B2u   Γ2-     1  -1    1  -1   -1   1   -1    1 
# B3u   Γ4-     1  -1   -1   1   -1   1    1   -1 
#REPRES_STOP
#
#  Plot representations of phonons
#
# Author: Eugene Roginskii
#

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

def cmpvec(vec1,vec2,accuracy):
    """ Compare vectors at a given accuracy
     -1 -- wrong vector size
     0 are not equivavent
     1 are equivavent 

    """

    if((len(vec1)<3) or (len(vec2)<3)):
        return(-1)
    else:
        k=0
        for i in range(len(vec1)):
            k+=abs(vec1[i]-vec2[i])
        if(k>accuracy):
            return(0)
        else:
            return(1)

import sys
import string
import re
import argparse
import numpy as np
import scipy as sc
from scipy import linalg

parser = argparse.ArgumentParser(description='Plot a table of Irreducible Representations using a given characters table. Read comments inside for instructions how to use.')
parser.add_argument('-i', '--input', required=True, type=str, dest="in_fn", help='input file name')
parser.add_argument('-o', '--output', required=True, type=str, dest="out_fn", help='output file name')

args = parser.parse_args()

if (args.in_fn == None):
    print('Error. No input file name was provided')
    sys.exit(1)

if (args.out_fn == None):
    print('Error. No output file name was provided')
    sys.exit(1)

try:
    in_fh = open(args.in_fn, 'r')
except IOError:
    print("ERROR Couldn't open abinit file, exiting...\n")
    sys.exit(1)

cdat=[]
for line in in_fh:
    if(re.match('^\s*#',line)):
        continue
    if(re.search('^\s*$', line)):
        continue
    if('REPRES_START' in line):
        break

for line in in_fh:
    if(re.match('^\s*#',line)):
        continue
    if(re.search('^\s*$', line)):
        continue
    if('REPRES_STOP' in line):
        break
    if(isNumeric(line.split()[2])):
        cdat.append([line.split()[0],[float(line.split()[i]) for i in range(2,len(line.split()))]])
for c in cdat:
    print(c)

try:
    out_fh = open(args.out_fn, 'w')
except IOError:
    print("ERROR Couldn't open abinit file, exiting...\n")
    sys.exit(1)

out_fh.write("  N  freq,cm-1  IR     Characters\n")

for line in in_fh:
    curir=""
    if(re.match('^\s*#',line)):
        continue
    if(re.search('^\s*$', line)):
        continue
    if(isNumeric(line.split()[3])):
        print("=== Mode %s ===" % line.split()[0])
        buf=[float(line.split()[i]) for i in range(3,len(line.split()))]
        print("%s %s" % (line.split()[1], "".join("%s " % buf[i] for i in range(len(buf)))))
        for i in range(len(cdat)):
            if(cmpvec(cdat[i][1],buf,0.01) == 1):
                print("Found chars: %s" % "".join("%3.1f " % cdat[i][1][k] for k in range(len(cdat[i][1]))))
                curir=cdat[i][0]
                break
        print("%10s  %4s" % (line.split()[1],curir))
        out_fh.write("%4d %8.2f %4s   %s\n" %( int(line.split()[0]), float(line.split()[1]),  curir, "".join("% 4.1f " % cdat[i][1][k] for k in range(len(cdat[i][1])))  ))
