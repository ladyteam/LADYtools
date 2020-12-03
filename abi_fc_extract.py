#!/usr/bin/env python3
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
# Extract Dynamical matrix from ABINIT output file and write  FORCE_CONSTANTS file in 
# PHONOPY software format. Aware to use asr variable in ABINIT input file.
#
# Author: Eugene Roginskii
#


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
from math import sqrt
import sys
import re
import argparse

# To produce Dynamic matrix with capabillity of PHONOPY code one need following constants
Ha2Ev=27.2107
Angstr2Bohr=1.8897261245650618

parser = argparse.ArgumentParser(description='''This program is to extract Force Constants from abinit output file in PHONOPY format file.
Please note, the input for the script is the output of abinit calculations, and the number of dataset is the one with rfphon activated.
For example if in input file  rfphon3=1 then the dataset is 3''')

parser.add_argument("-i", "--input", action="store", type=str, dest="abinit_fn",  help="Abinit output filename")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn",  default='FORCE_CONSTANTS', help="Output filename in pnonopy format. Default is FORCE_CONSTANTS")
parser.add_argument("-d", "--dset", action="store", type=int, dest="dset", default=0, help="Dataset number to extract dynmat from")

args = parser.parse_args()

if (args.dset==0):
    print('Error. No Data Set number provided')
    sys.exit(1)

if (args.abinit_fn == None):
    print('Error. No Abinit output filename was given.')
    sys.exit(1)

try:
    abinit_fh = open(args.abinit_fn, 'r')
except IOError:
    print('ERROR opening abinit output file, exiting...')
    sys.exit(1)

for line in abinit_fh:
    if ('Echo of variables that govern the present computation' in line):
        break
natom=0
for line in abinit_fh:
    if 'natom' in line:
        if(isNumeric(line.split()[1])):
            natom=int(line.split()[1])
            break
        else:
            print('Error getting number of atoms value')
            sys.exit(1)

if (natom==0):
    print('Error getting number of atoms value')
    sys.exit(1)

abinit_fh.seek(0)

# Read dynamical matrix from selected dataset
# If variable is provided loop through DATASETS
if(args.dset > 0):
    for line in abinit_fh:
        if ('Echo of variables that govern the present computation') in line:
            break

    for line in abinit_fh:
        if ('dataset' in line.lower()):
            if (re.match('^=+\s*DATASET\s+(\d+).*',line)):
                if (isNumeric(re.match('^=+\s*DATASET\s+(\d+).*',line).group(1))):
                    if(int(re.match('^=+\s*DATASET\s+(\d+).*',line).group(1)) == args.dset):
                        print('Extracting dynmat from DATASET %d' % args.dset)
                        break

for line in abinit_fh:
    if ('dataset' in line.lower()):
        if (re.match('^=+\s*DATASET\s+(\d+).*',line)):
            print('No dynmat in DATASET %d' % args.dset)
            sys.exit(1)
    if ('dynamical matrix, in cartesian coordinates' in line.lower()):
        break


col=0
fc=np.zeros((natom,natom,3,3),dtype='double')

for line in abinit_fh:
    if (len(line.split())<5):
        continue
    elif ((isNumeric(line.split()[4])!=True) or (isNumeric(line.split()[5]) != True) or
                    (isNumeric(line.split()[1])!=True) or (isNumeric(line.split()[3]) != True)):
        continue
    idx=[int(line.split()[i]) for i in range(4)]
    if ((idx[1]<=natom) and (idx[3]<=natom)):
        fc[idx[1]-1][idx[3]-1][idx[0]-1][idx[2]-1]=float(line.split()[4].replace('D', 'E'))*Ha2Ev*Angstr2Bohr
        col+=1
    if (col>=natom**2*9):
        break

try:
    out_fh = open(args.out_fn, 'w')
except IOError:
    print ("ERROR opening output file for writing failed, exiting...")
    sys.exit(1)

out_fh.write('%d\n' % natom)
for n1 in range(natom):
    for n2 in range(natom):
        out_fh.write("%d %d\n" %((n1+1),(n2+1)))
        for i in range(3):
            out_fh.write('  %s\n' % ' '.join('%13.9f' % fc[n1][n2][i][j] for j in range(3)))

