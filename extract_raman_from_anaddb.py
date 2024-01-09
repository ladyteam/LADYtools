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
# Extract Raman tensor from anaddb output file and calculate intencity for powder
#
# Author: Eugene Roginskii

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

import re
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Extract Raman tensor from anaddb output and plot table in output file')

parser.add_argument("-i", "--input", action="store", type=str, dest="in_fn",  help="Anaddb output filename")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn",   help="Output filename")

args = parser.parse_args()

if (args.in_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

try:
    in_fh = open(args.in_fn, 'r')
except IOError:
    print("ERROR Couldn't open input file, exiting...\n")
    sys.exit(1)

try:
    out_fh = open(args.out_fn, 'w')
except IOError:
    print("ERROR Couldn't open output file, exiting...\n")
    sys.exit(1)


for line in in_fh:
    if not line: break
    if (re.match('^\s*#',line)):
        continue
    if not line.strip():
        continue
    if ('#' in line):
        line=re.sub('\s*#.*','',line)
    if ('Raman susceptibilities of' in line):
        out_fh.write('# %s' % line)
        out_fh.write('# %s' % in_fh.readline())
        out_fh.write('#N    freq           xx            xy            xz            yx            yy            yz           zx            zy            zz          Alpha         Gamma2         Ipar         Iperp         Itot\n')
    if ('Raman susceptibility of' in line):
        out_fh.write('# %s' % line)
        out_fh.write('# %s' % in_fh.readline())
        out_fh.write('#N    freq           xx            xy            xz            yx            yy            yz           zx            zy            zz          Alpha         Gamma2         Ipar         Iperp         Itot\n')

    m=re.match('\s*Mod\s?\s?\s?(\d+)\s+\(\s*([+|-]?\d+\.\d+)\s*cm-1',line)
    if (m):
        raman=[]
        for j in range(3):
            buf=in_fh.readline()
            raman.append([float(buf.split()[i+1]) for i in range(3)])
        alpha=1.0/3*(raman[0][0]+raman[1][1]+raman[2][2])
        gamma2=1.0/2*((raman[0][0]-raman[1][1])**2+(raman[1][1]-raman[2][2])**2+(raman[0][0]-raman[2][2])**2)+\
                3.0/4*((raman[0][1]+raman[1][0])**2+(raman[0][2]+raman[2][0])**2+(raman[1][2]+raman[2][1])**2)
        Iperp=gamma2/15
        Iparal=(45*alpha**2+4*gamma2)/45
        Itot=Iparal+Iperp
        out_fh.write('%3d % 8.2f  ' % (int(m.group(1)),float(m.group(2))))
        for i in range(3):
            out_fh.write(''.join('% 12.10f ' % raman[i][j] for j in range(3)))
        out_fh.write('% 12.10f % 12.10f % 12.10f % 12.10f %12.9f\n'% (alpha**2,gamma2,Iparal,Iperp,Itot))

