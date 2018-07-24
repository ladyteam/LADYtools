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
# Exctract raman tensor from anaddb utility of ABINIT package
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

parser = argparse.ArgumentParser(description="Extract raman tensor and rotational invariants (Long's notation) from ABINIT anaddb output")
parser.add_argument('-i', '--input', required=True, type=argparse.FileType('r'), help='input file name')
parser.add_argument('-o', '--output', required=True, type=argparse.FileType('w'), help='output file name')
args = parser.parse_args()


raman=[]
freqs=[]


for line in args.input:
    if('Phonon frequencies in cm-1' in line):
        break

for line in args.input:
    if('Eigendisplacements' in line):
        break

    if(line[0] == '-'):
        for i in range(len(line.split())-1):
            freqs.append(float(line.split()[i+1]))
    else:
        break

for line in args.input:
    if ("Raman susceptibilities of transverse zone-center phonon modes" in line):
        break


for line in args.input:
    if('Electronic dielectric tensor' in line):
        break
    if('Mod' in line):
        freq=float(line.split('(')[1].split()[0])
        buf=[]
        for i in range(3):
            line=args.input.next()[1:]
            for j in range(3):
                buf.append(float(line.split()[j]))
        raman.append([freq,[buf[i] for i in range(9)]])


args.output.write("#Raman susceptibilities of transverse zone-center phonon modes, Intensities multiplied by 10^5\n")
args.output.write("#nu,cm-1   nu_asm       xx           xy           xz           yx           yy           yz           zx           zy           zz         alpha       gamma2     Ipar    Iperp    Itot\n")

for i in range(len(raman)):
    args.output.write("% 8.3f % 8.3f %s" % (raman[i][0], freqs[i], "".join("% 12.8f " % r for r in raman[i][1] )))
    alpha=(raman[i][1][0]+raman[i][1][4]+raman[i][1][8])/3
    gamma2=((raman[i][1][0]-raman[i][1][4])**2 + (raman[i][1][4]-raman[i][1][8])**2 +(raman[i][1][8]-raman[i][1][0])**2)/2
    gamma2+=3*(raman[i][1][1]**2+raman[i][1][2]**2+raman[i][1][5]**2)
    Iper=(45*alpha**2+4*gamma2)/45
    Ipar=gamma2/15
    args.output.write("% 10.8f % 10.8f % 8.5f % 8.5f % 8.5f\n" % (alpha,gamma2,(Iper*1e05),(Ipar*1e05),((Ipar+Iper)*1e05)))


# Next find LO intencity
args.input.seek(0)
LOexist=False

for line in args.input:
    if ("Raman susceptibility of zone-center phonons, with non-analyticity" in line):
        LOexist=True
        break

raman=[]
if (LOexist):
    args.output.write("#Raman susceptibilities of zone-center phonons, with non-analyticity in direction %s" % args.input.next().split(")")[1] )
    args.output.write("#nu,cm-1   nu_asm       xx           xy           xz           yx           yy           yz           zx           zy           zz         alpha       gamma2     Ipar    Iperp    Itot\n")
    for line in args.input:
        if('Electronic dielectric tensor' in line):
            break
        if('Mod' in line):
            freq=float(line.split('(')[1].split()[0])
            buf=[]
            for i in range(3):
                line=args.input.next()[1:]
                for j in range(3):
                    buf.append(float(line.split()[j]))
            raman.append([freq,[buf[i] for i in range(9)]])

    for i in range(len(raman)):
        args.output.write("% 8.3f % 8.3f %s" % (raman[i][0], freqs[i], "".join("% 12.8f " % r for r in raman[i][1] )))
        alpha=(raman[i][1][0]+raman[i][1][4]+raman[i][1][8])/3
        gamma2=((raman[i][1][0]-raman[i][1][4])**2 + (raman[i][1][4]-raman[i][1][8])**2 +(raman[i][1][8]-raman[i][1][0])**2)/2
        gamma2+=3*(raman[i][1][1]**2+raman[i][1][2]**2+raman[i][1][5]**2)
        Iper=(45*alpha**2+4*gamma2)/45
        Ipar=gamma2/15
        args.output.write("% 10.8f % 10.8f % 8.5f % 8.5f % 8.5f\n" % (alpha,gamma2,(Iper*1e05),(Ipar*1e05),((Ipar+Iper)*1e05)))



