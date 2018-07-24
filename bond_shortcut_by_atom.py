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
# Add atoms into the file with xyz format to shortcut open bonds.
# The additional file with parameters of newly added files is needed in 
# the following format each line contains four fields separated by spaces:
# N1 N2 AtomSymbol BondLength
# where N1 and N2 - (integer) is the number of first and second  atom in selected bond correspondingly
# AtomSymbol - (symbol) the symbol of atom to add
# BondLength - the length of newly generated bond between new atom and atom N2
#
# Author: Eugene Roginskii
#


import re
from math import sqrt
import sys
import os
import argparse

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



def vlength(a,b,c,t,dist):
    return ( sqrt((a*(t-1))**2 + (b*(t-1))**2 + (c*(t-1))**2) - dist )


def main():


    parser = argparse.ArgumentParser(description='Shortcut open bonds by atoms')
#    parser = OptionParser()

    parser.add_argument("-i", "--input", action="store", type=str, dest="data_fn",  help="Input data file name")
    parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn", default='out.dat', help="Output data file name")
    parser.add_argument("-p", "--param", action="store", type=str, dest="param_fn",  help="Filename with bonds to shurtcut data")

    args = parser.parse_args() 



    try:
        data_fh = open(args.data_fn, 'r')
    except IOError:
        print "ERROR Couldn't open data file, exiting...\n"
        sys.exit(1)

    try:
        out_fh = open(args.out_fn, 'w')
    except IOError:
        print "ERROR Couldn't open output file, exiting...\n"
        sys.exit(1)

    atoms=[]

    line=data_fh.readline()
    if (strType(line.split()[0])=='int'):
        out_fh.write("%6d\n" % int(line.split()[0]))

    else:
        print('Wrong format of xyz input file. First row should provide number of atoms')
        sys.exit(1)

    for line in data_fh:
        out_fh.write(line)
        print(line)
        if (len(line.split())<4):
            continue
        atoms.append([line.split()[0],float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])

    try:
        param_fh = open(args.param_fn, 'r')
    except IOError:
        print "ERROR Couldn't open param file, exiting...\n"
        sys.exit(1)

    nadat=0
    for line in param_fh:
        if (re.search('^\s*$',line)):
            continue
        if (re.search('^\s*#',line)):
            continue

        nadat+=1
        if (len(line.split())<4):
            print('Not enough data in parameter file for atom %d provided. Ignoring')
            continue
        pt1=[atoms[int(line.split()[0])-1][i] for i in range(1,4)]
        pt2=[atoms[int(line.split()[1])-1][i] for i in range(1,4)]
        label=line.split()[2]
        dist=float(line.split()[3])
        print(pt1,pt2,label,dist)
        a=pt2[0]-pt1[0]
        b=pt2[1]-pt1[1]
        c=pt2[2]-pt1[2]
#        m=1000
#        step=0
#        for t in range(1,300):
#            print(vlength(a,b,c,(t/100.0),dist),t)
#            if(abs(vlength(a,b,c,(t/100.0),dist))<m):
#                step=t
#                m=vlength(a,b,c,(t/100.0),dist)
        t=dist/sqrt(a**2+b**2+c**2)

        bl=sqrt(a**2+b**2+c**2) # Bond length
        print(a,b,c)
        out_fh.write("%s % 9.7f % 9.7f % 9.7f\n" % (label, pt2[0]+a*t, pt2[1]+b*t, pt2[2]+c*t ))

    out_fh.seek(0)
    out_fh.write("%6d\n" % (len(atoms)+nadat))

    #    x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function
#        x0=f.getparam()
#        finpars, ier = leastsq(objective, x0, args=(daty,datx)) #this is from scipy
#        print("".join("% E " %f for f in finpars))

if __name__ == '__main__':                                                                                                                                      
    main()
