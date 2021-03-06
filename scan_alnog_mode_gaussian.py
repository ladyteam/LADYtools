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
# The program to distort geometry along selected normal mode eigenvector. Abinit input format
#
#
# Author: Eugene Roginskii
#
# Some code was borrowed from very powerfool PHONOPY project
# http://phonopy.sourceforge.net/
# Copyright 2009, Atsushi Togo.

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



def genxyz(fn,modenum,cartshift,comment):
    out_fh=open(fn,'w')
    out_fh.write('Atom positions shifted along Mode %d %s\n\n0 1\n' % (modenum,comment))

    for cartat in cartshift:
        print(cartat)
        out_fh.write('%s    % 12.9f % 12.9f % 12.9f\n' % (cartat[0],cartat[1],cartat[2],cartat[3]))


import numpy as np
import re
from math import sqrt
import sys
from math import pi
from math import exp
from shutil import move
import os
import datetime
import time
import argparse
import xml.etree.cElementTree as etree


mau=1822.888485
AMU = 1.6605402e-27 # [kg]

# ABINIT freq to cm-1 factor
factorcm=716.851928469
# ABINIT freq to HZ factor
factorHz=21.49068e12

hbar=6.6260695729e-34 #J*s
#hbar=4.13566751691e−15 #eV*s
EV = 1.60217733e-19 # [J]
#
Angst2Bohr=1.889725989

atom_data = [ 
    [  0, "X", "X", 0], # 0
    [  1, "H", "Hydrogen", 1.00794], # 1
    [  2, "He", "Helium", 4.002602], # 2
    [  3, "Li", "Lithium", 6.941], # 3
    [  4, "Be", "Beryllium", 9.012182], # 4
    [  5, "B", "Boron", 10.811], # 5
    [  6, "C", "Carbon", 12.0107], # 6
    [  7, "N", "Nitrogen", 14.0067], # 7
    [  8, "O", "Oxygen", 15.9994], # 8
    [  9, "F", "Fluorine", 18.9984032], # 9
    [ 10, "Ne", "Neon", 20.1797], # 10
    [ 11, "Na", "Sodium", 22.98976928], # 11
    [ 12, "Mg", "Magnesium", 24.3050], # 12
    [ 13, "Al", "Aluminium", 26.9815386], # 13
    [ 14, "Si", "Silicon", 28.0855], # 14
    [ 15, "P", "Phosphorus", 30.973762], # 15
    [ 16, "S", "Sulfur", 32.065], # 16
    [ 17, "Cl", "Chlorine", 35.453], # 17
    [ 18, "Ar", "Argon", 39.948], # 18
    [ 19, "K", "Potassium", 39.0983], # 19
    [ 20, "Ca", "Calcium", 40.078], # 20
    [ 21, "Sc", "Scandium", 44.955912], # 21
    [ 22, "Ti", "Titanium", 47.867], # 22
    [ 23, "V", "Vanadium", 50.9415], # 23
    [ 24, "Cr", "Chromium", 51.9961], # 24
    [ 25, "Mn", "Manganese", 54.938045], # 25
    [ 26, "Fe", "Iron", 55.845], # 26
    [ 27, "Co", "Cobalt", 58.933195], # 27
    [ 28, "Ni", "Nickel", 58.6934], # 28
    [ 29, "Cu", "Copper", 63.546], # 29
    [ 30, "Zn", "Zinc", 65.38], # 30
    [ 31, "Ga", "Gallium", 69.723], # 31
    [ 32, "Ge", "Germanium", 72.64], # 32
    [ 33, "As", "Arsenic", 74.92160], # 33
    [ 34, "Se", "Selenium", 78.96], # 34
    [ 35, "Br", "Bromine", 79.904], # 35
    [ 36, "Kr", "Krypton", 83.798], # 36
    [ 37, "Rb", "Rubidium", 85.4678], # 37
    [ 38, "Sr", "Strontium", 87.62], # 38
    [ 39, "Y", "Yttrium", 88.90585], # 39
    [ 40, "Zr", "Zirconium", 91.224], # 40
    [ 41, "Nb", "Niobium", 92.90638], # 41
    [ 42, "Mo", "Molybdenum", 95.96], # 42
    [ 43, "Tc", "Technetium", 0], # 43
    [ 44, "Ru", "Ruthenium", 101.07], # 44
    [ 45, "Rh", "Rhodium", 102.90550], # 45
    [ 46, "Pd", "Palladium", 106.42], # 46
    [ 47, "Ag", "Silver", 107.8682], # 47
    [ 48, "Cd", "Cadmium", 112.411], # 48
    [ 49, "In", "Indium", 114.818], # 49
    [ 50, "Sn", "Tin", 118.710], # 50
    [ 51, "Sb", "Antimony", 121.760], # 51
    [ 52, "Te", "Tellurium", 127.60], # 52
    [ 53, "I", "Iodine", 126.90447], # 53
    [ 54, "Xe", "Xenon", 131.293], # 54
    [ 55, "Cs", "Caesium", 132.9054519], # 55
    [ 56, "Ba", "Barium", 137.327], # 56
    [ 57, "La", "Lanthanum", 138.90547], # 57
    [ 58, "Ce", "Cerium", 140.116], # 58
    [ 59, "Pr", "Praseodymium", 140.90765], # 59
    [ 60, "Nd", "Neodymium", 144.242], # 60
    [ 61, "Pm", "Promethium", 0], # 61
    [ 62, "Sm", "Samarium", 150.36], # 62
    [ 63, "Eu", "Europium", 151.964], # 63
    [ 64, "Gd", "Gadolinium", 157.25], # 64
    [ 65, "Tb", "Terbium", 158.92535], # 65
    [ 66, "Dy", "Dysprosium", 162.500], # 66
    [ 67, "Ho", "Holmium", 164.93032], # 67
    [ 68, "Er", "Erbium", 167.259], # 68
    [ 69, "Tm", "Thulium", 168.93421], # 69
    [ 70, "Yb", "Ytterbium", 173.054], # 70
    [ 71, "Lu", "Lutetium", 174.9668], # 71
    [ 72, "Hf", "Hafnium", 178.49], # 72
    [ 73, "Ta", "Tantalum", 180.94788], # 73
    [ 74, "W", "Tungsten", 183.84], # 74
    [ 75, "Re", "Rhenium", 186.207], # 75
    [ 76, "Os", "Osmium", 190.23], # 76
    [ 77, "Ir", "Iridium", 192.217], # 77
    [ 78, "Pt", "Platinum", 195.084], # 78
    [ 79, "Au", "Gold", 196.966569], # 79
    [ 80, "Hg", "Mercury", 200.59], # 80
    [ 81, "Tl", "Thallium", 204.3833], # 81
    [ 82, "Pb", "Lead", 207.2], # 82
    [ 83, "Bi", "Bismuth", 208.98040], # 83
    [ 84, "Po", "Polonium", 0], # 84
    [ 85, "At", "Astatine", 0], # 85
    [ 86, "Rn", "Radon", 0], # 86
    [ 87, "Fr", "Francium", 0], # 87
    [ 88, "Ra", "Radium", 0], # 88
    [ 89, "Ac", "Actinium", 0], # 89
    [ 90, "Th", "Thorium", 232.03806], # 90
    [ 91, "Pa", "Protactinium", 231.03588], # 91
    [ 92, "U", "Uranium", 238.02891], # 92
    [ 93, "Np", "Neptunium", 0], # 93
    [ 94, "Pu", "Plutonium", 0], # 94
    [ 95, "Am", "Americium", 0], # 95
    [ 96, "Cm", "Curium", 0], # 96
    [ 97, "Bk", "Berkelium", 0], # 97
    [ 98, "Cf", "Californium", 0], # 98
    [ 99, "Es", "Einsteinium", 0], # 99
    [100, "Fm", "Fermium", 0], # 100
    [101, "Md", "Mendelevium", 0], # 101
    [102, "No", "Nobelium", 0], # 102
    [103, "Lr", "Lawrencium", 0], # 103
    [104, "Rf", "Rutherfordium", 0], # 104
    [105, "Db", "Dubnium", 0], # 105
    [106, "Sg", "Seaborgium", 0], # 106
    [107, "Bh", "Bohrium", 0], # 107
    [108, "Hs", "Hassium", 0], # 108
    [109, "Mt", "Meitnerium", 0], # 109
    [110, "Ds", "Darmstadtium", 0], # 110
    [111, "Rg", "Roentgenium", 0], # 111
    [112, "Cn", "Copernicium", 0], # 112
    [113, "Uut", "Ununtrium", 0], # 113
    [114, "Uuq", "Ununquadium", 0], # 114
    [115, "Uup", "Ununpentium", 0], # 115
    [116, "Uuh", "Ununhexium", 0], # 116
    [117, "Uus", "Ununseptium", 0], # 117
    [118, "Uuo", "Ununoctium", 0], # 118
    ]



parser = argparse.ArgumentParser(description='The program to distort geometry along selected normal  mode eigenvector. Gaussian output format')


parser.add_argument("-i", "--input", action="store", type=str, dest="gaus_fn", help="Gaussian output filename")
parser.add_argument("-m", "--mode", action="store", type=int, dest="modenum", default=1, help="Number of mode")
parser.add_argument("-a", "--amplitudes", action="store", type=str, dest="amplstr", default='0.1 0.25 0.5 1', 
           help="Mode shift  amplitude, space separated")

args = parser.parse_args()

if (args.gaus_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

try:
    gaus_fh = open(args.gaus_fn, 'r')
except IOError:
    print("ERROR Couldn't open gaussian output file, exiting...")
    sys.exit(1)


ampl=[]
for i in range(len(args.amplstr.split())):
    if isNumeric(args.amplstr.split()[i]):
        ampl.append(float(args.amplstr.split()[i]))
    else:
        print('Error parsing amplitudes array')
        sys.exit(1)

modenum=args.modenum

natom=0

# Read Initial atomic positions
for line in gaus_fh:
    if("Standard orientation:" in line):
        break
    if("Input orientation:" in line):
        break

coor=[]
gaus_fh.readline()
gaus_fh.readline()
gaus_fh.readline()
gaus_fh.readline()

for line in gaus_fh:
    if(" ------------------" in line):
        break

    m=re.search('\s+(\d+)\s+(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)',line)
    if(m):
        coor.append([atom_data[int(m.group(2))][1],float(m.group(4)),float(m.group(5)),float(m.group(6))])

print(coor)

eigen=[]
# Read Eigenvalues and Eigenvectors
# skip till the line Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering 
for line in gaus_fh:
    if("Harmonic frequencies (cm**-1)," in line):
        break
stop=0
for j in range(len(coor)):
    if (stop==1):
        break
    freq=[]
    e=[]
    atnum=1
    for line in gaus_fh:
        if("Harmonic frequencies (cm**-1)," in line):
            stop=1
            break
        if ("Frequencies" in line):
            for i in range(len(line.split())-2):
                freq.append(float(line.split()[i+2]))
            continue
        m=re.search('\s+(\d+)\s+(\d+)\s+((\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s*)+)',line)
        if(m):
            buf=[]
            for k in range(int(len(m.group(3).split())/3)):
                buf.append([float(m.group(3).split()[l]) for l in range(3*k,3*(k+1))])

            e.append([buf[k] for k in range(len(buf))])
            atnum+=1
        if(atnum>len(coor)):
# iteration on each set of mode (3 max) and then each atoms DISPLACEMENTS
            for k in range(len(freq)):
                etmp=[]
                for l in range(len(coor)):
                    etmp.append(e[l][k])
                eigen.append([freq[k],etmp])
            break
#    for e in eigen:
#        print(e[0])
#        for ee in e[1]:
#            print(ee)


# GENERATION OF DISPLACEMENTS XYZ INPUT FILES
#j - mode number; i -atom number
j=args.modenum-1
for n in range(len(ampl)):
    cartshift=[]

    for i in range(len(coor)):
        shiftvec=[0.0e0,0.0e0,0.0e0]
#        print(eigen[j][1][i])
        for l in range(3):
            shiftvec[l]=ampl[n]*eigen[j][1][i][l]

        cartshift.append([coor[i][0],coor[i][1]-shiftvec[0], coor[i][2]-shiftvec[1], coor[i][3]-shiftvec[2]])

    outfn="shiftcell-%.3f.xyz" % ampl[n]

    genxyz(outfn,j+1,cartshift,'freq = %9.4f cm-1; Delta=%6.4f' % ((eigen[j][0]),ampl[n]))

