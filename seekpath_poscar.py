#!/usr/bin/env python
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
# The program to plot k-path for band structure using seekpath python library and input in POSCAR format
# 
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

# Compare vectors at a given accuracy
# -1 -- wrong vector size
# 0 are not equivavent
# 1 are equivavent
def cmpvec(vec1,vec2,accuracy):
    if((len(vec1)<3) or (len(vec2)<3)):
        return(-1)
    else:
        k=0
        for i in range(3):
            k+=abs(vec1[i]-vec2[i])
        if(k>accuracy):
            return(0)
        else:
            return(1)
def direct2cart(directpos,basis):
    cart=[]
    for atdirect in directpos:
        cart.append(np.dot(np.transpose(basis),atdirect))
    return np.array(cart)

def genposcar(poscarfn,modenum,freq,basis,natom,atom_types,species_atom,cartshiftdm):
    poscar_fh=open(poscarfn,'w')
    poscar_fh.write('Shifted geometry for %d mode, %f cm-1\n' % (modenum,freq))
    poscar_fh.write('1.0\n')
    for b in basis:
        poscar_fh.write('%20.16f %20.16f %20.16f\n' % (b[0], b[1], b[2]))
    poscar_fh.write(' '.join('%s' % s for s in atom_types) + '\n')
    poscar_fh.write(' '.join('%d' % s for s in species_atom)+ '\n')
    poscar_fh.write('Cartesian\n')

    for cartat in cartshiftdm:
        poscar_fh.write('%12.9f %12.9f %12.9f\n' % (cartat.tolist()[0], cartat.tolist()[1], cartat.tolist()[2]))
def chemelemnumber(chelm,atom_data):
    num=-1
    for el in atom_data:
        if(el[1]==chelm):
            num=el[0]
            break
    return num

import numpy as np
import re
import sys
from math import pi
from math import exp
import os
import argparse
import phonopy.structure.spglib as spl
import seekpath

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

parser = argparse.ArgumentParser(description='The program is to get primitive cell from POSCAR and plot k-path')


parser.add_argument("-i", "--input", action="store", type=str, dest="poscar_fn", default='POSCAR', help="POSCAR filename (default POSCAR)")

args = parser.parse_args()

try:
    poscar_fh=open(args.poscar_fn,'r')
except IOError:
    print('Error opening POSCAR file')
mult=1.0
line=poscar_fh.readline()
line=poscar_fh.readline()
if(isNumeric(line)):
    mult=float(line)
else:
    print('Error reading POSCAR file. Mult coefficient is not numeric')
    sys.exit(0)

b=[]
for i in range(3):
    line=poscar_fh.readline()
    try:
        b.append([float(line.split()[i])*mult for i in range(3)])
    except e:
        print('Error reading basis from POSCAR')
        sys.exit(1)
if (len(b) == 3):
    basis=np.array(b).reshape(3,3)
else:
    print('Error reading basis from POSCAR number of vectors not equival to 3')
    sys.exit(1)

# Read species and number of atoms
line=poscar_fh.readline()
species=[]
if (not isNumeric(line.split()[0])):
    try:
        for sp in line.split():
            species.append(sp)
    except:
        print('Error reading species array in POSCAR file')
    line=poscar_fh.readline() 

specnum=[]
if (not isNumeric(line.split()[0])):
    print('Error reading species from POSCAR file')
    sys.exit(1)
else:
    try:
        for sp in line.split():
            specnum.append(int(sp))
    except:
        print('Error reading number of species array in POSCAR file')

if (len(species)!=len(specnum)):
   print('Error in reading species, number of species and it\'s number not complies')
   sys.exit(1)
num_atom=0
for spn in specnum:
    num_atom+=spn
# Read coordinates. Only direct mode implemented. Do not use cartesian!
line=poscar_fh.readline()
if (not 'direct' in line.lower()):
    print('Error. Only direct positions of atom mode implemented. Make sure the keyword Direct exist in POSCAR')
    sys.exit(1)
dp=[]
try:
    for atnum in range(num_atom):
        line=poscar_fh.readline()
        dp.append([float(line.split()[i]) for i in range(3)]) 
except:
    print('Error reading atom positions from POSCAR file')
    sys.exit(1)
directpos=np.array(dp).ravel()
basis=np.array(b).ravel()

spec=[]
for i in range(len(specnum)):
    for j in range(specnum[i]):
        for a in atom_data:
            if(a[1]==species[i]):
                spec.append(a[0])
print("Input data:")
print("Cell:")
i=0
for v in b:
    print("a%d= % 9.7f % 9.7f % 9.7f" % (i,v[0],v[1],v[2]))
    i+=1

print(spec)
print('Atomic positions:')
print(dp)
structure=[b,dp,spec]

res=seekpath.getpaths.get_path(structure, with_time_reversal=True, recipe='hpkot', threshold=1e-07, symprec=1e-05, angle_tolerance=-1.0)
print("Output data:")
print("Space group: %s" %spl.get_spacegroup(cell=(b,dp,spec), symprec=1e-5))
print('Primitive lattice:')
for v in res['primitive_lattice']:
    print("a1= % 9.7f % 9.7f % 9.7f" % (v[0],v[1],v[2]))
print('Atomic positions:')
for i in range(len(res['primitive_positions'])):
    print("%2d  %s" % (res['primitive_types'][i], "".join("  % 12.9f" % val for val in res['primitive_positions'][i])))

print('Crystallographic lattice:')
for v in res['conv_lattice']:
    print("a1= % 9.7f % 9.7f % 9.7f" % (v[0],v[1],v[2]))
print('Points:')
for pt in res['point_coords'].items():
    print(pt)
print('PATH:')
print (res['path'])