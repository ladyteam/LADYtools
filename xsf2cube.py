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
# Data file in xsf format to CUBE converter (the case of CRYSTAL)
#
# Author: Eugene Roginskii
#
# xsf format description:
# http://www.xcrysden.org/doc/XSF.html
# CUBE format description
# http://www.xcrysden.org/doc/XSF.html
#
#  Instructions to use with Jmol:
# 1. Run Jmol
# 2. Open Jmol console (File-> console)
# 3. Run the command:  load "filename.cube"
# Additionaly run the command:  isosurface s1 sign cyan yellow cutoff 1 "filename.cube" color translucent 0.45
# where filename.cube is the filename given by parameter .cube


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

def direct2cart(directpos,basis):
    cart=[]
    for atdirect in directpos:
        cart.append(np.dot(np.transpose(basis),atdirect))
    return np.array(cart)

def cart2direct(cartpos,basis):
    direct=[]
    for atcart in cartpos:
        direct.append(np.dot(np.linalg.inv(np.transpose(basis)),atcart))
    return np.array(direct)

import numpy as np
import re
import sys
import argparse
import gzip

Angst2Bohr=1.889725989
HaTocm1=219474.631370515
factorcm=716.851928469
AMU = 1.6605402e-27 # [kg]
factorHz=21.49068e12
hbar=6.6260695729e-34 #J*s


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

parser = argparse.ArgumentParser(description='Data file in xsf format to CUBE converter (the case of CRYSTAL)')

parser.add_argument("-i", "--input", action="store", type=str, dest="xsf_fn",  help="Input filename")
parser.add_argument("-o", "--output", action="store", type=str, dest="cube_fn",   help="Output filename")
parser.add_argument("-c", "--comment", action="store", type=str, dest="comment",  default="", help="Comment to be added into output file")
parser.add_argument("-z", "--gzip", action="store_false", dest="tozip", help="Compress output. Warning unitcell data  won't be loaded in Jmol")

args = parser.parse_args()

if (args.xsf_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

try:
    xsf_fh = open(args.xsf_fn, 'r')
except IOError:
    print("ERROR Couldn't open input file, exiting...\n")
    sys.exit(1)

for line in xsf_fh:
    if not line: break
    if (re.match('^\s*#',line)):
        continue
    if not line.strip():
        continue
    if '#' in line:
        line=re.sub('\s*#.*','',line)
    if ('crystal' in line.lower()):
        break
    else:
        print('Moleclue mode is not yet implementet. Exit.')
        sys.exit(1)

# Read Primitive vectors

for line in xsf_fh:
    if not line: break
    if (re.match('^\s*#',line)):
        continue
    if not line.strip():
        continue
    if '#' in line:
        line=re.sub('\s*#.*','',line)
    if ('begin_block' in line.lower()):
        break
    if ('primvec' in line.lower()):
        break

primvec=[]
if ('primvec' in line.lower()):
    for i in range(3):
        line=xsf_fh.readline()
        primvec.append([float(l) for l in line.split()[0:3]])
else:
    print('Error. No primvec section found.')
    sys.exit(1)

# Read Primitive coordinates and atomic species

xsf_fh.seek(0)

for line in xsf_fh:
    if not line: break
    if (re.match('^\s*#',line)):
        continue
    if not line.strip():
        continue
    if '#' in line:
        line=re.sub('\s*#.*','',line)
    if ('begin_block' in line.lower()):
        break
    if ('primcoord' in line.lower()):
        break

primcoord=[]

if ('primcoord' in line.lower()):
    line=xsf_fh.readline()
    for i in range(int(line.split()[0])):
        line=xsf_fh.readline()
        atnum=-1
        for a in atom_data:
            if(line.split()[0].lower()==a[1].lower()):
                atnum=a[0]
        if (atnum==-1):
            print('Error reading atomic species. Exit.')
            sys.exit(1)
        primcoord.append([atnum, [float(c) for c in line.split()[1:4]]])
else:
    print('Error. No primcoord section found.')
    sys.exit(1)

# Read 3d DATAGRID data

for line in xsf_fh:
    if not line: break
    if (re.match('^\s*#',line)):
        continue
    if not line.strip():
        continue
    if '#' in line:
        line=re.sub('\s*#.*','',line)
    if ('begin_datagrid_3d' in line.lower()):
        break

if not ('begin_datagrid_3d' in line.lower()):
    print('Error. Block BEGIN_DATAGRID_3D not found.')
    sys.exit(1)

# Read number of points

line=xsf_fh.readline()

if(strType(line.split()[0]==int)):
    npntx=int(line.split()[0])
else:
    print('Error reading number of points.')
    sys.exit(1)

if(strType(line.split()[1]==int)):
    npnty=int(line.split()[1])
else:
    print('Error reading number of points.')
    sys.exit(1)

if(strType(line.split()[2]==int)):
    npntz=int(line.split()[2])
else:
    print('Error reading number of points.')
    sys.exit(1)

# Read origin
origin=[]

try:
    line=xsf_fh.readline()
    origin=[float(l) for l in line.split()]
except:
    print('Error reading origin.')
    sys.exit(1)

# Read spaning vector

svec=[]
try:
    for i in range(3):
        line=xsf_fh.readline()
        svec.append([float(l) for l in line.split()])
except:
    print('Error reading spanning vector')
    sys.exit(1)

print('header')
#buf='%s %s\n' % ('args.comment','jmols')
#print(buf.encode('latin-1'))
#Read data from file
value=[]

for line in xsf_fh:
    if not line: break
    if (re.match('^\s*#',line)):
        continue
    if not line.strip():
        continue
    if '#' in line:
        line=re.sub('\s*#.*','',line)
    if ('end_datagrid_3d' in line.lower()):
        break
    for l in line.split():
        if(isNumeric(l)):
            value.append(float(l))
        else:
            print('Error reading datagrid data')
            sys.exit(1)

# Compile sort data by x,y,z
data=[[[0 for z in range(npntz) ] for y in range(npnty)] for x in range(npntx)] 

n=0

for k in range(npntz):
    for j in range(npnty):
        for i in range(npntx):
            data[i][j][k]=value[n]
            n+=1

# Write to output file in CUBE format
#
# Non Zipped
if(args.tozip):
    if (args.cube_fn == None):
        print('Error. No output filename was given.')
        sys.exit(1)

    try:
        cube_fh = open(args.cube_fn, 'w')
    except IOError:
        print("ERROR Couldn't open output file, exiting...\n")
        sys.exit(1)
    # Generate jmolscript string
    basis=np.array(primvec).reshape(3,3)

    jmols='load "" {444 555 1} UNITCELL [' + ','.join('%8.7f' % b for b in basis.flatten()) + ']; isosurface s1 sign cyan yellow cutoff 1 "" color translucent 0.45'
    # Write comments
    cube_fh.write('%s %s\n' % (args.comment,jmols))
    cube_fh.write('%s %s\n' % ('jmolscript: ',jmols))

    #Angst2Bohr
    # Write number of atoms and origin
    cube_fh.write('%d %s\n' % (len(primcoord), ''.join(' % 10.6f' % (o*Angst2Bohr) for o in origin)))
    # Write number of voxels and axis vector
    cube_fh.write('%5d %s\n' % (npntx,''.join(' % 10.6f' % (s/npntx*Angst2Bohr) for s in svec[0])))
    cube_fh.write('%5d %s\n' % (npnty,''.join(' % 10.6f' % (s/npnty*Angst2Bohr) for s in svec[1])))
    cube_fh.write('%5d %s\n' % (npntz,''.join(' % 10.6f' % (s/npntz*Angst2Bohr) for s in svec[2])))
    # Write atoms

    for i in range(len(primcoord)):
        cube_fh.write('%3d %8.4f %s\n' % (primcoord[i][0],atom_data[primcoord[i][0]][3],''.join('% 10.6f' % (v*Angst2Bohr) for v in primcoord[i][1])))

    n=1
    for i in range(npntx):
        for j in range(npnty):
            for k in range(npntz):
                cube_fh.write('% 8.5E ' % data[i][j][k])
                if (n%6==0):
                    cube_fh.write('\n')
                n+=1
# Zipped
else:
    if (args.cube_fn == None):
        print('Error. No output filename was given.')
        sys.exit(1)

    try:
        cube_fh = gzip.open(args.cube_fn, 'wb')
    except IOError:
        print("ERROR Couldn't open output file, exiting...\n")
        sys.exit(1)
    # Generate jmolscript string
    basis=np.array(primvec).reshape(3,3)

    jmols='load "" {444 555 1} UNITCELL [' + ','.join('%8.7f' % b for b in basis.flatten()) + ']; isosurface s1 sign cyan yellow cutoff 1 "" color translucent 0.45'
    # Write comments
    buf='%s %s\n' % (args.comment,jmols)
    cube_fh.write(buf.encode('latin-1'))
    buf='%s %s\n' % ('jmolscript: ',jmols)
    cube_fh.write(buf.encode('latin-1'))

    #Angst2Bohr
    # Write number of atoms and origin
    buf='%d %s\n' % (len(primcoord), ''.join(' % 10.6f' % (o*Angst2Bohr) for o in origin))
    cube_fh.write(buf.encode('latin-1'))
    # Write number of voxels and axis vector
    buf='%5d %s\n' % (npntx,''.join(' % 10.6f' % (s/npntx*Angst2Bohr) for s in svec[0]))
    cube_fh.write(buf.encode('latin-1'))
    buf='%5d %s\n' % (npnty,''.join(' % 10.6f' % (s/npnty*Angst2Bohr) for s in svec[1]))
    cube_fh.write(buf.encode('latin-1'))
    buf='%5d %s\n' % (npntz,''.join(' % 10.6f' % (s/npntz*Angst2Bohr) for s in svec[2]))
    cube_fh.write(buf.encode('latin-1'))
    # Write atoms

    for i in range(len(primcoord)):
        buf='%3d %8.4f %s\n' % (primcoord[i][0],atom_data[primcoord[i][0]][3],''.join('% 10.6f' % (v*Angst2Bohr) for v in primcoord[i][1]))
        cube_fh.write(buf.encode('latin-1'))

    n=1
    for i in range(npntx):
        for j in range(npnty):
            for k in range(npntz):
                buf='% 8.5E ' % data[i][j][k]
                cube_fh.write(buf.encode('latin-1'))
                if (n%6==0):
                    cube_fh.write(b'\n')
                n+=1

