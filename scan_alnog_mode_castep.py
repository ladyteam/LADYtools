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
# The program to distort geometry along selected normal mode eigenvector. Castep output format
# (.phonon file)
#
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
def isDiag(M):
    i, j = np.nonzero(M)
    return np.all(i == j)

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

def gencastep(fn,modenum,basis,natom,typatsym,symprec,atpos):
    from phonopy import Phonopy
    import phonopy.structure.spglib as spg
    from phonopy.structure.atoms import PhonopyAtoms as Atoms
    from phonopy.structure.symmetry import Symmetry, find_primitive, get_pointgroup

    fh=open(fn,'w')
    unitcell = Atoms(symbols=typatsym, cell=basis, positions=atpos)
    pbasis=np.eye(3)
    for i in range(len(basis)):
        pbasis[i]=basis[i]/np.linalg.norm(basis[i])
    symmetry = Symmetry(unitcell, symprec)
    rotations = symmetry.get_symmetry_operations()['rotations']
    translations = symmetry.get_symmetry_operations()['translations']
    print('Space group International symbol: %s' % symmetry.get_international_table())

    fh.write('%BLOCK LATTICE_CART\n')
    for bl in basis:
        fh.write('%s\n' % ''.join(' %12.8f' % b for b in bl))
    fh.write('%ENDBLOCK LATTICE_CART\n\n')
    fh.write('%BLOCK POSITIONS_ABS\n')
    for i in range(len(typatsym)):
        fh.write("  %3s   " % typatsym[i])
        fh.write('%s\n' % ''.join(' %12.8f' % p for p in atpos[i].tolist()))
    fh.write('%ENDBLOCK POSITIONS_ABS\n\n')

    fh.write('SYMMETRY_TOL : %f ang\n' % symprec)
    fh.write('SYMMETRY_GENERATE \n')
    fh.write('#KPOINT_MP_GRID : 4 4 4\n#KPOINT_MP_OFFSET : 0.5 0.5 0.5\n')
#    fh.write('%BLOCK SYMMETRY_OPS\n')
#    for i in range(len(rotations)):
#        fh.write('# Sym. op. number %d\n' % (i+1))
#        if (isDiag(rotations[i])):
#            rot=rotations[i]
#        else:
#            rot=np.dot(np.dot(rotations[i],pbasis),rotations[i])

#        for symopl in rot:
#            fh.write('%s\n' % ''.join(' % 12.10f' % round(sop,6) for sop in symopl))
#        fh.write('%s\n' % ''.join(' % 12.10f' % round(tr) for tr in translations[i]))

#    fh.write('%ENDBLOCK SYMMETRY_OPS\n\n')

def genabinit(abinitfn,modenum,basis,natom,typat,znucl,cartshiftd,comment):
    abinit_fh=open(abinitfn,'w')
    abinit_fh.write('#Atom positions shifted along Mode %d\n' % modenum)
    abinit_fh.write('# %s\n' % comment)
    abinit_fh.write('acell 1.0 1.0 1.0 Angstrom\nrprim\n')
    for b in basis:
        abinit_fh.write('% 20.16f % 20.16f % 20.16f\n' % (b[0], b[1], b[2]))
    abinit_fh.write('typat ' + ' '.join('%d' % typat[i] for i in range(natom))+'\n')
    abinit_fh.write('znucl ' + ' '.join('%d' % znucl[i] for i in range(len(znucl)))+'\n')
    abinit_fh.write('natom %d\n' % natom)
    abinit_fh.write('ntypat %d\n' % len(znucl))
    abinit_fh.write('xangst\n')
    for cartat in cartshiftd:
        abinit_fh.write('%12.9f %12.9f %12.9f\n' % (cartat.tolist()[0], cartat.tolist()[1], cartat.tolist()[2]))

def get_epsilon(abinitfn):
    e=[]
    try:
        abinit_fh = open(abinitfn, 'r')
    except IOError:
        print ("ERROR Couldn't open abinit output file %s, exiting...\n" % abinitfn)
        sys.exit(1)

    while True:
        line=abinit_fh.readline()
        if not line: break
        if 'Dielectric tensor, in cartesian coordinates' in line:
            while True:
                sline=abinit_fh.readline()
                if not sline: break
                if 'Effective charges' in sline: break
                if (re.match('\s*\d+\s*\d+',sline)):
                    e.append(float(sline.split()[4]))
            break
    eps=np.array(e).reshape(3,3)
    return eps



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
#print sqrt(hbar/AMU/10e12)*10e10 #Angstrom
basedirname='epsilon'
Temp=300
ramanmult=0.0028

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




parser = argparse.ArgumentParser(description='The program to generate distortion structures along normal vibrational modes')


parser.add_argument("-i", "--input", action="store", type=str, dest="input_fn", help="Castep output .phonon filename")
parser.add_argument("-m", "--mode", action="store", type=int, dest="modenum", default=1, help="Number of mode")
parser.add_argument("-t", "--tolerance", action="store", type=float, dest="symprec", default=1e-04,
        help="Symmetry determine tolerance")
parser.add_argument("-a", "--amplitudes", action="store", type=str, dest="amplstr", default='0.1 0.25 0.5 1', 
           help="Mode shift  amplitude, space separated")

args = parser.parse_args()

if (args.input_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

try:
    input_fh = open(args.input_fn, 'r')
except IOError:
    print("ERROR Couldn't open qpoints file, exiting...")
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

for line in input_fh:
    if ( 'Number of ions' in line ):
        natom=int(line.split('ions')[1])
        break

for line in input_fh:
    if ( 'Number of branches' in line ):
        if (int(line.split('branches')[1])<modenum) or (modenum<1):
            print('Error. incorrect number of mode provided')
            sys.exit(1)
        break

for line in input_fh:
    if ( 'Unit cell vectors' in line ):
        break
buf=[]
rprim=[]
for line in input_fh:
    if ( 'Fractional Co-ordinates' in line ):
        break
    for coor in line.split():
        rprim.append(float(coor))

basis=np.array(rprim).reshape(3,3)
xred=[]
typat=[]
typatsym=[]

for line in input_fh:
    if ( 'END header' in line ):
        break
    for i in range(3):
        xred.append(float(line.split()[i+1]))
    for i in range(len(atom_data)):
        if(atom_data[i][1]==line.split()[4]):
            typat.append(atom_data[i][0])
    typatsym.append(line.split()[4])

for line in input_fh:
    if ( 'q-pt' in line ):
        qpt=line
        break

frequencies=[]
for line in input_fh:
    if ( 'Phonon Eigenvectors' in line ):
        break
    frequencies.append(float(line.split()[1]))

isMode=False
shiftvector=[]
for line in input_fh:
    if (isNumeric(line.split()[0])):
        if ((int(line.split()[0]) == args.modenum)):
            shiftvector.append(float(line.split()[2]))
            shiftvector.append(float(line.split()[4]))
            shiftvector.append(float(line.split()[6]))

znucl=[]
znucl.append(typat[0])
inArray=True
for i in range(natom):
    inArray=False

    for j in range(len(znucl)):
        if (znucl[j]==typat[i]):
            inArray=True
    if (not inArray):
        znucl.append(typat[i])
typata=[]
for i in range(natom):
    for j in range(len(znucl)):
        if(typat[i]==znucl[j]):
            typata.append(j+1)

print(basis)
print(frequencies)
print(xred)
print(typat)
print(typatsym)
print(shiftvector)
print(znucl)
print(natom)
print(len(shiftvector))
for v in np.array(shiftvector).reshape(natom,3):
    print(v)

directpos=np.array(xred).reshape(natom,3)
cartpos=direct2cart(directpos,basis)

genabinit('ideal.in',0,basis,natom,typata,znucl,cartpos,'ideal structure')
gencastep('ideal.cell',0,basis,natom,typatsym,args.symprec,cartpos)

j=args.modenum-1
for n in range(len(ampl)):
    cartshiftdm=[]
    cartshiftdp=[]

    print ('Process mode: %d  with freq: %8.5f  shifted by: %8.5f' % (j+1,frequencies[j],ampl[n] ))


    for i in range(natom):
        cartshiftdm.append(cartpos[i]-np.array(shiftvector).reshape(natom,3)[i]*ampl[n])

    abinitfn="shiftcell-%.3f.in" % ampl[n]
    castepfn="shiftcell-%.3f.cell" % ampl[n]


    genabinit(abinitfn,j+1,basis,natom,typata,znucl,cartshiftdm,'freq = %9.4f cm-1; Delta=%6.4f' % (frequencies[j],ampl[n]))
    gencastep(castepfn,j+1,basis,natom,typatsym,args.symprec,cartshiftdm)

print('Done.')
