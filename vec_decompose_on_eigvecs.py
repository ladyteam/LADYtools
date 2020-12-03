#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
The program for decomposition of the given vector in
cartesian coordinate on Dynamical matrix eigenvectors basis.
'''

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

def direct2cart(directpos,basis):
    cart=[]
    for atdirect in directpos:
        cart.append(np.dot(np.transpose(basis),atdirect))
    return np.array(cart)

def readvec(fn):
    lattice=np.zeros(9)

    xcart=[]
    xatom=[]
    atsym=[]
    try:
        in_fh = open(fn, 'r')
    except IOError:
        print("ERROR Couldn't open input files, exiting...")
        sys.exit(1)

    line=in_fh.readline()
    natom=int(line.split()[0])
    line=in_fh.readline()

    for i in range(9):
        lattice[i]=float(line.split('[')[1].split(',')[i].split(']')[0])

    n=0
    for line in in_fh:
        if(n>natom):
            break
        if(len(line.split())>=3):
            atsym.append(line.split()[0])
            xatom.append([float(line.split()[i+1]) for i in range(3)])
            xcart.append([float(line.split()[i+4]) for i in range(3)])
            n+=1
    if(len(xcart)!=natom):
        print('Error. Mismatch in size of atom coord and number of atoms')

    return([lattice.reshape(3,3), atsym, np.array(xatom).reshape(natom,3), np.array(xcart).reshape(natom,3)])


def chemelemnumber(chelm,atom_data):
    num=-1
    for el in atom_data:
        if(el[1]==chelm):
            num=el[0]
            break
    return num

import numpy as np
import re
from math import sqrt
from math import ceil
import sys
from math import pi
from shutil import move
import os
import datetime
import time
import argparse
import xml.etree.cElementTree as etree
import yaml

delta=1
AMU = 1.6605402e-27 # [kg]
Angst2Bohr=1.889725989
#factorcm=716.851928469
#VaspToTHz = sqrt(EV/AMU)/Angstrom/(2*pi)/1e12 # [THz] 15.633302 
hbar=6.6260695729e-34 #J*s
hbar2AMU=hbar/AMU
#hbar=4.13566751691e-15 #eV*s
EV = 1.60217733e-19 # [J]

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
    [ 23, "V", "Vanadium", 50.9415], # 23z
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


parser = argparse.ArgumentParser(description='Script to decompose LO modes by TO ones')

parser.add_argument("-v", "--vec", action="store", type=str, dest="vec_fn", help="filename with Vector to decompose (xyz format)")
parser.add_argument("-d", "--dynmat", action="store", type=str, dest="dynmat_fn", default='qpoints.yaml', help="Dynamic matrix in yaml format")
parser.add_argument("-c", "--cutoff", action="store", type=float, dest="cut", default=0.001, help="Minimal weight of mode eigenvector")
parser.add_argument("-f", "--factor", action="store", type=float, dest="factorcm", default=521.47083, help="Frequency factor. Default is " +
                                                             "521.47083 (convert to cm-1 when VASP calculator used). Use 716.851928469 for ABINIT ")
parser.add_argument("-a", "--ampl", action="store", type=float, dest="ampl", default=1, help="Amplitude for shift vector")
parser.add_argument("-o", "--out", action="store", type=str, dest="out_fn", help="filename with resulted distorted structure (xyz format)")
args = parser.parse_args()


lattice, atsym, xatom, xcart = readvec(args.vec_fn)

natom = len(xcart)

try:
    qdata = yaml.load(open(args.dynmat_fn),Loader=yaml.CSafeLoader)
    dynmat = []
    dynmat_data = qdata['phonon'][0]['dynamical_matrix']
    for row in dynmat_data:
        vals = np.reshape(row, (-1, 2))
        dynmat.append(vals[:, 0] + vals[:, 1] * 1j)

    dm = np.array(dynmat, dtype='double').transpose()

    eigvals, eigvecs, = np.linalg.eigh(dm)
except:
    print('Filed to read DYNAMICAL MATRIX')
    sys.exit(0)


#eigvals, evecs = np.linalg.eigh(dm)
#eigvals = eigvals.real
frequencies=[]
displ=[]
masses=[]

for el in atsym:
     elnum=chemelemnumber(el,atom_data)
     masses.append(atom_data[elnum][3])

#print('***eigvecs***')
#for e in eigvecs:
#    print(e)

#print('***masses***')
#print(masses)

frequencies=np.sqrt(np.abs(eigvals)) * np.sign(eigvals)

print("Mode frequencies in cm-1")
for freq in frequencies:
    nmod=frequencies.tolist().index(freq)+1
    print("%2d % 9.6f" %(nmod,freq*args.factorcm))

cvol=np.dot(lattice[0],np.cross(lattice[1],lattice[2]))

for j in range(natom*3):
    for i in range(natom):
        for l in range(3):
            displ.append(eigvecs[i*3+l][j]*sqrt(1/(masses[i]))/Angst2Bohr)

displ=np.array(displ).reshape(natom*3,natom*3)
#print("*** Mass weighted eigvec ***")
#for d in displ:
#    print(d)

c=np.zeros(natom*3)

print('Decomposition:')
for i in range(natom*3):
    c[i]=np.dot(displ[i],xcart.flatten())

cnorm=c.max()

modnums=[]

for i in range(natom*3):
    if (abs(c[i]/cnorm)>args.cut):
        modnums.append(i)
        print('Mode number %d with freq %f cm-1, weight=%f' % ((i+1),frequencies[i]*args.factorcm,c[i]/cnorm))

shiftvec=np.zeros(natom*3)
for i in range(len(modnums)):
    shiftvec+=displ[modnums[i]]*c[modnums[i]]

shiftvec.reshape(natom,3)
#for v in shiftvec.reshape(natom,3):
#    print(v)

xcart=shiftvec.reshape(natom,3)*args.ampl
xfin=xatom+xcart*args.ampl

try:
    out_fh = open(args.out_fn, 'w')
except IOError:
    print("ERROR Couldn't open output file, exiting...")
    sys.exit(1)

out_fh.write("%d\n" % len(xfin))
out_fh.write('load "" {1 1 1} UNITCELL ['+
                        ','.join('%8.7f' % b for b in lattice.flatten()) + ']\n')


for i in range(len(xfin)):
    out_fh.write("%2s " % atsym[i])
    out_fh.write(" ".join("% 12.9f " % x for x in xfin[i]))
    out_fh.write(" ".join("% 12.9f " % x for x in xcart[i]*args.ampl))
#    out_fh.write("#%12.9f" % np.linalg.norm(xcart[i])) 
    out_fh.write("\n")

