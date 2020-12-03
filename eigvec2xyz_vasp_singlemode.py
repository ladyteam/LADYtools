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
# Extract eigenvectoe from PHONOPY dynamical matrix format and construct file of XYZ format with displacements 
# (in cartezian coordinates).
#
# Author; Eugene Roginskii
#
# XYZ format:
# http://wiki.jmol.org/index.php/File_formats/Formats/XYZ


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

def _iterparse(fname, tag=None):
    for event, elem in etree.iterparse(fname):
        if tag is None or elem.tag == tag:
            yield event, elem

def _get_atom_types_from_vasprun_xml(root):
    atom_types = []
    masses = []
    num_atom = 0
    species_atom = []
    print(root.tag)
    for element in root:
      if element.tag == 'array':
        if 'name' in element.attrib:
            if element.attrib['name'] == 'atomtypes':
                for rc in element.findall('./set/rc'):
                    atom_info = [x.text for x in rc.findall('./c')]
                    num_atom += int(atom_info[0])
                    species_atom.append(int(atom_info[0]))
                    atom_types.append(atom_info[1].strip())
                    masses += ([float(atom_info[2])] * int(atom_info[0]))

    return atom_types, masses, num_atom, species_atom

def _get_basis(root):
    basis=[]
    for element in root:
        if element.tag == 'crystal':
            for rc in element:
                if rc.tag == 'varray':
                    if 'name' in rc.attrib:
                        if rc.attrib['name'] == 'basis':
                            for vectors in rc.findall('./v'):
                                basis.append([float(vec) for vec in vectors.text.split()])
    return [basis[i] for i in range(3)]

def _get_atpos(root):
    directpos=[]
    nat=0
    for element in root:
        if element.tag == 'varray':
            if 'name' in element.attrib:
                if element.attrib['name'] == 'positions':
                    for vectors in element.findall('./v'):
                        directpos.append([float(vec) for vec in vectors.text.split()])
                        nat+=1
    return [directpos[i] for i in range(nat)]

def direct2cart(directpos,basis):
    cart=[]
    for atdirect in directpos:
        cart.append(np.dot(np.transpose(basis),atdirect))
    return np.array(cart)

def genposcar(poscarfn,modenum,freq,basis,natom,atom_types,species_atom,cartshiftdm):
    poscar_fh=open(poscarfn,'w')
    poscar_fh.write('Epsilon calculation for %d mode, %f cm-1\n' % (modenum,freq))
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

delta=0.25
AMU = 1.6605402e-27 # [kg]

factorcm=521.47083
# VASPToTHz freq to HZ factor
factorTHz=15.633302
factorHz=15.633302e12
#VaspToTHz = sqrt(EV/AMU)/Angstrom/(2*pi)/1e12 # [THz] 15.633302 
hbar=6.6260695729e-34 #J*s
hbar2AMU=hbar/AMU
#hbar=4.13566751691e-15 #eV*s
EV = 1.60217733e-19 # [J]


parser = argparse.ArgumentParser(description='Script to generate XYZ file with atomic displacements using dynamical matrix calculated with VASP and PHONOPY')

parser.add_argument("-d", "--dynmat", action="store", type=str, dest="dynmat_fn", default='qpoints.yaml', help="Dynamic matrix in yaml format")
parser.add_argument("-i", "--poscar", action="store", type=str, dest="poscar_fn", default='POSCAR', help="vasp POSCAR file name")
parser.add_argument("-v", "--vasprun", action="store", type=str, dest="vasprun_fn",  default='vasprun.xml', help="vasp xml file name")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn", default='out.xyz', help="XYZ output filename file name")
parser.add_argument("-s", "--shift", action="store", type=str, dest="shift_fn",  help="Shift positions filename")
parser.add_argument("-m", "--modnum", action="store", type=int, dest="modnum", default=1, help="Mode number")

args = parser.parse_args()

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

try:
    qpoints_fh = open(args.dynmat_fn, 'r')
except IOError:
    print("ERROR Couldn't open qpoints file, exiting...\n")
    sys.exit(1)

natom=0

for line in qpoints_fh:
    if  'natom' in line:
        natom=int(line.split(':')[1])
        print(natom*3)
        break

qpoints_fh.close()

#  =========================== PHONOPY CODE =======================
try:
    qdata = yaml.load(open(args.dynmat_fn),Loader=yaml.CSafeLoader)
    dynmat = []
    dynmat_data = qdata['phonon'][0]['dynamical_matrix']
    for row in dynmat_data:
        vals = np.reshape(row, (-1, 2))
        dynmat.append(-vals[:, 0] - vals[:, 1] * 1j)

    dm = np.array(dynmat, dtype='double')

    eigvals, eigvecs, = np.linalg.eigh(dm)
except:
    print('Filed to read DYNAMICAL MATRIX')
    sys.exit(0)


#print('***dynmat***')
#for i in range(natom):
#    print (dm[i])
#print(eigvals)
#eigvals, evecs = np.linalg.eigh(dm)
#eigvals = eigvals.real
frequencies=[]

#print('***eigvecs***')
#print(eigvecs)

frequencies=np.sqrt(np.abs(eigvals)) * np.sign(eigvals) * -1
#print ("Mode frequencies in cm-1")
#for freq in frequencies:
#    print(freq*factorcm)



#  Read masses from vasprun.xml
for event, element in _iterparse(args.vasprun_fn,'atominfo'):
    (atom_types,  masses,  num_atom, species_atom) = _get_atom_types_from_vasprun_xml(element)
    break
if (natom != num_atom):
    print("Error parsing xml file, wrong num_atom value")
    sys.exit(1)

#print(masses)

if (os.path.isfile(args.poscar_fn)):
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

# Read species
    line=poscar_fh.readline()
    species=[]
    if (not isNumeric(line.split()[0])):
        try:
            for sp in line.split():
                species.append(sp)
        except:
            print('Error reading species array in POSCAR file')
        line=poscar_fh.readline() 
    else:
        print('Error reading species, please add text desription of atom speciec before numers of each atoms')
        sys.exit(1)
# Read number of species
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
# Construct array with atoms
    atype=[]
    for i in range(len(species)):
        for j in range(specnum[i]):
            atype.append(species[i])
#    print(atype)

# Check vasprun.xml and POSCAR data and fix masses values which are not very accurate in VASP for some reason
    for n in range(len(atype)):
            elnum=chemelemnumber(atype[n],atom_data)

            if (elnum<0):
                print('Error getting chemical element number')
                sys.exit(1)
            if (abs(atom_data[elnum][3] - masses[n])>0.08):
                print('Error.masses in periodic table and in vasprun.xml is not the same for element number %d: %s' % (n,atom_data[elnum][1]))
                print(
                atom_data[elnum][3])
                print('mass: %f' % masses[n])
                print('masses array is:')
                print(masses)
                sys.exit(0)
            else:
                masses[n]=atom_data[elnum][3]

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
    except e:
        print('Error reading atom positions from POSCAR file')
        sys.exit(1)
    directpos=np.array(dp).reshape(num_atom,3)
else:
    print('Warning! POSCAR was not found, will use vasp.xml to read crystallographic data.')
#  Read vectors from vasprun.xml 
    atype=[]
    for i in range(len(atom_types)):
        if(strType(species_atom[i])=='int'):
            for k in range(int(species_atom[i])):
                atype.append(atom_types[i])
        else:
            print('Error constructing atom type array')
            sys.exit(1)

    for event, element in _iterparse(args.vasprun_fn,'structure'):
        if 'name' in element.attrib:
            if element.attrib['name'] == 'initialpos':
                basis=np.array(_get_basis(element)).reshape(3,3)
                directpos=np.array(_get_atpos(element)).reshape(num_atom,3)

#cvol=np.dot(basis[0],np.cross(basis[1],basis[2]))

#print (atom_types)
#print(basis)
print ('Jmol command to plot unitcell (Angstroms):\n load "" {1 1 1} UNITCELL ['+
                        ','.join('%8.7f' % b for b in basis.flatten()) + ']')


cartpos=direct2cart(directpos,basis)

# All this shifting stuff
if (args.shift_fn == None):
    print('No shifts data was provded. Using input strcuture')
else:
    try:
        shift_fh = open(args.shift_fn, 'r')
    except IOError:
        print("ERROR Couldn't open shift file, exiting...\n")
        sys.exit(1)
    shifts=[]
    l=0
    for line in shift_fh:
        if l>natom:
            break
        shifts.append([float(line.split()[i]) for i in range(3)])
        l=l+1
    if (len(shifts)!=natom):
        print('The shift array have a wrong dimention')
        sys.exit(1)
    shiftv=np.array(shifts).reshape(natom,3)
    print('shift vector')
    print(shiftv)
    for i in range(natom):
        directpos[i]=directpos[i]+shiftv[i]
    cartpos=direct2cart(directpos,basis)

try:
    out_fh = open(args.out_fn, 'w')
except IOError:
    print("ERROR Couldn't open output file for writing, exiting...\n")
    sys.exit(1)

j=natom*3-args.modnum
out_fh.write('%d\n' % natom)
out_fh.write('Mode %d %fcm-1 ' % ((args.modnum), frequencies[j]*factorcm ))
out_fh.write ('load "" {1 1 1} UNITCELL ['+
                        ','.join('%8.7f' % b for b in basis.flatten()) + ']\n')

for i in range(natom):
    shiftvec=[0.0e0,0.0e0,0.0e0]
    for l in range(3):
        shiftvec[l]=-eigvecs[i*3+l,j]*sqrt(1/(masses[i]))
    out_fh.write('%s '% atype[i])
    out_fh.write(' '.join(' % 11.8f' % cartpos[i][l] for l in range(3)))
    out_fh.write(' '.join(' % 11.8f' % shiftvec[l] for l in range(3)))
    out_fh.write('\n')

