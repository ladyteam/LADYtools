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
# Extract eigenvectors from PHONOPY dynamical matrix format and construct file of XYZ format with displacements.
# 
# Abinit input file format required
#
# Author: Eugene Roginskii
#
# XYZ format:
# http://wiki.jmol.org/index.php/File_formats/Formats/XYZ

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
import phonopy
import re
import yaml
from math import sqrt
import sys
from math import pi
from math import exp
from optparse import OptionParser
import xml.etree.cElementTree as etree
import argparse

Angst2Bohr=1.889725989
HaTocm1=219474.631370515
factorcm=716.851928469
AMU = 1.6605402e-27 # [kg]
factorHz=21.49068e12
hbar=6.6260695729e-34 #J*s

parser = argparse.ArgumentParser(description='Script to generate rprimd and reciprocal base vectros from abinit input file')

parser.add_argument("-i", "--input", action="store", type=str, dest="abinit_fn",  help="Abinit filename")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn",  default='abinit.xyz', help="Output filename")
parser.add_argument("-f", "--fc", action="store", type=str, dest="fc_fn", default='FORCE_CONSTANTS', help="Force constants filename")
parser.add_argument("-s", "--shift", action="store", type=str, dest="shift_fn",  help="Shift positions filename")
parser.add_argument("-m", "--mult", action="store", type=int, dest="mult", default=10, help="Displacement multiplyer")
parser.add_argument('--nac', dest='nac', action='store_false')
parser.add_argument("--factor", action="store", type=float, dest="factor", 
                     default=716.8519, help="Frequency factor. Default for cm-1 521.47083 (vasp), 716.8519 (abinit)")

args = parser.parse_args()

if (args.abinit_fn == None):
    print('Error. No Abinit input filename was given.')
    sys.exit(1)

try:
    abinit_fh = open(args.abinit_fn, 'r')
except IOError:
    print("ERROR Couldn't open abinit file, exiting...\n")
    sys.exit(1)
sc=np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)
pm=np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)

ph = phonopy.load(supercell_matrix=sc,
                  primitive_matrix=pm,
                  unitcell_filename=args.abinit_fn,
                  is_nac=args.nac,  calculator="abinit", factor=args.factor,
                  force_constants_filename=args.fc_fn)

#print(ph.primitive.get_masses())
basis=ph.get_unitcell().get_cell()
xred=ph.get_unitcell().get_scaled_positions()

natom=len(xred)
masses=ph.get_unitcell().get_masses()
chemel=ph.get_unitcell().get_chemical_symbols()
print('basis, Bohr')
print(basis)
print ('Jmol command to plot unitcell (Angstroms):\n load "" {1 1 1} UNITCELL ['+
                        ','.join('%8.7f' % (b/Angst2Bohr) for b in basis.flatten()) + ']')


# All this shifting stuff
if (args.shift_fn == None):
    print('No shifts data was provded. Using input strcuture')
else:
    try:
        shift_fh = open(args.shift_fn, 'r')
    except IOError:
        print("ERROR Couldn't open abinit file, exiting...\n")
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
#    print('shift vector')
#    print(shiftv)
    for i in range(natom):
        xred[i]=xred[i]+shiftv[i]

cartpos=direct2cart(xred,basis)

#print ('cartpos (Bohr)')
#print (cartpos)
#print ('directpos')
#print (directpos)

try:
    out_fh = open(args.out_fn, 'w')
except IOError:
    print("ERROR Couldn't open output file for writing, exiting...\n")
    sys.exit(1)

dynmat_data = ph.get_dynamical_matrix_at_q([0,0,0])

dynmat = []
#i = 0
try:
    for row in dynmat_data:
        dynmat.append(row.real + row.imag)
    dm = np.array(dynmat, dtype='double')
    eigvals, eigvecs, = np.linalg.eigh(dm)
except IOError as e:
    print('Filed to read DYNAMICAL MATRIX %s %s' % (e.errno, e.strerror))
    sys.exit(0)

#print('***dynmat***')
#print (dm)
frequencies=[]

#print('***eigvecs***')
#print(eigvecs)

frequencies=np.sqrt(np.abs(eigvals)) * np.sign(eigvals)
print ("Mode frequencies in cm-1")
for freq in frequencies:
    print(freq*args.factor)


for j in range(natom*3):
    out_fh.write('%d\n' % natom)
    out_fh.write('Mode %d %fcm-1\n' % ((j+1), frequencies[j]*factorcm ))

    for i in range(natom):
        shiftvec=[0.0e0,0.0e0,0.0e0]
        for l in range(3):
# just to escape zero div for accoustic modes
            if (abs(frequencies[j])<0.1):
                frequencies[j]=0.1
#            shiftvec[l]=eigvecs[i*3+l,natom*3-1-j]*sqrt(hbar/(AMU*masses[i]*abs(frequencies[natom*3-1-j])*factorHz))*1e9*Angst2Bohr*args.mult
            shiftvec[l]=eigvecs[i*3+l,j]*sqrt(1/masses[i])/Angst2Bohr*args.mult
#        print(["%10.7f" % shiftvec[l] for l in range(3)])

        out_fh.write('%s '% chemel[i])
        out_fh.write(' '.join(' % 11.8f' % (cartpos[i][l]/Angst2Bohr) for l in range(3)))
        out_fh.write(' '.join(' % 11.8f' % -shiftvec[l] for l in range(3)))
        out_fh.write('\n')


out_fh.write ('Jmol command to plot unitcell (Angstroms):\n load "" {1 1 1} UNITCELL ['+
                        ','.join('%8.7f' % (b/Angst2Bohr) for b in basis.flatten()) + ']')


