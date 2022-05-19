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
# Extract eigenvectors from PHONOPY FORCE_SET matrix format and construct file of XYZ format with displacements.
# 
# Author: Eugene Roginskii
#
# XYZ format:
# http://wiki.jmol.org/index.php/File_formats/Formats/XYZ


def direct2cart(directpos,basis):
    cart=[]
    for atdirect in directpos:
        cart.append(np.dot(np.transpose(basis),atdirect))
    return np.array(cart)

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



parser = argparse.ArgumentParser(description='Script to generate atomic displacements ini xyz format')

parser.add_argument("-i", "--input", action="store", type=str, dest="in_fn",  help="Input filename")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn",  default='modes.xyz', help="Output filename")
parser.add_argument("-c", "--calc", action="store", type=str, dest="calc",  help="Calculator (abinit,vasp,qe,castep,cp2k)")
parser.add_argument("-f", "--fc", action="store", type=str, dest="fc_fn", default='FORCE_CONSTANTS', help="Force constants filename")
parser.add_argument("-s", "--shift", action="store", type=str, dest="shift_fn",  help="Shift positions filename")
parser.add_argument("-m", "--mult", action="store", type=int, dest="mult", default=10, help="Displacement multiplyer")
parser.add_argument('--nac', dest='nac', action='store_true')
parser.add_argument("--factor", action="store", type=float, dest="factor", 
                     default=716.8519, help="Frequency factor. Default for cm-1 521.47083 (vasp), 716.8519 (abinit), 3739.4256800756 (cp2k)")
parser.add_argument("--q-direction", action="store", type=str, dest="nacqdir", 
                      help="Direction of q-vector for non-analytical term")
parser.add_argument("-d", "--dim", action="store", type=str, dest="dim", default="",
                                                help="Supercell Transformation matrix")
parser.add_argument("--irreps", dest="irreps", action="store_true",
                                                    default=False, help="Find Irreducible representations and print in input files")

Angst2Bohr=1.889725989

args = parser.parse_args()

if (args.in_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

try:
    abinit_fh = open(args.in_fn, 'r')
except IOError:
    print("ERROR Couldn't open input file, exiting...\n")
    sys.exit(1)
if ( len(args.dim.split()) > 0 ):
    if (len(args.dim.split()) > 3):
        sc=np.array(float(d) for d in args.dim).reshape(3,3)
    else:
        sc=np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)
        for i in range(3):
            sc[i][i]=float(args.dim.split()[i])
else:
    sc=np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)
pm=np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)

ph = phonopy.load(supercell_matrix=sc,
                  primitive_matrix=pm,
                  unitcell_filename=args.in_fn,
                  is_nac=args.nac,  calculator=args.calc, factor=args.factor,
                  force_constants_filename=args.fc_fn)

basis = ph.get_unitcell().get_cell()
xred = ph.get_unitcell().get_scaled_positions()

natom = len(xred)
masses = ph.get_unitcell().get_masses()
chemel = ph.get_unitcell().get_chemical_symbols()
ph.set_irreps([0.0,0.0,0.0])

if (args.irreps):
    ir_labels=[]
    for ir in ph.get_irreps()._ir_labels:
        if ('T' in ir):
            for i in range(3):
                ir_labels.append(ir)
            continue
        elif  ('E' in ir):
            for i in range(2):
                ir_labels.append(ir)
            continue
        else:
            ir_labels.append(ir)
else:
    ir_labels=["" for x in range(natom*3)]

print(ir_labels)
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
    if (len(shifts) != natom):
        print('The shift array have a wrong dimention')
        sys.exit(1)
    shiftv = np.array(shifts).reshape(natom, 3)

    for i in range(natom):
        xred[i] = xred[i] + shiftv[i]
        print(xred[i])

cartpos=direct2cart(xred, basis)


try:
    out_fh = open(args.out_fn, 'w')
except IOError:
    print("ERROR Couldn't open output file for writing, exiting...\n")
    sys.exit(1)

if(args.nacqdir):
    ph.run_qpoints([[0, 0, 0]], nac_q_direction=[int(l) for l in args.nacqdir.split()],with_dynamical_matrices=True)
    dynmat_data=ph.get_qpoints_dict()['dynamical_matrices'][0]
else:
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

frequencies=[]


frequencies=np.sqrt(np.abs(eigvals)) * np.sign(eigvals)
print ("Mode frequencies in cm-1")
for freq in frequencies:
    print(freq*args.factor)

# Abinit use atomic units
if (args.calc=='abinit'):
    lUnits=1.889725989
else:
    lUnits=1

for j in range(natom*3):
    out_fh.write('%d\n' % natom)
    out_fh.write('Mode %d %fcm-1 IR: %s\n' % ((j+1), frequencies[j]*args.factor, ir_labels[j] ))

    for i in range(natom):
        shiftvec=[0.0e0,0.0e0,0.0e0]
        for l in range(3):
            shiftvec[l]=eigvecs[i*3+l,j]*sqrt(1/masses[i])/lUnits*args.mult

        out_fh.write('%s '% chemel[i])
        out_fh.write(' '.join(' % 11.8f' % (cartpos[i][l]/lUnits) for l in range(3)))
        out_fh.write(' '.join(' % 11.8f' % -shiftvec[l] for l in range(3)))
        out_fh.write('\n')

out_fh.write ('Jmol command to plot unitcell (Angstroms):\n load "" {1 1 1} UNITCELL ['+
                        ','.join('%8.7f' % (b/lUnits) for b in basis.flatten()) + ']')


