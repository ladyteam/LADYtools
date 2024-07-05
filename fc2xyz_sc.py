#!/usr/bin/env python3.9
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

import numpy as np

def direct2cart(directpos,basis):
    cart=[]
    for atdirect in directpos:
        cart.append(np.dot(np.transpose(basis),atdirect))
    return np.array(cart)

def coordination_sphere(supercell, atnums, uc_center):
    dists = []
    for n in atnums:
        dist = 0
        for i in range(3):
            dist += (supercell.positions[n][i] - uc_center[i])**2
        dists.append([n, dist**0.5])

    distsnp = np.array(dists)
    return([dists[i] for i in np.argsort(distsnp[:,1])])

import phonopy
from phonopy.api_phonopy import Phonopy
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
                                                help="Supercell Transformation matrix (calculated)")
parser.add_argument("-t", "--tmat", action="store", type=str, dest="tmat", default="",
                                                help="Supercell Transformation matrix (for visualization)")
parser.add_argument("--irreps", dest="irreps", action="store_true",
                                                    default=False, help="Find Irreducible representations and print in input files")
parser.add_argument("--fc-symmetry", dest="fcsymmetry", action="store_true",
                                                    default=False, help="Apply accoustic sume rule")

Angst2Bohr=1.889725989

args = parser.parse_args()
phase = np.array([1, 1, 1, 1, 1])
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
        sc=np.array([float(d) for d in args.dim.split()]).reshape(3,3)
    else:
        sc=np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)
        for i in range(3):
            sc[i][i]=float(args.dim.split()[i])
else:
    sc=np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)

if ( len(args.tmat.split()) > 0 ):
    if (len(args.tmat.split()) > 3):
        print(args.tmat.split())
        tdim=np.array([float(d) for d in args.tmat.split()]).reshape(3,3)
    else:
        tdim=np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)
        for i in range(3):
            tdim[i][i]=float(args.tmat.split()[i])
else:
    tdim=np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)
print("tdim: ", tdim)
print("tdim determinant = %d" % int(np.linalg.det(tdim)))
pm = np.array([1,0,0,0,1,0,0,0,1]).reshape(3,3)

ph = phonopy.load(supercell_matrix=sc,
                  primitive_matrix=pm,
                  unitcell_filename=args.in_fn,
                  is_nac=args.nac,  calculator=args.calc, factor=args.factor,
                  force_constants_filename=args.fc_fn)

if (args.fcsymmetry):
    ph.symmetrize_force_constants()

basis = ph.get_unitcell().get_cell()
xred = ph.get_unitcell().get_scaled_positions()

natom = len(xred)
masses = ph.get_unitcell().get_masses()
chemel = ph.get_unitcell().get_chemical_symbols()


if (args.irreps):
# Hack to do irrep analysis. Spin order will be restored later
    ph.primitive.set_magnetic_moments(None)
# Set IR for Gamma point
    if(args.nacqdir):
        ph.set_irreps([0.0,0.0,0.0],
                      nac_q_direction=[int(d) for d in args.nacqdir.split()])
    else:
        ph.set_irreps([0.0,0.0,0.0])
    ir_labels = ['N' for i in range(natom*3)]
#    print(ph.get_irreps()._get_degenerate_sets()[0])
#    print(ir_labelstmp)
    for deg_set, ir in zip(ph.get_irreps()._get_degenerate_sets(),
                           ph.get_irreps()._get_ir_labels()):
        print(deg_set, ir)
        for d in deg_set:
            if(ir is None):
                ir_labels[d] = 'Non'
            else:
                ir_labels[d] = ir
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
        print("ERROR Couldn't open shift file, exiting...\n")
        sys.exit(1)
    shifts=[]
    l=0
    for line in shift_fh:
        if (l > natom*np.linalg.det(tdim)):
            break
        shifts.append([float(line.split()[i]) for i in range(3)])
        l=l+1
    if (len(shifts) != int(natom*np.linalg.det(tdim))):
        print('The shift array have a wrong dimention %d' % len(shifts))
        sys.exit(1)
    shiftv = np.array(shifts).reshape(int(natom*np.linalg.det(tdim)), 3)

#    for i in range(natom):
#        xred[i] = xred[i] + shiftv[i]
#        print(xred[i])

# Generate supercell for visualization
ph_sc = Phonopy(ph.get_unitcell(), supercell_matrix = tdim)
print("ph_sc_species: ", ph_sc.supercell.symbols)
print(ph_sc.symmetry.get_independent_atoms())
print(ph_sc.symmetry.get_map_atoms())
print(ph_sc.primitive.symbols)
print("uc2sc_map")
print(ph_sc.unitcell.symbols)
uc_center = np.array([0.0, 0.0, 0.0])
for coor in ph_sc.unitcell.positions:
    uc_center += coor
# Calculate center of unitcell
uc_center = uc_center/len(ph_sc.unitcell.positions)

u2s_map = {}

for i in range(len(ph_sc.supercell.u2s_map)):
    print("atom in unitcell: %d %s" % (i, ph_sc.unitcell.symbols[i]))
    print("Maps to atoms in supercell:")
    atnums = np.where(ph_sc.supercell.s2u_map == ph_sc.supercell.u2s_map[i])[0]
    print(atnums)
    print(coordination_sphere(ph_sc.supercell, atnums, uc_center))
    u2s_map.update({i: coordination_sphere(ph_sc.supercell, atnums, uc_center)})
print(u2s_map)
#print(ph_sc.supercell.s2u_map)
#print(ph_sc.primitive.p2s_map)
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


frequencies = np.sqrt(np.abs(eigvals)) * np.sign(eigvals)
#print ("Mode frequencies in cm-1")
#for freq in frequencies:
#    print(freq*args.factor)

# Abinit use atomic units
if (args.calc=='abinit'):
    lUnits=1.889725989
else:
    lUnits=1

if (args.shift_fn == None):
    shiftv_cart = np.zeros(len(ph_sc.supercell.symbols)*3).reshape(len(ph_sc.supercell.symbols),3)
else:
    shiftv_cart = direct2cart(shiftv, ph_sc.supercell.cell)

for j in range(natom*3):
    print("Mode %d" % j)
    out_fh.write('%d\n' % len(ph_sc.supercell.positions))
    out_fh.write('Mode %d %fcm-1 IR: %s\n' % ((j+1), frequencies[j]*args.factor, ir_labels[j] ))
    s = 0
    for i in range(natom):
        shiftvec=[0.0e0,0.0e0,0.0e0]
        for l in range(3):
            shiftvec[l]=eigvecs[i*3+l,j]*sqrt(1/masses[i])/lUnits*args.mult
# Write eigenvector for each atom in supercell
        for k in range(len(u2s_map[i])):
        #    print(k, " : ", u2s_map[i][k])
#            print("Shift atom %d by vector " % k, shiftv_cart[k])
            out_fh.write('%2s '% chemel[i])
            out_fh.write(' '.join(' % 12.8f' % (ph_sc.supercell.positions[u2s_map[i][k][0]][l]/lUnits + 
                                                shiftv_cart[s][l]) for l in range(3)))
            out_fh.write(' '.join(' % 12.8f' % shiftvec[l]*phase[k] for l in range(3)))
            out_fh.write('\n')
            s+=1

# Change to supercell basis
out_fh.write ('Jmol command to plot unitcell (Angstroms):\n load "" {1 1 1} UNITCELL ['+
                        ','.join('%8.7f' % (b/lUnits) for b in ph_sc.supercell.cell.flatten()) + ']')


