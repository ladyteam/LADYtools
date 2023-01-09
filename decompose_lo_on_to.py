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
# Decompose LO modes eigenvectors on TO modes eigenvectors basis set
# 
# Author: Eugene Roginskii

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

parser = argparse.ArgumentParser(description='Decompose LO mode on TO basis')

parser.add_argument("-i", "--input", action="store", type=str, dest="in_fn",  help="Input filename")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn",  default='modes.xyz', help="Output filename")
parser.add_argument("-c", "--calc", action="store", type=str, dest="calc",  help="Calculator (abinit,vasp,qe,castep,cp2k)")
parser.add_argument("-f", "--fc", action="store", type=str, dest="fc_fn", default='FORCE_CONSTANTS', help="Force constants filename")
parser.add_argument("-n", "--lonum", action="store", type=int, dest="modenum", default=1, help="Number of LO mode, starting from 1")
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
                  is_nac=True,  calculator=args.calc, factor=args.factor,
                  force_constants_filename=args.fc_fn)
#                  born_filename='BORN')

basis = ph.get_unitcell().get_cell()
xred = ph.get_unitcell().get_scaled_positions()

natom = len(xred)
masses = ph.get_unitcell().get_masses()
chemel = ph.get_unitcell().get_chemical_symbols()


if (args.irreps):
# Hack to do irrep analysis. Spin order will be restored later
    ph.primitive.set_magnetic_moments(None)
# Set IR for Gamma point
    ph.set_irreps([0.0,0.0,0.0])
    ir_labels = ['N' for i in range(natom*3)]
#    print(ph.get_irreps()._get_degenerate_sets()[0])
#    print(ir_labelstmp)
    for deg_set, ir in zip(ph.get_irreps()._get_degenerate_sets(),
                           ph.get_irreps()._get_ir_labels()):
#        print(deg_set, ir)
        for d in deg_set:
            if(ir is None):
                ir_labels[d] = 'Non'
            else:
                ir_labels[d] = ir
else:
    ir_labels=["" for x in range(natom*3)]

#print(ir_labels)

try:
    out_fh = open(args.out_fn, 'w')
except IOError:
    print("ERROR Couldn't open output file for writing, exiting...\n")
    sys.exit(1)

ph.run_qpoints([[0, 0, 0]], with_eigenvectors=True)
eigens_to = ph.get_qpoints_dict()['eigenvectors'][0].real
freqs_to = ph.get_qpoints_dict()['frequencies'][0]

#Debug TO eigens and orthogonality of TO basis
#
#for i in range(natom*3):
#    print("Freq % 9.4f" % freqs_to[i])
#    for j in range(natom):
#        print(''.join('% 9.6f ' % eigens_to[:,i][j*3+k] for k in range(3)))
#    weights=np.dot(eigens_to[:,i], eigens_to)
#    for j in range(len(weights)):
#        print("% 4.2f % 9.5f %s" % (weights[j], freqs_to[j], ir_labels[j]))

ph.run_qpoints([[0, 0, 0]], nac_q_direction=[float(l) for l in args.nacqdir.split()], with_eigenvectors=True, with_dynamical_matrices=True)
eigens_lo = ph.get_qpoints_dict()['eigenvectors'][0].real
freqs_lo = ph.get_qpoints_dict()['frequencies'][0]
#for i in range(natom*3):
#    print(eigens_lo[:,i])

#Debug treat all LO modes
#for i in range(natom*3):
#    print("Freq % 9.4f" % freqs_lo[i])
#    weights=np.dot(eigens_lo[:,i], eigens_to)
#    for j in range(natom):
#        print(''.join('% 8.6f ' % eigens_lo[:,i][3*j+k] for k in range(3)))
#    for j in range(len(weights)):
#        print("% 4.2f % 9.5f %s" % (weights[j], freqs_to[j], ir_labels[j]))

print("Mode number %d Freq % 9.4f" % (args.modenum, freqs_lo[args.modenum-1]))
weights=np.dot(eigens_lo[:,args.modenum-1], eigens_to)
for j in range(natom):
    print(''.join('% 8.6f ' % eigens_lo[:,args.modenum-1][3*j+k] for k in range(3)))
for j in range(len(weights)):
    print("% 4.2f % 10.4f %s" % (weights[j], freqs_to[j], ir_labels[j]))

