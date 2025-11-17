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
# Distort initial structure by eigenvector of a given mode by amplitude
# 
# Author: Eugene Roginskii
#


def direct2cart(directpos,basis):
    cart=[]
    for atdirect in directpos:
        cart.append(np.dot(np.transpose(basis),atdirect))
    return np.array(cart)

def cart2direct(cartpos,basis):
    direct=[]
    for atcart in cartpos:
        direct.append(np.dot(np.linalg.inv(np.transpose(basis)),atcart))
    return np.array(direct).reshape(len(cartpos),3)

def savestruct(calc, amp, basis, xred, species, conv_numbers, writeposcar=True):
    from phonopy.structure.atoms import PhonopyAtoms as Atoms
    from phonopy.interface.crystal import write_crystal
    from phonopy.interface.castep import write_castep
    from phonopy.interface.vasp import write_vasp
    from phonopy.interface.abinit import write_abinit
    cell = Atoms(cell=basis, symbols=species, 
            scaled_positions=xred)
    if (calc == 'crystal'):
        fn = 'shiftcell-%5.3f' % amp
        write_crystal(fn, cell, conv_numbers)
        if (writeposcar):
            fn='POSCAR-%5.3f' % amp
            write_vasp(fn, cell)
    elif (calc == 'vasp'):
#        fn = ''.join('%s' % bfn)
        fn='POSCAR-%5.3f' % amp
        write_vasp(fn, cell)
    elif (calc == 'abinit'):
#        fn = ''.join('%s' % bfn)
        fn='shiftcell-%5.3f' % amp
        write_abinit(fn, cell)
    elif (calc == 'cp2k'):
        fn='shiftcell-%5.3f.cell' % amp
        write_castep(fn, cell)

    else:
        print('Writting structures for %s is not implemented yet' % calc)

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
parser.add_argument("-m", "--modnum", action="store", type=str, dest="modnum",  help="Modes number")
parser.add_argument("-c", "--calc", action="store", type=str, dest="calc",  help="Calculator (abinit,vasp,qe,castep,cp2k)")
# read force constants instead FORCE_SET. Fix to true-false
parser.add_argument("-f", "--fc", action="store", type=str, dest="fc_fn", default='FORCE_CONSTANTS', help="Force constants filename")
parser.add_argument("-a", "--ampl", action="store", type=str, dest="ampl",  help="Displacement multiplyer amplitudes")
parser.add_argument("-w", "--weight", action="store", type=str, dest="weight",  help="Weight for each mode")
#parser.add_argument('--nac', dest='nac', action='store_true')
parser.add_argument("--factor", action="store", type=float, dest="factor", 
                     default=716.8519, help="Frequency factor. Default for cm-1 521.47083 (vasp), 716.8519 (abinit), 3739.4256800756 (cp2k)")
#parser.add_argument("--q-direction", action="store", type=str, dest="nacqdir", 
#                      help="Direction of q-vector for non-analytical term")
parser.add_argument("-d", "--dim", action="store", type=str, dest="dim", default="",
                                                help="Supercell Transformation matrix")


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

weights = [1.0]*len(args.modnum.split())
if(args.weight):
    for i in range(len(args.weight.split())):
        weights[i]=float(args.weight.split()[i])
#print(weights)
if ( len(args.dim.split()) > 0 ):
    if (len(args.dim.split()) > 3):
        sc=np.array([float(d) for d in args.dim.split()]).reshape(3,3)
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
                  is_nac=False,  calculator=args.calc, factor=args.factor,
                  force_constants_filename=args.fc_fn)

basis = ph.get_unitcell().get_cell()
xred = ph.get_unitcell().get_scaled_positions()

natom = len(xred)
masses = ph.get_unitcell().get_masses()
chemel = ph.get_unitcell().get_chemical_symbols()
#chemel = ph.get_unitcell().symbols()
chemnum = ph.get_unitcell().get_atomic_numbers()

ph.run_qpoints([[0, 0, 0]], with_eigenvectors=True)
eigens_to = ph.get_qpoints_dict()['eigenvectors'][0].real
freqs_to = ph.get_qpoints_dict()['frequencies'][0]

cartpos=direct2cart(xred, basis)

dynmat_data = ph.get_dynamical_matrix_at_q([0,0,0])
eigens = ph.get_qpoints_dict()['eigenvectors'][0].real
freqs = ph.get_qpoints_dict()['frequencies'][0]

# Abinit use atomic units
if (args.calc=='abinit'):
    lUnits=1.889725989
else:
    lUnits=1

displcart=np.zeros(natom*3).reshape(natom,3)

for i in range(natom):
    for n in range(len(args.modnum.split())):
        for l in range(3):
            displcart[i][l]+=eigens[i*3+l, int(args.modnum.split()[n])]*sqrt(1/masses[i])/lUnits*weights[n]


print("Displacement vector: ")
for i in range(len(displcart)):
    print("%2d %s" % ((i+1), ''.join("%9.5f " % c for c in displcart[i])))

for a in args.ampl.split():
    displ=np.zeros(natom*3).reshape(natom,3)
    for i in range(natom):
        for l in range(3):
            displ[i][l]=displcart[i][l]*float(a)
    displred=cart2direct(displ, basis)
    savestruct(args.calc, float(a), basis, (xred+displred),chemel,chemnum)
    print('Struct saved.  Ampl %f Freqs: %s' % (float(a),
                                                ''.join('%9.5f ' % freqs[int(mm)] for mm in  args.modnum.split() )))

try:
    out_fh = open('mode_scan.xyz', 'w')
except IOError:
    print("ERROR Couldn't open output file for writing, exiting...\n")
    sys.exit(1)

cartpos=direct2cart(xred, basis)

out_fh.write('%d\n' % natom)

out_fh.write ('jmolscript: load "" {1 1 1} UNITCELL ['+
                        ','.join('%8.7f' % (b/lUnits) for b in basis.flatten()) + ']; vector on; vector 5; vector scale 15\n')

for i in range(natom):
    out_fh.write('%s '% chemel[i])
    out_fh.write(' '.join(' % 11.8f' % (cartpos[i][l]/lUnits) for l in range(3)))
    out_fh.write(' '.join(' % 11.8f' % displcart[i][l] for l in range(3)))
    out_fh.write('\n')

out_fh.write('Modes %s\n' % ''.join('%d ' % int(m) for m in args.modnum.split()) )