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


import phonopy
from phonopy.interface.castep import write_castep
from phonopy.interface.calculator import read_crystal_structure
from phonopy.structure.atoms import PhonopyAtoms as Atoms
import numpy as np
import re
import sys
from math import pi
from math import exp
import os
import argparse
import spglib as spl
import seekpath

parser = argparse.ArgumentParser(description='The program is to get primitive cell from POSCAR and plot k-path')


parser.add_argument("-i", "--input", action="store", type=str, dest="in_fn", default='POSCAR', help="Structure filename (default POSCAR)")
parser.add_argument("--calc", action="store", type=str, dest="calc", default='vasp', help="Calculator (vasp,abinit,castep)")
parser.add_argument("--tolerance", action="store", type=float, dest="tol", default=1e-05, help="Symmetry tolerance")

args = parser.parse_args()

cell, optional_structure_info = read_crystal_structure(
        filename=args.in_fn,
        interface_mode=args.calc)

print("Cell:")
for v in cell.cell:
    print("% 9.7f % 9.7f % 9.7f" % (v[0],v[1],v[2]))
   

print("Species:")
print(cell.symbols)
print('Atomic positions:')
print(cell.scaled_positions)

structure=[cell.cell,cell.scaled_positions,cell.numbers]

res=seekpath.getpaths.get_path(structure, with_time_reversal=True, recipe='hpkot', threshold=1e-07, symprec=1e-05, angle_tolerance=-1.0)
print("Output data:")
print("Space group: %s" %spl.get_spacegroup(cell=(cell.cell,cell.scaled_positions,cell.numbers), symprec=1e-3))
print('Primitive lattice:')

b=[]
for v in res['primitive_lattice']:
    print("a1= % 9.7f % 9.7f % 9.7f" % (v[0],v[1],v[2]))
    b.append([v[0],v[1],v[2]])

basis=np.array(b).reshape(3,3)
cvol=np.dot(basis[0],np.cross(basis[1],basis[2]))
b1=np.cross(basis[1],basis[2])/cvol*pi
b2=np.cross(basis[2],basis[0])/cvol*pi
b3=np.cross(basis[0],basis[1])/cvol*pi

print('unit cell volume=%7.4f' % cvol)
print('Reciprocal basis pi units:')
print(b1)
print(b2)
print(b3)
print('Reciprocal basis reciprocal units:')
print(b1/pi)
print(b2/pi)
print(b3/pi)


print('Atomic positions:')
for i in range(len(res['primitive_positions'])):
    print("%2d  %s" % (res['primitive_types'][i], "".join("  % 12.9f" % val for val in res['primitive_positions'][i])))

print('Crystallographic lattice:')
for v in res['conv_lattice']:
    print("a1= % 9.7f % 9.7f % 9.7f" % (v[0],v[1],v[2]))
print('Points:')
for pt in res['point_coords'].items():
    if("gamma" in pt[0].lower()):
        continue
    print("%s  %3s" % (" ".join("% 9.6f" % p for p in pt[1]),pt[0]))

print('PATH IN LETTERS:')
print (res['path'])
print('PATH IN COORDINATES:\n        From                     To')
for seg in res['path']:
    print("%s   %s" % ("".join("%5.4f " % c for c in res['point_coords'][seg[0]]), 
                       "".join("%5.4f " % c for c in res['point_coords'][seg[1]])))

