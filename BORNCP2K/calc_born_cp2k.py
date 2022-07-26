#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
#
# The code is based on PHONOPY API https://phonopy.github.io/phonopy/
#
# Author: Eugene Roginskii

def genxyz(fn, basis, species, cart, comment=''):
#    print("Generating file %s" % fn)
    try:
        fh = open(fn, 'w')
    except OSError as err:
        print("ERROR Couldn't open output file %s for writing: %err" % (fn, err))
        return(-1)
    fh.write('%d\n' % len(species))
    fh.write('%s\n' % comment)
    for i in range(len(species)):
        fh.write("%s  %s\n" % (species[i],
                             "".join("    % 15.10f" % c for c in cart[i])))
def get_cp2kver(fn):
    try:
        fh = open(fn, 'r')
        for line in fh:
            if ('CP2K| version string' in line):
                ver = float(line.split('CP2K version')[1])
                return(int(ver))
                break
    except OSError as err:
        print('Error opening file %s: %s' % (fn, err))


def get_dipole(fn):
    dipole=[]
    Debye2au = 2.54174
    xyz=['X=', 'Y=', 'Z=']
    try:
        fh = open(fn, 'r')
        for line in fh:
            if ('Dipole moment [Debye]' in line):
                break
        line = fh.readline()
        for i in range(3):
            dipole.append(float(line.split(xyz[i])[1].split()[0])/Debye2au)
    except OSError as err:
        print('Error opening file %s: %s' % (fn, err))
    return(dipole)

def get_epsilon_cp2k(fn, ucvol):
    from numpy import pi
    epsilon=np.zeros(9)
    try:
        fh = open(fn, 'r')
    except IOError:
        print("ERROR Couldn't open output file %s for reading" % fn)
        return(-1)
    for line in fh:
        if ('Polarizability tensor [a.u.]' in line):
            break
# Components sequence in CP2K output
# xx,yy,zz
# xy,xz,yz
# yx,zx,zy
# Components sequence in epsilon array:
#                     0    1    2    3    4    5    6    7    8
#      9-components:  xx   xy   xz   yx   yy   yz   zx   zy   zz
    line=fh.readline()
    for i in range(3):
        epsilon[i*4]=float(line.split()[i+2])/ucvol*4*pi + 1

    line=fh.readline()
    for i in range(2):
        epsilon[1+i]=float(line.split()[i+2])/ucvol*4*pi
    epsilon[5]=float(line.split()[4])/ucvol*4*pi

    line=fh.readline()
    epsilon[3]=float(line.split()[2])/ucvol*4*pi
    for i in range(2):
        epsilon[6+i]=float(line.split()[i+3])/ucvol*4*pi
# Symmitrize
    epsilon[1]=(epsilon[1]+epsilon[3])/2
    epsilon[3]=epsilon[1]
    epsilon[2]=(epsilon[2]+epsilon[6])/2
    epsilon[6]=epsilon[2]
    epsilon[5]=(epsilon[5]+epsilon[7])/2
    epsilon[7]=epsilon[5]
#    print(epsilon)

    return(epsilon)

def get_epsilon_cp2kv6(fn, ucvol):
    from numpy import pi
    epsilon=np.zeros(9)
    print('cp2kv6')
    try:
        fh = open(fn, 'r')
    except IOError:
        print("ERROR Couldn't open output file %s for reading" % fn)
        return(-1)
    for line in fh:
        if ('POLARIZABILITY TENSOR (atomic units)' in line):
            print(line)
            break

    line=fh.readline()
    for i in range(3):
        epsilon[i*4]=float(line.split()[i+1])/ucvol*4*pi

    line=fh.readline()
    for i in range(2):
        epsilon[1+i]=float(line.split()[i+1])/ucvol*4*pi
    epsilon[5]=float(line.split()[3])/ucvol*4*pi

    line=fh.readline()
    epsilon[3]=float(line.split()[1])/ucvol*4*pi
    for i in range(2):
        epsilon[6+i]=float(line.split()[i+2])/ucvol*4*pi
# Symmitrize
    epsilon[1]=(epsilon[1]+epsilon[3])/2
    epsilon[3]=epsilon[1]
    epsilon[2]=(epsilon[2]+epsilon[6])/2
    epsilon[6]=epsilon[2]
    epsilon[5]=(epsilon[5]+epsilon[7])/2
    epsilon[7]=epsilon[5]
#    print(epsilon)

    return(epsilon)

from ase import Atoms
from ase.data import (covalent_radii,
                      atomic_numbers)
from ase import io
from ase.neighborlist import NeighborList
from phonopy.interface.cp2k import read_cp2k
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.structure.symmetry import Symmetry
import numpy as np
import argparse
from sys import exit
from copy import deepcopy
from os import path

parser = argparse.ArgumentParser(description='The code generates input files and compute Born charges and dielectric tensor')

parser.add_argument("-i", "--input", action="store", type=str, dest="str_fn", help="Input filename with structure")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn", help="Output filename prefix")

parser.add_argument("-p", "--policy", action="store", type=str, dest="policy", default='displ',
  help="Script modes. 'displ' -- generate input files; 'calc' -- calculate Born charges and Dielectric tensor")
parser.add_argument("-D", "--delta", action="store", type=float, dest="delta", default=0.01, help="Atomic displacement in Angstrom")

args = parser.parse_args()

xyz=['x', 'y', 'z']

b2a = 0.52917721092
atoms = read_cp2k(args.str_fn)
basis = atoms[0].get_cell()
cart = atoms[0].get_positions()
ucvol = atoms[0].get_volume()/(b2a**3)
chemsyms = atoms[0].get_chemical_symbols()
symmetry = Symmetry(atoms[0])


if (args.policy == 'displ'):
    print('Generating file with ideal structure for calculations of dielectric tensor')
    fn = ''.join('%s_ideal.xyz' % args.out_fn)
    genxyz(fn, basis, chemsyms, cart, 
                   comment='Ideal structure')
    for i in symmetry.get_independent_atoms():
        for d in range(3):
            fn = ''.join('%s_%d_%s_plus.xyz' % (args.out_fn, i, xyz[d]))
            print('Generating file with atom number %d shifted along direction %s' % (i, xyz[d]))
            cartdist = deepcopy(cart)
            cartdist[i][d] += args.delta
            genxyz(fn, basis, chemsyms, cartdist, 
                   comment=''.join('Atom number %d shifted forward along %s' %
                                   (i, xyz[d])))
            fn = ''.join('%s_%d_%s_minus.xyz' % (args.out_fn, i, xyz[d]))
            print('Generating file with atom number %d shifted along direction %s' % (i, xyz[d]))
            cartdist = deepcopy(cart)
            cartdist[i][d] -= args.delta
            genxyz(fn, basis, chemsyms, cartdist, 
                   comment=''.join('Atom number %d shifted backward along %s' %
                                   (i, xyz[d])))
else:
    born = []
    for i in symmetry.get_independent_atoms():
        borntmp = np.zeros(9).reshape(3,3)
        for d in range(3):
            fn = path.join('%s_%d_%s_plus' % (args.out_fn, i, xyz[d]),
                           'dipole.out')
            print(fn)
            dipolep = get_dipole(fn)
            print(dipolep)
            fn = path.join('%s_%d_%s_minus' % (args.out_fn, i, xyz[d]),
                           'dipole.out')
            print(fn)
            dipolem = get_dipole(fn)
            print(dipolem)
            for j in range(3):
                borntmp[j][d]=(dipolep[j]-dipolem[j])/2/(args.delta/b2a)
        print(borntmp)
        born.append(borntmp)
# epsilon
    fn = path.join('%s_ideal' % args.out_fn, 'polar.out')
    if (get_cp2kver(fn) <= 6):
        epsilon = get_epsilon_cp2kv6(fn, ucvol)
    else:
        epsilon = get_epsilon_cp2k(fn, ucvol)
    print('epsilon')
    print(epsilon)
