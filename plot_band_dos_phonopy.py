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
# This script plots phonon band structure and density of states
#
# The code is based on PHONOPY API https://phonopy.github.io/phonopy/
#
# Author: Eugene Roginskii

import matplotlib.pyplot as plt
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from phonopy.phonon.band_structure import get_band_qpoints_by_seekpath
import argparse
import seekpath
import numpy as np
import sys

parser = argparse.ArgumentParser(description='''This script plots phonon\
 band structure along with projected density of states obtained with\
 VASP/CASTEP/CP2K/ABINIT/CRYSTAL codes''')
parser.add_argument("-i", "--input", action = "store", type = str,
                    dest = "str_fn", help = '''Input filename with
                    structure (POSCAR/castep.cell/input.inp/input.abi)''')
parser.add_argument("-o", "--output", action="store",
                    type = str, dest = "out_fn", help = "Output filename")
parser.add_argument("--dim", action = "store", type = str,
                    dest = "dim", help = "Supercell matrix")
parser.add_argument("-f", action = "store", type = str, dest = "fsetfn",
                    default = 'FORCE_SETS', help = "Force_sets filename")
parser.add_argument("--readfc", dest = "read_force_constants",
                    action = "store_true", default = False,
                    help = "Read FORCE_CONSTANTS")
parser.add_argument("--nac", dest = "nac", action = "store_true",
                    default = False, help = "Include non-analytical term")
parser.add_argument("-s", "--soft", dest = "calc", action = "store",
                    help="Calculator software vasp/castep/cp2k/ABINIT/crystal")
parser.add_argument("-n", "--numpoints", dest = "npts",
                    action = "store", type = int, default = 51,
                    help = "Number of points for each segment")
parser.add_argument("--path", action = "store", type = str, dest = "path",
                    help = "Path in BZ as sequence of points. If skipped then\
 will be generated automatically with seekpath")
parser.add_argument("--pathsym", action = "store", type = str, dest = "pathsym",
                    help = "Symbols for path in BZ. If skipped then\
 will be generated automatically with seekpath")
parser.add_argument("-r", "--range", action = "store",
                    type = str, dest = "range",
                    help = "Range of frequencies")

args = parser.parse_args()

params = {
'mathtext.default': 'regular',
'axes.linewidth': 1.2,
'axes.edgecolor': 'Black',
'figure.dpi' : 150,
'figure.figsize' : [8.0, 4.0]  #  A4 -- 8.27 x 11.69
}

plt.rcParams.update(params)

if(args.calc == None):
    print('Error. Calculator name is missed')
    sys.exit(1)
calc=""
if(args.calc=="castep"):
    factorcm=521.47083
    calc='castep'
elif(args.calc=="vasp"):
    factorcm=521.47083
    calc='vasp'
elif(args.calc=="cp2k"):
    factorcm=3739.4256800756
    calc='cp2k'
elif(args.calc=="cp2kv6"):
    factorcm=3739.4256800756
    calc='cp2k'
elif(args.calc=="abinit"):
    factorcm=716.85192105135115965589
    calc='abinit'
elif(args.calc=="crystal"):
    factorcm=521.47083
    calc='crystal'
else:
    print('Wrong calculator name %s' % args.calc)
    sys.exit(1)

dim = []
if (args.dim):
    for d in args.dim.split():
        dim.append(int(d))
else:
    dim = [1, 0, 0, 0, 1, 0, 0, 0, 1]

if (len(dim) > 3):
    dim = np.array(dim).reshape(3, 3)

if (args.read_force_constants):
    print("Read force contants from FORCE_CONSTANTS file")
    ph = phonopy.load(supercell_matrix = dim,
                  primitive_matrix = [1, 0, 0, 0, 1, 0, 0, 0, 1],
                  unitcell_filename = args.str_fn,
                  calculator = calc, factor = factorcm,
                  force_constants_filename = 'FORCE_CONSTANTS',
                  symmetrize_fc = True, is_nac = args.nac,
                  born_filename = 'BORN')
else:
    ph = phonopy.load(supercell_matrix = dim,
                  primitive_matrix = [1, 0, 0, 0, 1, 0, 0, 0, 1],
                  unitcell_filename = args.str_fn,
                  calculator = calc, factor = factorcm,
                  force_sets_filename = args.fsetfn,
                  symmetrize_fc = True, is_nac = args.nac,
                  born_filename = 'BORN')

# species=ph.primitive.get_chemical_symbols()
numbers=ph.primitive.get_atomic_numbers()
basis = ph.primitive.get_cell()
xred = ph.primitive.get_scaled_positions()
natom = ph.primitive.get_number_of_atoms()
cell = [basis, xred, numbers] 
path = []
labels=[]

if (args.path):
    if (len(args.pathsym.split()) != (len(args.path.split()) / 3)):
        print('Error. Number of points and symbols in path not equivalent.')
        sys.exit(1)
    for i in range(len(args.pathsym.split())):
        if ('gamma' in args.pathsym.split()[i].lower()):
            labels.append(r'$\Gamma$')
        else:
            labels.append(args.pathsym.split()[i])
        path.append([float(eval(pt)) for pt in args.path.split()[i*3:i*3+3]])

# Construct path using seekpath. Only first segments without interruptions
else:    
    res = seekpath.getpaths.get_path(cell, with_time_reversal = True,
                                   recipe = 'hpkot', threshold = 1e-07,
                                   symprec = 1e-05, angle_tolerance = -1.0)

    for i in range(len(res['path']) - 1):
        if (res['path'][i][1] != res['path'][i+1][0]):
            if ('gamma' in res['path'][i][0].lower()):
                labels.append(r'$\Gamma$')
            else:
                labels.append(res['path'][i][1])
            path.append(res['point_coords'][res['path'][i][1]])
            
            break
        if ('gamma' in res['path'][i][0].lower()):
            labels.append(r'$\Gamma$')
        else:
            labels.append(res['path'][i][0])
        path.append(res['point_coords'][res['path'][i][0]])

print(labels)
print(path)

# Alternative automatic way
# (qpoints,
#  labels,
#  connections) = get_band_qpoints_by_seekpath(ph.primitive, 
#                                                  npts,
#                                                  is_const_interval = False) 

qpoints, connections = get_band_qpoints_and_path_connections([path], npoints=args.npts)
ph.run_band_structure(qpoints, with_eigenvectors = True,
                      is_band_connection = True,
                      path_connections = connections,
                      labels = labels)

band_dict = ph.get_band_structure_dict()
ph.run_mesh()
ph.run_total_dos(sigma = 1.5)
dos = ph.get_total_DOS()

dist = band_dict['distances']  # dist is an array of segments array. Corresponds to x-axis
freq = band_dict['frequencies']
kpt = qpoints
eigs = band_dict['eigenvectors']

plt.rcParams.update(params)
fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]},
                       sharey = True,
                       constrained_layout = True)

# Plot Bandstructure Iterate over each segment
for i in range(len(dist)):
# Iteraton over each band
    for nbnd in range(len(freq[i][0])):
        x=[]
        y=[]
        for j in range(len(dist[i])):
            x.append(dist[i][j]) # j-point in the i-th segment 
            y.append(freq[i][j][nbnd]) # frequency at j-point in the i-th segment for band nbnd
        ax[0].plot(x, y, c='black', lw=1.0, alpha=0.7)

ax[0].set_ylabel(r'Frequency, cm$^{-1}$', fontsize=12)
ax[0].set_xlim([dist[0][0], dist[len(dist)-1][-1]])
xticks=[dist[i][0] for i in range(len(dist))]
xticks.append(dist[len(dist)-1][-1])
ax[0].set_xticks(xticks)
ax[0].set_xticklabels(labels)
ax[0].grid(which='major', axis='both', c='gray', linestyle='-', alpha=0.8)
ax[1].grid(which='major', axis='both', c='gray', linestyle='-', alpha=0.8)

ax[1].plot(dos[1], dos[0], c='black', lw=1.0, alpha=0.7)

if (args.range):
  if (len(args.range.split())==2):
    ax[0].set_ylim(float(args.range.split()[0]),
                   float(args.range.split()[1]))

plt.savefig(args.out_fn, bbox_inches='tight')

