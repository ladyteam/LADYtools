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
#
# Example: ./plot_phbands_and_dos.py -i rf2_5.out_PHBST.nc -d rf2_5.out_PHDOS.nc -o band_dos.pdf --nac --range="0 250" -s 2
#
# If the following error arise:
# ValueError: Non analytical contribution has not been calculated for
# reduced direction [0.   0.05 0.  ]
# And for example the crystal is Face-centere:
# Real(R)+Recip(G) space primitive vectors, cartesian coordinates:
# R(1)=4.2736156 -7.2385958  0.0000000  G(1)=0.1169970 -0.0690742  0.0000000
# R(2)=4.2736156  7.2385958  0.0000000  G(2)=0.1169970  0.0690742  0.0000000
# R(3)=0.0000000  0.0000000  8.5491209  G(3)=0.0000000  0.0000000  0.1169711
# Then you simple have to add G(2) vector in qph2l array:
# qph2l = 0.1169970  0.0690742  0.0000000 0.0
# Evaluate in cartesian coordinate using numpy:
# b=np.array([[0.1169970, -0.0690742,  0.0000000],[ 0.1169970, 0.0690742,  0.0000000],[0.0000000,  0.0000000,  0.1169711]])
# p=np.array([-0.0263,  0.0263,  0.0263])
# np.dot(p,b)
# [0.         0.0036333  0.00307634]
r"""
Phonon Band structures
======================

Script to plot the phonon band structure.
"""
from abipy.abilab import abiopen
from abipy.dfpt.phonons import PhononDos
from scipy.signal import savgol_filter
import abipy.data as abidata
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='The program to plot\
 Phonon Bandstructure and DOS')


parser.add_argument("-i", "--input", action = "store", type = str,
                    dest="iphon_fn",
                    help="Anaddb output filename *PHBST.nc with\
 phonon bandstructure data")
parser.add_argument("-d", "--dos", action = "store",
                    type = str, dest = "dos_fn",
                    help = "Anaddb output filename *PHDOS.nc with\
 phonon DOS data")
parser.add_argument("-o", "--out", action = "store", type = str,
                    dest = "out_fn",
                    help = "output filename")
parser.add_argument("-t", "--title", action = "store",
                    type = str, dest = "title",
                    default = "Phonon Bandstructure in cm-1",
                    help = "Title")
parser.add_argument("--nac", action = "store_true", dest = "nac",
                    default = "Phonon Bandstructure in cm-1",
                    help = "Take non-analitical term into account\
 anaddb.nc file should exist and contain NLO data (see help on\
 qph2l variable of anaddb")
parser.add_argument("-r", "--range", action = "store",
                    type = str, dest = "range",
                    help = "Range of frequencies")
parser.add_argument("-s", "--smooth", dest = 'smp', default = 0, type = int,
                    help="Window Savitzky-Golay filter length for smoothing")

args = parser.parse_args()

# Open the PHBST file produced by anaddb and get the phonon bands.
with abiopen(abidata.ref_file(os.path.join(os.getcwd(),
                                           args.iphon_fn))) as ncfile:
    phbands = ncfile.phbands

# Read the Phonon DOS from the netcd file produced by anaddb (prtdos 2)
with abiopen(abidata.ref_file(os.path.join(os.getcwd(),
                                           args.dos_fn))) as ncfile:
    phdos = ncfile.phdos

if (args.nac == True):
    phbands.read_non_anal_from_file(abidata.ref_file(os.path.join(os.getcwd(),
                                                     "anaddb.nc")))
params = {
'mathtext.default': 'regular',
'axes.linewidth': 1.2,
'axes.edgecolor': 'Black',
'figure.dpi' : 150,
'figure.figsize' : [8.0, 4.0]  #  A4 -- 8.27 x 11.69
}

plt.rcParams.update(params)

fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]},
                       sharey = True,
                       constrained_layout = True)


branches =  phbands.plot_ax(ax[0], branch=None, units = 'cm-1',
                            match_bands=False,
                            c = 'black', lw = 1.0, alpha = 0.7)
# Spline 
if (args.smp > 0):
    for b in branches:
#   Do not spline accoustic branches
        if (np.min(b.get_ydata()) > 5):
            y = savgol_filter(b.get_ydata(), (args.smp*2+1), 3)
            b.set_ydata(y)

ticks, labels = phbands._make_ticks_and_labels(None)

if ticks:
    # Don't show label if previous k-point is the same.
    for il in range(1, len(labels)):
        if labels[il] == labels[il-1]: labels[il] = ""
    ax[0].set_xticks(ticks, minor = False)
    ax[0].set_xticklabels(labels, fontdict = None, minor = False, size = 12)
    #print("ticks", len(ticks), ticks)
    ax[0].set_xlim(ticks[0], ticks[-1])

phdos.plot_dos_idos(ax[1], what = "d", units = 'cm-1', exchange_xy=True, 
                    c = 'black', lw = 1.0, alpha = 0.7)
if (len(args.range.split())==2):
    ax[0].set_ylim(float(args.range.split()[0]),
                   float(args.range.split()[1]))

ax[0].grid(which='major', axis='both', c='gray', linestyle='-', alpha=0.8)
ax[0].set_ylabel(r'Frequency, cm$^{-1}$', fontsize=12)

ax[1].grid(which='major', axis='both', c='gray', linestyle='-', alpha=0.8)

if (args.out_fn):
    plt.savefig(args.out_fn, bbox_inches='tight')
else:
    plt.show()
