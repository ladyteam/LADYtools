#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Test with ./plot_band_part.py -i init/displ.in --dim "2 2 2" --band "1/3 1/3 0.0   1/2 0.0 0.0  0.0 0.0 0.0  0.0 0.0 0.5   1/3 1/3 1/2" --band-labels 'K  M $\Gamma$ A  H' --atomset "1 2 3 4 -5 -6 -7 -8" --nac --calculator "abinit" --factor=716.8519 --cmap-labels "AlN GaN" -o band_partial.pd
# ./plot_band_partial.py -i init/displ.in --dim "2 2 2" --band "0.0 0.0 0.0  0.0 0.0 0.5" --band-labels '$\Gamma$ A ' --atomset "1 2 3 4 -5 -6 -7 -8" --nac --calculator "abinit" --factor=716.8519 --cmap-labels "AlN GaN" -o band_partial_G-A.pdf --color-smooth=10 --np=301 --ylim="0 950"

import h5py
import matplotlib.pyplot as plt 
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
import argparse
import warnings
import matplotlib
import numpy as np
import sys
from math import sqrt
from scipy.signal import savgol_filter

warnings.filterwarnings("ignore") # Ignore all warnings

def displstrength(eig,atomsetp,atomsetn,masses):
    e=[]
    ep=[]
    en=[]
    for i in range(len(masses)):
        for j in range(3):
            e.append(eig[i*3+j]/sqrt(masses[i]))
# Check if current atom number persist in list atomsetp (positive list of atoms)
        for a in atomsetp:
            if(a==(i+1)):
                for j in range(3):
                    ep.append(eig[i*3+j]/sqrt(masses[i]))
# Check if current atom number persist in list atomsetn
        if(len(atomsetn)>0):
          for a in atomsetn:
            if(a==(i+1)):
                for j in range(3):
                    en.append(eig[i*3+j]/sqrt(masses[i]))

    enorm=np.linalg.norm(e)
    epnorm=np.linalg.norm(ep)
    if(len(atomsetn)==0):
        ennorm=0
    else:
        ennorm=np.linalg.norm(en)
    return(epnorm/enorm-ennorm/enorm)

parser = argparse.ArgumentParser(description='The program plots partial bandstructure')
parser.add_argument("-i", "--input", action="store", type=str, dest="in_fn", help="Input filename")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn", help="Output filename")
parser.add_argument("--calculator", action="store", type=str, dest="calc", help="Ab initio program (vasp, abinit, crystal)")
parser.add_argument("-f", "--forcesets", action="store", type=str, dest="fsetfn", default='FORCE_SETS', help="FORCE_SETS filename")
parser.add_argument("--factor", action="store", type=float, dest="factor", default=521.47083, help="Frequency factor. Default for cm-1 521.47083 (vasp), 716.8519 (abinit)")
parser.add_argument("--pa", action="store", type=str, dest="prim_mat", default="1 0 0 0 1 0 0 0 1", help="Transformation matrix for primitive cell")
parser.add_argument("--dim", action="store", type=str, dest="sc_mat", default="1 0 0 0 1 0 0 0 1", help="Trasformation matrix for supercell")
parser.add_argument("--band", action="store", type=str, dest="path", help="Set of k-points to plot dispertion")
parser.add_argument("--band-labels", action="store", type=str, dest="labels", help="Labels for k-points")
parser.add_argument("--cmap-labels", action="store", type=str, dest="cmlabels", default="", help="Labels for color map")
parser.add_argument("--atomset", action="store", type=str, dest="atset", 
                      help="Atomic numbers. Could be negative for negative intensity. Do not use only negative nubers. At least one number should be positive!")
parser.add_argument("--np", action="store", type=int, dest="npoints", default=51, help="Number of points for interpalation. Default is 51")
parser.add_argument('--nac', dest='nac', action='store_true')
parser.add_argument('--color-smooth', dest='cmapintpnp', default=0, type=int, help="Window length of Savitzky-Golay filter for smoothing of the colormap artifacts")
parser.add_argument("--ylim", action="store", type=str, dest="ylim", help="Limits for y axis")

args = parser.parse_args()

if (args.in_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

# Primitive and supercell_matrix
if(len(args.prim_mat.split())==9):
    pm=np.array([eval("".join("1.0*%s") % p) for p in args.prim_mat.split()]).reshape(3,3)
else:
    print('Error. Number of components in primitive matrix less then 9')
    sys.exit(1)
if(len(args.sc_mat.split())==9):
    sc=np.array([int(s) for s in args.sc_mat.split()]).reshape(3,3)
elif(len(args.sc_mat.split())==3): 
    sc=np.array([int(s) for s in args.sc_mat.split()])
else:
    print('Error. Number of components in supercell matrix less then 9')
    sys.exit(1)

path=np.array([eval("".join("1.0*%s") % p) for p in args.path.split()]).reshape(int(len(args.path.split())/3),3)

atomsetp=[]
atomsetn=[]
for a in args.atset.split():
    try:
        if(int(a)>0):
            atomsetp.append(int(a))
        else:
            atomsetn.append(-1*int(a))
    except:
        print('Error reading numbers of atoms')
        sys.exit(1)

if (len(atomsetp)==0):
    print('Error. The set of atomic numbers (positive) for partial Bandstructure is empty')
    sys.exit(1)

labels=[]
for t in args.labels.split():
    labels.append(t)

print(pm)
print(sc)

ph = phonopy.load(supercell_matrix=sc,
                  primitive_matrix=pm,
                  unitcell_filename=args.in_fn,
                  is_nac=args.nac,  calculator=args.calc, factor=args.factor,
                  force_sets_filename=args.fsetfn)

print(ph.primitive)

qpoints, connections = get_band_qpoints_and_path_connections([path], npoints=args.npoints)
ph.run_band_structure(qpoints,with_eigenvectors=True,is_band_connection=True, path_connections=connections, labels=labels)
band_dict=ph.get_band_structure_dict()

dist=band_dict['distances']  # dist is an array of segments array. Corresponds to x-axis
freq=band_dict['frequencies']
kpt=qpoints
eigs=band_dict['eigenvectors']

# Plotting routine
params = {
'mathtext.default': 'regular',
'axes.linewidth': 1.2,
'axes.edgecolor': 'Black',
'figure.dpi' : 70
}

plt.rcParams.update(params)
fig, ax = plt.subplots() 

# Find min and max 
mini=0
maxi=0
intens=[[[0.0 for i in range(args.npoints)] for j in range(len(freq[0][0]))] for k in range(len(dist))]

# Plotting the band
# Iterate over each segment
for i in range(len(dist)):
# Iteraton over each band
    for nbnd in range(len(freq[i][0])):
        for j in range(len(dist[i])):
            intens[i][nbnd][j]=displstrength(eigs[i][j].T[nbnd].real,atomsetp,atomsetn,ph.primitive.get_masses())
            if(intens[i][nbnd][j]>maxi):
               maxi=intens[i][nbnd][j]
            if(intens[i][nbnd][j]<mini):
               mini=intens[i][nbnd][j]

# If only positive set
if(len(atomsetn)==0):
    mini=0
# If two sets negative and positive mixed states should be displayed in balanced color map
else:
    if(abs(maxi)>abs(mini)):
        mini=-abs(maxi)
    else:
        maxi=abs(mini)

print('mini=%f, maxi=%f' % (mini,maxi))
# Set manually minimal amd maximal color scale belo if needed by uncommenting the following 2 lines
#mini=-1
#maxi=1
norm = plt.Normalize(mini, maxi)

if (abs(mini)<0.1):
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#111111", "#0FFEF9"])
else:
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#FF000D", "#111111", "#0FFEF9"])

# Iterate over each segment
for i in range(len(dist)):
# Iteraton over each band

    for nbnd in range(len(freq[i][0])):
        x=[]
        y=[]
        if(args.cmapintpnp>2):
            # Use 3-d order polynom
            z=savgol_filter(intens[i][nbnd],(args.cmapintpnp*2+1), 3)
        else:
            z=intens[i][nbnd]

        for j in range(len(dist[i])):
            x.append(dist[i][j]) # j-point in the i-th segment 
            y.append(freq[i][j][nbnd]) # frequency at j-point in the i-th segment for band nbnd

        plt.scatter(x, y, c=z, cmap=cmap, marker='s', s=1, norm=norm)

# Labels
ax.set_ylabel(r'$\nu$, cm$^{-1}$', fontsize=12)
ax.set_xlim([dist[0][0],dist[len(dist)-1][-1]])
# Set y limit if set
if (args.ylim):
    ax.set_ylim(float(args.ylim.split()[0]),float(args.ylim.split()[1]))
xticks=[dist[i][0] for i in range(len(dist))]
xticks.append(dist[len(dist)-1][-1])
ax.set_xticks(xticks)
print(labels)
ax.set_xticklabels(labels)
ax.grid(which='major', axis='x', c='gray', linestyle='--', alpha=0.8)

cbar=plt.colorbar(cax = fig.add_axes([0.92, 0.4, 0.02, 0.48]))
cbar.ax.get_yaxis().set_ticks([])

if(len(args.cmlabels)<1):
    if(len(atomsetn)<1):
        cbar.ax.text(1.2, 0.0, '0', ha='left', va='center')
        cbar.ax.text(1.2, 1.0, ''.join('%.1f' % maxi), ha='left', va='center')
    else:
        cbar.ax.text(1.2, 0.0, ''.join('%.1f' % mini), ha='left', va='center')
        cbar.ax.text(1.2, 1.0, ''.join('%.1f' % maxi), ha='left', va='center')
else:
    if(len(atomsetn)<1):
        cbar.ax.text(1.2, 0.0, '0', ha='left', va='center')
        cbar.ax.text(1.2, 1.0, args.cmlabels.split()[0], ha='left', va='center')
    else:
        if(len(args.cmlabels.split())<2):
            cbar.ax.text(1.2, 0.0, ''.join('%.1f' % mini), ha='left', va='center')
            cbar.ax.text(1.2, 1.0, ''.join('%.1f' % maxi), ha='left', va='center')
        else:
            cbar.ax.text(1.2, 0.0, args.cmlabels.split()[0], ha='left', va='center')
            cbar.ax.text(1.2, 1.0, args.cmlabels.split()[1], ha='left', va='center')

plt.savefig(args.out_fn, bbox_inches='tight')