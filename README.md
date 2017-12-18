# LADYtools
LAttice DYnamic tools

The set of helpful tools for ab initio calculation software packages like ABINIT VASP CRYSTAL etc.

1. scan_alnog_mode_vasp.py - the program to distort geometry along selected normal mode eigenvector. VASP POSCAR format.
1a. make_scan_dirs.sh - script to copy geometry and other VASP files in new directories
1b. To plot table Energy vs shift one may use bash command:
VASP:
for i in SCAN*; do cd $i; echo -n "${i##SCAN-} "; grep TOTEN OUTCAR | tail -n1 | awk '{print ($5)}'; cd ..; done
ABINIT:
for i in displ*; do cd $i; echo -n "${i##displ-} "; grep etotal displ.out | tail -n1 | awk '{print ($2)}'; cd ..; done
2. Simple script to extract FORCE_CONSTANTS from abinit output file in PHONOPY software format.
3. Two utilites make_repr.py and make_repr_table.py to plot characters table and Irreps using ABINIT anaddb utility.
4. abinit-patterns.tar.bz2 a set of patterns for ABINIT.
5. extract_raman_abinit.py  Extract raman tensor and rotational invariants (Long's notation) from ABINIT anaddb output
6. plot_ramanspectrum.py   Plot raman spectrum with given intensities and line halfwidths.