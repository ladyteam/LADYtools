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