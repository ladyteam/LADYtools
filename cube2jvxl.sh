#!/bin/bash
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
# Isosurface CUBE to JVXL file format converter using JMOL (www.jmol.org)
# Will convert all files with .cube suffix. Warning! This version is for zipped
# CUBE files. Replace all zcat by cat for unzipped CUBE files.
#
# Author: Eugene Roginskii

atom_data=(X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge \
As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy \
Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm \
Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Uut Uuq Uup Uuh Uus Uuo)

#temporary jmolscript file for converting isourface in JVXL format
sptfn=jmolscript_$$".spt"

for i in *.cube; do
# isosurface jvxl filename
    o=${i%cube}"jvxl"
# structural data filename
    xyz=${i%cube}"xyz"
# persistent jmolscript filename
    scriptfn=${i%cube}"spt"
# write script and do isosurface conversion
    echo "load \"$i\"" > $sptfn
    echo "write isosurface \"$o\"" >> $sptfn

    echo "Converting file $i"
    /opt/jmol/jmol.sh -ions $sptfn
    rm $sptfn
# Read atomic data from zipped cube file and write in xyz format
    nat=`zcat wannier90_00002.cube | head -n3 | tail -n1 | awk '{print($1)}'`
    ((nl=nat+6))
    echo $nat > $xyz
    zcat wannier90_00002.cube | head -n2 | tail -n1 | awk -F";" '{print $1}' >> $xyz
    zcat wannier90_00002.cube | head -n $nl | tail -n $nat | \
    awk -v var="${atom_data[*]}" 'BEGIN{split(var,list," ");a2b=0.52917721092}{printf("%3s % 10.6f % 10.6f % 10.6f\n", list[$1+1],$3*a2b,$4*a2b,$5*a2b)}' >> $xyz
# Write jmolscript file
    echo "load \"$xyz\"" > $scriptfn
    echo "isosurface s1 sign cyan yellow cutoff 1 \"$o\" color translucent 0.45" >> $scriptfn
done
