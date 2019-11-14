#!/bin/bash
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
# Conver CUBE isosurface data into JVXL format and pack all data into ZIP
# using JMOL software  (http://jmol.sourceforge.net)
#
# Note: In order to save isosurface data one have to put specific jmolscript
# command inside CUBE input file (second line) For example: 
# jmolscript: isosurface s1 sign cyan yellow cutoff 1.000000 "" color translucent 0.45
# or provide such command  with -c parameter. For example:
# cube2jvxl.sh -i input.cube -c 'isosurface s1 sign cyan yellow cutoff 1.0 "" color translucent 0.45'
# Notice on single and double quote. If you whant to use only double quotes don't forget to slash it:
# cube2jvxl.sh -i input.cube -c "isosurface s1 sign cyan yellow cutoff 1.0 \"\" color translucent 0.45"
# 
# First string after comments is the path to Jmol executable. Set it according to your system installation.
#
# Example how to convert all .cube files in the directory:
# $ for i in *.cube; do ./cube2jvxl.sh -i $i; done
#
# Author: Eugene Roginskii

#Do not modify here, see below
if [[ -n $JMOLEXE ]]; then
    Jmolexe=$JMOLEXE
    echo "Jmol path is set according to JMOLEXE variable"
else
# Modify line below to provide full path to Jmol executable
    Jmolexe=/opt/jmol/jmol.sh
fi



if [[ ! -x ${Jmolexe} ]]; then
    echo ${Jmolexe} not found
    if [[ -n `which jmol.sh 2>&1 | xargs file  | grep "No such file or directory"` ]]; then
        echo "Jmol not found. Provide full path to Jmol executable"
        echo " by setting JMOLEXE global variable in your .bashrc"
        exit 1
    else
        Jmolexe=`which jmol.sh` || exit 1
    fi
fi

atom_data=(X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge \
As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy \
Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm \
Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Uut Uuq Uup Uuh Uus Uuo)

usage()
{
cat << EOF >&2
usage: $0 options

Converter CUBE isosurface file format to jvxl one
Note: In order to save isosurface data one have to put specific jmolscript
command inside CUBE input file (second line) For example: 
jmolscript: isosurface s1 sign cyan yellow cutoff 1.000000 "" color translucent 0.45
or provide such command  with -c parameter. For example:
cube2jvxl.sh -i input.cube -c 'isosurface s1 sign cyan yellow cutoff 1.0 "" color translucent 0.45'
Notice on single and double quote. If you whant to use only double quotes don't forget to slash it:
cube2jvxl.sh -i input.cube -c "isosurface s1 sign cyan yellow cutoff 1.0 \"\" color translucent 0.45"
See isosurface help of Jmol documentation.


OPTIONS:
   -h      Show this message
   -i      Input filename.
   -c      Jmolscript command
   -d      Debug
EOF
}

while getopts "hdi:c:" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         d)
             DEBUG="TRUE"
             ;;
         i)
             FN=$OPTARG
             ;;
         c)
             CMD=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z "$FN" ]]; then
    echo "Error: No input file name was given"
    usage
    exit 1
fi

if [[ ! -e "$FN" ]]; then
    echo "Error: Input file does not exist"
    usage
    exit 1
fi


#temporary jmolscript file for converting isourface in JVXL format
sptfn=jmolscript_$(((RANDOM%10000)+1000))".spt"

# isosurface jvxl filename
o=${FN%.*}".jvxl"
# structural data filename
xyz=${FN%.*}".xyz"
# write script and do isosurface conversion
echo "load \"$FN\"" > $sptfn
if [[ -n $CMD ]]; then
    echo $CMD >> $sptfn
fi
echo "write isosurface \"$o\"" >> $sptfn

echo "Converting file $FN"
sync
$Jmolexe -ions $sptfn

if [[ $DEBUG != "TRUE" ]]; then
    rm $sptfn
fi

sptfn=jmolscript_$(((RANDOM%10000)+1000))".spt"

echo "load \"$FN\"" > $sptfn
# Check if compressed
type=`file ${FN} | grep zip`
if [[ -z "$type" ]]; then
    echo "Proceed noncompressed data"
    CAT=cat
else 
    echo "Proceed compressed data"
    CAT=zcat
fi
# Read atomic data from zipped cube file and write in xyz format
nat=`$CAT $FN | head -n3 | tail -n1 | awk '{print($1)}'`
((nl=nat+6))
echo $nat > $xyz
$CAT $FN | head -n2 | tail -n1 | awk -F";" '{print $1}' >> $xyz
$CAT $FN | head -n $nl | tail -n $nat | \
awk -v var="${atom_data[*]}" 'BEGIN{split(var,list," ");a2b=0.52917721092}{printf("%3s % 10.6f % 10.6f % 10.6f\n", list[$1+1],$3*a2b,$4*a2b,$5*a2b)}' >> $xyz

#Do zip
echo "load \"$xyz\"" > $sptfn
echo "background [120,180,180]" >> $sptfn
echo "isosurface s1 sign cyan yellow \"$o\" color translucent 0.45" >> $sptfn
echo "write ZIP  \"${FN%.*}.zip\"" >> $sptfn
sync
$Jmolexe -ions $sptfn
if [[ $DEBUG != "TRUE" ]]; then
    rm  $sptfn $xyz $o
fi
