#!/bin/bash


for i in POSCAR-*; do
    num=${i##*-}
    dir="SCAN-"$num
    if [[ -d $dir ]]; then
	echo "Error. Directory $dir already exist."
	exit 1
    fi
    
    mkdir $dir
    cp KPOINTS INCAR  POTCAR $dir
    cp $i $dir/POSCAR
done
