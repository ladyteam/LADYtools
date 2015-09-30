#!/bin/bash


for i in POSCAR-*; do
    num=${i##*-}
    dir="SCAN-"$num
    if [[ -d $dir ]]; then
	echo "Error. Directory $dir already exist."
    else
	mkdir $dir
	cp KPOINTS INCAR  POTCAR $dir
	cp $i $dir/POSCAR
    fi
done
