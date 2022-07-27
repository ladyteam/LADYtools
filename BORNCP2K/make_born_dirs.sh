#!/bin/bash 

for i in born*.xyz; do
    dir=${i%.*}
    if [[ -d $dir ]]; then
        echo "Error. Directory $dir already exist."
        exit 1
    fi

    mkdir $dir

    cp $i $dir/born.xyz
    if [[ $dir == "born_ideal" ]]; then
        cp polar.inp $dir
    else
    cp dipole.inp $dir
    fi
done
