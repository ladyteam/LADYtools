#!/bin/bash

INITDIR=$1
TAIL=$2
CWD=`pwd`

if [[ ! -d $INITDIR ]]; then
    echo "Error init directory with abinit input files does not exist."
    echo "Usage $0 init_dir abinit_tail.in"
    exit 1
fi

if [[ ! -f $TAIL ]]; then
    echo "Error abinit_Tail file does not exist."
    echo "Usage $0 init_dir abinit_tail.in"
    exit 1
fi

if [[ ! -f $INITDIR/displ.files ]]; then
    echo "Error. You should prepare init directory carefully. Make shure displ.files file is exist"
    exit 1
fi

for file in shiftcell-*; do
	fn=${file%.in}
	num=${fn##*-}
	if [[ -d "displ-$num" ]]; then
	    echo "Error directory displ-$num already exist"
	    exit 1
	fi
	echo "Making displ-$num directory"
        mkdir displ-$num
	cp $CWD/$INITDIR/* displ-$num/
	ain=`ls -1 displ-$num/*.in | grep -m1 .in`
	ain=`basename $ain`
	mv "displ-${num}/${ain}" "displ-${num}/${ain}.init"
	rm displ-${num}/*.out
	echo "Generating new input using ${TAIL} pattern"
	mv $file "displ-$num"/displ.in
	echo -e "\n" >> "displ-$num"/displ.in
	cat $TAIL >> "displ-$num"/displ.in
done
