#!/bin/bash

# this code allows you to get the good edges from TTMC quickly and then
# run FASTME on them as a dissimililarity matrix

TTMC=/fs/cbcb-lab/ekmolloy/trubel/repos/TREE-TMC/TREE-TMC
FASTME=/fs/cbcb-lab/ekmolloy/trubel/software/fastme-2.1.5/binaries/fastme-2.1.5-linux64

INPUT=$1

OUTPUT=$2

R=$RANDOM

$TTMC -i $INPUT -m -o $OUTPUT-$R.tmp

$FASTME -i $OUTPUT-$R.tmp -o $OUTPUT

rm $OUTPUT-$R.tmp
