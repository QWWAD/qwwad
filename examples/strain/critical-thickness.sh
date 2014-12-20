#! /bin/sh
set -e

outfile=critical-thickness.dat

# SiGe on silicon
critical-thickness
mv hc.r $outfile

# AlGaAs on GaAs
printf "\n" >> $outfile
critical-thickness --C110 118.8 --C111 120.2 --C120 53.8 --C121 57.0 --a0 5.65325 --a1 5.6611
cat hc.r >> $outfile

# InGaAs on GaAs
printf "\n" >> $outfile
critical-thickness --C110 118.8 --C111 83.4 --C120 53.8 --C121 45.4 --a0 5.65325 --a1 6.0583
cat hc.r >> $outfile
