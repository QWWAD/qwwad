#!/bin/sh
set -e

# Initialise files
outfile=uncertainty-principle-finite-well.dat
rm -f $outfile

# Loop for different well widths
for LW in 20 30 40 50 60 70 80 90 100 120 140 160 180 200; do

 # Calculate ground state energy and wave function as a function 
 # of well width for GaAs
 efsqw --well-width $LW

 # Search for line in standard output from hup and write to file 
 data=`hup | awk '/Delta_z.Delta_p/{printf("%8.3f\n",$2)}'`

 printf "%d\t%s\n" $LW $data >> $outfile
done
