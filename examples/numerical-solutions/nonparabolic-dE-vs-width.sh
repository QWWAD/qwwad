#!/bin/sh
set -e

outfile=nonparabolic-dE-vs-width.dat
rm -f $outfile

# Loop over barrier concentrations
for X in 0.2 0.4 0.6 0.8 1.0; do

 # Loop over well width
 for LW in 20 40 60 80 100 120 140 160 180 200; do

     # First generate structure definition `s.r' file
     echo 200 $X 0.0 > s.r
     echo $LW 0.0 0.0 >> s.r
     echo 200 $X 0.0 >> s.r

     find_heterostructure # generate alloy concentration as a function of z
     efxv 		  # generate potential data, and bandgap

     # Calculate ground state energy with band non-parabolicity
     efss --solver shooting-nonparabolic --nst-max 1

     # Write energy to output file
     E1_np=`awk '{printf("\t%f",$2)}' Ee.r`

     # Calculate ground state energy without band non-parabolicity
     efss --solver shooting-variable-mass

     # Write energy to output file
     E1_parab=`awk '{printf("\t%f",$2)}' Ee.r`

     # Now calculate difference between `with' and `without' band non-parabolicity
     dE=`echo $E1_parab $E1_np | awk '{print $2-$1}'`

     echo $LW $dE >> $outfile
 done # LW

 printf "\n" >> $outfile
done # X
