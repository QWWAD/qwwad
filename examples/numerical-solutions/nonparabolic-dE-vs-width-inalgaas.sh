#!/bin/sh
set -e

# Define alloy components X and Y in In(1-x-y)Al(x)Ga(y)As for both well
# and barrier

# Quantum wells, In(0.53)Ga(0.47)As
Yw=0.47
Xw=0.0

# Barriers, In(0.52)Al(0.48)As
Yb=0.0
Xb=0.48

# Define output file
outfile=nonparabolic-dE-vs-width-inalgaas.dat
rm -f $outfile

# Loop over well width
for LW in 20 30 40 50 60 80 100 120 140 160 180 200; do
 # First generate structure definition `s.r' file
 echo 200 $Xb $Yb 0 > s.r
 echo $LW $Xw $Yw 0 >> s.r
 echo 200 $Xb $Yb 0 >> s.r
 
 # Remember to switch material system
 find_heterostructure	  # generate alloy concentration as a function of z
 efxv --material inalgaas # generate potential data, and bandgap
 
 # Calculate ground state energy with band non-parabolicity
 efss --solver shooting-nonparabolic

 # Write energy to output file
 E_np=`awk '{printf("\t%f",$2)}' Ee.r`

 # Calculate ground state energy without band non-parabolicity
 efss --solver shooting-variable-mass

 # Write energy to output file
 E_parab=`awk '{printf("\t%f",$2)}' Ee.r`

 # Now calculate difference between `with' and `without' band non-parabolicity
 dE=`echo $E_parab - $E_np | awk '{print $1 - $2}'`

 echo $LW $dE >> $outfile
done # LW
