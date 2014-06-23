#!/bin/sh
set -e

# Define output file
outfile=double-quantum-well-E-vs-LB.dat

# Initialise files
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=167.0985

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.0836

# Define well width here
LW=60

# Loop over barrier width, execute out of order to retain 40 Angstrom data
for LB in 10 20 30 40 50 60 70 80 90 100 120 140 160 180 200; do

 # perform numerical solution
 #
 # First generate structure definition `s.r' file
 echo 200 0.2 0.0 > s.r
 echo $LW 0.0 0.0 >> s.r
 echo $LB 0.2 0.0 >> s.r
 echo $LW 0.0 0.0 >> s.r
 echo 200 0.2 0.0 >> s.r

 # Work out how many points we need for a 1-angstrom sampling period
 nz=`echo $LW $LB | awk '{print $1*2 + $2 + 401}'`

 find_heterostructure --nz $nz
 efxv # generate potential data

 efss --nst-max 2 # calculate 2 lowest energy levels

 E1_numerical=`awk '/^1/{print $2}' Ee.r`
 E2_numerical=`awk '/^2/{print $2}' Ee.r`

 if [ x$E1_numerical = "x" ]; then
     E1_numerical="--"
 fi

 if [ x$E2_numerical = "x" ]; then
     E2_numerical="--"
 fi

 # Save files
 mv Ee.r Ee-$LB.r
 mv wf_e1.r wf_e1-$LB.r
 mv wf_e2.r wf_e2-$LB.r

 printf "%e\t%s\t%s\n" $LB $E1_numerical $E2_numerical >> $outfile
done

wfplot --energy-input Ee-40.r --wf-input-ext "-40.r"
