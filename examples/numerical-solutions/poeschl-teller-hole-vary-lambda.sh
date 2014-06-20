#!/bin/sh
set -e

# Define output file
OUT=E-LAMBDA.r

# Initialise files

rm -f Ee-LAMBDA*
rm -f v-LAMBDA*
rm -f $OUT

# Define length of Poschl-Teller potential
L=300

# Loop over depth parameter lambda
for LAMBDA in 0.01 0.02 0.05 0.1 0.2 0.5 0.75 1 1.5 2.0 5.0 10.0
do
{
 pth -a 0.05 -l $LAMBDA -L $L
 echo -n $LAMBDA >> $OUT

 # Now perform numerical solution
 efss --nst-max 2 --mass 0.067 --solver matrix-constant-mass	# calculate 2 lowest energy levels
 
 # Save lowest two energies to file
 awk '/^1/||/^2/{printf("\t%f",$2)}' Ee.r >> $OUT	
 cp Ee.r Ee-LAMBDA=$LAMBDA.r
 
 echo -n -e "\n" >> $OUT
 cp v.r v-LAMBDA=$LAMBDA.r
}
done
