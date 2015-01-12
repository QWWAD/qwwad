#! /bin/sh
set -ev

# Calculation of the mean e-e scattering rate over two subband populations
# as a function of temperature

# Initialise files
OUT=cc-intra-T.dat
rm -f $OUT

# Define the well width
LW=300

# Generate infinitely deep well solutions
efiw -L $LW -N 300 -s 2

# Define the required e-e rate
echo 2 2 2 2 > rr.r

# Write carrier density to file
cat > N.r << EOF
1 10
2 10
EOF

# Loop over different temperatures
for T in 10 20 40 77 120 160 220 300
do
 # Write temperature to output file
 printf "%f " $T >> $OUT

 # Calculate the distribution functions
 sbp --Te $T
 
 # Calculate carrier-carrier (e-e) scattering rate
 srcc -T $T

 # Sort and store in output file
 awk '/2 2 2 2/{printf(" %e\n",$5)}' ccABCD.r >> $OUT
done
