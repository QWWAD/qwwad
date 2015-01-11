#! /bin/sh
set -ve

# Calculation of the mean e-e scattering rate over two subband populations
# as a function of temperature
# Define output file

OUT=cc-T.r

# Initialise files

rm -f $OUT

# Define the well width
LW=300

# Generate infinitely deep well solutions
efiw -L $LW -N 300 -s 2

# Define the required e-e rate
echo 2 2 1 1 > rr.r

# Write carrier density to file
cat > N.r << EOF
1 10
2 10
EOF

# Loop over different temperatures
for T in 4 10 15 20 30 40 77 120 160 220 300 500
do
{
 # Calculate the distribution functions
 sbp --Te $T
 
 # Calculate carrier-carrier (e-e) scattering rate
 srcc -T $T

 # Sort and store in output file
 W=`awk '/2 2 1 1/{printf(" %e\n",$5)}' ccABCD.r`
 echo $T $W >> $OUT
}
done

