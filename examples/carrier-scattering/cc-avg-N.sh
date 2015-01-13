#! /bin/sh
set -e

# Calculation of the mean e-e scattering rate over two subband populations
# as a function of carrier density

# Initialise files
OUT=cc-avg-N.dat
rm -f $OUT

# Define the well width
LW=300

# Generate infinitely deep well solutions
efiw -L $LW -N 300 -s 2

# Define the required e-e rate
echo 2 2 1 1 > rr.r

T=300

# Loop over subband populations
for N in 1 2 5 10 20 50 100 
do
{
 # Write carrier density to file
cat > N.r << EOF
1 $N
2 $N
EOF

 # and save in output file
 printf "%f " $N >> $OUT

 # Calculate the distribution functions
 sbp --Te $T
 
 # Calculate carrier-carrier (e-e) scattering rate
 srcc -T $T --Ecutoff 800

 # Sort and store in output file
 awk '/2 2 1 1/{printf(" %e",$5)}' ccABCD.r >> $OUT

# End line in output file
printf "\n" >> $OUT
}
done
