#! /bin/sh
set -ev

# Calculation of the mean e-e scattering rate over two subband populations
# as a function of well width

# Initialise files
OUT=cc-intra-N.r
rm -f $OUT

# Define the well width
LW=300

# Generate infinitely deep well solutions
efiw -L $LW -N 301 -s 2

# Define the required e-e rate
echo 2 2 2 2 > rr.r

# Loop over subband populations
for N in 1 2 5 10 20 50 100 
do
 # Write carrier density to file
cat > N.r << EOF
1 $N
2 $N
EOF

 # and save in output file
 printf "%f " $N >> $OUT

 # Loop over different temperatures
 for T in 4 77 300
 do
     # Calculate the distribution functions
     sbp --Te $T
 
     # Calculate carrier-carrier (e-e) scattering rate
     srcc -T $T

     # Sort and store in output file
     awk '/2 2 2 2/{printf(" %e",$5)}' ccABCD.r >> $OUT
 done

# End line in output file
printf "\n" >> $OUT
done
