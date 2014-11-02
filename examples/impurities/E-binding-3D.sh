#!/bin/sh
# 3D trial wave fucntion donor binding energy calculationas a function of
# well width for different barrier heights
# Define the barrier alloy concentration

# Initialise files
outfile_E=E-binding-3D.dat
outfile_ED0=E-binding-3D-ED0.dat
rm -f $outfile_E $outfile_ED0

# Loop over well width
for LW in 10 20 30 40 50 60 80 100 140 180 220 260 300 400 500 600 800 1000; do

# Define structure
cat > s.r << EOF
200 0.1 0.0
$LW 0.0 0.0
200 0.1 0.0
EOF

# Generate alloy profile
find_heterostructure

# Generate potential profile for Ga(1-x)Al(x)As, can use defaults
efxv --mass 0.067 

# Create r_d.r with single entry at centre of well
echo $LW | awk '{print (200+$1/2)/1e10}' > r_d.r

# Start donor binding energy calculation
d02D --lambdastart 10 --lambdastop 300 --symmetry 3D > garbage.r

# Calculate electron energy for same quantum well but without donor
efss

# Energy with donor present
E=`awk '{printf(" %e",$2)}' e.r`

# Energy without donor present
E0=`awk '{printf(" %20.17e\n",$2)}' Ee.r`

# Store data to file, i.e. energy with donor (from e.r), energy
# without donor (from Ee.r) versus well width (lw)
echo $LW $E $E0 >> $outfile_E

awk '{print $1, $3 - $2, 11.7}' < $outfile_E > $outfile_ED0
done

