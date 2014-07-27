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

# Generate structure definition file
cat > s.r << EOF
200 0.2 0.0
$LW 0.0 0.0
$LB 0.2 0.0
$LW 0.0 0.0
200 0.2 0.0
EOF

find_heterostructure --dz-max 0.25
efxv # generate potential data

efss --nst-max 2 # calculate 2 lowest energy levels

E1_numerical=`awk '/^1/{print $2}' Ee.r`
E2_numerical=`awk '/^2/{print $2}' Ee.r`

# Save files
mv Ee.r    Ee-$LB.r
mv v.r     v-$LB.r
mv wf_e1.r wf_e1-$LB.r
mv wf_e2.r wf_e2-$LB.r

printf "%e\t%s\t%s\n" $LB $E1_numerical $E2_numerical >> $outfile
done

wfplot --plot-wf --energy-input Ee-40.r --wf-input-ext "-40.r" --potential-input "v-40.r"
