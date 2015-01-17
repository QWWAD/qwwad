#! /bin/sh
set -e

# Initialise files
OUT=E-F.r
OUTsr=cc-auger-F.dat
rm -f $OUT $OUTsr

# Define temperature
T=300

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=167.0985

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.0836

# Define well width here
LW=200

# perform numerical solution
#
# First generate structure definition `s.r' file
cat > s.r << EOF
100 0.2 0.0
$LW 0.0 0.0
100 0.2 0.0
EOF

find_heterostructure  # generate alloy concentration as a function of z
efxv		      # generate potential data
mv v.r vcb.r # Save conduction-band potential for later

# Define required e-e scattering rates
cat > rr.r << EOF
2 2 2 1
2 2 1 1
2 1 1 1
EOF

# Define subband populations
cat > N.r << EOF
1 1
2 1
EOF

# Loop over electric field 
for F in 0 1 2 5 10 20 50
do
 printf "%f " $F >> $OUT
 printf "%f " $F >> $OUTsr

 # Add electric field to potential
 find_poisson_potential --field $F --uncharged --centred --Vbasefile vcb.r

 # Need a small energy difference to split nearly degenerate levels
 efss --nst-max 2

 # Write energy to output file
 awk '{printf("\t%f",$2)}' Ee.r >> $OUT	
 printf "\n" >> $OUT

 # Calculate subband populations
 sbp --Te $T

 # Implement e-e scattering rate calculation
 srcc -T $T

 # Collate results
 awk '{printf(" %e",$5)}' ccABCD.r >> $OUTsr

 # end line
 printf "\n" >> $OUTsr
done
