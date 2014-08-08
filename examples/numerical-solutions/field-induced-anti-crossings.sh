#!/bin/sh
set -e

# Define output file

OUT=E-F.r

# Initialise files

rm -f $OUT

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2

V=167.0985

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2

MB=0.0836

# Define well widths here

LW1=60
LW2=50

# perform numerical solution
#
# First generate structure definition `s.r' file
echo 200 0.2 0.0 > s.r
echo $LW1 0.0 0.0 >> s.r
echo 60 0.2 0.0 >> s.r
echo $LW2 0.0 0.0 >> s.r
echo 200 0.2 0.0 >> s.r

find_heterostructure 	# generate alloy concentration as a function of z
efxv			# generate potential data

# Loop over electric field
for F in 0 1 2 3 4 5 6 7 8 9 10 11 12 15 20 25 30 40
do
{
 # Add electric field to potential
 find_poisson_potential --centred --field $F --uncharged
 paste v.r v_p.r | awk '{print $1, $2+$4}' > v_t.r

 efss --nst-max 2 --v-file v_t.r # calculate ground and first excited states

 # Write energy to output file
 E1=`awk '/^1/{print $2}' Ee.r`
 E2=`awk '/^2/{print $2}' Ee.r`

 echo $F $E1 $E2 >> $OUT
}
done
