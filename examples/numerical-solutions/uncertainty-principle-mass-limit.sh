#!/bin/sh
set -e

# Initialise files
outfile=uncertainty-principle-mass-limit.dat
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.4
V=334.196

# Loop over well width
for LW in 20 50 100; do

    # Loop for different barrier heights
    for MB in 0.01 0.02 0.03 0.04 0.05 0.06 0.067 0.07 0.08 0.09 0.10 0.20 0.50 1.0; do
 
        # Calculate ground state energy for different well and barrier masses
        efsqw --well-width $LW --well-mass 0.067 --barrier-mass $MB --potential $V

        # Search for line in standard output from hup and write to file 
        data=`hup | awk '/Delta_z.Delta_p/{printf("%8.3f\n",$2)}'`

        echo $MB $data >> $outfile
    done

    printf "\n" >> $outfile
done

