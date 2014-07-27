#!/bin/sh
set -e

# Define output file
outfile=variable-mass-0.75-alloy-dz.dat

# Initialise files
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=626.6175

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.12925

# Define a set well width
LW=20

# Loop over spatial resolution [points-per-angstrom]
for N in 2 4 6 8 10 12; do
    # Calculate lowest 2 levels with analytical form
    efsqw --well-width $LW --barrier-mass $MB --nst 2 --potential $V

    E1_analytical=`awk '/^1/{print $2}' Ee.r`

    # Fill in the blanks in the table if there is no solution
    if [ x$E1_analytical = "x" ]; then
        E1_analytical="--"
    fi

    # Now perform numerical solution

    # First generate structure definition `s.r' file
    echo 200 0.75 0.0 > s.r
    echo $LW 0.0 0.0  >> s.r
    echo 200 0.75 0.0 >> s.r

    # Work out how many points we need for the desired sampling period
    nz=`echo $LW $N | awk '{print ($1 + 400) * $2 + 1}'`

    find_heterostructure --nz-1per $nz # generate alloy concentration as a function of z
    efxv			  # generate potential data

    efss --nst-max 1 --solver matrix-variable-mass # calculate lowest energy level

    E1_numerical=`awk '/^1/{print $2}' Ee.r`

    if [ x$E1_numerical = "x" ]; then
        E1_numerical="--"
    fi

    printf "%e\t%s\t%s\n" $N $E1_analytical $E1_numerical >> $outfile
done
