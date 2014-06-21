#!/bin/sh
set -e

# Define output file
outfile=variable-mass-0.75-alloy-well-width.dat

# Initialise files
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=626.6175

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.12925

# Loop over well width
for LW in 20 40 60 80 100 120 160 200; do
    # Calculate lowest 2 levels with analytical form
    efsqw --well-width $LW --barrier-mass $MB --nst 2 --potential $V

    E1_analytical=`awk '/^1/{print $2}' Ee.r`
    E2_analytical=`awk '/^2/{print $2}' Ee.r`

    # Fill in the blanks in the table if there is no solution
    if [ x$E1_analytical = "x" ]; then
        E1_analytical="--"
    fi

    if [ x$E2_analytical = "x" ]; then
        E2_analytical="--"
    fi

    # Save lowest two energies to file

    # Now perform numerical solution
    #
    # First generate structure definition `s.r' file
    echo 200 0.75 0.0 > s.r
    echo $LW 0.0 0.0 >> s.r
    echo 200 0.75 0.0 >> s.r

    # Work out how many points we need for a 1-angstrom sampling period
    nz=`echo $LW | awk '{print $1 + 401}'`

    find_heterostructure --nz $nz # generate alloy concentration as a function of z
    efxv			  # generate potential data

    efss --nst-max 2 --solver matrix-variable-mass # calculate 2 lowest energy levels

    E1_numerical=`awk '/^1/{print $2}' Ee.r`
    E2_numerical=`awk '/^2/{print $2}' Ee.r`

    if [ x$E1_numerical = "x" ]; then
        E1_numerical="--"
    fi

    if [ x$E2_numerical = "x" ]; then
        E2_numerical="--"
    fi

    printf "%e\t%s\t%s\t%s\t%s\n" $LW $E1_analytical $E2_analytical $E1_numerical $E2_numerical >> $outfile
done
