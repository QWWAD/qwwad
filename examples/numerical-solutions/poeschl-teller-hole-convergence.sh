#!/bin/sh
set -e

# Define output file
outfile=poeschl-teller-hole-convergence.dat

# Initialise files
rm -f $outfile
    
# Define length of Poschl-Teller potential [angstrom]
L=300

# Loop over number of points per Angstrom
for N in 1 2 5 10 20 50 100; do

    # Compute the total number of points needed
    nz=`echo $N $L | awk '{print $1 * $2 + 1}'`

    pth --alpha 0.05 --lambda 5 --length $L --nz $nz

    E1_analytical=`awk '/^1/{print $2}' Ee.r`
    E2_analytical=`awk '/^2/{print $2}' Ee.r`

    # Fill in the blanks in the table if there is no solution
    if [ x$E1_analytical = "x" ]; then
        E1_analytical="--"
    fi

    if [ x$E2_analytical = "x" ]; then
        E2_analytical="--"
    fi

    # Now perform numerical solution
    efss --nst-max 2 --mass 0.067 --solver matrix-constant-mass	# calculate 2 lowest energy levels

    E1_numerical=`awk '/^1/{print $2}' Ee.r`
    E2_numerical=`awk '/^2/{print $2}' Ee.r`

    if [ x$E1_numerical = "x" ]; then
        E1_numerical="--"
    fi

    if [ x$E2_numerical = "x" ]; then
        E2_numerical="--"
    fi

    printf "%e\t%s\t%s\t%s\t%s\n" $N $E1_analytical $E2_analytical $E1_numerical $E2_numerical >> $outfile
done
