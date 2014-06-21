#!/bin/sh
set -e

# Define output file
outfile_E=poeschl-teller-hole-E-lambda.dat
outfile_V=poeschl-teller-hole-V.dat

# Initialise files
rm -f $outfile_E $outfile_V

# Define length of Poschl-Teller potential
L=300

# Loop over depth parameter lambda
for LAMBDA in 0.75 1 1.5 2.0 5.0 10.0; do

    # Find potential profile and get analytical solution
    pth --alpha 0.05 --lambda $LAMBDA --length $L

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

    printf "%e\t%s\t%s\t%s\t%s\n" $LAMBDA $E1_analytical $E2_analytical $E1_numerical $E2_numerical >> $outfile_E
    cp v.r v-$LAMBDA.r
done

awk '{print $1*1e10, $2*1000/1.6e-19}' v-0.75.r >> $outfile_V
printf "\n" >> $outfile_V
awk '{print $1*1e10, $2*1000/1.6e-19}' v-1.r >> $outfile_V
printf "\n" >> $outfile_V
awk '{print $1*1e10, $2*1000/1.6e-19}' v-1.5.r >> $outfile_V
printf "\n" >> $outfile_V
awk '{print $1*1e10, $2*1000/1.6e-19}' v-2.0.r >> $outfile_V
printf "\n" >> $outfile_V
