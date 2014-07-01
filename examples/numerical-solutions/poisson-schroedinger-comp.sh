#!/bin/sh
set -e

# Define output file
outfile=poisson-schroedinger-comp-E-I.dat

# Initialise files
rm -f $outfile wf_e1*.r

# First generate structure definition `s.r' file
echo 200 0.2 0.0 0.0  > s.r
echo 100 0.0 0.0 4e18 >> s.r
echo 200 0.2 0.0 0.0  >> s.r
 
# Loop over number of points along z-axis
for N in 10 5 1; do	# reverse loop and hence retain diverging wave functions
    nz=`echo $N | awk '{print 500*$1 + 1}'`

    find_heterostructure --nz $nz	# generate alloy concentration as a function of z
    efxv			# generate potential data
    cp v.r vcb.r # Save conduction-band energy
  
    for I in 0 1 2 3 4 5 6 7; do

        # Calculate ground state Schroedinger solution
        efss --nst-max 1

        densityinput # Generate an estimate of the population density
        chargedensity # Compute charge density profile
  
        # save wave function is separate file
        cp wf_e1.r wf_e1-I=$I.r

        # Write energy to output file
        E1=`awk '{printf("\t%20.17e\n",$2)}' Ee.r`

        echo $I $E1 >> $outfile

        # Implement self consistent Poisson calculation
        find_poisson_potential
        paste vcb.r v_p.r | awk '{print $1, $2+$4}' > v.r
    done # X

   printf "\n" >> $outfile 
done # N
