#!/bin/sh
set -e

# Define output file
outfile=multiple-qw-40-40-E-N.dat

# Initialise files
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=167.0985

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.0836

# Define well and barrier widths here
LW=40
LB=40

# To compare with infinite superlattice
efkpsl --well-width $LW --barrier-width $LB --barrier-mass $MB --potential $V
E1_analytical=`awk '/^1/{print $2}' Ee.r`

for N in 1 2 3 4 5 6 7 8 9 10; do
 
    # Write first barrier and initiate file
    echo 200 0.2 0.0 > s.r	

 # Could only think of an awk script as a replacement `for(i=1;i<N-1;i++)'
 # loop---I'm sure you must be able to do this within `/bin/sh'
 echo $N $LW $LB | awk '{
 Nwells=$1; LW=$2; LB=$3;
 while(Nwells-- > 1)
     {
         printf("%i 0.0 0.0\n",LW)	
         printf("%i 0.2 0.0\n",LB)	
     }
 }' >> s.r
		 
 # Write last well and barrier
 echo $LW 0.0 0.0 >> s.r
 echo 200 0.2 0.0 >> s.r

 find_heterostructure --dz-max 0.25 # generate alloy concentration as a function of z
 efxv			# generate potential data

 efss --nst-max 1

 # Write energy to output file
 E1_numerical=`awk '/^1/{print $2}' Ee.r`
 printf "%d\t%e\t%e\n" $N $E1_analytical $E1_numerical >> $outfile
done

wfplot --plot-wf --plot-file multiple-qw-40-40-wf.dat
