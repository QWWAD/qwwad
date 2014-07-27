#!/bin/sh
set -e

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=334.196

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.1002

# Define well and barrier widths here
LW=100
LB=100
N=4

 # Write first barrier and initiate file

 echo 100 0.4 0.0 > s.r		

 # Could only think of an awk script as a replacement `for(i=1;i<N-1;i++)'
 # loop---I'm sure you must be able to do this within `/bin/sh'
 echo $N $LW $LB | awk '{
 Nwells=$1; LW=$2; LB=$3;
 while(Nwells-->1)
     {
         printf("%i 0.0 0.0\n",LW)
         printf("%i 0.4 0.0\n",LB)
     }
 }' >> s.r
		 
 # Write last well and barrier
 echo $LW 0.0 0.0 >> s.r
 echo 100 0.4 0.0 >> s.r

 find_heterostructure --dz-max 0.25 # generate alloy concentration as a function of z
 efxv			# generate potential data

 efss --nst-max 1
 wfplot --plot-wf --plot-file multiple-qw-100-100-wf.dat
