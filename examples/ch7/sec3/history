#!/bin/sh
# Infinitely deep square cross-sectional quantum wire calculation
# Generate 4 lowest energy solutions and collate energy eigenvalues as a
# function of well width.  Store the cross-sectional charge densities for
# one well width only.
# Define and initialise output file

OUTPUT=Ee-L.r

rm -f $OUTPUT

# Loop over wire side

for L in 10 20 30 40 50 60 70 80 90 100 120 140 160 200 240 300
do
{
 efiwire -y $L -z $L -s 2
 
 echo -n $L >> $OUTPUT
 nawk '{printf(" %e",$2)}' Ee.r >> $OUTPUT
 echo -n -e "\n" >> $OUTPUT
}
done


