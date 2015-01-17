#! /bin/sh
set -e
# Calculation of the mean e-e scattering rate over two subband populations
# as a function of well width

# Initialise files
OUT=cc-avg-DE.dat
rm -f $OUT

# Define subband populations in file `N.r'
 
echo 1 1 > N.r
echo 2 1 >> N.r

# Define the required e-e rate
echo 2 2 1 1 > rr.r

# Loop over infinite well width

for LW in 100 200 300 400 500 600
do
{
 # Generate infinitely deep well solutions
 efiw -L $LW -N 300 -s 2

 # Save lowest two energies to file
 E1=`sed -n 1p < Ee.r | awk '{print $2}'`
 E2=`sed -n 2p < Ee.r | awk '{print $2}'`
 DE=`echo $E1 $E2 | awk '{print $2-$1}'`

 printf "%f " $DE >> $OUT

 # Loop over different temperatures
 for T in 4 300
 do
 {
 # Calculate the distribution functions
 sbp --Te $T
 
 # Calculate carrier-carrier (e-e) scattering rate
 srcc -T $T

 # Sort and store in output file
 awk '/2 2 1 1/{printf("%e ",$5)}' ccABCD.r >> $OUT
 }
 done

# End line in output file
printf "\n" >> $OUT

}
done
