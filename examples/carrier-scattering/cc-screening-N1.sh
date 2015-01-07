#!/bin/sh
# Calculation of the e-e scattering rate over two subband populations
# with and without screening
# Define temperature

T=77

# Initialise files

rm -f $OUT

# Scattering rate codes (at the moment) are angled towards numerical
# solutions and require a potential barrier height as an upper limit
# for integration, so define an artificial structure
# making sure it contains the same number of points as below

echo 100 1.0 0.0 > s.r
echo 100 0.0 0.0 >> s.r
echo 100 1.0 0.0 >> s.r

# Now convert structure into potential data
find_heterostructure
efxv

# Define width of infinite well 

LW=400

# Generate infinitely deep well solutions

efiw -L $LW -N 300 -s 2

# Define subband populations in file `N.r'
 
echo 1 1 > N.r
echo 2 1 >> N.r

# Calculate the distribution functions
sbp --Te $T
 
# Define required e-e rate

echo 2 2 1 1 > rr.r

# Calculate e-e rate WITH screening

srcc -T $T

# Save data
 
mv cc2211.r cc2211S.r
 
# Now calculate WITHOUT screening
 
srcc -T $T -S
 
