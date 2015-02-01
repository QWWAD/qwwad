#!/bin/sh
# Scattering rate codes (at the moment) are angled towards numerical
# solutions and require a potential barrier height as an upper limit
# for integration, so define an artificial structure
# making sure it contains the same number of points as below, just for
# this purpose

echo 100 1.0 0.0 > s.r
echo 100 0.0 0.0 >> s.r
echo 100 1.0 0.0 >> s.r

# Now convert structure into potential data
find_heterostructure
efxv

# Generate solutions to 100 A infinite well

efiw -N 300 -s 3

# Perform scattering rate calculations
# Define subband populations in file `N.r'
 
echo 1 1 > N.r
echo 2 1 >> N.r
echo 3 1 >> N.r

# Calculate the distribution functions

sbp 

# Define required e-LO phonon rates and calculate 

echo 1 1  > rrp.r
echo 2 2 >> rrp.r
echo 3 3 >> rrp.r
echo 2 1 >> rrp.r
echo 3 1 >> rrp.r

# Perform scattering rate calculation, need the `-a' option to output form
# factors

srelo -a

# Tidy up, remove all output data files except form factors

rm -f s.r x.r v.r Ee.r wf* Ef.r LO* 

