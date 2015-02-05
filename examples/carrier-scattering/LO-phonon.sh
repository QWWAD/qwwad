#! /bin/sh
set -e

# Generate solutions to 100 A infinite well
efiw -N 300 -s 3

# Perform scattering rate calculations
# Define subband populations in file `N.r'
cat > N.r << EOF
1e14
1e14
1e14
EOF

# Calculate the distribution functions
sbp 

# Define required e-LO phonon rates and calculate 
cat > rrp.r << EOF
1 1
2 2
3 3
2 1
3 1
EOF

# Perform scattering rate calculation, need the `-a' option to output form
# factors
srelo -a

# Tidy up, remove all output data files except form factors
rm -f s.r x.r v.r Ee.r wf* Ef.r LO*
