#! /bin/sh
set -e
# Spin-flip Raman scattering calculation

# Define structure
cat > s.r << EOF
200 0.15 0.0
60  0.0  0.0
200 0.15 0.0
EOF

# Generate alloy profile
find_heterostructure

# Generate potential profile, note the use of the paramagnetic
# Cd(1-x)Mn(x)Te
efxv -M cdmnte

# Add an 8 Tesla magnetic field to the potential
efmfv -B 8 	# add Zeeman splitting due to 8 T, default electron, spin +

seq 0 10e-10 230e-10 > r_d.r

# Calculate electron-donor energy with cdmnte parameters, let's just use
# the 3D (spherical) trial wave function---it's very quick
qwwad_find_donor_state --symmetry 3D -m 0.096 -e 10.6 --lambdastart 40 --lambdastop 300 > output 

# Save all data for this the `+' spin state

mv e.r e.r+	
mv l.r l.r+
mv v.r v.r+

# Now repeat for the `-' spin state

efmfv -B 8 -s -		# Generate the potential profile

qwwad_find_donor_state --symmetry 3D -m 0.096 -e 10.6 --lambdastart 40 --lambdastop 300 >> output		# donor calculation

mv e.r e.r-		# Save all data
mv l.r l.r-
mv v.r v.r-

# Now produce the energy difference between the states

nawk '{Eplus=$2;getline<"e.r-";printf("%e %e\n",$1,Eplus-$2)}' e.r+ > e_sf.r-raw

# With only a few donor points, calculating the spin-flip spectra produces
# a very spiky Intensity-energy curve.  So take a spline of the spin-flip
# energies in e_sf.r-raw to simulate a continuous donor distribution and
# save in file `e_sf.r'.

sfr -l 0.5 -s 13 -t 0.1 -u 23
mv I.r Il=0.5.r

sfr -s 13 -t 0.1 -u 23
mv I.r I1=1.0.r

sfr -l 2 -s 13 -t 0.1 -u 23
mv I.r Il=2.0.r
