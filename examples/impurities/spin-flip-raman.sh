#! /bin/sh
set -e

# Spin-flip Raman scattering calculation
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.5
#
# (c) Copyright 1996-2016
#     Paul Harrison  <p.harrison@shu.ac.uk>
#     Alex Valavanis <a.valavanis@leeds.ac.uk>
#
# QWWAD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# QWWAD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with QWWAD.  If not, see <http://www.gnu.org/licenses/>.

# Initialise files


# Define structure
cat > s.r << EOF
200 0.15 0.0
60  0.0  0.0
200 0.15 0.0
EOF

# Generate alloy profile
qwwad_mesh

# Generate potential profile, note the use of the paramagnetic
# Cd(1-x)Mn(x)Te
qwwad_ef_band_edge --material cdmnte

# Add an 8 Tesla magnetic field to the potential
qwwad_ef_zeeman --magneticfield 8 --spinup # add Zeeman splitting due to 8 T, default electron, spin up

exit

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

# <Edit accordingly, to describe each output file>
cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - <Whatever>
  COLUMN 2 - <Something>
  <...>

  <Remove this chunk if only 1 data set is in the file>
  The file contains <x> data sets, each set being separated
  by a blank line, representing <whatever>:

  SET 1 - <Description of the 1st set>
  SET 2 - <Description of the 2nd set>
  <...>

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# <Delete all temporary files you created>
rm -f *.r
