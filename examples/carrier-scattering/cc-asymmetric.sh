#! /bin/sh
set -e

# Calculation of the intrasubband Auger-type scattering rates
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2014
#     Paul Harrison <p.harrison@shu.ac.uk>
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
OUTsr=cc-asymmetric-F.dat
rm -f $OUTsr

# Define temperature
T=300

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=167.0985

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.0836

# Define well width here
LW=200

# perform numerical solution
#
# First generate structure definition `s.r' file
cat > s.r << EOF
100 0.2 0.0
$LW 0.0 0.0
100 0.2 0.0
EOF

find_heterostructure	# generate alloy concentration as a function of z
efxv			# generate potential data
mv v.r vcb.r # Save conduction band plot for later

# Define required e-e scattering rates
cat > rr.r << EOF
1 1 1 1
2 1 2 1
2 1 1 2
EOF

# Define subband populations [m^{-2}]
cat > N.r << EOF
1e14
1e14
EOF

# Loop over electric field 
for F in 0 1 2 5 10 20 50
do
 printf "%f " $F >> $OUTsr

 # Add electric field to potential
 find_poisson_potential --field $F --uncharged --Vbasefile vcb.r

 # Need a small energy difference to split nearly degenerate levels
 efss --nst-max 2

 # Calculate subband populations
 sbp --Te $T

 # Implement e-e scattering rate calculation
 srcc -T $T

 # Collate results
 awk '{printf(" %e",$5)}' ccABCD.r >> $OUTsr

 # end line
 printf "\n" >> $OUTsr
done

cat << EOF
Results have been written to $OUTsr in the format:

  COLUMN 1 - Applied field [kV/cm]
  COLUMN 2 - Average |1,1> -> |1,1> scattering rate [s^{-1}]
  COLUMN 3 - Average |2,1> -> |2,1> scattering rate [s^{-1}]
  COLUMN 4 - Average |2,1> -> |1,2> scattering rate [s^{-1}]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
#rm -f *.r
