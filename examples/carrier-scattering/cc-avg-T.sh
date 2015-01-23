#! /bin/sh
set -e

# Calculation of the mean e-e scattering rate over two subband populations
# as a function of temperature
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
OUT=cc-avg-T.dat
rm -f $OUT

# Define the well width
LW=300

# Generate infinitely deep well solutions
efiw -L $LW -N 300 -s 2

# Define the required e-e rate
echo 2 2 1 1 > rr.r

# Write carrier density to file
cat > N.r << EOF
10e14
10e14
EOF

# Loop over different temperatures
for T in 4 10 15 20 30 40 77 120 160 220 300 500
do
    # Calculate the distribution functions
    sbp --Te $T

    # Calculate carrier-carrier (e-e) scattering rate
    srcc -T $T

    # Sort and store in output file
    W=`awk '/2 2 1 1/{printf(" %e\n",$5)}' ccABCD.r`
    echo $T $W >> $OUT
done

cat << EOF
Results have been written to $OUT in the format:

  COLUMN 1 - Temperature of electron distribution [K]
  COLUMN 2 - Average |2,2> -> |1,1> scattering rate [s^{-1}]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
#rm -f *.r
