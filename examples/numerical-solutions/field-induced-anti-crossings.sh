#!/bin/sh
set -e

# Calculate the field-induced anti-crossing between first 2 states in double
# quantum well
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2014
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
outfile=field-induced-anti-crossings-E-F.dat
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=167.0985

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.0836

# Define well widths here
LW1=60
LW2=50

# perform numerical solution

# First generate structure definition `s.r' file
echo 200 0.2 0.0 > s.r
echo $LW1 0.0 0.0 >> s.r
echo 60 0.2 0.0 >> s.r
echo $LW2 0.0 0.0 >> s.r
echo 200 0.2 0.0 >> s.r

find_heterostructure 	# generate alloy concentration as a function of z
efxv			# generate potential data

# Loop over electric field
for F in 0 1 2 3 4 5 6 7 8 9 10 11 12 15 20 25 30 40; do
 # Add electric field to potential
 find_poisson_potential --centred --field $F --uncharged --Vbasefile v.r --potential-file v_t.r

 efss --nst-max 2 --v-file v_t.r # calculate ground and first excited states

 # Write energy to output file
 E1=`awk '/^1/{print $2}' Ee.r`
 E2=`awk '/^2/{print $2}' Ee.r`

 echo $F $E1 $E2 >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Electric field [kV/cm]
  COLUMN 2 - Energy of state |1> [meV]
  COLUMN 3 - Energy of state |2> [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
