#! /bin/sh
set -e

# Find eigenstates of an infinite rectangular quantum wire of varying widths
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2014
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

# Infinitely deep square cross-sectional quantum wire calculation
# Generate 4 lowest energy solutions and collate energy eigenvalues as a
# function of well width.  Store the cross-sectional charge densities for
# one well width only.
# Define and initialise output file

# Initialise files
outfile=infinite-rectangular-wire.dat
rm -f $outfile

# Loop over wire side
for L in 10 20 30 40 50 60 70 80 90 100 120 140 160 200 240 300
do
{
 efiwire -y $L -z $L -s 2
 
 E=`awk '{printf(" %e",$3)}' Ee.r`

 echo $L $E >> $outfile
}
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Width of wire [angstrom]
  COLUMN 2 - Energy of state 1,1
  COLUMN 3 - Energy of state 1,2
  COLUMN 4 - Energy of state 2,1
  COLUMN 5 - Energy of state 2,2

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# rm -f *.r
