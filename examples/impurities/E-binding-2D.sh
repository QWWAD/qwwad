#!/bin/sh
set -e

# Find binding-energy for a 2D donor in a QW of variable width
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

# Initialise files
outfile_E=E-binding-2D.dat
outfile_ED0=E-binding-2D-ED0.dat
rm -f $outfile_E $outfile_ED0

# Loop over well width
for LW in 20 30 40 60 80 100 140 180 220 300 400; do

# Define structure
cat > s.r << EOF
200 0.1 0.0
$LW 0.0 0.0
200 0.1 0.0
EOF

# Generate alloy profile
find_heterostructure 

# Generate potential profile for Cd(1-x)Mn(x)Te
efxv --material cdmnte --mass 0.096

# Create r_d.r with single entry at centre of well
echo $LW | awk '{print (200+$1/2)/1e10}' > r_d.r

# Start donor binding energy calculation
qwwad_find_donor_state --epsilon 10.6 --mass 0.096 --lambdastart 25 --lambdastop 150 > garbage.r

# Calculate electron energy for same quantum well but without donor
efss

# Energy with donor present
E=`awk '{printf(" %e",$2)}' e.r`

# Energy without donor present
E0=`awk '{printf(" %20.17e\n",$2)}' Ee.r`

# Store data to file, i.e. energy with donor (from e.r), energy
# without donor (from Ee.r) versus well width (lw)
echo $LW $E $E0 >> $outfile_E

awk '{print $1, $3 - $2, 11.7}' < $outfile_E > $outfile_ED0
done

cat << EOF
Results have been written to $outfile_E and ${outfile_ED0}.

$outfile_E is in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Total carrier energy [meV]
  COLUMN 3 - Carrier energy without donor [meV]

$outfile_ED0 is in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Binding energy [meV]
  COLUMN 3 - Bulk binding energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF


# Clean up workspace
rm -f *.r
