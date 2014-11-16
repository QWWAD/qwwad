#! /bin/sh
set -e

# Compare total energy for 2D, 3D and variable symmetry donors at various locations in a QW
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

outfile=E-donor-comparison.dat
rm -f $outfile

# Define structure
cat > s.r << EOF
200 0.1 0.0
60  0.0 0.0
200 0.1 0.0
EOF

# Generate alloy profile
find_heterostructure

# Generate potential profile
efxv --material cdmnte --mass 0.096

# Define donor positions
seq 150e-10 10e-10 230e-10 > r_d.r

# Perform 2D donor calculation and save results to output file
qwwad_find_donor_state --mass 0.096 --epsilon 10.6 --lambdastart 25 --lambdastop 300 --symmetry 2D > garbage.r
mv e.r e-2D.r

# Perform 3D donor calculation and save results to file
qwwad_find_donor_state --mass 0.096 --epsilon 10.6 --lambdastart 25 --lambdastop 300 --symmetry 3D > garbage.r
mv e.r e-3D.r

# Perform variable symmetry calculation
qwwad_find_donor_state --mass 0.096 --epsilon 10.6 --lambdastart 25 --zetastart 0.5 --symmetry variable > garbage.r

paste e-2D.r e-3D.r e.r | awk '{print $1, $2, $4, $6}' > $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Total carrier energy (2D symmetry) [meV]
  COLUMN 3 - Total carrier energy (3D symmetry) [meV]
  COLUMN 4 - Total carrier energy (variable symmetry) [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
