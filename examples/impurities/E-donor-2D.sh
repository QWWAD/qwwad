#!/bin/sh
set -e

# Plot binding-energy for a 2D donor wave function at various positions in a QW
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
outfile_E=E-donor-2D-total.dat
outfile_ED=E-donor-2D-binding.dat
outfile_lambda=E-donor-2D-lambda.dat
rm -f $outfile_E $outfile_ED

# Define structure
cat > s.r << EOF
200 0.1 0.0
60  0.0 0.0
200 0.1 0.0
EOF

find_heterostructure    # generates x.r file
efxv --material cdmnte --mass 0.096	# converts x.r into v.r appropriate to Cd(1-x)Mn(x)Te

# Just calculate energy of same structure without a donor, for deduction
# of donor binding energy
efss # calculate electron energy without donor, m=constant
E1=`awk '/^1/{print $2}' Ee.r`

# Define donor positions
cat > r_d.r << EOF
0.0e-10
2.0e-9
4.0e-9
6.0e-9
8.0e-9
1.0e-8
1.2e-8
1.4e-8 
1.6e-8 
1.8e-8 
2.0e-8 
2.2e-8 
2.3e-8
EOF

# Initiate donor binding energy calculation
# Direct output to garbage file... we don't need it for this example
qwwad_find_donor_state --mass 0.096 --epsilon 10.6 --lambdastart 25 --lambdastop 300 > garbage.r

# Output the total energy and binding energy to files
awk '{print $1, $2, E1}' E1=$E1 < e.r > $outfile_E
awk '{print $1, E1-$2}' E1=$E1 < $outfile_E > $outfile_ED

mv l.r $outfile_lambda

cat << EOF
Results have been written to $outfile_E, ${outfile_ED} and ${outfile_lambda}.

$outfile_E is in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Total carrier energy [meV]
  COLUMN 3 - Carrier energy without donor [meV]

$outfile_ED is in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Binding energy [meV]

$outfile_lambda is in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Bohr radius [Angstrom]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
