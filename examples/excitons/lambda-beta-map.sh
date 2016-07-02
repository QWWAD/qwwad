#! /bin/sh
set -e

# Generate a table of exciton binding-energies as a function of variational parameters
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.6
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
outfile=lambda-beta-map.dat

# Infinite well electron and hole ground states
export QWWAD_WELLWIDTH=100
qwwad_ef_infinite_well --mass 0.096 --particle e
qwwad_ef_infinite_well --mass 0.6   --particle h

# Find donor binding energies
qwwad_ef_exciton --dcpermittivity 10.6   \
	         --electronmass    0.096 \
		 --holemass        0.6   \
		 --lambdastart    50     \
		 --lambdastop     70     \
		 --betastart       0.70  \
		 --betastep        0.01  \
		 --betastop        0.91  \
		 --searchlogfile $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Bohr radius [Angstrom]
  COLUMN 2 - Variational parameter, beta
  COLUMN 3 - Binding energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
