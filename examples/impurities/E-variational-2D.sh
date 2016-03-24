#! /bin/sh
set -e

# Mapping of lambda-energy space for a 2D donor wave function in a quantum well
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
outfile=E-variational-2D.dat
rm -f $outfile

# Generate structure file
cat > s.r << EOF
200 0.1 0.0
60  0.0 0.0
200 0.1 0.0
EOF

# Convert structure file into alloy profile
qwwad_mesh

# Generate potential profile
qwwad_ef_band_edge --material cdmnte --bandedgepotentialfile v.r

# Place donor at centre of well
echo "230e-10" > r_d.r

# Start donor binding energy calculation,
# force output of energy E versus lambda data for lambda=40 to 90A in 2A steps
# in order to illustrate the variational principle
qwwad_ef_donor_specific --lambdastart 40 --lambdastep 2 --lambdastop 90 --dcpermittivity 10.6 --mass 0.096 --searchmethod linear > output.tmp

# Filter `output' file to give energy versus lambda data
awk '{printf("%e %e\n",$4,$6)}' output.tmp > $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Bohr radius [Angstrom]
  COLUMN 2 - Electron energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
#rm -f *.r
