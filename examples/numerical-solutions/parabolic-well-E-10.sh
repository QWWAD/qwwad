#!/bin/sh
set -e

# Computes the first ten states in a parabolic well
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.3
#
# (c) Copyright 1996-2016
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
outfile=parabolic-well-E-10.dat
rm -f $outfile

# Calculate alloy profile along z-axis, choose artificially high potential
qwwad_ef_parabolic_well --xmax 10 --nz 1000

# Convert alloy concentration `x' into effective mass `m', specify constant
# mass
qwwad_ef_band_edge --mass 0.067 --bandedgepotentialfile v.r

# Implement shooting technique
qwwad_ef_generic --nstmax 10

mv Ee.r $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - State index
  COLUMN 2 - Energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
