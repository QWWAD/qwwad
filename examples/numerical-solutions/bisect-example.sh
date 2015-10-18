#!/bin/sh
set -e

# Calculates the matching equation for a 200-angstrom QW
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2015
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
outfile=bisect-example.dat
rm -f $outfile

qwwad_ef_square_well --wellwidth 200 --barrierpotential 100 --outputequations --nst 1

paste lhs.r rhs_1.r | awk '{if ($1 < 1.5) print $1, $2 - $4}' > bisect-example.dat

# <Edit accordingly, to describe each output file>
cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well wave-vector (normalised)
  COLUMN 2 - Matching equation

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f lhs.r rhs_?.r Ee.r wf_?.r
