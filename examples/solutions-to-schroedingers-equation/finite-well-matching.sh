#! /bin/sh
set -e

# Generates a plot of the matching equations for a finite quantum well
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
outfile=finite-well-matching.dat
rm -f $outfile

efsqw --well-width 200 --potential 100 --output-equations

cat lhs.r   >> $outfile
printf "\n" >> $outfile
cat rhs_1.r >> $outfile
printf "\n" >> $outfile
cat rhs_2.r >> $outfile
printf "\n" >> $outfile
cat rhs_3.r >> $outfile
printf "\n" >> $outfile
cat rhs_4.r >> $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well wave vector (normalised)
  COLUMN 2 - Barrier decay constant (normalised)

  The file contains 5 data sets, each set being separated
  by a blank line, representing part of the matching
  equations for the system:

  SET 1 - The left-hand-side function for this specific system
  SET 2 - The right-hand-side function for the 1st state of a generic finite well
  SET 3 - The right-hand-side function for the 2nd state of a generic finite well
  SET 4 - The right-hand-side function for the 3rd state of a generic finite well
  SET 5 - The right-hand-side function for the 4th state of a generic finite well

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f lhs.r rhs_?.r wf*.r Ee.r
