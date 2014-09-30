#! /bin/sh
set -e

# Calculates the first 3 energy level in a finite GaAs quantum well
# and plots the wavefunctions and potential.
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
outfile=finite-well-wf.dat

# Find wave functions
efsqw --well-width 200 --potential 100 --nst 3 --output-potential

# Generate plot file
wfplot --plot-wf --plot-file $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Position [angstrom]
  COLUMN 2 - Plot data [wavefunction in a.u. or potential in meV]

  The file contains 4 data sets, each being separated by a blank line

  SET 1 - The confining potential [meV]
  SET 2 - Wavefunction for state 1 [arb. units]
  SET 3 - Wavefunction for state 2 [arb. units]
  SET 4 - Wavefunction for state 3 [arb. units]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f v.r Ee.* wf*
