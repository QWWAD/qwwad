#! /bin/sh
set -e

# Calculates the first 3 energy level in a 100-angstrom-wide quantum
# well and outputs the wavefunctions
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2015
#     Paul Harrison <p.harrison@shu.ac.uk>
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

# Initialise output file
outfile=infinite-well-wavefunctions.dat
rm -f $outfile

# Solve Schroedinger equation
qwwad_ef_infinite_well --width 100 --nz 1000 --nst 3

# Shift wavefunctions and scale positions to angstrom
awk '{printf("%3.5f    %.12e\n", $1*1e10, $2)}' wf_e1.r >> $outfile
printf "\n" >> $outfile
awk '{printf("%3.5f    %.12e\n", $1*1e10, $2 + 30000)}' wf_e2.r >> $outfile
printf "\n" >> $outfile
awk '{printf("%3.5f    %.12e\n", $1*1e10, $2 + 60000)}' wf_e3.r >> $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Spatial position [angstrom]
  COLUMN 2 - Wavefunction amplitude [a.u.]

The file contains three data sets, separated by a blank line,
which represent each of the first three states in an infinite
well.

Note that an arbitrary offset has been added onto the amplitude
of the second and third wavefunctions such that they can be
plotted neatly on the same axis.

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

rm wf_e?.r Ee.r
