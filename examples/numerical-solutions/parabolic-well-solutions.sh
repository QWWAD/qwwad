#!/bin/sh
set -e

# Computes the eigenstates of a parabolic quantum well
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
outfile=parabolic-well-solutions.dat
outfile_E=parabolic-well-energies.dat
rm -f $outfile $outfile_E

# Calculate alloy profile along z-axis
efpqw --xmax 1

# Convert alloy concentration `x' into effective mass `m'
efxv -m 0.067

# Implement shooting technique
efss --nst-max 3

wfplot --plot-wf --plot-file $outfile
mv Ee.r $outfile_E

cat << EOF
Results have been written to $outfile and $outfile_E, which
contain the wavefunctions and energies, respectively, for a 
parabolic GaAs quantum well.

$outfile is in the format:

  COLUMN 1 - Spatial location [angstrom]
  COLUMN 2 - Confining potential [meV] or wavefunction amplitude [a.u.]

  The file contains 4 data sets, each set being separated
  by a blank line, representing the confining potential and the
  wavefunctions of the first three eigenstates, which are plotted
  at their respective energies:

  SET 1 - Confining potential
  SET 2 - State |1>
  SET 3 - State |2>
  SET 4 - State |3>

$outfile_E is in the format:

  COLUMN 1 - State index
  COLUMN 2 - Energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f files.r
