#! /bin/sh
set -e

# Calculates the wavefunctions for the first two states in a quantum well
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
outfile_odd=shooting-method-wavefunction-odd.dat
outfile_even=shooting-method-wavefunction-even.dat

# Describe square well structure (thickness, alloy, doping)
cat > s.r << EOF
150 0.2 0.0
100 0.0 0.0
150 0.2 0.0
EOF

# Create alloy concentration file
qwwad_mesh

# Find band-edge parameters
efxv

# Find lowest two states
efss --nst-max 2

# Rescale to angstrom units and centre at zero
awk '{print $1*1e10 - 200, $2}' wf_e1.r > $outfile_even
awk '{print $1*1e10 - 200, $2}' wf_e2.r > $outfile_odd

cat << EOF
Results have been written to $outfile_even and $outfile_odd
in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Wavefunction amplitude [m^{-0.5}]

$outfile_even contains the (even parity) ground state and
$outfile_odd contains the (odd parity) first excited state

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# <Delete all temporary files you created>
rm -f *.r
