#!/bin/sh
set -e

# Find ground-state of a 100/100 4-quantum-well system.
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2015
#     Alex Valavanis <a.valavanis@leeds.ac.uk>
#     Paul Harrison  <p.harrison@leeds.ac.uk>
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
outfile=multiple-qw-100-100-wf.dat
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
export QWWAD_BARRIERPOTENTIAL=334.196

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
export QWWAD_BARRIERMASS=0.1002

# Define well and barrier widths here
export QWWAD_WELLWIDTH=100
export QWWAD_BARRIERWIDTH=100

# Write first barrier and initiate file
cat > s.r << EOF
$QWWAD_BARRIERWIDTH 0.4 0.0
EOF

# Now append a set of wells and barriers
for i in `seq 1 4`; do

cat >> s.r << EOF
$QWWAD_WELLWIDTH    0.0 0.0
$QWWAD_BARRIERWIDTH 0.4 0.0
EOF

done

# Generate a mesh and find band parameters
qwwad_mesh
qwwad_ef_band_edge --bandedgepotentialfile v.r

# Find the ground state numerically and output a plottable file
qwwad_ef_generic --nstmax 1
qwwad_ef_plot --style wf --plotfile $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Spatial location [angstrom]
  COLUMN 2 - Potential [meV] or wavefunction amplitude [a.u.]

  The file contains 2 data sets, each set being separated
  by a blank line, representing the potential profile and the
  wavefunctions for the first two states:

  SET 1 - Potential profile
  SET 2 - Wavefunction for state |1> offset by state energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
