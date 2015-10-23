#!/bin/sh
set -e

# Calculates the shooting method wavefunction at either side of correct eigenvalue
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
outfile=shooting-convergence.dat
rm -f $outfile

# Describe square well structure (thickness, alloy, doping)
cat > s.r << EOF
150 0.2 0
100 0.0 0
150 0.2 0
EOF

# Always use a fixed effective mass and shooting method solver
export QWWAD_MASS=0.067
export QWWAD_SOLVER=shooting

# Create alloy concentration file
qwwad_mesh

# Find band-edge parameters
qwwad_ef_band_edge --bandedgepotentialfile v.r

# Find correct wavefunction
qwwad_ef_generic --nstmax 1

# Rescale to angstrom and send to output file
awk '{print $1*1e10, $2}' wf_e1.r > $outfile

# Find trial wavefunction just below state
qwwad_ef_generic --tryenergy 29.3

printf "\n" >> $outfile
awk '{print $1*1e10, $2}' wf_eE.r >> $outfile

# Find trial wavefunction just above state
qwwad_ef_generic --tryenergy 29.5

printf "\n" >> $outfile
awk '{print $1*1e10, $2}' wf_eE.r >> $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Spatial location [angstrom]
  COLUMN 2 - Wavefunction amplitude [m^{-0.5}]

  The file contains 3 data sets, each set being separated
  by a blank line, representing the trial wavefunction at
  a different energy:

  SET 1 - Energy = E1 (i.e., the correct value)
  SET 2 - Energy 1 meV below eigenvalue
  SET 3 - Energy 1 meV above eigenvalue

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
