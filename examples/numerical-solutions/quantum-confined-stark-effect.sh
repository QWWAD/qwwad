#!/bin/sh
set -e

# Find ground state of a quantum well in a variable electric field
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.2
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
outfile=quantum-confined-stark-effect-E-F.dat
rm -f $outfile

# First generate structure definition `s.r' file
cat > s.r << EOF
200 0.2 0.0
 60 0.0 0.0
200 0.2 0.0
EOF

export QWWAD_NSTMAX=1 # Only consider ground state

qwwad_mesh         # generate alloy concentration as a function of z
qwwad_ef_band_edge # generate potential data

# Loop over electric field 
for F in 0 1 2 3 4 5 6 7 8 9 10 12 15 20 25 30 40 50 60 70; do
 # Add electric field to potential
 qwwad_poisson --field $F --centred

 qwwad_ef_generic

 # Write energy to output file
 E1=`awk '/^1/{print $2}' Ee.r`
 echo $F $E1 >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Electric field [kV/cm]
  COLUMN 2 - Energy of ground state [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
