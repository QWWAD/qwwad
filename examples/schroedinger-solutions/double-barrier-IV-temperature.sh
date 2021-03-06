#!/bin/sh
set -e

# Calculates the current-voltage relation through a double barrier
# at a range of temperatures
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

# Initialise files
outfile=double-barrier-IV-temperature.dat
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
X=0.2

# Use V=0.67*1247*x
export QWWAD_BARRIERPOTENTIAL=`echo $X | awk '{print 0.67*1247*$1}'`

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x
export QWWAD_BARRIERMASS=`echo $X | awk '{print 0.067+0.083*$1}'`

# Keep well width constant
export QWWAD_LEFTBARRIERWIDTH=100
export QWWAD_WELLWIDTH=50
export QWWAD_RIGHTBARRIERWIDTH=100

# Loop over carrier temperature
for T in 2 77 300; do
    qwwad_tx_double_barrier_iv --Te $T
    cat IV.r >> $outfile
    printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Voltage drop across structure [V]
  COLUMN 2 - Current [a.u.]

The file contains 3 data sets representing characteristics at different
temperatures:

  SET 1 - 2 K
  SET 2 - 77 K
  SET 3 - 300 K

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm T.r IV.r
