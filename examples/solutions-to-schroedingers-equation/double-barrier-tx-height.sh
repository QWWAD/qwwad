#!/bin/sh
set -e

# Calculates the transmission coefficient through a double barrier
# with a range of heights
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2014
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

outfile=double-barrier-tx-height.dat
rm -f $outfile

# Keep well width constant
L2=50

for X in 0.1 0.2 0.3 0.4; do
    # Use V=0.67*1247*x
    V=`echo $X | awk '{print 0.67*1247*$1}'`

    # Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
    # Use MB=0.067+0.083*x
    MB=`echo $X | awk '{print 0.067+0.083*$1}'`

    tdb -b $L2 --potential $V -n $MB
    cat T.r >> $outfile
    printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Energy of incident electron [meV]
  COLUMN 2 - Transmission coefficient

The file contains 4 data sets representing barriers with different alloy
concentrations.  Each set is separated by a blank line.

  SET 1 - 10% alloy
  SET 2 - 20% alloy
  SET 3 - 30% alloy
  SET 4 - 40% alloy

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm T.r
