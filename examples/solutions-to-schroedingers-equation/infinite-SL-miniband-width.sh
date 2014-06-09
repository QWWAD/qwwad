#!/bin/sh
set -e

# Calculates the lowest miniband-width in a Kronig-Penney superlattice
# and compares it with a single quantum well
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

# Initialise files
outfile=infinite-SL-miniband-width.dat
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.4
V=334.1965

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x
MB=0.1002

# Define number of states
S=1

# Loop for well and barrier widths

for LW in `seq 20 5 200`; do
    printf "%d\t" "$LW" >> $outfile # write well width to file
 
    # Loop for carrier momentum at centre and edge of Brillouin zone only
    for K in 0.0 1.0; do
        # Calculate energies for different wave vectors

        efkpsl -a $LW -b $LW -m 0.067 -n $MB --potential $V -k $K -s $S
        awk '{printf("%9.3f",$2)}' Ee.r >> $outfile	# send data to file
    done	# done loop over k

    # Compare with energy of single quantum well
    efsqw --well-width $LW --well-mass 0.067 --barrier-mass $MB --potential $V --nst $S
    awk '{printf("%9.3f\n",$2)}' Ee.r >> $outfile	# send data to file
done	# loop over LW

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well width [angstrom]
  COLUMN 2 - Energy of lowest state in miniband [meV]
  COLUMN 3 - Energy of highest state in miniband [meV]
  COLUMN 4 - Energy of ground state in single quantum well [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

rm Ee.r
