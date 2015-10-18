#!/bin/sh
set -e

# Calculates the Landau levels derived from the first two subbands in an
# quantum well, as a function of magnetic field
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
outfile=landau-quantisation.dat
rm -f $outfile

# Find first two states in an infinite well
qwwad_ef_square_well --wellwidth 200 --barrierpotential 100 --nst 2
E1=`sed -n 1p < Ee.r | awk '{print $2}'`
E2=`sed -n 2p < Ee.r | awk '{print $2}'`

# Loop over the first 4 Landau levels
for j in `seq 0 3`; do

    # Loop over magnetic field
    for B in `seq 0 8`; do
        echo $B $E1 $E2 $j | awk '{
        B=$1;
        E1=$2;
        E2=$3;
        j=$4;
        hBar=1.055e-34;

        # Find cyclotron frequency and Landau levels
        omega_c=1000*B/(0.067*9.11e-31);
        E1j = E1 + (j + 0.5) * hBar * omega_c;
        E2j = E2 + (j + 0.5) * hBar * omega_c;
        print B, E1j, E2j;
                               }' >> $outfile
    done
    printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Magnetic field [T]
  COLUMN 2 - Energy of Landau level in subband 1 [meV]
  COLUMN 3 - Energy of Landau level in subband 2 [meV]

The file contains 4 data sets, representing the first
4 Landau levels (j) associated with the first two subbands
of a quantum well.  Each set is separated by a blank line:

  SET 1 - Landau levels 0
  SET 2 - Landau levels 1
  SET 3 - Landau levels 2
  SET 4 - Landau levels 3

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm Ee.r wf_e1.r wf_e2.r
