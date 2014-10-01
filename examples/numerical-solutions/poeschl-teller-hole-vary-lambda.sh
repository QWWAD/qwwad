#!/bin/sh
set -e

# Computes the eigenstates of Poeschl-Teller holes of various depths
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
outfile_E=poeschl-teller-hole-E-lambda.dat
outfile_V=poeschl-teller-hole-V.dat
rm -f $outfile_E $outfile_V

# Define length of Poschl-Teller potential
L=300

# Loop over depth parameter lambda
for LAMBDA in 0.75 1 1.5 2.0 5.0 10.0; do

    # Find potential profile and get analytical solution
    pth --alpha 0.05 --lambda $LAMBDA --length $L

    E1_analytical=`awk '/^1/{print $2}' Ee.r`
    E2_analytical=`awk '/^2/{print $2}' Ee.r`

    # Fill in the blanks in the table if there is no solution
    if [ x$E1_analytical = "x" ]; then
        E1_analytical="--"
    fi

    if [ x$E2_analytical = "x" ]; then
        E2_analytical="--"
    fi

    # Now perform numerical solution
    efss --nst-max 2 --mass 0.067 # calculate 2 lowest energy levels

    E1_numerical=`awk '/^1/{print $2}' Ee.r`
    E2_numerical=`awk '/^2/{print $2}' Ee.r`

    if [ x$E1_numerical = "x" ]; then
        E1_numerical="--"
    fi

    if [ x$E2_numerical = "x" ]; then
        E2_numerical="--"
    fi

    printf "%e\t%s\t%s\t%s\t%s\n" $LAMBDA $E1_analytical $E2_analytical $E1_numerical $E2_numerical >> $outfile_E
    cp v.r v-$LAMBDA.r
done

awk '{print $1*1e10, $2*1000/1.6e-19}' v-0.75.r >> $outfile_V
printf "\n" >> $outfile_V
awk '{print $1*1e10, $2*1000/1.6e-19}' v-1.r >> $outfile_V
printf "\n" >> $outfile_V
awk '{print $1*1e10, $2*1000/1.6e-19}' v-1.5.r >> $outfile_V
printf "\n" >> $outfile_V
awk '{print $1*1e10, $2*1000/1.6e-19}' v-2.0.r >> $outfile_V
printf "\n" >> $outfile_V

cat << EOF
Results have been written to $outfile_V and
$outfile_E, containing the potential profiles and the energies
respectively.

$outfile_V in the format:

  COLUMN 1 - Spatial location [angstrom]
  COLUMN 2 - Potential [meV]

  The file contains 4 data sets, each set being separated
  by a blank line, representing different depth parameters:

  SET 1 - lambda=0.75
  SET 2 - lambda=1.00
  SET 3 - lambda=1.50
  SET 4 - lambda=2.00

$outfile_E in the format:

  COLUMN 1 - Depth parameter
  COLUMN 2 - Energy of state |1> (analytical) [meV]
  COLUMN 3 - Energy of state |2> (analytical) [meV]
  COLUMN 4 - Energy of state |1> (numerical) [meV]
  COLUMN 5 - Energy of state |2> (numerical) [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
