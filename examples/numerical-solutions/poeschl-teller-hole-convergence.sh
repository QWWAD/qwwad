#!/bin/sh
set -e

# Computes the convergence of Poeschl-Teller hole Schroedinger solution
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.3
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

# Define output file
outfile=poeschl-teller-hole-convergence.dat

# Initialise files
rm -f $outfile
    
# Define a fixed geometry Poschl-Teller potential [angstrom]
export QWWAD_LENGTH=300
export QWWAD_WIDTHPARAMETER=0.05
export QWWAD_DEPTHPARAMETER=5

# Loop over number of points per Angstrom
for N in 1 2 5 10 20 50 100; do

    # Compute the total number of points needed
    nz=`echo $N $QWWAD_LENGTH | awk '{print $1 * $2 + 1}'`

    qwwad_ef_poeschl_teller --nz $nz

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
    qwwad_ef_generic --nstmax 2 --mass 0.067 # calculate 2 lowest energy levels

    E1_numerical=`awk '/^1/{print $2}' Ee.r`
    E2_numerical=`awk '/^2/{print $2}' Ee.r`

    if [ x$E1_numerical = "x" ]; then
        E1_numerical="--"
    fi

    if [ x$E2_numerical = "x" ]; then
        E2_numerical="--"
    fi

    printf "%e\t%s\t%s\t%s\t%s\n" $N $E1_analytical $E2_analytical $E1_numerical $E2_numerical >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Number of points per angstrom
  COLUMN 2 - Analytical solution for state |1> [meV]
  COLUMN 3 - Analytical solution for state |2> [meV]
  COLUMN 4 - Numerical solution for state |1> [meV]
  COLUMN 5 - Numerical solution for state |2> [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

rm -f *.r
