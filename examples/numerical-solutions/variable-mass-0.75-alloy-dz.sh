#!/bin/sh
set -e

# Computes the states in a square well with spatially-variable effective mass
# as a function of spatial resolution
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
outfile=variable-mass-0.75-alloy-dz.dat
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=626.6175

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.12925

# Define a set well width
LW=20

# Loop over spatial resolution [points-per-angstrom]
for N in 2 4 6 8 10 12; do
    # Calculate lowest 2 levels with analytical form
    efsqw --well-width $LW --barrier-mass $MB --nst 2 --potential $V

    E1_analytical=`awk '/^1/{print $2}' Ee.r`

    # Fill in the blanks in the table if there is no solution
    if [ x$E1_analytical = "x" ]; then
        E1_analytical="--"
    fi

    # Now perform numerical solution

    # First generate structure definition `s.r' file
    echo 200 0.75 0.0 > s.r
    echo $LW 0.0 0.0  >> s.r
    echo 200 0.75 0.0 >> s.r

    # Work out how many points we need for the desired sampling period
    nz=`echo $LW $N | awk '{print ($1 + 400) * $2 + 1}'`

    find_heterostructure --nz-1per $nz # generate alloy concentration as a function of z
    efxv			  # generate potential data

    efss --nst-max 1 --solver matrix-variable-mass # calculate lowest energy level

    E1_numerical=`awk '/^1/{print $2}' Ee.r`

    if [ x$E1_numerical = "x" ]; then
        E1_numerical="--"
    fi

    printf "%e\t%s\t%s\n" $N $E1_analytical $E1_numerical >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well width [angstrom]
  COLUMN 2 - Analytical solution for |1> [meV]
  COLUMN 3 - Numerical solution for |1> [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
