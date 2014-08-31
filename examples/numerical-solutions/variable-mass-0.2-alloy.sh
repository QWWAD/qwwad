#!/bin/sh
set -e

# Computes the states in a square well with spatially-variable effective mass
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
outfile=variable-mass-0.2-alloy.dat
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
V=167.0985

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.0836

# Loop over well width
for LW in 20 40 60 80 100 120 160 200; do
    # Calculate lowest 2 levels with analytical form
    efsqw --well-width $LW --barrier-mass $MB --nst 2 --potential $V

    # Save lowest two energies to file
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
    # First generate structure definition `s.r' file
    echo 200 0.2 0.0 > s.r
    echo $LW 0.0 0.0 >> s.r
    echo 200 0.2 0.0 >> s.r

    find_heterostructure --dz-max 1 # generate alloy concentration as a function of z
    efxv			    # generate potential data

    efss --nst-max 2 --solver matrix-variable-mass # calculate 2 lowest energy levels

    E1_numerical=`awk '/^1/{print $2}' Ee.r`
    E2_numerical=`awk '/^2/{print $2}' Ee.r`

    if [ x$E1_numerical = "x" ]; then
        E1_numerical="--"
    fi

    if [ x$E2_numerical = "x" ]; then
        E2_numerical="--"
    fi

    printf "%e\t%s\t%s\t%s\t%s\n" $LW $E1_analytical $E2_analytical $E1_numerical $E2_numerical >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well width [angstrom]
  COLUMN 2 - Analytical solution for |1> [meV]
  COLUMN 3 - Analytical solution for |2> [meV]
  COLUMN 4 - Numerical solution for |1> [meV]
  COLUMN 5 - Numerical solution for |2> [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
