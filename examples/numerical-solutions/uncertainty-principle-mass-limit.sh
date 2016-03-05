#!/bin/sh
set -e

# Find uncertainty (position x momentum) for wells of varying width and
# barrier mass
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

# Initialise files
outfile=uncertainty-principle-mass-limit.dat
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.4
export QWWAD_BARRIERPOTENTIAL=334.196
export QWWAD_WELLMASS=0.067

# Loop over well width
for LW in 20 50 100; do

    # Loop for different barrier heights
    for MB in 0.01 0.02 0.03 0.04 0.05 0.06 0.067 0.07 0.08 0.09 0.10 0.12 0.14 0.16 0.18 0.20 0.3 0.4 0.50 0.6 0.7 0.8 0.9 1.0; do
 
        # Calculate ground state energy for different well and barrier masses
        qwwad_ef_square_well --wellwidth $LW --barriermass $MB

        # Search for line in standard output from qwwad_uncertainty and write to file 
        data=`qwwad_uncertainty | awk '/Delta_z.Delta_p/{printf("%8.3f\n",$2)}'`

        echo $MB $data >> $outfile
    done

    printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well width [angstrom]
  COLUMN 2 - Uncertainty [x hbar/(2 pi)]
  
  The file contains 3 data sets, each set being separated
  by a blank line, representing different well widths:

  SET 1 - 20-angstrom well
  SET 2 - 50-angstrom well
  SET 3 - 100-angstrom well

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
