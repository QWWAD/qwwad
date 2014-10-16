#!/bin/sh -e

# Computes a diffusion profile for QW with constant diff. coeff. and varying time
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
outfile_D=constant-D.dat
outfile_E=constant-D-E-t.dat
rm -f $outfile_D $outfile_E

# Create structure definition file
cat > s.r << EOF
200 0.1 0.0
200 0.0 0.0
200 0.1 0.0
EOF

find_heterostructure --dz-max 1

# Run diffusion `simulation' for various times
for t in 0 10 20 50 100 200 500 1000; do

    # Generate diffusion profile
    gde --coeff 10 --time $t

    # Store diffusion profiles
    awk '{print $1*1e10, $2}' X.r >> $outfile_D
    printf "\n" >> $outfile_D

    # Now solve for the electron energy
    efxv --alloyfile X.r    # Find potential profile
    efss                    # calculate ground state electron energy

    # Save electron energy in file
    E1=`awk '/^1/{printf("\t%e\n",$2)}' Ee.r`
    echo $t $E1 >> $outfile_E
done

cat << EOF
Results have been written to $outfile_D and $outfile_E.

$outfile_D contains the diffusion profiles in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Alloy concentration

  The file contains 8 data sets, each set being separated
  by a blank line, representing different diffusion times:

  SET 1 - t = 0 s
  SET 2 - t = 10 s
  SET 3 - t = 20 s
  SET 4 - t = 50 s
  SET 5 - t = 100 s
  SET 6 - t = 200 s
  SET 7 - t = 500 s
  SET 8 - t = 1000 s

$outfile_E contains the variation in ground-state energy with
respect to diffusion time, in the format:

  COLUMN 1 - Diffusion time [s]
  COLUMN 2 - Energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
