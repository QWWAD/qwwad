#!/bin/sh -e

# Computes a diffusion profile for QW with constant diff. coeff.
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.4
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
outfile=constant-Dt.dat
rm -f $outfile

# Create structure definition file
cat > s.r << EOF
200 0.1 0.0
200 0.0 0.0
200 0.1 0.0
EOF

# Generate alloy concentration (diffusant) profile
qwwad_mesh --dzmax 1

# Run diffusion `simulation' for various times
for t in 100 1000; do

    # Find diffuse profile
    qwwad_diffuse --coeff 10 --time $t
 
    # Store diffusion profiles
    awk '{print $1*1e10, $2}' X.r >> $outfile
    printf "\n" >> $outfile
done

# Now repeat calculations for different D but equivalent Dt
for D in 100 1000; do

    # Need to reduce time-step to ensure stability
    qwwad_diffuse --coeff $D --time 10 --dt 0.0001

    # Store diffusion profiles
    awk '{print $1*1e10, $2}' X.r >> $outfile
    printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Alloy concentration

  The file contains 4 data sets, each set being separated
  by a blank line, representing different combinations of
  diffusion time and diffusion coefficient:

  SET 1 - t = 100 s,  D = 10 Angstrom^2/s
  SET 2 - t = 1000 s, D = 10 Angstrom^2/s
  SET 3 - t = 10 s,   D = 100 Angstrom^2/s
  SET 4 - t = 10 s,   D = 1000 Angstrom^2/s

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
