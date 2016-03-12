#!/bin/sh -e

# Computes a simple step diffusion profile with constant diff. coeff.
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.4
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
outfile=constant-step.dat
rm -f $outfile

# Create structure definition file
cat > s.r << EOF
200 0.1 0.0
200 0.0 0.0
EOF

# Generate alloy concentration (diffusant) profile
qwwad_mesh --dzmax 1

# Run diffusion `simulation' for various times
for t in 0 1 10 100 1000 10000
do
    # Use constant diffusion coefficient
    qwwad_diffuse --coeff 10 --time $t

    # Store diffusion profiles
    awk '{print $1*1e10, $2}' X.r >> $outfile
    printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Alloy concentration

  The file contains 6 data sets, each set being separated
  by a blank line, representing different diffusion times:

  SET 1 - t = 0 s
  SET 2 - t = 1 s
  SET 3 - t = 10 s
  SET 4 - t = 100 s
  SET 5 - t = 1000 s
  SET 6 - t = 10000 s

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
