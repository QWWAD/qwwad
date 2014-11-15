#! /bin/sh
set -e

# Variable symmetry donor calculation versus donor position across a quantum well
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2014
#     Paul Harrison  <p.harrison@shu.ac.uk>
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
outfile_E=E-donor-variable.dat
outfile_zl=E-donor-variable-zeta-lambda.dat
rm -f $outfile_E $outfile_zl

# Define structure
cat > s.r << EOF
200 0.1 0.0
100 0.0 0.0
200 0.1 0.0
EOF

# Generate alloy profile
find_heterostructure --dz-max 1		

# Generate potential profile
efxv		

seq 0 10e-10 250e-10 > r_d.r

d02D --symmetry variable --lambdastart 50 --zetastart 0.65 > garbage.r

mv e.r $outfile_E
mv l.r $outfile_zl
printf "\n" >> $outfile_zl
cat zeta.r >> $outfile_zl

cat << EOF
Results have been written to $outfile_E and $outfile_zl

$outfile_E is in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Total carrier energy [meV]

$outfile_zl is in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Hydrogenic wavefunction parameter

  The file contains 2 data sets, each set being separated
  by a blank line, representing the paramaters needed to
  describe the hydrogenic waveform:

  SET 1 - Bohr radius [Angstrom]
  SET 2 - Symmetry parameter

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
