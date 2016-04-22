#! /bin/sh
set -e

# Map total energy of electrons in QW with variable-symmetry donors
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.5
#
# (c) Copyright 1996-2016
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
outfile=E-zeta.dat
rm -f $outfile

# Define structure
cat > s.r << EOF
200 0.1 0.0
100 0.0 0.0
200 0.1 0.0
EOF

# Create alloy profile
qwwad_mesh --dzmax 1

# Create potential profile
qwwad_ef_band_edge --bandedgepotentialfile v.r

# Initiate variable symmetry donor calculation, note specifying the final
# lambda and zeta values forces this domain to be covered
qwwad_ef_donor_specific --donorposition 250      \
                        --lambdastart   89       \
                        --lambdastop    91       \
                        --zetastart     0.6      \
                        --zetastep      0.02     \
                        --zetastop      0.8      \
                        --symmetry      variable \
                        --searchmethod  linear

# Now collate data from the search log into file suitable for plotting
# Do each zeta value at a time
for lambda in 89 90 91; do
	awk '{if($1==lambda*1e-10) print $2, $3/1.6e-19*1000}' lambda=$lambda < searchlog.r >> $outfile
	printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile.

$outfile is in the format:

  COLUMN 1 - Symmetry parameter (zeta)
  COLUMN 2 - Total carrier energy [meV]

  The file contains 3 data sets, each set being separated
  by a blank line, representing different Bohr radii.

  SET 1  - lambda = 89 Angstrom
  SET 2  - lambda = 90 Angstrom
  SET 3  - lambda = 91 Angstrom

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
#rm -f *.r
