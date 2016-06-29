#! /bin/sh
set -e

# Calculate binding energy for 3D donors in QWs with varying width and depth
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.5
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
outfile_E=E-binding-alternative.dat
rm -f $outfile_E

# Loop over alloy concentration
for x in 0.1 0.2 0.3; do

# Loop over well width
for LW in 10 20 30 40 50 60 80 100 140 180 220 250; do

# Define structure
cat > s.r << EOF
200 $x  0.0
$LW 0.0 0.0
200 $x  0.0
EOF

# Generate alloy profile
qwwad_mesh --dzmax 1

# Generate potential profile for Ga(1-x)Al(x)As, can use defaults
qwwad_ef_band_edge --mass 0.067 --bandedgepotentialfile v.r 

# Calculate electron energy for same quantum well but without donor
qwwad_ef_generic --nst 1

# Energy without donor present
E0=`awk '{printf(" %20.17e\n",$2)}' Ee.r`

# Start donor binding energy calculation
qwwad_ef_donor_specific --lambdastart 10 --lambdastop 300 --symmetry 3D

# Energy with donor present
E=`awk '{printf(" %e",$2)}' Ee.r`

# Store the Bohr radius to file
lambda=`awk '{print $2}' l.r`

# Repeat donor calculation using new method
qwwad_ef_generic
qwwad_ef_donor_generic

# Energy with donor present
Enew=`awk '{printf(" %e",$2)}' Ee.r`

# Store data to file, i.e. energy with donor (from e.r), energy
# without donor (from Ee.r) versus well width (lw)
echo $LW $E $Enew $E0 | awk '{print $1, $4 - $2, $4 - $3}' >> $outfile_E

done # End loop over well width

printf "\n" >> $outfile_E
done # End loop over alloy concentration

cat << EOF
Results have been written to $outfile_E in the format:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Binding energy [meV] (semi-analytical)
  COLUMN 3 - Binding energy [meV] (direct integration)

  The file contains 3 data sets, each set being separated
  by a blank line, representing wells with different barrier
  compositions.

  SET 1  - 10% aluminium barriers
  SET 2  - 20% aluminium barriers
  SET 3  - 30% aluminium barriers 

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
