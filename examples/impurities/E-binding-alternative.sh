#! /bin/sh
set -e

# Calculate binding energy for 3D donors in QWs with varying width and depth
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
outfile_E=E-binding-3D.dat
outfile_lambda=E-binding-3D-lambda.dat
rm -f $outfile_E $outfile_lambda

# Loop over alloy concentration
for x in 0.1; do

# Loop over well width
for LW in 10 20 30 40 50 60 80 100 140 180 220 250; do

# Define structure
cat > s.r << EOF
200 $x  0.0
$LW 0.0 0.0
200 $x  0.0
EOF

# Generate alloy profile
find_heterostructure --dz-max 1

# Generate potential profile for Ga(1-x)Al(x)As, can use defaults
efxv --mass 0.067 

# Create r_d.r with single entry at centre of well
echo $LW | awk '{print (200+$1/2)/1e10}' > r_d.r

# Calculate electron energy for same quantum well but without donor
efss

# Start donor binding energy calculation
qwwad_find_donor_state --lambdastart 10 --lambdastop 300 --symmetry 3D > garbage.r

# Energy with donor present
E=`awk '{printf(" %e",$2)}' e.r`

# Store the Bohr radius to file
lambda=`awk '{print $2}' l.r`

# Energy without donor present
E0=`awk '{printf(" %20.17e\n",$2)}' Ee.r`

# Repeat donor calculation using new method
i0

# Energy with donor present
Enew=`awk '{printf(" %e",$2)}' e.r`

# Store data to file, i.e. energy with donor (from e.r), energy
# without donor (from Ee.r) versus well width (lw)
echo $LW $E $Enew $E0 | awk '{print $1, $4 - $2, $4 - $3, 5.3}' >> $outfile_E

# Store the Bohr radius to file
lambda_new=`awk '{print $2}' l.r`

echo $LW $lambda $lambda_new >> $outfile_lambda

done # End loop over well width

printf "\n" >> $outfile_E
printf "\n" >> $outfile_lambda
done # End loop over alloy concentration

# Output the bulk Bohr radius
cat >> $outfile_lambda << EOF
10   103
1000 103
EOF

cat << EOF
Results have been written to $outfile_E and $outfile_lambda.

$outfile_E is in the format:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Binding energy [meV]

$outfile_lambda is in the format:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Bohr radius [Angstrom]

  Each file contains 4 data sets, each set being separated
  by a blank line, representing wells with different barrier
  compositions, and the donor in bulk GaAs for comparison.

  SET 1  - 10% aluminium barriers
  SET 2  - 20% aluminium barriers
  SET 3  - 30% aluminium barriers 
  SET 4  - Bulk material

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
