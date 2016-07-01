#! /bin/sh
set -e

# Calculate binding energy for 1s and 2s states in QWs of varying width
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
outfile=E-1s-2s.dat
outfile_ratio=E-1s-2s-ratio.dat

rm -f $outfile $outfile_ratio

# Loop over well width
for LW in 10 20 30 40 50 60 80 100 140 180 220 250 400 500 600 800 1000; do

# Define structure
cat > s.r << EOF
200 0.1 0.0
$LW 0.0 0.0
200 0.1 0.0
EOF

# Generate alloy profile
qwwad_mesh --dzmax 1

# Generate potential profile for Ga(1-x)Al(x)As, can use defaults
qwwad_ef_band_edge --mass 0.067 --bandedgepotentialfile v.r 

# Calculate electron energy for same quantum well but without donor
qwwad_ef_generic --nst 1

# Energy without donor present
E0=`awk '{printf(" %20.17e\n",$2)}' Ee.r`

# Donor calculation using new method
qwwad_ef_donor_generic --impuritystate 1s

# Energy with donor present
E1s=`awk '{printf(" %e",$2)}' Ee.r`

qwwad_ef_donor_generic --impuritystate 2s
E2s=`awk '{printf(" %e",$2)}' Ee.r`

# Calculate binding energies
Eb_1s=`echo $E1s $E0 | awk '{print $2 - $1}'`
Eb_2s=`echo $E2s $E0 | awk '{print $2 - $1}'`

Eb_ratio=`echo $Eb_1s $Eb_2s | awk '{print $1/$2}'`

# Store data to file, i.e. energy with donor (from e.r), energy
# without donor (from Ee.r) versus well width (lw)
echo $LW $Eb_1s $Eb_2s >> $outfile
echo $LW $Eb_ratio >> $outfile_ratio

done # End loop over well width

cat << EOF
Results have been written to $outfile and $outfile_ratio.

$outfile contains the binding energies for the 1s and 2s states
in the format:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Binding energy [meV] (1s)
  COLUMN 3 - Binding energy [meV] (2s)

$outfile_ratio contains the ratio of 1s:2s binding energies:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Ratio of 1s:2s binding energies

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
