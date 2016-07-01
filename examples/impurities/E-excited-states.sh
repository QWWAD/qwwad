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
outfile_E_1s_2s=E-1s-2s.dat
outfile_E_1s_2s_ratio=E-1s-2s-ratio.dat
outfile_E_2s_2px=E-2s-2px.dat
outfile_E_2pz=E-2pz.dat
outfile_lambda_2px_1s=lambda-2px-1s.dat
outfile_lambda_2px_1s_ratio=lambda-2px-1s-ratio.dat

rm -f $outfile_E_1s_2s       \
      $outfile_E_1s_2s_ratio \
      $outfile_E_2s_2px      \
      $outfile_E_2pz         \
      $outfile_lambda_2px_1s \
      $outfile_lambda_2px_1s_ratio

# Loop over well width
for LW in 10 20 30 40 50 60 80 100 140 180 220 250 400 500 600 800 1000; do

echo "Solving for $LW"

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

echo "...bulk"
# Calculate electron energy for same quantum well but without donor
qwwad_ef_generic --nst 1

# Energy without donor present
E0=`awk '{printf(" %20.17e\n",$2)}' Ee.r`

echo "...1s"
# Donor calculation using new method
qwwad_ef_donor_generic --impuritystate 1s

# Energy and Bohr radius for 1s state
E1s=`awk '{printf(" %e",$2)}' Ee.r`
lambda_1s=`awk '{print $2}' l.r`

echo "...2s"
# Now find higher-order states
qwwad_ef_donor_generic --impuritystate 2s
E2s=`awk '{printf(" %e",$2)}' Ee.r`

echo "...2px"
qwwad_ef_donor_generic --impuritystate 2px
E2px=`awk '{printf(" %e",$2)}' Ee.r`
lambda_2px=`awk '{print $2}' l.r`

echo "...2pz"
qwwad_ef_donor_generic --impuritystate 2pz
E2pz=`awk '{printf(" %e",$2)}' Ee.r`
lambda_2pz=`awk '{print $2}' l.r`

# Calculate binding energies
Eb_1s=`echo $E1s $E0 | awk '{print $2 - $1}'`
Eb_2s=`echo $E2s $E0 | awk '{print $2 - $1}'`
Eb_2px=`echo $E2px $E0 | awk '{print $2 - $1}'`
Eb_2pz=`echo $E2pz $E0 | awk '{print $2 - $1}'`

# Calculate useful ratios
Eb_1s_2s_ratio=`echo $Eb_1s $Eb_2s | awk '{print $1/$2}'`
lambda_2px_1s_ratio=`echo $lambda_2px $lambda_1s | awk '{print $1/$2}'`

# Store data to file
echo $LW $Eb_1s $Eb_2s   >> $outfile_E_1s_2s
echo $LW $Eb_2s $Eb_2px  >> $outfile_E_2s_2px
echo $LW $Eb_2pz         >> $outfile_E_2pz

echo $LW $Eb_1s_2s_ratio >> $outfile_E_1s_2s_ratio

echo $LW $lambda_2px_1s_ratio   >> $outfile_lambda_2px_1s_ratio
echo $LW $lambda_2px $lambda_1s >> $outfile_lambda_2px_1s


done # End loop over well width

cat << EOF
Results have been written to:
$outfile_E_1s_2s
$outfile_E_1s_2s_ratio
$outfile_E_2s_2px
$outfile_lambda_2px_1s
$outfile_lambda_2px_1s_ratio
and
$outfile_E_2pz

$outfile_E_1s_2s contains the binding energies for the
1s and 2s states in the format:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Binding energy [meV] (1s)
  COLUMN 3 - Binding energy [meV] (2s)

$outfile_E_1s_2s_ratio contains the ratio of 1s:2s
binding energies:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Ratio of 1s:2s binding energies

$outfile_E_2s_2px contains the binding energies for the
2s and 2px states in the format:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Binding energy [meV] (2s)
  COLUMN 3 - Binding energy [meV] (2px)

$outfile_lambda_2px_1s contains the Bohr radii for the
2px and 1s states in the format:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Bohr radius [Angstrom] (2px)
  COLUMN 3 - Bohr radius [Angstrom] (1s)

$outfile_lambda_2px_1s_ratio contains the ratio of 2px:1s
Bohr radii in the format:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Ratio of 2px:1s Bohr radii

$outfile_E_2pz contains the binding energies for the
2px state in the format:

  COLUMN 1 - Quantum well width [Angstrom]
  COLUMN 2 - Bohr radius [Angstrom] (2pz)

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
