#! /bin/sh
set -e

# Generate a table of exciton binding parameters for finite wells
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.6
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
outfile_EX0=finite-well-EX0-lw.r
outfile_lambda=finite-well-lambda-lw.r
outfile_beta=finite-well-beta-lw.r
outfile_zeta=finite-well-zeta-lw.r
rm -f $outfile_EX0 $outfile_lambda $outfile_beta $outfile_zeta

export QWWAD_DZMAX=1

# Loop over barrier alloy
for xb in 0.1 0.2 0.3; do

# Loop for different well widths
for LW in 20 30 40 50 60 70 80 90 100 120 140 160 180 200; do

echo Calculating for alloy: $xb, width: $LW Angstrom

cat > s.r << EOF
200 $xb 0.0
$LW 0.0 0.0
200 $xb 0.0
EOF

qwwad_mesh

# Solve Schroedinger equation for electrons
qwwad_ef_band_edge --mass 0.067 --bandedgepotentialfile v.r
qwwad_ef_generic

# Solve Schroedinger equation for holes
qwwad_ef_band_edge --mass 0.62  --bandedgepotentialfile v.r --particle h
qwwad_ef_generic --particle h

# Find exciton binding energy
qwwad_ef_exciton --lambdastart 85    \
	         --lambdastop  120   \
	         --betastart   0.75  \
		 --betastop    0.95  \
		 --betastep    0.01

EX0=`awk '{print $1}' EX0.r`
lambda=`awk '{print $2}' EX0.r`
beta=`awk '{print $3}' EX0.r`
zeta=`awk '{print sqrt(1.0 - $3*$3)}' EX0.r`

echo $LW $EX0 >> $outfile_EX0
echo $LW $lambda >> $outfile_lambda
echo $LW $beta >> $outfile_beta
echo $LW $zeta >> $outfile_zeta

done # End loop over well width

printf "\n" >> $outfile_EX0
printf "\n" >> $outfile_lambda
printf "\n" >> $outfile_beta
printf "\n" >> $outfile_zeta
done # End loop over barrier alloy concentration

cat << EOF
Results have been written to $outfile_EX0, $outfile_lambda, $outfile_beta
and $outfile_zeta.

$outfile_EX0 is in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Binding energy [meV]

$outfile_lambda is in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Bohr radius [Angstrom]

$outfile_beta is in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Symmetry parameter, beta

$outfile_zeta is in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Symmetry parameter, zeta

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
