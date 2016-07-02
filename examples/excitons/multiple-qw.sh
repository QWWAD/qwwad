#! /bin/sh
set -e

# Calculate exciton binding-energies in a multiple quantum well system
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
outfile=multiple-qw-EX0-lw.dat
rm -f $outfile

# Loop for different alloys
for xb in 0.1 0.2; do

# Loop for different well widths
for LW in 20 30 40 50 60 70 80; do

echo Calculating for alloy: $xb, width: $LW Angstrom...

cat > s.r << EOF
200 $xb 0.0
$LW 0.0 0.0
$LW $xb 0.0
$LW 0.0 0.0
$LW $xb 0.0
$LW 0.0 0.0
$LW $xb 0.0
$LW 0.0 0.0
$LW $xb 0.0
$LW 0.0 0.0
200 $xb 0.0
EOF

qwwad_mesh --dzmax 1
qwwad_ef_band_edge --mass 0.067 --bandedgepotentialfile v.r
qwwad_ef_generic

cp wf_e1.r wf_${LW}_e1.r # Save wave functions for later analysis

qwwad_ef_band_edge --mass 0.62 --particle h --bandedgepotentialfile v.r
qwwad_ef_generic --particle h
 
cp wf_h1.r wf_${LW}_h1.r	

# Find exciton binding energies
qwwad_ef_exciton --lambdastart 80    \
	         --betastart   0.001 \
		 --betastep    0.01

EX0=`awk '{print $1}' EX0.r`

echo $LW $EX0 >> $outfile
done

printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well and barrier width [Angstrom]
  COLUMN 2 - Binding energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
