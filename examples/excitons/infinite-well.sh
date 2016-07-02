#! /bin/sh
set -e

# Find binding energy for an exciton in a QW of variable width
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
outfile_EX0=infinite-well-EX0-lw.r
outfile_lambda=infinite-well-lambda-lw.r
rm -f $outfile_EX0 $outfile_lambda

# Loop for different well widths
export QWWAD_WELLWIDTH

for QWWAD_WELLWIDTH in 0.1 0.5 1 2 5 10 20 30 40 50 60 70 80 90 100 120 140 160 180 200 300 400 500 600 800 1000;
do
	export QWWAD_NZ=1000
	# Find ground states for electron and hole
	qwwad_ef_infinite_well --mass 0.096 --particle e # infinite well electron ground state
	qwwad_ef_infinite_well --mass 0.6   --particle h # infinite well hole ground state

	qwwad_ef_exciton --dcpermittivity 10.6 --electronmass 0.096 --holemass 0.6 --lambdastart 30 	   # start variational calculation

	# Save data to file
	EX0=`awk '{print $1}' EX0.r`
	lambda=`awk '{print $2}' EX0.r`

       	# Write data to file
	echo $QWWAD_WELLWIDTH $EX0 >> $outfile_EX0
	echo $QWWAD_WELLWIDTH $lambda >> $outfile_lambda
done

cat << EOF
Results have been written to $outfile_EX0 and ${outfile_lambda}.

$outfile_EX0 is in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Binding energy [meV]

$outfile_lambda is in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Bohr radius [Angstrom]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
