#!/bin/sh
set -e

# Compute eigenstates of a double well as a function of width
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2015
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
outfile=double-quantum-well-E-vs-LB.dat
outfile_wf=double-quantum-well-wf.dat
rm -f $outfile

# Define well width here
LW=60

# Loop over barrier width, execute out of order to retain 40 Angstrom data
for LB in `seq 10 10 200`; do

# Generate structure definition file
cat > s.r << EOF
200 0.2 0.0
$LW 0.0 0.0
$LB 0.2 0.0
$LW 0.0 0.0
200 0.2 0.0
EOF

qwwad_mesh --dzmax 0.25
qwwad_ef_band_edge --bandedgepotentialfile v.r # generate potential data

qwwad_ef_generic --nstmax 2 # calculate 2 lowest energy levels

E1_numerical=`awk '/^1/{print $2}' Ee.r`
E2_numerical=`awk '/^2/{print $2}' Ee.r`

# Save files
mv Ee.r    Ee-$LB.r
mv v.r     v-$LB.r
mv wf_e1.r wf_e1-$LB.r
mv wf_e2.r wf_e2-$LB.r

printf "%e\t%s\t%s\n" $LB $E1_numerical $E2_numerical >> $outfile
done

# Generate a 'pretty' plot of bandstructure for the 40-angstrom wells
qwwad_ef_plot --style wf --energyfile Ee-40.r --wffileext "-40.r" --totalpotentialfile "v-40.r" --plotfile $outfile_wf

cat << EOF
Results have been written to $outfile and $outfile_wf, containing
the energies of states as a function of barrier width, and a plot of
the wavefunctions respectively.

$outfile is in the format:

  COLUMN 1 - Barrier width [angstrom]
  COLUMN 2 - Energy of state |1> [meV]
  COLUMN 3 - Energy of state |2> [meV]

$outfile_wf is in the format:

  COLUMN 1 - Spatial location [angstrom]
  COLUMN 2 - Potential [meV] or wavefunction amplitude [a.u.]

  The file contains 3 data sets, each set being separated
  by a blank line, representing the potential profile and the
  wavefunctions for the first two states:

  SET 1 - Potential profile
  SET 2 - Wavefunction for state |1> offset by state energy [meV]
  SET 2 - Wavefunction for state |2> offset by state energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
