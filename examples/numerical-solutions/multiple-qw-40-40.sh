#!/bin/sh
set -e

# Find ground-state of a multi-quantum-well system with varying number
# of wells.
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
outfile=multiple-qw-40-40-E-N.dat
outfile_wf=multiple-qw-40-40-wf.dat
rm -f $outfile $outfile_wf

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.2
export QWWAD_BARRIERPOTENTIAL=167.0985

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
export QWWAD_BARRIERMASS=0.0836

# Define well and barrier widths here
export QWWAD_WELLWIDTH=40
export QWWAD_BARRIERWIDTH=40

# Specify padding thickness at the start and end of the structure,
# such that startpadding = finalpadding + barrierwidth
startpadding=200
endpadding=160

# To compare with infinite superlattice
qwwad_ef_superlattice
E1_analytical=`awk '/^1/{print $2}' Ee.r`

# Start a temporary layer definition file that we'll use to "build"
# a stack of quantum wells inside the following loop
# First, we'll just define the initial "padding" layer
cat > s.r.tmp << EOF
$startpadding 0.2 0.0
EOF

# Loop over number of wells
for N in `seq 1 10`; do

# Now, add a quantum well and a barrier to the temporary layer file
cat >> s.r.tmp << EOF
$QWWAD_WELLWIDTH    0.0 0.0
$QWWAD_BARRIERWIDTH 0.2 0.0
EOF

# Finally, copy the temporary layer file and add a final bit of padding
# to give a symmetrical structure
cp s.r.tmp s.r

cat >> s.r << EOF
$endpadding 0.2 0.0 
EOF
		 
# Now generate mesh and bandstructure for the MQW system
qwwad_mesh
qwwad_ef_band_edge --bandedgepotentialfile v.r

# Solve the Schroedinger equation numerically
qwwad_ef_generic --nstmax 1

# Write energy to output file
E1_numerical=`awk '/^1/{print $2}' Ee.r`
printf "%d\t%e\t%e\n" $N $E1_analytical $E1_numerical >> $outfile
done

qwwad_ef_plot --style wf --plotfile $outfile_wf

cat << EOF
Results have been written to $outfile and
$outfile_wf, containing the energy of the ground-state
and the wavefunctions respectively.

$outfile is in the format:

  COLUMN 1 - Number of quantum wells
  COLUMN 2 - Energy of ground state in superlattice [meV]
  COLUMN 3 - Energy of ground state in MQW system [meV]

$outfile_wf is in the format:

  COLUMN 1 - Spatial location [angstrom]
  COLUMN 2 - Potential [meV] or wavefunction amplitude [a.u.]

  The file contains 2 data sets, each set being separated
  by a blank line, representing the potential profile and the
  wavefunctions for the first two states:

  SET 1 - Potential profile
  SET 2 - Wavefunction for state |1> offset by state energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
