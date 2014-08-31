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
# (c) Copyright 1996-2014
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
V=167.0985

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x, keep x=0.2
MB=0.0836

# Define well and barrier widths here
LW=40
LB=40

# To compare with infinite superlattice
efkpsl --well-width $LW --barrier-width $LB --barrier-mass $MB --potential $V
E1_analytical=`awk '/^1/{print $2}' Ee.r`

# Loop over number of wells
for N in 1 2 3 4 5 6 7 8 9 10; do
 
 # Write first barrier and initiate file
 echo 200 0.2 0.0 > s.r	

 # Could only think of an awk script as a replacement `for(i=1;i<N-1;i++)'
 # loop---I'm sure you must be able to do this within `/bin/sh'
 echo $N $LW $LB | awk '{
 Nwells=$1; LW=$2; LB=$3;
 while(Nwells-- > 1)
     {
         printf("%i 0.0 0.0\n",LW)	
         printf("%i 0.2 0.0\n",LB)	
     }
 }' >> s.r
		 
 # Write last well and barrier
 echo $LW 0.0 0.0 >> s.r
 echo 200 0.2 0.0 >> s.r

 find_heterostructure --dz-max 0.25 # generate alloy concentration as a function of z
 efxv			# generate potential data

 efss --nst-max 1

 # Write energy to output file
 E1_numerical=`awk '/^1/{print $2}' Ee.r`
 printf "%d\t%e\t%e\n" $N $E1_analytical $E1_numerical >> $outfile
done

wfplot --plot-wf --plot-file $outfile_wf

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

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
