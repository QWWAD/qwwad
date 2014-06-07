#!/bin/sh
set -e

# Calculates the ground state energy in a finite GaAs quantum well
# with a range of barrier heights
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2014
#     Paul Harrison <p.harrison@shu.ac.uk>
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
outfile=infinite-barrier-height-limit.dat
rm -f $outfile

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x for x=0.4---keep constant
MB=0.1002
MW=0.067 # Well mass (GaAs)

# Keep well width constant at 100 Angstom
LW=100

# Calculate the infinite well solution and store ground
# state energy as Einf
efiw --width $LW --mass $MW
Einf=`awk '{printf("%8.3f", $2)}' Ee.r`

# Loop for different barrier heights
for iV in `seq 2 0.1 8`
do
{
 # Generate potentials exponentially so we get a smooth curve at small
 # widths
 V=`echo $iV | awk '{print 10^$1}'`

 # Write potential and infinite well solution to file
 printf "%f\t%f\t" "$V" "$Einf" >> $outfile

 # Calculate ground state energy for barrier mass (MB)=well mass (0.067)
 efsqw -a $LW -m $MW -n $MW --potential $V
 awk '{printf("%8.3f",$2)}' Ee.r >> $outfile # send data to file

 # Calculate ground state energy for different well and barrier masses
 efsqw -a $LW -m $MW -n $MB --potential $V
 awk '{printf("%8.3f\n",$2)}' Ee.r >> $outfile # send data to file
}
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Barrier height [meV]
  COLUMN 2 - Ground state in infinite well [meV]
  COLUMN 3 - Ground state in finite well (constant mass) [meV]
  COLUMN 4 - Ground state in finite well (variable mass) [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up the workspace
rm wf_e?.r Ee.r
