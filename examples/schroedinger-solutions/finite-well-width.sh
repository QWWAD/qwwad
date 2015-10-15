#! /bin/sh
set -e

# Calculates the first 3 energy level in a finite GaAs quantum well
# with a range of well widths assuming parabolic dispersion.
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
outfile=finite-well-E-L.dat
rm -f Ee1-lw.r Ee2-lw.r Ee3-lw.r Ee4-lw.r $outfile

# Append data to the file Eei-lw.r if the state exists
# Param 1: Number of state
output_E_if_real()
{
 # Extract energies from output file
 Ee=`awk '{print $2}' Ee.r|sed -n "$1p"`

 # Check that the energy was read successfully.
 # If it was, then add it to the end of the output file
 if test "x$Ee" != "x"; then
     printf "%f\t%.17e\n" $LW $Ee >> Ee$1-lw.r
 fi
}

# Set fixed parameters
export QWWAD_NST=3 # Number of states to find

# Loop for different well widths
for i in `seq 1 0.05 2.3`
do
 # Generate well-widths exponentially so we get a smooth curve at small
 # widths
 LW=`echo $i | awk '{print 10^$1}'`

 # Calculate energy levels
 qwwad_ef_square_well --wellwidth $LW

 # Output energies to file
 for state in `seq 1 $QWWAD_NST`; do
	 output_E_if_real $state
 done
done

cat Ee1-lw.r >> $outfile
printf "\n"   >> $outfile
cat Ee2-lw.r >> $outfile
printf "\n"   >> $outfile
cat Ee3-lw.r >> $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Energy of state [meV]

  The file contains 3 data sets, each set being separated
  by a blank line, representing each state in the system:

  SET 1 - Results for state 1
  SET 2 - Results for state 2
  SET 3 - Results for state 3

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f wf*.r Ee.r Ee?-lw.r
