#!/bin/sh
set -e

# Computes error in ground-state energy when nonparabolicity is ignored
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.3
#
# (c) Copyright 1996-2016
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
outfile=nonparabolic-dE-vs-width.dat
rm -f $outfile

export QWWAD_NSTMAX=1 # Only study ground states

# Loop over barrier concentrations
for X in 0.4 0.6 0.8 1.0; do

 # Loop over well width
 for LW in 20 22 24 26 28 30 35 40 50 60 80 100 120 140 160 180 200; do
     printf "\rCalculating barrier = %f; width = %d" $X $LW

     # First generate structure definition `s.r' file
     echo 200 $X  0.0  > s.r
     echo $LW 0.0 0.0 >> s.r
     echo 200 $X  0.0 >> s.r

     qwwad_mesh --dzmax 0.25 # generate alloy concentration as a function of z
     qwwad_ef_band_edge --bandedgepotentialfile v.r    # generate potential data, and bandgap

     # Calculate ground state energy with band non-parabolicity
     qwwad_ef_generic --solver shooting-nonparabolic

     # Get energy from file
     E1_np=`awk '/^1/{printf $2}' Ee.r`

     # Calculate ground state energy without band non-parabolicity
     qwwad_ef_generic

     # Get energy from file
     E1_parab=`awk '/^1/{printf $2}' Ee.r`

     # Now calculate difference between `with' and `without' band non-parabolicity
     dE=`echo $E1_parab $E1_np | awk '{print $2-$1}'`

     echo $LW $dE >> $outfile
 done # LW

 printf "\n" >> $outfile
done # X

printf "\n"

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well width [angstrom]
  COLUMN 2 - Error in ground-state energy [meV]

  The file contains 4 data sets, each set being separated
  by a blank line, representing <whatever>:

  SET 1 - 40% AlAs barriers
  SET 2 - 60% AlAs barriers
  SET 3 - 80% AlAs barriers
  SET 4 - 100% AlAs barriers

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
