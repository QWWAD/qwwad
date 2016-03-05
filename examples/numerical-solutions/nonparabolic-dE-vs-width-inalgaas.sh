#!/bin/sh
set -e

# Computes error in ground-state energy when nonparabolicity is ignored
# in an InGaAs/InAlAs well
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
outfile=nonparabolic-dE-vs-width-inalgaas.dat
rm -f $outfile

# Define alloy components X and Y in In(1-x-y)Al(x)Ga(y)As for both well
# and barrier

# Quantum wells, In(0.53)Ga(0.47)As
Yw=0.47
Xw=0.0

# Barriers, In(0.52)Al(0.48)As
Yb=0.0
Xb=0.48

export NSTMAX=1 # Only study ground states

# Loop over well width
for LW in 20 30 40 50 60 80 100 120 140 160 180 200; do
 # First generate structure definition `s.r' file
 echo 200 $Xb $Yb 0 > s.r
 echo $LW $Xw $Yw 0 >> s.r
 echo 200 $Xb $Yb 0 >> s.r
 
 # Remember to switch material system
 qwwad_mesh --dzmax 0.25  # generate alloy concentration as a function of z
 qwwad_ef_band_edge --material inalgaas --bandedgepotentialfile v.r # generate potential data, and bandgap
 
 # Calculate ground state energy with band non-parabolicity
 qwwad_ef_generic --solver shooting-nonparabolic

 # Get energy from file
 E_np=`awk '/^1/{printf $2}' Ee.r`

 # Calculate ground state energy without band non-parabolicity
 qwwad_ef_generic

 # Get energy from file
 E_parab=`awk '/^1/{printf $2}' Ee.r`

 # Now calculate difference between `with' and `without' band non-parabolicity
 dE=`echo $E_parab $E_np | awk '{print $2 - $1}'`

 echo $LW $dE >> $outfile
done # LW

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well width [angstrom]
  COLUMN 2 - Error in ground-state energy [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
