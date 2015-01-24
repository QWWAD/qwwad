#! /bin/sh
set -e

# Alloy disorder scattering for a finite QW with varying alloy composition in the barriers
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
outfile=ado-x-barrier.dat
rm -f $outfile

# Define well width
LB=400
LW=400
nz=801

# Define required rate
cat > rrp.r << EOF
2 2
2 1
EOF

# Electron temperature
T=4
for x in `seq 0.05 0.05 1`; do

# Generate array of doping
# Set volume doping to give sheet density of 1e10 cm^{-2} in each level
d=5e15
cat > s.r << EOF
$LB $x 0
$LW 0  $d
$LB $x 0
EOF

find_heterostructure --nz-1per $nz
efxv

# Solve well
efss --nst 2

# Define subband populations in file `N.r'
densityinput --type even

# Calculate distribution function
sbp --Te $T

# Find alloy disorder scattering
ado --temperature $T
rate22=`awk '/^2 2/{print $3}' ado-avg.dat`
rate21=`awk '/^2 1/{print $3}' ado-avg.dat`
printf "%f %e %e\n" $x $rate22 $rate21 >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Alloy fraction in well
  COLUMN 2 - |2> -> |2> scattering rate [s^{-1}]
  COLUMN 3 - |2> -> |1> scattering rate [s^{-1}]
  COLUMN 4 - |1> -> |2> scattering rate [s^{-1}]
  COLUMN 5 - |1> -> |1> scattering rate [s^{-1}]

The file contains 2 data sets:

  SET 1 - T = 4 K
  SET 2 - T = 300 K

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# rm -f *.r
