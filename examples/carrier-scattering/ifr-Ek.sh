#! /bin/sh
set -e

# Interface roughness scattering for a finite QW as a function of initial KE
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
outfile=ifr-Ek.dat
rm -f $outfile

# Define well width
LB=400
LW=400
nz=801

# Define required rate
cat > rrp.r << EOF
2 2
2 1
1 2
1 1
EOF

# Electron temperature
T=300
x=0.15

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

# Find interface-roughness scattering
ifr --temperature $T --Ecutoff 800
cat ifr22.r >> $outfile
printf "\n" >> $outfile
cat ifr21.r >> $outfile
printf "\n" >> $outfile
cat ifr12.r >> $outfile
printf "\n" >> $outfile
cat ifr11.r >> $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Alloy fraction in barrier
  COLUMN 2 - |2> -> |2> scattering rate [s^{-1}]
  COLUMN 3 - |2> -> |1> scattering rate [s^{-1}]

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
