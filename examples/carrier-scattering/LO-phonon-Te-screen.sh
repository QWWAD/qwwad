#! /bin/sh
set -e

# LO-phonon scattering for a finite QW with varying carrier temperature, with and without screening.
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
outfile=LO-phonon-Te-screen.dat
rm -f $outfile

# Define well width
LW=220
LB=60
nz=421

# Define required rate
cat > rrp.r << EOF
2 1
EOF

# Electron temperature
Tl=15
x=0.3

cat > s.r << EOF
$LB $x 0
$LW 0  0
$LB $x 0
EOF

find_heterostructure
efxv

# Solve well
efss --nst 2

for Te in `seq 15 5 300`; do 

# Generate array of doping [1e11 cm^{-2} in each]
cat > N.r << EOF
1e15
1e15
EOF

# Calculate distribution function
sbp --Te $Te

# Find scattering with screening
srelo --Te $Te --Tl $Tl
rate21_11_screened=`awk '/^2\t1/{print $3}' LOe-if.r`

srelo --Te $Te --Tl $Tl --noscreening
rate21_11_unscreened=`awk '/^2\t1/{print $3}' LOe-if.r`

# Generate array of doping [1e12 cm^{-2} in each]
cat > N.r << EOF
1e16
1e16
EOF

# Calculate distribution function
sbp --Te $Te

# Find scattering with screening
srelo --Te $Te --Tl $Tl
rate21_12_screened=`awk '/^2\t1/{print $3}' LOe-if.r`

srelo --Te $Te --Tl $Tl --noscreening
rate21_12_unscreened=`awk '/^2\t1/{print $3}' LOe-if.r`

printf "%d %e %e %e %e\n" $Te $rate21_11_screened $rate21_11_unscreened $rate21_12_screened $rate21_12_unscreened >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well width
  COLUMN 2 - |2> -> |1> scattering rate [s^{-1}] (with screening)
  COLUMN 3 - |2> -> |1> scattering rate [s^{-1}] (without screening)

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# rm -f *.r
