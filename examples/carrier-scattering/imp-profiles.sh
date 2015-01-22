#! /bin/sh
set -e

# Impurity scattering with various doping profiles 
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
outfile=imp-profiles.dat
rm -f $outfile

# Generate square well with delta doping in middle of well
# Set volume doping to give sheet density of 1e10 cm^{-2} in each level
cat > s.r << EOF
200 0.15 0.0
195 0.0  0.0
10  0.0  1e17
195 0.0  0.0
200 0.15 0.0
EOF

nz=801

find_heterostructure --nz-1per $nz
efxv

# Solve Schroedinger equation
efss --nst 2

# Define subband populations in file `N.r'
densityinput --type even

# Define required rate
echo "2 1" > rrp.r

for T in `seq 10 10 300`; do
    # Calculate distribution function
    sbp --Te $T

    # Find impurity scattering WITH screening
    imp --temperature $T

    rate=`awk '{print $3}' imp-avg.dat`

    printf "%d %e\n" $T $rate >> $outfile
done

printf "\n" >> $outfile

# Generate square well with constant doping in entire well
# Set volume doping to give sheet density of 1e10 cm^{-2} in each level
cat > s.r << EOF
200 0.15 0.0
400 0.0  2.5e15
200 0.15 0.0
EOF

nz=801

find_heterostructure --nz-1per $nz
efxv

# Solve Schroedinger equation
efss --nst 2

# Define subband populations in file `N.r'
densityinput --type even

# Define required rate
echo "2 1" > rrp.r

for T in `seq 10 10 300`; do
    # Calculate distribution function
    sbp --Te $T

    # Find impurity scattering WITH screening
    imp --temperature $T

    rate=`awk '{print $3}' imp-avg.dat`

    printf "%d %e\n" $T $rate >> $outfile
done

printf "\n" >> $outfile

# Generate square well with modulation doping in barriers
# Set volume doping to give sheet density of 1e10 cm^{-2} in each level
cat > s.r << EOF
200 0.15 2.5e15
400 0.0  0.0
200 0.15 2.5e15
EOF

nz=801

find_heterostructure --nz-1per $nz
efxv

# Solve Schroedinger equation
efss --nst 2

# Define subband populations in file `N.r'
densityinput --type even

# Define required rate
echo "2 1" > rrp.r

for T in `seq 10 10 300`; do
    # Calculate distribution function
    sbp --Te $T

    # Find impurity scattering WITH screening
    imp --temperature $T

    rate=`awk '{print $3}' imp-avg.dat`

    printf "%d %e\n" $T $rate >> $outfile
done
cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Temperature [K]
  COLUMN 2 - |2> -> |1> scattering rate [s^{-1}]

The file contains 3 data sets:

  SET 1 - Delta-doping in centre of well
  SET 2 - Constant doping throughout well
  SET 3 - Constant doping in barriers

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# rm -f *.r
