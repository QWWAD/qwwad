#! /bin/sh
set -ve

# Impurity scattering with and without screening for an infinite QW 
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
outfile=imp-screening.dat
rm -f $outfile

# Define well width
LW=400

# Electron temperature
T=77

# Generate array of doping
d=1e18
cat > s.r << EOF
$LW 0.0 $d
EOF

nz=301

find_heterostructure --nz-1per $nz

# Solve infinite well
efiw --width $LW --nz $nz --nst 2

# Define subband populations in file `N.r'
densityinput --type even

# Calculate distribution function
sbp --Te $T

# Define required rate
echo "2 1" > rrp.r

# Find impurity scattering WITH screening
imp --temperature $T --Ecutoff 800

cat imp21.r > $outfile
printf "\n" >> $outfile

# Find impurity scattering WITHOUT screening
imp --temperature $T --Ecutoff 800 --noscreening
cat imp21.r >> $outfile

# <Edit accordingly, to describe each output file>
cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - <Whatever>
  COLUMN 2 - <Something>
  <...>

  <Remove this chunk if only 1 data set is in the file>
  The file contains <x> data sets, each set being separated
  by a blank line, representing <whatever>:

  SET 1 - <Description of the 1st set>
  SET 2 - <Description of the 2nd set>
  <...>

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# rm -f *.r
