#!/bin/sh
set -e

# Calculates the ground state wavefunction in a finite GaAs quantum well
# with a range of barrier masses
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
outfile=large-barrier-mass-wf.dat
rm -f Ee-mb.r $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.4
V=334.196

# Keep well width constant at 100 Angstom
LW=100

# Loop for different barrier heights
for MB in 0.1002 1 10 100
do
    # Calculate ground state energy for different well and barrier masses
    efsqw --well-width $LW --well-mass 0.067 --barrier-mass $MB --potential $V
    mv wf_e1.r wf_e1-${MB}.r
done

awk '{print $1*1e10, $2}' wf_e1-0.1002.r >> $outfile
printf "\n" >> $outfile
awk '{print $1*1e10, $2 + 20000}' wf_e1-1.r >> $outfile
printf "\n" >> $outfile
awk '{print $1*1e10, $2 + 40000}' wf_e1-10.r >> $outfile
printf "\n" >> $outfile
awk '{print $1*1e10, $2 + 60000}' wf_e1-100.r >> $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Spatial location [angstrom]
  COLUMN 2 - Wavefunction amplitude [normalised]

  The values in column 2 are offset so that plots
  appear stacked above each other

The file contains 4 data sets for various barrier masses:

  SET 1 - Mass = 0.1002 m0
  SET 2 - Mass = 1 m0
  SET 3 - Mass = 10 m0
  SET 4 - Mass = 100 m0

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

rm Ee.r wf_e*.r
