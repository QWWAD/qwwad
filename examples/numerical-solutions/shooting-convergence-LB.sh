#!/bin/sh
set -e

# Calculates shooting method solutions as function of barrier width
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
outfile_E=shooting-convergence-LB-E.dat
outfile_wf=shooting-convergence-LB-wf.dat
rm -f $outfile_E $outfile_wf

# Loop over barrier width
for LB in `seq 60 20 300` 600; do

# Describe square well structure (thickness, alloy, doping)
cat > s.r << EOF
$LB 0.2 0
100 0.0 0
$LB 0.2 0
EOF

nz=`echo $LB | awk '{print 5*$1*2 + 501}'`

# Create alloy concentration file
find_heterostructure --nz-1per $nz

# Find band-edge parameters
efxv

# Find ground-state solution
efss --solver shooting-constant-mass --nst-max 1

E=`cut -f2 Ee.r`
echo $LB $E >> $outfile_E

mv wf_e1.r wf_$LB.r

done

# Copy selected wavefunctions to output file
# rescale lengths to Angstrom and offset
awk '{print $1*1e10, $2}' wf_60.r > $outfile_wf
printf "\n" >> $outfile_wf
awk '{print $1*1e10, $2 + 15000}' wf_200.r >> $outfile_wf
printf "\n" >> $outfile_wf
awk '{print $1*1e10, $2 + 30000}' wf_600.r >> $outfile_wf

cat << EOF
Results have been written to $outfile_E and ${outfile_wf}.

$outfile_E contains the energy of the ground state as a function
of barrier thickness, and is in the format:

  COLUMN 1 - Barrier thickness [angstrom]
  COLUMN 2 - Energy of ground state [meV]

$outfile_wf contains the wavefunction for the ground state and
is in the format:

  COLUMN 1 - Spatial location [angstrom]
  COLUMN 2 - Wavefunction amplitude [m^{-0.5}]

  The file contains 3 data sets, each set being separated
  by a blank line, representing the wavefunctions for different
  barrier thicknesses:

  SET 1 - Thickness = 60 Angstrom
  SET 2 - Thickness = 200 Angstrom
  SET 3 - Thickness = 600 Angstrom

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# <Delete all temporary files you created>
rm -f *.r
