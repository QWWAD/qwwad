#! /bin/sh
set -e

# Calculates density of states in bulk and 2D system as a function of energy
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2015
#     Paul Harrison  <p.harrison@shu.ac.uk>
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

# Initialise output files
outfile_bulk=density-of-states-bulk.dat
outfile_2D=density-of-states-2D.dat
rm -f $outfile_bulk $outfile_2D

# Set fixed parameters
L=200 # Well-width (angstrom)
nst=3 # Number of states to find

# Find first three states in an infinite well
efiw --width $L --nst $nst

# Find density of states for bulk and 2D system
dos

awk '{print $1, $2*(1.6e-19*1e-27)}' rho.r >> $outfile_bulk
awk '{print $1, $3*(1.6e-19*1e-18)}' rho.r >> $outfile_2D

# Repeat for a nonparabolic well
alpha=0.7
efiw --width $L --nst $nst --alpha $alpha
dos --alpha $alpha

printf "\n" >> $outfile_bulk
printf "\n" >> $outfile_2D

awk '{print $1, $2*(1.6e-19*1e-27)}' rho.r >> $outfile_bulk
awk '{print $1, $3*(1.6e-19*1e-18)}' rho.r >> $outfile_2D

cat << EOF
Results have been written to $outfile_bulk and $outfile_2D.

$outfile_bulk contains the density-of-states for bulk GaAs in the format:

  COLUMN 1 - Energy above band edge [meV]
  COLUMN 2 - Density of states [eV^{-1} nm^{-3}]

$outfile_2D contains the density-of-states for an infinite quantum well:

  COLUMN 1 - Energy above band edge [meV]
  COLUMN 2 - Density of states [eV^{-1} nm^{-2}]
  
Each output file contains two data sets, which are separated by a blank line.

  SET 1 - Results for parabolic dispersion
  SET 2 - Results for nonparabolic dispersion (alpha=0.7 eV^{-1})

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

rm -f rho.r Ee.r wf_*.r
