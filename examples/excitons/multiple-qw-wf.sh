#! /bin/sh
set -e

# Calculate wavefunctions and overlap probability in a multiple quantum well system
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.6
#
# (c) Copyright 1996-2016
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

# Initialise files
outfile_wf=multiple-qw-wf.dat
outfile_p=multiple-qw-p-vs-a.dat
rm -f $outfile_wf $outfile_p

xb=0.1
LW=50

cat > s.r << EOF
200 $xb 0.0
$LW 0.0 0.0
$LW $xb 0.0
$LW 0.0 0.0
$LW $xb 0.0
$LW 0.0 0.0
$LW $xb 0.0
$LW 0.0 0.0
$LW $xb 0.0
$LW 0.0 0.0
200 $xb 0.0
EOF

qwwad_mesh --dzmax 1
qwwad_ef_band_edge --mass 0.067 --bandedgepotentialfile v.r
qwwad_ef_generic

qwwad_ef_band_edge --mass 0.62 --particle h --bandedgepotentialfile v.r
qwwad_ef_generic --particle h
 
# Find exciton binding energies
qwwad_ef_exciton --lambdastart 80    \
	         --betastart   0.001 \
		 --betastep    0.01

# Output wavefunction data to file
awk '{print $1*1e10, $2}' wf_e1.r >> $outfile_wf
printf "\n" >> $outfile_wf
awk '{print $1*1e10, $2}' wf_h1.r >> $outfile_wf
printf "\n" >> $outfile_wf
awk '{print $1*1e10, $2*70000}' x.r >> $outfile_wf

mv p.r $outfile_p

cat << EOF
Results have been written to $outfile_wf and $outfile_p.

$outfile_wf contains the wavefunctions and barrier potential
in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Wavefunction or barrier potential

  The file contains 3 data sets, each set being separated
  by a blank line, representing either the wavefunction or
  barrier potential in arbitrary units:

  SET 1 - Electron wavefunction
  SET 2 - Hole wavefunction
  SET 3 - Barrier potential

$outfile_p contains the separation probability in the format:

  COLUMN 1 - Electron/hole separation [Angstrom]
  COLUMN 2 - Probability [1/Angstrom]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
#rm -f *.r
