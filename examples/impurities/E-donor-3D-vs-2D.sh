#! /bin/sh
set -e

# Compare total energy for 2D and 3D donors at various locations in a QW
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2014
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

outfile=E-donor-3D-vs-2D.dat
outfile_wf=E-donor-3D-wf.dat
outfile_wf_diff=E-donor-3D-wf-vs-2D.dat
rm -f $outfile $outfile_wf $outfile_wf_diff

# 3D trial wave function donor binding energy calculation

# Define structure
cat > s.r << EOF
200 0.1 0.0
60  0.0 0.0
200 0.1 0.0
EOF

# Define donor positions
cat > r_d.r << EOF
0.0e-10
2.0e-9
4.0e-9
6.0e-9
8.0e-9
1.0e-8
1.2e-8
1.4e-8 
1.6e-8 
1.8e-8 
2.0e-8 
2.2e-8
2.3e-8
EOF

# Generate alloy profile
find_heterostructure --dz-max 1

# Generate potential profile
efxv --material cdmnte --mass 0.096

# Perform 2D donor calculation and save results to output file
d02D --mass 0.096 --epsilon 10.6 --lambdastart 25 --lambdastop 300 --symmetry 2D > garbage.r
mv e.r e-2D.r

# Perform 3D donor calculation and save results to file
d02D --mass 0.096 --epsilon 10.6 --lambdastart 25 --lambdastop 300 --symmetry 3D > garbage.r

paste e-2D.r e.r | awk '{print $1, $4 - $2}' > $outfile

# Pack all wavefunction files into a single output file
for iwf in `seq 0 12`; do
    # Figure out the filename
    wf_file=wf${iwf}.r

    awk '{print $1*1e10, $2 + iwf*5000}' iwf=$iwf $wf_file >> $outfile_wf
    printf "\n" >> $outfile_wf
done

# Copy desired wavefunctions to output file
awk '{print $1*1e10, $2}' wf0.r >> $outfile_wf_diff
printf "\n" >> $outfile_wf_diff
awk '{print $1*1e10, $3}' wf0.r >> $outfile_wf_diff
printf "\n" >> $outfile_wf_diff
awk '{print $1*1e10, $2+10000}' wf10.r >> $outfile_wf_diff
printf "\n" >> $outfile_wf_diff
awk '{print $1*1e10, $3+10000}' wf10.r >> $outfile_wf_diff
printf "\n" >> $outfile_wf_diff
awk '{print $1*1e10, $2+20000}' wf12.r >> $outfile_wf_diff
printf "\n" >> $outfile_wf_diff
awk '{print $1*1e10, $3+20000}' wf12.r >> $outfile_wf_diff
printf "\n" >> $outfile_wf_diff

cat << EOF
Results have been written to $outfile, $outfile_wf and $outfile_wf_diff

$outfile is in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Difference in total carrier energy (3D - 2D) [meV]

$outfile_wf is in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Wavefunction amplitude [m^{-1/2}] (offset)

  The file contains 13 data sets, each set being separated
  by a blank line, representing donors at different locations
  relative to the left-hand-side of the structure:

  SET 1  - Donor location = 0 Angstrom
  SET 2  - Donor location = 20 Angstrom
  SET 3  - Donor location = 40 Angstrom
  SET 4  - Donor location = 60 Angstrom
  SET 5  - Donor location = 80 Angstrom
  SET 6  - Donor location = 100 Angstrom
  SET 7  - Donor location = 120 Angstrom
  SET 8  - Donor location = 140 Angstrom
  SET 9  - Donor location = 160 Angstrom
  SET 10 - Donor location = 180 Angstrom
  SET 11 - Donor location = 200 Angstrom
  SET 12 - Donor location = 220 Angstrom
  SET 13 - Donor location = 230 Angstrom

$outfile_wf_diff is in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Wavefunction amplitude [m^{-1/2}] (offset)

  The file contains 6 data sets, each set being separated by a blank line,
  representing the wavefunctions (with or without hydrogenic contributions)
  of donors at different locations relative to the left-hand-side of the
  structure:

  SET 1  - Donor location = 0 Angstrom (full wavefunction)
  SET 2  - Donor location = 0 Angstrom (no hyd. component)
  SET 3  - Donor location = 200 Angstrom (full wavefunction)
  SET 4  - Donor location = 200 Angstrom (no hyd. component)
  SET 5  - Donor location = 230 Angstrom (full wavefunction)
  SET 6  - Donor location = 230 Angstrom (no hyd. component)

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
