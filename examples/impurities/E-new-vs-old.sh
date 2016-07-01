#! /bin/sh
set -e

# Calculate 3D wave function donor energy as a function of well width for
# different barrier heights
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.5
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
outfile=E-new-vs-old.dat
rm -f $outfile

# Define structure
cat > s.r << EOF
200 0.1 0.0
100 0.0 0.0
200 0.1 0.0
EOF

# Generate alloy profile
qwwad_mesh --dzmax 1

# Generate potential profile for Ga(1-x)Al(x)As, can use defaults
qwwad_ef_band_edge --mass 0.067 --bandedgepotentialfile v.r 

# Calculate electron energy for same quantum well but without donor
qwwad_ef_generic --nst 1

# Energy without donor present
E0=`awk '{printf(" %20.17e\n",$2)}' Ee.r`

export QWWAD_DONORPOSITION

for QWWAD_DONORPOSITION in `seq 0 30 210` 240 250 260 290 320;
do
    qwwad_ef_generic --nst 1
    qwwad_ef_donor_generic

    # Energy with donor present
    Enew=`awk '{printf(" %e",$2)}' Ee.r`

    # Start donor binding energy calculation
    qwwad_ef_donor_specific --lambdastart 10 --lambdastop 300 --symmetry 3D
    Eold=`awk '{printf(" %e",$2)}' Ee.r`

    printf "%f %f %f %f\n" $QWWAD_DONORPOSITION $Eold $Enew $E0 >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Carrier energy [meV] (semi-analytical)
  COLUMN 3 - Carrier energy [meV] (direct integration)
  COLUMN 4 - Carrier energy [meV] (no doping)

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
