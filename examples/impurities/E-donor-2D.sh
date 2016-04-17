#!/bin/sh
set -e

# Plot binding-energy for a 2D donor wave function at various positions in a QW
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
outfile_E=E-donor-2D-total.dat
outfile_ED=E-donor-2D-binding.dat
outfile_lambda=E-donor-2D-lambda.dat
rm -f $outfile_E $outfile_ED $outfile_lambda

# Define structure
cat > s.r << EOF
200 0.1 0.0
60  0.0 0.0
200 0.1 0.0
EOF

# Set a fixed effective mass throughout structure
export QWWAD_MASS=0.096

# Generate CdMnTe band profile
qwwad_mesh
qwwad_ef_band_edge --material cdmnte --bandedgepotentialfile v.r

# Calculate energy of structure without a donor, for deduction
# of donor binding energy
qwwad_ef_generic

# Store the ground-state for later
E1=`awk '/^1/{print $2}' Ee.r`

# Loop over donor positions
for r_d in `seq 0.0 20 220` 230; do
    # Donor state calculation
    qwwad_ef_donor_specific --dcpermittivity 10.6 --lambdastart 25 --lambdastop 300 --donorposition $r_d

    # Append the total energy and Bohr radius to output files
    awk '{print r_d, $2, E1}' r_d=$r_d E1=$E1 < Ee.r >> $outfile_E
    awk '{print r_d, $2}'     r_d=$r_d < l.r >> $outfile_lambda
done

# Process the energy file to get the binding energy
awk '{print $1, $3-$2}' < $outfile_E > $outfile_ED

cat << EOF
Results have been written to $outfile_E, ${outfile_ED} and ${outfile_lambda}.

$outfile_E is in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Total carrier energy [meV]
  COLUMN 3 - Carrier energy without donor [meV]

$outfile_ED is in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Binding energy [meV]

$outfile_lambda is in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Bohr radius [Angstrom]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
