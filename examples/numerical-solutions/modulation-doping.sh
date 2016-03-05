#!/bin/sh
set -e

# Calculates the bandstructure for a modulation-doped quantum well
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.3
#
# (c) Copyright 1996-2016
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
outfile=modulation-doping-V.dat
rm -f $outfile wf_e1*.r

# First generate structure definition `s.r' file
cat > s.r << EOF
100 0.2 2e17
100 0.0 0.0
100 0.2 2e17
EOF
 
qwwad_mesh --dzmax 1	# generate alloy concentration as a function of z
qwwad_ef_band_edge      # generate potential data

cp v_b.r v.r # Copy conduction-band potential to total potential for first iteration
  
for I in 0 1 2 3; do
 # Calculate ground state Schroedinger solution
 qwwad_ef_generic --nstmax 1

 qwwad_population_init # Generate an estimate of the population density
 qwwad_charge_density  # Compute charge density profile

 # Implement self consistent Poisson calculation
 qwwad_poisson
done

# Save data to output file [meV]
awk '{print $1*1e10, $2/1.6e-19 * 1000}' v_b.r  > $outfile # band-edge potential
printf "\n"                                    >> $outfile # blank line
awk '{print $1*1e10, $2/1.6e-19 * 1000}' v.r   >> $outfile # total potential

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Potential [meV]

  The file contains 2 data sets, each set being separated
  by a blank line, representing:

  SET 1 - The conduction band edge
  SET 2 - The total potential, including charge effects

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
