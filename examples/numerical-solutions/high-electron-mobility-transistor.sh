#!/bin/sh
set -e

# Find eigenstates of a HEMT structure self-consistently
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
outfile_V=high-electron-mobility-transistor-V.dat
outfile_wf=high-electron-mobility-transistor-wf.dat
rm -f $outfile_V $outfile_wf wf_e1*.r v*.r

# First generate structure definition `s.r' file
cat > s.r << EOF
200 0.2 2e17
200 0.0 0.0
EOF
 
qwwad_mesh          # generate alloy concentration as a function of z
qwwad_ef_band_edge  # generate potential data

cp v_b.r v.r # Copy conduction-band potential to total potential for first iteration
  
for I in `seq 0 6`; do
 # Calculate ground state energy
 qwwad_ef_generic --nstmax 1

 qwwad_population_init # Generate an estimate of the population density
 qwwad_charge_density  # Compute charge density profile

 # save wave function and potential in separate files
 cp wf_e1.r wf_e1-I=$I.r
 cp v.r v-I=$I.r

 # Implement self consistent Poisson calculation
 qwwad_poisson --mixed
done

# Save data to output files
awk '{print $1*1e10, $2}' wf_e1-I=0.r >> $outfile_wf
printf "\n"                           >> $outfile_wf
awk '{print $1*1e10, $2}' wf_e1-I=2.r >> $outfile_wf
printf "\n"                           >> $outfile_wf
awk '{print $1*1e10, $2}' wf_e1-I=6.r >> $outfile_wf

awk '{print $1*1e10, $2*1000/1.6e-19}' v-I=0.r >> $outfile_V
printf "\n"                                    >> $outfile_V
awk '{print $1*1e10, $2*1000/1.6e-19}' v-I=2.r >> $outfile_V
printf "\n"                                    >> $outfile_V
awk '{print $1*1e10, $2*1000/1.6e-19}' v-I=6.r >> $outfile_V

cat << EOF
Results have been written to $outfile_V and
$outfile_wf.

$outfile_V contains the total potential profile in the format:

  COLUMN 1 - spatial location [Angstrom]
  COLUMN 2 - potential [meV]

$outfile_wf contains the wave function of the ground state
in the format:

  COLUMN 1 - spatial location [Angstrom]
  COLUMN 2 - wave function amplitude [m^{-0.5}]

Each file contains 3 data sets representing different
iterations of the self-consistent solution:

  SET 1 - Iteration 0
  SET 2 - Iteration 2
  SET 3 - Iteration 6

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
