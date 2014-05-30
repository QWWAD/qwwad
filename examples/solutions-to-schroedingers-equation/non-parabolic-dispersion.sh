#! /bin/sh
set -e

# Calculates subband dispersion in infinite quantum wells with various
# nonparabolicity parameters
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 2014
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

# Initialise output file
outfile=non-parabolic-dispersion.dat
rm -f $outfile

# Set fixed parameters (200 angstrom, GaAs well)
mass=0.067
LW=200
nst=2  # Number of subbands to find

# Find parabolic solution first
efiw --width $LW --mass $mass --states $nst
awk '{print $1, mass*9.11e-31}' mass=$mass < wf_e1.r > massd.dat

# A little hack needed to provide the subband class with the necessary
# data for construction (just fills some arrays with placeholder data)
awk '{print 1}' < Ee.r > populations.dat
cp Ee.r Ef.dat

# Compute the dispersion relation for all states in the system
dispersion_relation --parabolic --disp-ext "_0.dat"

# Output a zero potential profile to file (again, just a hack to
# make the subband class play nicely)
awk '{print $1, 0}' < wf_e1.r > Vtotal.dat

# Repeat the calculation using nonparabolicity
for alpha in 0.7 5; do
    efiw --width $LW --mass $mass --states $nst --alpha $alpha

    # Rescale energies to Joules for use with dispersion-relation code
    awk '{print $1, alpha/1.6e-19}' alpha=$alpha < wf_e1.r > alphad.dat

    dispersion_relation --disp-ext "_$alpha.dat"
done

# Now, glue all our output files together into one convenient data file
awk '{print $1/1e9, $2}' dr_e1_0.dat >> $outfile
printf "\n" >> $outfile
awk '{print $1/1e9, $2}' dr_e2_0.dat >> $outfile
printf "\n" >> $outfile
awk '{print $1/1e9, $2}' dr_e1_0.7.dat >> $outfile
printf "\n" >> $outfile
awk '{print $1/1e9, $2}' dr_e2_0.7.dat >> $outfile
printf "\n" >> $outfile
awk '{print $1/1e9, $2}' dr_e1_5.dat >> $outfile
printf "\n" >> $outfile
awk '{print $1/1e9, $2}' dr_e2_5.dat >> $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Wave vector [nm^{-1}]
  COLUMN 2 - Wavefunction amplitude [a.u.]

The file contains six data sets, separated by a blank line,
which represent the following dispersion relations in an infinite
well:

  SET 1 - State 1, parabolic dispersion
  SET 2 - State 2, parabolic dispersion
  SET 3 - State 1, nonparabolic (alpha = 0.7)
  SET 4 - State 2, nonparabolic (alpha = 0.7)
  SET 5 - State 1, nonparabolic (alpha = 5)
  SET 6 - State 2, nonparabolic (alpha = 5)

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f dr_*.dat Vtotal.dat alphad.dat Ee.* wf_* populations.dat massd.dat Ef.dat
