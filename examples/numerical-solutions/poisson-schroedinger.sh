#!/bin/sh
set -e

# Find eigenstates of a doped well self-consistently
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

# Initialise files
outfile=poisson-schroedinger-E-I.dat
outfile_sigma=poisson-schroedinger-sigma.dat
outfile_field=poisson-schroedinger-field.dat
outfile_Vp=poisson-schroedinger-Vp.dat
outfile_Vt=poisson-schroedinger-Vt.dat
rm -f $outfile wf_e1*.r

# First generate structure definition `s.r' file
cat > s.r << EOF
200 0.2 0.0
100 0.0 2e18
200 0.2 0.0
EOF
 
find_heterostructure --dz-max 1	# generate alloy concentration as a function of z
efxv                            # generate potential data

cp v.r vcb.r     # Save conduction-band potential for use as a baseline

for I in `seq 0 7`; do
 # Calculate ground state energy
 efss --nst-max 1

 densityinput  # Estimate the population density in each state
 chargedensity # Compute charge density profile
 
 # save wave function in separate file
 cp wf_e1.r wf_e1-I=$I.r

 # Write energy to output file
 E1=`awk '{printf("\t%20.17e\n",$2)}' Ee.r`

 echo $I $E1 >> $outfile

 # Implement self consistent Poisson calculation
 find_poisson_potential --Vbasefile vcb.r
done # X

# Convert volume density to sheet density within slices
awk '{print $1*1e10, $2/(1e14*1.6e-19) * 1e-10}' rho.r > $outfile_sigma

awk '{print $1*1e10, $2/1e6}' field.r > $outfile_field
awk '{print $1*1e10, $2*1000/1.6e-19}' v_p.r > $outfile_Vp
awk '{print $1*1e10, $2*1000/1.6e-19}' v.r > $outfile_Vt

cat << EOF
Results have been written to $outfile,
$outfile_sigma,
$outfile_field,
$outfile_Vp and
$outfile_Vt.

$outfile contains the convergence of energy in the format:

  COLUMN 1 - iteration number
  COLUMN 2 - energy of ground state [meV]

$outfile_sigma contains the charge density profile in the format:

  COLUMN 1 - spatial location [Angstrom]
  COLUMN 2 - charge density [10^{14} e m^{-2}]

$outfile_field contains the electric field profile in the format:

  COLUMN 1 - spatial location [Angstrom]
  COLUMN 2 - electric field [10^6 V/m]

$outfile_Vp contains the Poisson potential profile in the format:

  COLUMN 1 - spatial location [Angstrom]
  COLUMN 2 - potential [meV]

$outfile_Vt contains the total potential profile in the format:

  COLUMN 1 - spatial location [Angstrom]
  COLUMN 2 - potential [meV]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
