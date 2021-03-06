#!/bin/sh
# Define output file
set -e

# Calculates the effective mass at the bottom edge of the lowest miniband
# in a Kronig-Penney superlattice with a range of well/barrier widths
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2015
#     Paul Harrison <p.harrison@shu.ac.uk>
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
tempfile=Ee-k-LW.r
outfile=infinite-SL-mass.dat
rm -f $tempfile $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.4
export QWWAD_BARRIERPOTENTIAL=334.1965

# Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
# Use MB=0.067+0.083*x
export QWWAD_BARRIERMASS=0.1002

# Define number of states
export QWWAD_NST=1

# Loop for well and barrier widths
for LW in `seq 1 1 60`
do
    printf "%f\t" "$LW" >> $tempfile # write well width to file [angstrom]

    # calculate E(k) at k=0 and k=dk for 2nd derivative near Gamma
    # Note that we only need to evaluate one side of the dispersion curve
    # since it is even-symmetric about k=0
    for K in 0.0 0.01
    do
        qwwad_ef_superlattice --wellwidth $LW --barrierwidth $LW --wavevector $K

        # Send E(K) to file, in meV
        awk '{E = $2 * 1.60219e-22;
              printf("%22.17e\t", E)}' Ee.r >> $tempfile
    done # done loop over k

    printf "\n" >> $tempfile	# end line
done	# loop over LW

# Now calculate the effective mass along growth (z-) axis
awk '{
        hBar=1.054588757e-34;
        m0=9.109534e-31;
        pi = 3.14159;

        # Compute second derivative using E(0) and E(dk)
        E_0  = $2;
        E_dk = $3;

        # Step size in wave-vector [1/m]
        dk = 0.01*pi / ($1 * 1e-10)

        D2 = 2*(E_dk - E_0) / (dk*dk);
        mstar = hBar*hBar/D2/m0;
        printf("%f %f\n",$1,mstar)
     }' $tempfile >> $outfile
     
cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Well (and barrier) width [angstrom]
  COLUMN 2 - Effective mass [relative to m0]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm Ee.r $tempfile
