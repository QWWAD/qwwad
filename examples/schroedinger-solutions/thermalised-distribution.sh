#!/bin/sh
set -e

# Calculates the distribution of carriers between the first 10 subbands
# in a quantum well, as a function of energy
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2015
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
outfile=thermalised-distribution.dat
rm -f Ef-T*.r N.r FD*.r N?-T.r $outfile

# Set global options
export QWWAD_ALPHA=0.701 # Nonparabolicity parameter for GaAs [1/eV]

# Solve single quantum well
qwwad_ef_infinite_well --wellwidth 200 --nst 6

# Set the global population
N=10

# Loop for different temperatures
for T in `seq 10 10 1000`
do
    # Calculate the global Fermi energy for the system
    qwwad_fermi_distribution --global-population $N --Te $T

    # Get the populations of each state from the data file 
    N1=`awk '/1\t/{printf("%8.3f\n", $2)}' N-out.r`
    N2=`awk '/2\t/{printf("%8.3f\n", $2)}' N-out.r`
    N3=`awk '/3\t/{printf("%8.3f\n", $2)}' N-out.r`
    N4=`awk '/4\t/{printf("%8.3f\n", $2)}' N-out.r`
    N5=`awk '/5\t/{printf("%8.3f\n", $2)}' N-out.r`
    N6=`awk '/6\t/{printf("%8.3f\n", $2)}' N-out.r`

    # Find cumulative populations up to, and including each subband
    N2all=`echo "$N1 + $N2"|bc`
    N3all=`echo "$N2all + $N3"|bc`
    N4all=`echo "$N3all + $N4"|bc`
    N5all=`echo "$N4all + $N5"|bc`
    N6all=`echo "$N5all + $N6"|bc`

    # Dump the cumulative populations to the output file
    printf "%d\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n" $T $N1 $N2all $N3all $N4all $N5all $N6all >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Temperature [K]
  COLUMN 2 - Population of subband 1 [x10^{10} cm^{-2}]
  COLUMN 3 - Cumulative population of subbands 1 & 2 [x10^{10} cm^{-2}]
  COLUMN 4 - Cumulative population of subbands 1 - 3 [x10^{10} cm^{-2}]
  COLUMN 5 - Cumulative population of subbands 1 - 4 [x10^{10} cm^{-2}]
  COLUMN 6 - Cumulative population of subbands 1 - 5 [x10^{10} cm^{-2}]
  COLUMN 7 - Cumulative population of subbands 1 - 6 [x10^{10} cm^{-2}]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f wf*.r v.r N-*.r Ee.r Ef.r
