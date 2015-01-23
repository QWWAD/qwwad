#! /bin/sh
set -e

# Calculation of the mean e-e scattering rate over two subband populations
# as a function of well width
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2014
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
outfile=cc-avg-DE.dat
rm -f $outfile

# Define subband populations in file `N.r' [1e10 cm^{-2} in each]
N=1e14 # Define in m^{-2}
cat > N.r << EOF
$N
$N
EOF

# Define the required e-e rate
echo 2 2 1 1 > rr.r

# Loop over infinite well width
for LW in 100 200 300 400 500 600
do
 # Generate infinitely deep well solutions
 efiw -L $LW -N 300 -s 2

 # Save lowest two energies to file
 E1=`sed -n 1p < Ee.r | awk '{print $2}'`
 E2=`sed -n 2p < Ee.r | awk '{print $2}'`
 DE=`echo $E1 $E2 | awk '{print $2-$1}'`

 printf "%f " $DE >> $outfile

 # Loop over different temperatures
 for T in 4 300
 do
     # Calculate the distribution functions
     sbp --Te $T

     # Calculate carrier-carrier (e-e) scattering rate
     srcc -T $T

     # Sort and store in output file
     awk '/2 2 1 1/{printf("%e ",$5)}' ccABCD.r >> $outfile
 done

# End line in output file
printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Separation between |2> and |1> subband minima [meV]
  COLUMN 2 - Average |2,2> -> |1,1> scattering rate at 4 K [s^{-1}]
  COLUMN 3 - Average |2,2> -> |1,1> scattering rate at 300 K [s^{-1}]

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
#rm -f *.r
