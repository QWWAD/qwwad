#! /bin/sh
set -e

# Calculation of the mean e-LO scattering rate over two subband populations
# as a function of well width
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2014
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
OUT=sr-LW.r
rm -f $OUT

# Define carrier and lattice temperatures as being identical
Tl=77
Te=77

# Define required rate
echo "2 1" > rrp.r

# Define subband populations in file `N.r' [1e10 cm^{-2} in each]
N=1e14 # Define in m^{-2}
cat > N.r << EOF
$N
$N
EOF

nz=601

# Loop over infinite well width
for LW in `seq 70 5 360`; do
 # Generate infinitely deep well solutions
 efiw --width $LW --nz $nz --nst 2

 # Save lowest two energies to file
 E1=`sed -n 1p < Ee.r | awk '{print $2}'`
 E2=`sed -n 2p < Ee.r | awk '{print $2}'`
 DE=`echo $E1 $E2 | awk '{print $2-$1}'`

 printf "%f " $DE >> $OUT

 # Calculate the distribution functions at 77K
 sbp --Te $Te
 
 # Calculate scattering rate
 srelo --Te $Te --Tl $Tl --Ecutoff 800 --noblocking

 # Sort and store in output file
 awk '{print $3}' LOe-if.r >> $OUT
done
