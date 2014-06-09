#!/bin/sh
set -e

# Calculates the ground state energy in a finite GaAs quantum well
# with a range of barrier masses
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
outfile=infinite-barrier-mass-limit.dat
rm -f $outfile

# Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
# Use V=0.67*1247*x, keep x=0.4
V=334.1965

LW=100 # Well width [angstrom]
MW=0.067 # Well mass (GaAs)

calculate_plot_using_well_width()
{
    LW=$1

    tempfile=Ee-mb-${LW}.dat

    # Loop for different barrier masses
    for iMB in `seq -2 0.1 3`
    do
        {
            # Generate barrier mass exponentially so we get a smooth curve
            MB=`echo $iMB | awk '{print 10^$1}'`
            printf "%8.3e\t" "$MB" >> $tempfile	# write potential to file

            # Calculate ground state energy for different well and barrier masses
            efsqw --well-width $LW --well-mass $MW --barrier-mass $MB --potential $V
            awk '{printf("%8.3f\n",$2)}' Ee.r >> $tempfile   # send data to file
        }
    done
}

for width in 20 50 100 200; do
    calculate_plot_using_well_width $width

    cat Ee-mb-${width}.dat >> $outfile
    printf "\n" >> $outfile
done

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Barrier effective mass [relative to free electron mass]
  COLUMN 2 - Ground state in finite well [meV]

The file contains 4 data sets for various well widths:

  SET 1 - Width = 20 angstrom
  SET 2 - Width = 50 angstrom
  SET 3 - Width = 100 angstrom
  SET 4 - Width = 200 angstrom

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

rm Ee-mb-*.dat wf_*.r Ee.r
