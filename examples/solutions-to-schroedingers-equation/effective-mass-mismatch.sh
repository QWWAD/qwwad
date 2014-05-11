#!/bin/sh

# Calculates the ground state energy in a finite GaAs quantum well
# with a range of barrier heights; both with and without effective mass
# mismatch being included in the model.
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
outfile=effective-mass-mismatch-E-L.dat
outfile_DE=effective-mass-mismatch-DE-L.dat
rm $outfile $outfile_DE

# Sweeps through a range of well widths and calculates
# ground state energy.
#
# Argument 1: Alloy composition to be used
sweep_well_width_using_alloy()
{
    rm -f Ee-lw.r

    # Read the barrier alloy from command-line argument
    if [ "z$1" = "z" ]; then
        echo "You need to specify the barrier alloy as an argument to this script"
        echo "For example, to use 10% alloy in the barrier, enter:"
        echo "  $0 0.1"
        exit 1
    else
        x=$1
    fi

    # Calculate conduction band barrier height for GaAs/Ga(1-x)Al(x)As
    # Use V=0.67*1247*x
    V=`echo 0.67*1247*$x|bc`

    # Calculate bulk effective mass of electron in Ga(1-x)Al(x)As
    # Use MB=0.067+0.083*x
    MB=`echo 0.067 + 0.083*$x|bc`
    MW=0.067 # Well mass (GaAs)

    # Loop for different well widths
    for i in `seq 1 0.05 2.3`
    do
        # Generate well-widths exponentially so we get a smooth curve at small
        # widths
        LW=`echo $i | awk '{print 10^$1}'`

        printf "%e\t" "$LW" >> Ee-lw.r	# write well width to file

        # Calculate ground state energy for barrier mass (MB)=well mass (0.067)
        efsqw -a $LW -m $MW -n $MW --potential $V
        awk '{printf("%8.3f",$2)}' Ee.r >> Ee-lw.r	# send data to file

        # Calculate ground state energy for correct barrier mass 
        efsqw -a $LW -m $MW -n $MB --potential $V
        awk '{printf("%8.3f\n",$2)}' Ee.r >> Ee-lw.r  # send data to file
    done

    cp Ee-lw.r x\=$x-Ee-lw.r
}

# Loop over alloy fractions
for lw in 0.1 0.2 0.3 0.4; do
    sweep_well_width_using_alloy $lw

    cat x\=$x-Ee-lw.r >> $outfile
    printf "\n"       >> $outfile

    awk '{print $1, ($3 > $2 ? $3-$2 : $2 - $3)}' < x\=$x-Ee-lw.r >> $outfile_DE
    printf "\n"       >> $outfile_DE
done

cat << EOF
Results have been written to $outfile and $outfile_DE.

$outfile contains data in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Energy of ground state with constant mass [meV]
  COLUMN 3 - Energy of ground state with correct barrier mass [meV]

$outfile_DE contains data in the format:

  COLUMN 1 - Well width [Angstrom]
  COLUMN 2 - Difference in energy between constant mass and correct mass
             models [meV]

Both files contain 4 data sets, with each set being separated
by a blank line, representing different barrier compositions:

  SET 1 - 10% alloy in barriers
  SET 2 - 20% alloy in barriers
  SET 3 - 30% alloy in barriers
  SET 4 - 40% alloy in barriers

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# clean up workspace
rm -f *Ee-lw* wf_* Ee.r
