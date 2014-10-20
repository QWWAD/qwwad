#!/bin/sh
set -e

# Computes self-consistent solution for a QW with a diffuse dopant profile
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
outfile_d=delta-doping-d.dat  # Charge density profile
outfile_V=delta-doping-vp.dat # Poisson potential
outfile_E=delta-doping-E.dat  # Ground-state energy
rm -f $outfile_d $outfile_V $outfile_E
rm -f wf_*.r Ee.r Eh.r d.r D.r

# Define band edge profile of single quantum well with a
# delta-doped layer
#
# TODO: The doping profile in Fig. 4.14 appears to use 1e18
#       doping, but the self-consistent solution appears to
#       use 2e18. Some investigation needed.
cat > s.r << EOF
50.0000 0.2 0.0
48.5875 0.0 0.0
 2.8250 0.0 2.0e18
48.5875 0.0 0.0
50.0000 0.2 0.0
EOF

# Generate quantum well profile and initial dopant profile
find_heterostructure --dz-max 0.25

# Run diffusion `simulation' for various times
for t in 0 10 20 50 100 200; do
    # Generate diffuse dopant profile using constant diffusion coefficient
    gde --coeff 1 --time $t --infile d.r --outfile D.r

    # Save doping profile to file
    awk '{print $1*1e10, $2/1e24}' D.r >> $outfile_d
    printf '\n' >> $outfile_d

    # Find valence band edge
    efxv --particle h

    cp v.r vvb.r # Save valence-band potential for use as a baseline

    # Perform iterative solution
    for I in `seq 0 7`; do
        # Calculate ground state Schroedinger solution
        efss --particle h --nst-max 1

        # Estimate population density and charge density profile
        densityinput --energyfile Eh.r --dopingfile D.r
        chargedensity --ptype --dopingfile D.r --wf-input-prefix wf_h --energy-input Eh.r

        # Implement self consistent Poisson calculation
        find_poisson_potential --Vbasefile vvb.r --ptype --field 0
    done

    # Write energy to output file
    E1=`awk '{printf("\t%20.17e\n",$2)}' Eh.r`
    echo $t $E1 >> $outfile_E

    # Write Poisson potential to output file
    awk '{print $1*1e10, $2*1000/1.6e-19}' v_p.r >> $outfile_V
    printf '\n' >> $outfile_V
done

cat << EOF
Results have been written to $outfile_d, $outfile_V and $outfile_E.

$outfile_d contains the dopant profiles in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Doping concentration [1e{18} cm^{-3}]

$outfile_V contains the Poisson potential in the format:

  COLUMN 1 - Spatial location [Angstrom]
  COLUMN 2 - Potential [meV]

$outfile_E contains the ground-state energy in the format:

  COLUMN 1 - Diffusion time
  COLUMN 2 - Energy [meV]

Each file contains 6 data sets, each set being separated
by a blank line, representing different diffusion times:

  SET 1 - t = 0 s
  SET 2 - t = 10 s
  SET 3 - t = 20 s
  SET 4 - t = 50 s
  SET 5 - t = 100 s
  SET 6 - t = 200 s

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
