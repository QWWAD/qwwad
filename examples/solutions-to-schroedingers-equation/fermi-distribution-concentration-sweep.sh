#! /bin/sh
set -e

# Calculates quasi-Fermi distribution for lowest subband in a quantum well
# at a range of different electron concentrations.
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
outfile_f_E=fermi-distribution-concentration-sweep.dat
outfile_Ef_N=fermi-energy-concentration-sweep.dat
rm -f $outfile_f_E $outfile_Ef_N
rm -f v.r N.r Ef-N.r Ef-N-np.r

# Solve single quantum well
efiw -L 200 -s 1

# Need to artificially create a file v.r containing the barrier height,
# this used later in numerical solutions of Schrodinger's equation
echo 0.00000 1.60219e-19 > v.r	# i.e., 1eV

# Loop for different concentrations
for N in 0.1 0.2 0.5 1 2 5 10 20 50 100
do
    # Write populations file N.r first
    echo 1 $N > N.r

    # Calculate Fermi energies and population distribution
    sbp -f -T 77
    mv FD1.r FD1N=$N.r

    # Write Fermi energy to file
    echo -n "$N\t" >> $outfile_Ef_N
    awk '{printf("%8.3f\n",$2)}' Ef.r >> $outfile_Ef_N
done

# Output Fermi function vs. energy to plottable file
# for various concentrations
cat FD1N=1.r  >> $outfile_f_E
printf "\n"   >> $outfile_f_E
cat FD1N=2.r  >> $outfile_f_E
printf "\n"   >> $outfile_f_E
cat FD1N=5.r  >> $outfile_f_E
printf "\n"   >> $outfile_f_E
cat FD1N=10.r >> $outfile_f_E

cat << EOF
Results have been written to $outfile_f_E and $outfile_Ef_N.

$outfile_f_E contains the Fermi occupation number for states as a function of
energy in the format:

  COLUMN 1 - Energy (relative to band edge) [meV]
  COLUMN 2 - Fermi occupation

  The file contains four data sets for various electron concentrations, each
  set being separated by a blank line:

  SET 1 - Results at 1x10^{10} cm^{-2}
  SET 2 - Results at 2x10^{10} cm^{-2}
  SET 3 - Results at 5x10^{10} cm^{-2}
  SET 4 - Results at 10x10^{10} cm^{-2}

$outfile_Ef_N contains the quasi-Fermi energy for the lowest subband as a
function of electron density:

  COLUMN 1 - Electron concentration [x10^{10} cm^{-2}]
  COLUMN 2 - Quasi-Fermi energy (relative to band edge) [meV]
  
This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f Ee.r Ef.r FD1N* N.r v.r wf*
