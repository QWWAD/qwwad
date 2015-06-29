#!/bin/sh
set -e

# Calculates quasi-Fermi distribution for first 3 subbands in a quantum well
# at a range of different temperatures, with and without nonparabolicity
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2015, ch.2
#
# (c) Copyright 1996-2015
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


# Find the Fermi-Dirac functions and subband energies for a given
# non-parabolicity value
#
# parameter 1: The nonparabolicty [1/eV]
solve_for_alpha()
{
    # Solve single quantum well
    qwwad_ef_infinite_well --alpha $1

    # Calculate Fermi energies and population distributions
    qwwad_fermi_distribution --fd --Te $T --alpha $1
    mv FD1.r FD1T=$T-np-alpha$1.r
    mv FD2.r FD2T=$T-np-alpha$1.r
    mv FD3.r FD3T=$T-np-alpha$1.r

    # Write Fermi energy of states to file
    printf "%d\t" $T >> Ef1-T-alpha$1.r
    printf "%d\t" $T >> Ef2-T-alpha$1.r
    printf "%d\t" $T >> Ef3-T-alpha$1.r
    awk '/^1/{printf("%8.3f\n",$2)}' Ef.r >> Ef1-T-alpha$1.r
    awk '/^2/{printf("%8.3f\n",$2)}' Ef.r >> Ef2-T-alpha$1.r
    awk '/^3/{printf("%8.3f\n",$2)}' Ef.r >> Ef3-T-alpha$1.r
}

# Initialise files
outfile_T=fermi-distribution-T-sweep.dat
outfile_np=fermi-energy-T-sweep.dat
rm -f Ef?-T*.r N.r FD*.r $outfile_np $outfile_T

# set fixed parameters
export QWWAD_WELLWIDTH=200 # Well width
export QWWAD_NST=3  # Number of states

# Generate table of populations
N=1e14
cat > N.r << EOF
1 $N
2 $N
3 $N
EOF

# Loop for different temperatures
for T in 2 20 40 60 77 100 140 180 220 260 300
do
    # Find Fermi distribution at this temperature using various degrees
    # of nonparabolicity
    solve_for_alpha 0
    solve_for_alpha 0.7
    solve_for_alpha 5
done

# Dump the variation in Fermi energy vs. temperature to output file
cat Ef1-T-alpha0.r >> $outfile_np
printf "\n" >> $outfile_np
cat Ef2-T-alpha0.r >> $outfile_np
printf "\n" >> $outfile_np
cat Ef3-T-alpha0.r >> $outfile_np
printf "\n" >> $outfile_np
cat Ef1-T-alpha0.7.r >> $outfile_np
printf "\n" >> $outfile_np
cat Ef2-T-alpha0.7.r >> $outfile_np
printf "\n" >> $outfile_np
cat Ef3-T-alpha0.7.r >> $outfile_np
printf "\n" >> $outfile_np
cat Ef1-T-alpha5.r >> $outfile_np
printf "\n" >> $outfile_np
cat Ef2-T-alpha5.r >> $outfile_np
printf "\n" >> $outfile_np
cat Ef3-T-alpha5.r >> $outfile_np

# Dump the Fermi distribution to output file
cat FD1T\=2-np-alpha0.r >> $outfile_T
printf "\n" >> $outfile_T
cat FD2T\=2-np-alpha0.r >> $outfile_T
printf "\n" >> $outfile_T
cat FD3T\=2-np-alpha0.r >> $outfile_T
printf "\n" >> $outfile_T
cat FD1T\=77-np-alpha0.r >> $outfile_T
printf "\n" >> $outfile_T
cat FD2T\=77-np-alpha0.r >> $outfile_T
printf "\n" >> $outfile_T
cat FD3T\=77-np-alpha0.r >> $outfile_T
printf "\n" >> $outfile_T
cat FD1T\=300-np-alpha0.r >> $outfile_T
printf "\n" >> $outfile_T
cat FD2T\=300-np-alpha0.r >> $outfile_T
printf "\n" >> $outfile_T
cat FD3T\=300-np-alpha0.r >> $outfile_T
printf "\n" >> $outfile_T

cat << EOF
Results have been written to $outfile_np and $outfile_T.

$outfile_T contains the Fermi occupation number for states as a function of
energy in the format:

  COLUMN 1 - Energy (relative to band edge) [meV]
  COLUMN 2 - Fermi occupation

  The file contains nine data sets for various temperatures and states, each
  set being separated by a blank line:

  SET 1 - Results for subband 1 at 2 K
  SET 2 - Results for subband 2 at 2 K
  SET 3 - Results for subband 3 at 2 K
  SET 4 - Results for subband 1 at 77 K
  SET 5 - Results for subband 2 at 77 K
  SET 6 - Results for subband 3 at 77 K
  SET 7 - Results for subband 1 at 300 K
  SET 8 - Results for subband 2 at 300 K
  SET 9 - Results for subband 3 at 300 K

$outfile_np contains the quasi-Fermi energy for each subband as a function of
temperature, both with and without nonparabolicity effects:

  COLUMN 1 - Temperature [K]
  COLUMN 2 - Quasi-Fermi energy (relative to band edge) [meV]
  
  The file contains nine data sets, which are separated by a blank line.

  SET 1 - Results for subband 1 (parabolic dispersion)
  SET 2 - Results for subband 2 (parabolic dispersion)
  SET 3 - Results for subband 3 (parabolic dispersion)
  SET 4 - Results for subband 1 (alpha = 0.7 eV^{-1})
  SET 5 - Results for subband 2 (alpha = 0.7 eV^{-1})
  SET 6 - Results for subband 3 (alpha = 0.7 eV^{-1})
  SET 7 - Results for subband 1 (alpha = 5 eV^{-1})
  SET 8 - Results for subband 2 (alpha = 5 eV^{-1})
  SET 9 - Results for subband 3 (alpha = 5 eV^{-1})

This script is part of the QWWAD software suite.

(c) Copyright 1996-2015
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f FD* Ef* Ee.r wf_* N.r v.r
