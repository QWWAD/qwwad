#! /bin/sh
set -e

# Calculates the first 3 energy level in a 100-angstrom-wide quantum
# well and outputs the wavefunctions
#
# Column 1: Position along growth axis (z) [angstrom]
# Column 2: Wavefunction for |1> [m^{-1/2}]
# Column 3: Wavefunction for |1> [m^{-1/2}]
# Column 4: Wavefunction for |1> [m^{-1/2}]
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

# Initialise output file
outfile=infinite-well-wavefunctions.sh
rm -f $outfile

# Solve Schroedinger equation
efiw -L 100 -N 1000 -s 3

# Shift wavefunctions and scale positions to angstrom
awk '{print $1*1e10, $2}' wf_e1.r >> $outfile
printf "\n" >> $outfile
awk '{print $1*1e10, $2+30000}' wf_e2.r >> $outfile
printf "\n" >> $outfile
awk '{print $1*1e10, $2+60000}' wf_e3.r >> $outfile

rm wf_e?.r Ee.r
