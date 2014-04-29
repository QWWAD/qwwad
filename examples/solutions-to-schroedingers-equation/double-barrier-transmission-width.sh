#!/bin/sh
set -e

# Calculates the transmission coefficient as a function of the energy
# through a double barrier of width 100 angstrom and height 100 meV.
# Data is generated three different well widths: 20, 50 and 100 angstrom.
#
# The output file contains three data sets (for each well width in
# ascending order). Each data set is separated by a blank line and contains
# two columns of data:
#   Column 1: Energy [meV]
#   Column 2: Transmission coefficient
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
outfile=double-barrier-transmission-width.r
rm -f $outfile

# Loop for well widths, use default parameters otherwise
for L2 in 20 50 100
do
{
 tdb -b $L2 
 #mv T.r TL2=$L2.r
 cat T.r >> $outfile
 printf "\n" >> $outfile
} done

# Clean up the workspace
rm T.r
