#!/bin/sh -e

# Computes exact and approximate specific heat capacity for GaAs and Au
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
outfile=specific-heat.dat
rm -f $outfile

# GaAs (exact)
specheat --debye 360 --molarmass 0.1446 --natoms 2
cat c.r    >> $outfile
printf "\n" >> $outfile

# GaAs (approx)
specheat --debye 360 --molarmass 0.1446 --natoms 2 --approx
cat c.r    >> $outfile
printf "\n" >> $outfile

# Gold (exact)
specheat --debye 165 --molarmass 0.1970 --natoms 1
cat c.r    >> $outfile
printf "\n" >> $outfile

# Gold (exact)
specheat --debye 165 --molarmass 0.1970 --natoms 1 --approx
cat c.r    >> $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Temperature [K]
  COLUMN 2 - Specific heat capacity [J/kg/K]

  The file contains 4 data sets, each set being separated
  by a blank line, representing:

  SET 1  - GaAs (exact)
  SET 2  - GaAs (approximate)
  SET 3  - Au (exact)
  SET 4  - Au (approximate)

This script is part of the QWWAD software suite.

(c) Copyright 1996-2014
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
