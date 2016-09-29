#! /bin/sh
set -e

# Compute electric-field effect in structure with variable permittivity
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.3
#
# (c) Copyright 1996-2016
#     Alex Valavanis <a.valavanis@leeds.ac.uk>
#     Paul Harrison  <p.harrison@leeds.ac.uk>
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
outfile=field-effect-variable-epsilon.dat
rm -f $outfile

cat > s.r << EOF
200 1 0
100 0 0
200 1 0
EOF

qwwad_mesh --dzmax 1
qwwad_ef_band_edge

# Create a permittivity profile using 
# 5 eps0 in barrier; 10 eps0 in well
awk '{
eps0=8.854e-12;

if ($2 > 0.5)
    print $1, 5*eps0;
else
    print $1, 10*eps0;
}' < x.r > eps_dc.r

qwwad_poisson --centred --uncharged --field 10 --poissonpotentialfile v_10_5.r

# 10 eps0 in barrier; 10 eps0 in well
awk '{
eps0=8.854e-12;

if ($2 > 0.5)
    print $1, 10*eps0;
else
    print $1, 10*eps0;
}' < x.r > eps_dc.r

qwwad_poisson --centred --uncharged --field 10 --poissonpotentialfile v_10_10.r

# 15 eps0 in barrier; 10 eps0 in well
awk '{
eps0=8.854e-12;

if ($2 > 0.5)
    print $1, 15*eps0;
else
    print $1, 10*eps0;
}' < x.r > eps_dc.r

qwwad_poisson --centred --uncharged --field 10 --poissonpotentialfile v_10_15.r

awk '{print $1*1e10, $2*1000/1.6021766e-19}' v_10_5.r >> $outfile
printf "\n" >> $outfile
awk '{print $1*1e10, $2*1000/1.6021766e-19}' v_10_10.r >> $outfile
printf "\n" >> $outfile
awk '{print $1*1e10, $2*1000/1.6021766e-19}' v_10_15.r >> $outfile
printf "\n" >> $outfile

cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - Spatial location [angstrom]
  COLUMN 2 - Potential [meV]

  The file contains 3 data sets, each set being separated
  by a blank line, representing the solution with different
  permittivity in the barrier layers:

  SET 1 - Barrier permittivity = 5 eps0
  SET 2 - Barrier permittivity = 10 eps0
  SET 3 - Barrier permittivity = 15 eps0

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
rm -f *.r
