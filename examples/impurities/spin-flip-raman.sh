#! /bin/sh
set -e

# Spin-flip Raman scattering calculation
#
# This script is part of the QWWAD software suite. Any use of this code
# or its derivatives in published work must be accompanied by a citation
# of:
#   P. Harrison and A. Valavanis, Quantum Wells, Wires and Dots, 4th ed.
#    Chichester, U.K.: J. Wiley, 2016, ch.5
#
# (c) Copyright 1996-2016
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
outfile_E_pm="spin-flip-E.dat"
outfile_dE="spin-flip-dE.dat"
rm -f $outfile_E_pm $outfile_dE

# Define structure
cat > s.r << EOF
200 0.15 0.0
60  0.0  0.0
200 0.15 0.0
EOF

# Set some fixed material system parameters
export QWWAD_MATERIAL=cdmnte
export QWWAD_MASS=0.096
export QWWAD_DCPERMITTIVITY=10.6

export QWWAD_MAGNETICFIELD=8

# Generate alloy profile and band-edge parameters
qwwad_mesh
qwwad_ef_band_edge

# Find Zeeman splitting to conduction-band due to magnetic field for both spin-up and spin-down cases
qwwad_ef_zeeman --spinup --totalpotentialfile v_up.r
qwwad_ef_zeeman --totalpotentialfile v_down.r

export QWWAD_DONORPOSITION
for QWWAD_DONORPOSITION in `seq 0 10 230`; do

	# Find impurity states for spin-up and spin-down cases
	qwwad_ef_donor_specific --symmetry 3D --lambdastart 10 --lambdastop 1000 --totalpotentialfile v_up.r

	# Save all data for the spin-up state
	Eplus=`awk '{print $2}' Ee.r`

	# Now repeat for the spin-down state
	qwwad_ef_donor_specific --symmetry 3D  --lambdastart 10 --lambdastop 1000 --totalpotentialfile v_down.r
	Eminus=`awk '{print $2}' Ee.r`

	printf "%d\t%e\t%e\n" $QWWAD_DONORPOSITION $Eplus $Eminus >> $outfile_E_pm
done

# Now produce the energy difference between the states
awk '{print $1, $2-$3}' $outfile_E_pm > $outfile_dE

# With only a few donor points, calculating the spin-flip spectra produces
# a very spiky Intensity-energy curve.  So take a spline of the spin-flip
# energies in e_sf.r-raw to simulate a continuous donor distribution and
# save in file `e_sf.r'.

export QWWAD_WAVENUMBERMIN=13
export QWWAD_WAVENUMBERMAX=23
export QWWAD_WAVENUMBERSTEP=0.1
export QWWAD_SPINFLIPFILE=$outfile_dE

qwwad_spin_flip_raman -l 0.5
mv I.r Il-0.5.dat

qwwad_spin_flip_raman
mv I.r Il-1.0.dat

qwwad_spin_flip_raman -l 2
mv I.r Il-2.0.dat

# <Edit accordingly, to describe each output file>
cat << EOF
Results have been written to $outfile in the format:

  COLUMN 1 - <Whatever>
  COLUMN 2 - <Something>
  <...>

  <Remove this chunk if only 1 data set is in the file>
  The file contains <x> data sets, each set being separated
  by a blank line, representing <whatever>:

  SET 1 - <Description of the 1st set>
  SET 2 - <Description of the 2nd set>
  <...>

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# <Delete all temporary files you created>
rm -f *.r
