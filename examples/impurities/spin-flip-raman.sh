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
export QWWAD_SPINFLIPFILE="spin-flip-dE.dat"
outfile_spectra="spin-flip-spectra.dat"
rm -f $outfile_E_pm $outfile_dE $outfile_spectra

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
for QWWAD_DONORPOSITION in `seq 0 5 230`; do

	# Find impurity states for spin-up and spin-down cases
	qwwad_ef_donor_specific --symmetry 3D --lambdastart 10 --lambdastop 1000 --totalpotentialfile v_up.r

	# Save all data for the spin-up state
	Eplus=`awk '{print $2}' Ee.r`

	# Now repeat for the spin-down state
	qwwad_ef_donor_specific --symmetry 3D  --lambdastart 10 --lambdastop 1000 --totalpotentialfile v_down.r
	Eminus=`awk '{print $2}' Ee.r`

	printf "%d\t%e\t%e\n" $QWWAD_DONORPOSITION $Eplus $Eminus >> E_pm.tmp
done

# Now store the spin-flip energy between the states
# as a function of the donor position
awk '{print $1, $2-$3}' E_pm.tmp > $QWWAD_SPINFLIPFILE

# Run the spectral calculation for a range of linewidths
export QWWAD_WAVENUMBERMIN=13
export QWWAD_WAVENUMBERMAX=23
export QWWAD_WAVENUMBERSTEP=0.1

export QWWAD_LINEWIDTH
for QWWAD_LINEWIDTH in 0.5 1.0 2.0; do
	qwwad_spin_flip_raman

	cat I.r     >> $outfile_spectra
	printf "\n" >> $outfile_spectra
done

# Save the energies and the band potentials to file
awk '{print $1, $2}' E_pm.tmp > $outfile_E_pm
printf "\n"  >> $outfile_E_pm
awk '{print $1, $3}' E_pm.tmp >> $outfile_E_pm
printf "\n"  >> $outfile_E_pm
awk '{print $1*1e10, $2*1000/1.6e-19}' v_up.r   >> $outfile_E_pm
printf "\n"  >> $outfile_E_pm
awk '{print $1*1e10, $2*1000/1.6e-19}' v_down.r >> $outfile_E_pm

cat << EOF
Results have been written to $outfile_E_pm, $QWWAD_SPINFLIPFILE and $outfile_spectra.

$outfile_E_pm is in the format:

  COLUMN 1 - Donor position or spatial coordinate [Angstrom]
  COLUMN 2 - Energy or potential [meV]

  The file contains 4 data sets, each set being separated
  by a blank line, representing the Zeeman-split potentials and states:

  SET 1 - Spin-up state as function of donor position
  SET 2 - Spin-down state as function of donor position
  SET 3 - Zeeman up-state potential as a function of position
  SET 4 - Zeeman down-state potential as a function of position

$QWWAD_SPINFLIPFILE is in the format:

  COLUMN 1 - Donor position [Angstrom]
  COLUMN 2 - Spin-flip energy [meV]

$outfile_spectra is in the forma:

  COLUMN 1 - Raman shift [1/cm]
  COLUMN 2 - Intensity [arb. units]

  The file contains 3 data sets, each set being separated
  by a blank line, representing the spectra with different linewidths:

  SET 1 - Linewidth = 0.5 /cm
  SET 2 - Linewidth = 1.0 /cm
  SET 3 - Linewidth = 2.0 /cm

This script is part of the QWWAD software suite.

(c) Copyright 1996-2016
    Alex Valavanis <a.valavanis@leeds.ac.uk>
    Paul Harrison  <p.harrison@leeds.ac.uk>

Report bugs to https://bugs.launchpad.net/qwwad
EOF

# Clean up workspace
# <Delete all temporary files you created>
rm -f *.r *.tmp
