#!/bin/sh
set -e

# Calculates the first 3 energy level in an infinite GaAs quantum well
# with a range of well widths. The energies are output to the file
# "Ee-lw.r" (assuming parabolic dispersion) and to "Ee-lw-alpha.r" for
# the nonparabolic case.  In both cases, the files are laid out as:
#
# Column 1: Layer width [angstrom]
# Column 2: Energy of 1st state [meV]
# Column 3: Energy of 2nd state [meV]
# Column 4: Energy of 3rd state [meV]
# 
# (c) 1996-2014
#     Paul Harrison <p.harrison@shu.ac.uk>
#     Alex Valavanis <a.valavanis@leeds.ac.uk>

# Initialise files
outfile=infinite-well-energy-vs-width.r
rm -f $outfile

# Set fixed parameters
mass=0.067 # Effective mass relative to a free electron
nst=3      # Number of states

# Loop for different well widths
for i in `seq 1 0.1 2.3`
do
{
 # Generate well-widths exponentially so we get a smooth curve at small
 # widths
 LW=`echo $i | awk '{print 10^$1}'`

 # Calculate first 3 energy levels as a function of well width for GaAs
 efiw -L $LW -m $mass -s $nst

 printf "%f\t" "$LW" >> $outfile	# write well width to file

 energies=`awk '{print $2}' < Ee.r`
 echo $energies >> $outfile
}
done

rm wf_e?.r Ee.r
