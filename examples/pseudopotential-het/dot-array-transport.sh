#! /bin/sh
set -ve

# Transport through dot array calculation
# Initialise files
rm -f ank*.r Ek*.r

# Define the unit cell extent

NX=2
NY=2
NZ=2

# and the lattice constant in Angstrom
# Use value for Si
A0=5.431

# Create simple cube unit cell of (NX,NY,NZ) lattice constants
cszb -a SI -c SI -x $NX -y $NY -z $NZ -A $A0

# Create quantum dot by substituting atoms explicitly

sed '3,10s/SI/GE/' atoms.xyz > atoms.tmp
mv atoms.tmp atoms.xyz

# Convert also to pdb format for visual check
xyz2pdb atoms

# Specify k-points for calculation, to calculate dispersion curves need more
# points than zone center, note with n_z=2, zone edge at 1/4 (2*pi/A0)
cat > k.r << EOF
0.0 0.0 0.000
0.0 0.0 0.025
0.0 0.0 0.050
0.0 0.0 0.075
0.0 0.0 0.100
0.0 0.0 0.125
0.0 0.0 0.150
0.0 0.0 0.175
0.0 0.0 0.200
0.0 0.0 0.225
0.0 0.0 0.250
EOF

# Generate reciprocal lattice vectors for this simple cube, and sort
rlv-sc -x $NX -y $NY -z $NZ

ppsg

# Implement pseudopotential band structure calculation, use default lattice
# constant, output all valence band eigenfunctions.  Note the number of
# valence band states is twice the number of atoms in the basis

# Calculate all VB and lowest CB state
N=`echo $NX $NY $NZ | awk '{print 2*$1*$2*$3*8}'`
M=`echo $NX $NY $NZ | awk '{print 2*$1*$2*$3*8 + 1}'`
echo $N $M
pplb -n $N -m $M -A $A0

# Collate top VB and bottom CB states

# Create dummy file with single line in to allow use of awk for the division

echo dummy line > dummy.file

for K in 0 1 2 3 4 5 6 7 8 9 10
do
{
 nawk "{printf(\"%f\",$K/40)}" dummy.file >> Ek.r
 nawk '{printf(" %f",$1)}' Ek$K.r >> Ek.r
 echo -n -e "\n" >> Ek.r
}
done

# Tidy up

#rm -f ank*.r
