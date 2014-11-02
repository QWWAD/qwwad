#! /bin/sh
set -e

outfile=E-donor-3D-vs-2D.dat
rm -f $outfile

# 3D trial wave function donor binding energy calculation

# Define structure
cat > s.r << EOF
200 0.1 0.0
60  0.0 0.0
200 0.1 0.0
EOF

# Define donor positions
echo '0.0e-10
2.0e-9
4.0e-9
6.0e-9
8.0e-9
1.0e-8
1.2e-8
1.4e-8 
1.6e-8 
1.8e-8 
2.0e-8 
2.2e-8
2.3e-8' > r_d.r

# Generate alloy profile
find_heterostructure

# Generate potential profile
efxv --material cdmnte --mass 0.096

# Perform 2D donor calculation
d02D --mass 0.096 --epsilon 10.6 --lambdastart 25 --lambdastop 300 --symmetry 2D > garbage.r
mv e.r e-2D.r
d02D --mass 0.096 --epsilon 10.6 --lambdastart 25 --lambdastop 300 --symmetry 3D > garbage.r

paste e-2D.r e.r | awk '{print $1, $2 - $4}' > $outfile
