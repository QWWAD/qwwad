#! /bin/sh
set -e

# Mapping of lambda-energy space for a two-dimensional donor wave function
# in a quantum well

# Generate structure file
cat > s.r << EOF
200 0.1 0.0
60  0.0 0.0
200 0.1 0.0
EOF

# Convert structure file into alloy profile
find_heterostructure

# Generate potential profile
efxv --material cdmnte

# Place donor at centre of well
echo "230e-10" > r_d.r

# Start donor binding energy calculation,
# force output of energy E versus lambda data for lambda=40 to 90A in 2A steps
# in order to illustrate the variational principle
d02D -s 40 -t 2 -u 90 -e 10.6 -m 0.096 > output

# Filter `output' file to give energy versus lambda data
awk '{printf("%e %e\n",$4,$6)}' output > e-lambda.r
