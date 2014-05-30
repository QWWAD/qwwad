#! /bin/sh
set -e

cat > s.r << EOF
150 0.2 0.0
100 0.0 0.0
150 0.2 0.0
EOF

find_heterostructure

# Create alloy concentration file
efxv

# Find lowest two states
efss --nst-max 2

awk '{print $1*1e10 - 200, $2}' wf_e1.r > shooting-method-wavefunction-1.dat
awk '{print $1*1e10 - 200, $2}' wf_e2.r > shooting-method-wavefunction-2.dat
