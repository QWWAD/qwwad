#! /bin/sh
set -e

cat > s.r << EOF
150 0.2 0.0
100 0.0 0.0
150 0.2 0.0
EOF

find_heterostructure --alloy-file x.r
