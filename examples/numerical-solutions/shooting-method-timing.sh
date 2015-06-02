#! /bin/sh

outfile=shooting-method-timing.dat
rm -f $outfile

for L in `seq 100 1000 10000`; do

cat > s.r << EOF
100 0.2 0.0
$L  0.0 0.0
100 0.2 0.0
EOF

NZ=`echo $L| awk '{print $1+200+1}'`

find_heterostructure --nz-1per $NZ
efxv

# Solve analytically
V=`head -n1 v.r | awk '{print $2/1.6e-19 * 1000}'`
/usr/bin/time -f "%e" efsqw --well-width $L --potential $V 2> T-analytic.log
T_analytic=`cat T-analytic.log`

# Determine the minimum sized search grid
dE=`head -n1 Ee.r | awk '{print $2/2}'`

/usr/bin/time -f "%e" efss --nst-max 1000 --solver shooting --dE $dE 2> T-shooting.log
T_shooting=`cat T-shooting.log`

/usr/bin/time -f "%e" efss --nst-max 1000 --solver matrix 2> T-matrix.log
T_matrix=`cat T-matrix.log`

echo $L $T_analytic $T_shooting $T_matrix >> $outfile
done
