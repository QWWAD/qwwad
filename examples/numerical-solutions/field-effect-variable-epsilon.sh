#! /bin/sh
set -e

outfile=field-effect-variable-epsilon.dat
rm -f $outfile

cat > s.r << EOF
200 1 0
100 0 0
200 1 0
EOF

find_heterostructure
efxv

mv eps-dc.r eps-orig.r

# 5 eps0 in barrier; 10 eps0 in well
awk '{
eps0=8.854e-12;

if ($2 < 12*eps0)
    print $1, 5*eps0;
else
    print $1, 10*eps0;
}' < eps-orig.r > eps-dc.r

find_poisson_potential --centred --uncharged --field 10 --potential-file v_10_5.r

# 10 eps0 in barrier; 10 eps0 in well
awk '{
eps0=8.854e-12;

if ($2 < 12*eps0)
    print $1, 10*eps0;
else
    print $1, 10*eps0;
}' < eps-orig.r > eps-dc.r

find_poisson_potential --centred --uncharged --field 10 --potential-file v_10_10.r

# 15 eps0 in barrier; 10 eps0 in well
awk '{
eps0=8.854e-12;

if ($2 < 12*eps0)
    print $1, 15*eps0;
else
    print $1, 10*eps0;
}' < eps-orig.r > eps-dc.r

find_poisson_potential --centred --uncharged --field 10 --potential-file v_10_15.r

awk '{print $1*1e10, $2*1000/1.6e-19}' v_10_5.r >> $outfile
printf "\n" >> $outfile
awk '{print $1*1e10, $2*1000/1.6e-19}' v_10_10.r >> $outfile
printf "\n" >> $outfile
awk '{print $1*1e10, $2*1000/1.6e-19}' v_10_15.r >> $outfile
printf "\n" >> $outfile

