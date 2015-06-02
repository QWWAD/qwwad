#! /bin/sh
set -e

outfileEF=E-F.dat
outfileLO=LO.dat
outfilecc=cc.dat
outfilepop=pop-ratio.dat
rm -f $outfile $outfileLO $outfilecc $outfilepop

cat > s.r << EOF
100   0.2 0
56.6  0.0 0
56.5  0.2 0
96.1  0.0 0
28.25 0.2 0
84.8  0.0 0
100   0.2 0
EOF

find_heterostructure
efxv
mv v.r vcb.r

Tl=77
Te=77

# Define required scattering rates
cat > rrp.r << EOF
2 1
3 2
EOF

# Carrier-carrier rates
cat > rr.r << EOF
3 3 2 3
3 3 2 2
3 3 2 1
3 2 2 2
3 2 2 1
3 1 2 3
3 1 2 2
3 1 2 1
2 3 1 3
2 3 1 2
2 3 1 1
2 2 1 3
2 2 1 2
2 2 1 1
2 1 1 3
2 1 1 1
EOF

# Define subband populations
cat > N.r << EOF
1e14
1e14
1e14
EOF

for field in `seq 0 20`; do
    find_poisson_potential --centred --field $field --uncharged --Vbasefile vcb.r
    efss --nst-max 3

    E1=`cut -f2 Ee.r | sed -n "1p"`
    E2=`cut -f2 Ee.r | sed -n "2p"`
    E3=`cut -f2 Ee.r | sed -n "3p"`

    # Generate Fermi distribution for each subband
    sbp --Te $Te

    # Compute LO phonon rates
    srelo --Te $Te --Tl $Tl
    LO21=`cut -f3 LOe-if.r | sed -n "1p"`
    LO32=`cut -f3 LOe-if.r | sed -n "2p"`

    # Compute carrier-carrier rates
    srcc  --temperature $Te
    cp ccABCD.r ccABCD_${field}.r
    cc3323=`cut -f5 -d' ' ccABCD.r | sed -n "1p"`
    cc3322=`cut -f5 -d' ' ccABCD.r | sed -n "2p"`
    cc3321=`cut -f5 -d' ' ccABCD.r | sed -n "3p"`
    cc3222=`cut -f5 -d' ' ccABCD.r | sed -n "4p"`
    cc3221=`cut -f5 -d' ' ccABCD.r | sed -n "5p"`
    cc3123=`cut -f5 -d' ' ccABCD.r | sed -n "6p"`
    cc3122=`cut -f5 -d' ' ccABCD.r | sed -n "7p"`
    cc3121=`cut -f5 -d' ' ccABCD.r | sed -n "8p"`
    cc2313=`cut -f5 -d' ' ccABCD.r | sed -n "9p"`
    cc2312=`cut -f5 -d' ' ccABCD.r | sed -n "10p"`
    cc2311=`cut -f5 -d' ' ccABCD.r | sed -n "11p"`
    cc2213=`cut -f5 -d' ' ccABCD.r | sed -n "12p"`
    cc2212=`cut -f5 -d' ' ccABCD.r | sed -n "13p"`
    cc2211=`cut -f5 -d' ' ccABCD.r | sed -n "14p"`
    cc2113=`cut -f5 -d' ' ccABCD.r | sed -n "15p"`
    cc2111=`cut -f5 -d' ' ccABCD.r | sed -n "16p"`

    CC32=`echo $cc3323 $cc3322 $cc3321 $cc3222 $cc3221 $cc3123 $cc3122 $cc3121 | awk '{print $1+2*$2+$3+$4+$5+$6+$7+$8}'`
    CC21=`echo $cc2313 $cc2312 $cc2311 $cc2213 $cc2212 $cc2211 $cc2113 $cc2111 | awk '{print $1+$2+$3+$4+$5+2*$6+$7+$8}'`

    echo $field $E1 $E2 $E3 >> $outfile
    echo $field $LO21 $LO32 >> $outfileLO
    echo $field $CC21 $CC32 >> $outfilecc
done

# Now find lifetime ratio
cat $outfileLO $outfilecc | awk '{print $1, ($2 + $5)/($3 + $6)}' > $outfilepop

# Finally, regenerate wavefunctions for a pretty picture
find_poisson_potential --centred --field 10 --uncharged --Vbasefile vcb.r
efss --nst-max 3
wfplot --plot-wf

rm *.r
