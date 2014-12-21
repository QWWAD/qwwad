#! /usr/bin/gnuplot

# Plot each of the charge density plots
set term postscript eps enhanced size 3.3,3.3

unset xtics
unset ytics
unset colorbox
unset key
set palette grey

set output 'infinite-rectangular-wire-wf11.eps'
plot 'cd11.r' using 1:2:($3*$3) with image
unset output

set output 'infinite-rectangular-wire-wf12.eps'
plot 'cd12.r' using 1:2:($3*$3) with image
unset output

set output 'infinite-rectangular-wire-wf21.eps'
plot 'cd21.r' using 1:2:($3*$3) with image
unset output

set output 'infinite-rectangular-wire-wf22.eps'
plot 'cd22.r' using 1:2:($3*$3) with image
unset output
