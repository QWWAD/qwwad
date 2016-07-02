#! /bin/sh
set -e

# Initialise files

rm -f EX0-lw.r

# Loop for different well widths
export QWWAD_WELLWIDTH

for QWWAD_WELLWIDTH in 0.1 0.5 1 2 5 10 20 30 40 50 60 70 80 90 100 120 140 160 180 200 300 400 500 600 800 1000;
do
	# Find ground states for electron and hole
	qwwad_ef_infinite_well --mass 0.096 --particle e # infinite well electron ground state
	qwwad_ef_infinite_well --mass 0.6   --particle h # infinite well hole ground state

	qwwad_ef_exciton --dcpermittivity 10.6 --electronmass 0.096 --holemass 0.6 --lambdastart 30 	   # start variational calculation
	printf "%f\t" $QWWAD_WELLWIDTH >> EX0-lw.r # write well width to file
	cat EX0.r >> EX0-lw.r			   # send data to file
done
