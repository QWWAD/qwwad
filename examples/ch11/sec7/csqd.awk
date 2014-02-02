#=====================================================================
#   csqd     Crystal Structure Quantum Dot
#=====================================================================
#
# This awk script places a cubic GE dot at the centre of a unit cell.
# It would be executed with a command of the form
#
# nawk -v BARRIER=1 -v SIDE=5 -v A0=5.65 -f csqd.awk atoms.xyz > qd.atoms.xyz
#
# The side of the dot is represented by `WELL' and the thickness of the
# surrounding host material by `BARRIER', both in lattice constants `A0'
# which is in Angstroms

BEGIN{A=(BARRIER-0.125)*A0;B=SIDE*A0}
	
	{X=$2;Y=$3;Z=$4
	 if((X>A)&&(X<A+B)&&(Y>A)&&(Y<A+B)&&(Z>A)&&(Z<A+B))
	 printf("GE %9.3f %9.3f %9.3f\n",$2,$3,$4);
         else print $0
	}

