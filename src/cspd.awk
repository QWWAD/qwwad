#=====================================================================
#   cspd     Crystal Structure Pyramidal Dot
#=====================================================================
#
# This awk script places a pyramidal GE dot at the centre of a unit cell.
# It would be executed with a command of the form
#
# nawk -v BARRIER=1 -v BASE=5 -v HEIGHT=5 -v A0=5.65 -f cspd.awk atoms.xyz > qd.atoms.xyz
#
# The side of the base is represented by `BASE', the height of the apex by
# `HEIGHT' and the thickness of the surrounding host material by `BARRIER', 
# all in lattice constants `A0' which is in Angstroms

BEGIN{A=(BARRIER-0.125)*A0;B=BASE*A0;H=HEIGHT*A0
	print 8*(2*BARRIER+BASE)*(2*BARRIER+BASE)*(2*BARRIER+HEIGHT)}
	
	{X=$2;Y=$3;Z=$4
	 if((X>A)&&(X<A+B)&&(Y>A)&&(Y<A+B)&&(Z>A)&&(Z<2*H/B*(X-A)+A)&&(Z<-2*H/B*(X-(A+B/2))+A+H)&&(Z<2*H/B*(Y-A)+A)&&(Z<-2*H/B*(Y-(A+B/2))+A+H))
	 printf("GE %9.3f %9.3f %9.3f\n",$2,$3,$4);
	}

# Note this version only prints the dot itself, not the surrounding lattice
