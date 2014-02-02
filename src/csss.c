/*=================================================================
       csss     Crystal Structure Single Spiral
  =================================================================

   This program generates the atomic positions of a single spiral
   along the z-axis of a zinc blende crystal and writes them in 
   XYZ format to the file zb.xyz

   Paul Harrison, July 1998		                         */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include "struct.h"

main(int argc,char *argv[])
{
void	write_ap();	/* writes the atomic positions/species to file	*/

double	A0;		/* the lattice constant				*/
int	n_z;		/* number of lattice points along z-axis of cell*/
char	cation[12];	/* cation species				*/
char	anion[12];	/* anion species				*/
vector	a;		/* lattice vectors 				*/
vector	T[4];		/* basis vectors				*/

/* default values	*/

n_z=1;
A0=5.65;		/* break all the rules and keep in Angstrom	*/
sprintf(cation,"GA");
sprintf(anion,"AS");

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'A':
           A0=atof(argv[2]);
           break;
  case 'a':
           sprintf(anion,argv[2]);
           break;
  case 'c':
           sprintf(cation,argv[2]);
           break;
  case 'z':
           n_z=atoi(argv[2]);
           break;
  default :
	   printf("Usage:  csss [-a anion \033[1mAS\033[0m][-c cation \033[1mGA\033[0m]\n");
	   printf("             [-z # cells \033[1m1\033[0m][-A lattice constant (\033[1m5.65\033[0mA)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

a.x=0;a.y=0;a.z=1.0;          /* in order to construct a */

T[0].x=-1.0/8;T[0].y=-1.0/8;T[0].z=-1.0/8;
T[1].x=+1.0/8;T[1].y=+1.0/8;T[1].z=+1.0/8;
T[2].x=-1.0/8;T[2].y=+3.0/8;T[2].z=+3.0/8;
T[3].x=-3.0/8;T[3].y=+1.0/8;T[3].z=+5.0/8;

write_ap(A0,n_z,a,T,anion,cation);

}/* end main */





void 
write_ap(A0,n_z,a,T,anion,cation)

double	A0;		/* lattice constant				*/
int	n_z;		/* number of lattice points along z-axis of cell */
vector	a;		/* lattice vectors (a1, a2, a3 plus null vector) */
vector	T[];     	/* basis vector   */
char	anion[];	/* anion species				*/
char	cation[];	/* cation species				*/
{
 int    i_n_z;  /* index to n_z */
 int    i_T;    /* index across basis vectors `+T' and `-T' */
 vector t;      /* general vertor representing atom within cell  */
 FILE	*Fap;	/* pointer to output file	*/

 Fap=fopen("atoms.xyz","w");

 /* Write number of atoms to first line of file, then leave blank line	*/

 fprintf(Fap,"%i\n\n",4*n_z);

 for(i_n_z=0;i_n_z<n_z;i_n_z++)
 {
  t.x=T[0].x;
  t.y=T[0].y;
  t.z=i_n_z*a.z+T[0].z;
  fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",cation,t.x*A0,t.y*A0,t.z*A0);
  t.x=T[1].x;
  t.y=T[1].y;
  t.z=i_n_z*a.z+T[1].z;
  fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",anion,t.x*A0,t.y*A0,t.z*A0);
  t.x=T[2].x;
  t.y=T[2].y;
  t.z=i_n_z*a.z+T[2].z;
  fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",cation,t.x*A0,t.y*A0,t.z*A0);
  t.x=T[3].x;
  t.y=T[3].y;
  t.z=i_n_z*a.z+T[3].z;
  fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",anion,t.x*A0,t.y*A0,t.z*A0);
 }

 fclose(Fap);
}

