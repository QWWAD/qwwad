/*=================================================================
       cszb     Crystal Structure Zinc Blende 
  =================================================================

   This program generates the atomic positions of a zinc blende 
   crystal and writes them in XYZ format to the file atoms.xyz

   Paul Harrison, July 1998		                         */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

typedef struct
{
 double x;
 double y;
 double z;
}vector;

main(int argc,char *argv[])
{
int	*read_as();	/* function to read in atomic species		*/
void	write_ap();	/* writes the atomic positions/species to file	*/

double	A0;		/* the lattice constant				*/
int	n_x;		/* number of lattice points along x-axis of cell*/
int	n_y;		/* number of lattice points along y-axis of cell*/
int	n_z;		/* number of lattice points along z-axis of cell*/
char	cation[12];	/* cation species				*/
char	anion[12];	/* anion species				*/
vector	a[4];		/* lattice vectors (a1, a2, a3 plus null vector)*/
vector	T[2];		/* basis vector					*/

/* default values	*/

n_x=1;n_y=1;n_z=1;
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
  case 'x':
           n_x=atoi(argv[2]);
           break;
  case 'y':
           n_y=atoi(argv[2]);
           break;
  case 'z':
           n_z=atoi(argv[2]);
           break;
  default :
	   printf("Usage:  cszb [-a anion \033[1mGA\033[0m][-c cation \033[1mAS\033[0m]\n");
	   printf("             [-x # cells along x-axis \033[1m1\033[0m][-y # cells \033[1m1\033[0m][-z # cells \033[1m1\033[0m]\n");
	   printf("             [-A lattice constant (\033[1m5.65\033[0mA)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

a[0].x=0;a[0].y=0;a[0].z=0;          /* in order to construct a */
a[1].x=1.0/2;a[1].y=1.0/2;a[1].z=0;    /* cubic crystal, this alg-*/
a[2].x=0;a[2].y=1.0/2;a[2].z=1.0/2;    /* orithm requires null ve-*/
a[3].x=1.0/2;a[3].y=0;a[3].z=1.0/2;    /* ctor-called a[0]        */

T[0].x=-1.0/8;T[0].y=-1.0/8;T[0].z=-1.0/8;
T[1].x=+1.0/8;T[1].y=+1.0/8;T[1].z=+1.0/8;

write_ap(A0,n_x,n_y,n_z,a,T,anion,cation);

}/* end main */





void 
write_ap(A0,n_x,n_y,n_z,a,T,anion,cation)

double	A0;	/* lattice constant				*/
int	n_x;	/* number of lattice points along x-axis of cell */
int	n_y;	/* number of lattice points along y-axis of cell */
int	n_z;	/* number of lattice points along z-axis of cell */
vector a[];     /* lattice vectors (a1, a2, a3 plus null vector) */
vector T[];     /* basis vector   */
char	anion[4];	/* anion species				*/
char	cation[4];	/* cation species				*/
{
 int    i_a;    /* index across primitive vectors `a' */
 int    i_n_x;  /* index to n_x */
 int    i_n_y;  /* index to n_y */
 int    i_n_z;  /* index to n_z */
 int    i_T;    /* index across basis vectors `+T' and `-T' */
 vector t;      /* general vertor representing atom within cell  */
 FILE	*Fap;	/* pointer to output file	*/

 Fap=fopen("atoms.xyz","w");

 /* Write number of atoms to first line of file, then leave blank line	*/

 fprintf(Fap,"%i\n\n",8*n_x*n_y*n_z);

 for(i_n_x=0;i_n_x<n_x;i_n_x++)
  for(i_n_y=0;i_n_y<n_y;i_n_y++)
   for(i_n_z=0;i_n_z<n_z;i_n_z++)
    for(i_a=0;i_a<=3;i_a++)
     {
      t.x=i_n_x+a[i_a].x+T[0].x;
      t.y=i_n_y+a[i_a].y+T[0].y;
      t.z=i_n_z+a[i_a].z+T[0].z;
      fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",cation,t.x*A0,t.y*A0,t.z*A0);
      t.x=i_n_x+a[i_a].x+T[1].x;
      t.y=i_n_y+a[i_a].y+T[1].y;
      t.z=i_n_z+a[i_a].z+T[1].z;
      fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",anion,t.x*A0,t.y*A0,t.z*A0);
     }

 fclose(Fap);
}

