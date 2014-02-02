/*=================================================================
                    ppsg PseudoPotential Sort G             
  =================================================================

   This program sorts the reciprocal lattice or G vectors into
   ascending magnitude.

   Input files:
                  G.r       reciprocal lattice vectors

   Output files:

                  G.r~      original reciprocal lattice vectors
                  G.r       sorted reciprocal lattice vectors

   Paul Harrison, October 1996                                */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"struct.h"
#include"maths.h"
#include"const.h"

main()
{
vector	*read_rlv();	/* function to read reciprocal lattice vectors	*/
int	i;		/* index					*/
int	iG;		/* index over G vectors				*/
int	N;		/* number of reciprocal lattice vectors		*/
vector	*G;		/* reciprocal lattice vectors			*/
vector	Gt;		/* temporary storage of G			*/
FILE	*FG;		/* file pointer to wavefunction file		*/

G=read_rlv(&N);

system("cp G.r G.r~");

/* Need to repeat N times, to allow G[0] to propogate to G[N] */
for(i=0;i<N;i++)    
{
 for(iG=0;iG<N-1;iG++)                    
 {
  if(vmod(*(G+iG))>vmod(*(G+iG+1)))
  {
   Gt=*(G+iG+1);
   *(G+iG+1)=*(G+iG);
   *(G+iG)=Gt;
  }
 }
}

FG=fopen("G.r","w");

for(iG=0;iG<N;iG++)                    
 fprintf(FG,"%f %f %f\n",(G+iG)->x,(G+iG)->y,(G+iG)->z);

fclose(FG);

free(G);

}/* end main */


/******************************************************************************/



vector
*read_rlv(N)

/* This function reads the reciprocal lattice vectors (defined in 
   the file G.r) into the array G[] and then converts into SI units */

int	*N;
{
 int	i=0;
 vector	*G;
 FILE 	*FG;           /* file pointer to wavefunction file */

 if((FG=fopen("G.r","r"))==0)
 {
  fprintf(stderr,"Error: Cannot open input file 'G.r'!\n");
  exit(0);
 }

 *N=0;
 while(fscanf(FG,"%*lf %*lf %*lf")!=EOF)
  (*N)++;
 rewind(FG);

 G=(vector *)calloc(*N,sizeof(vector));
 if (G==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FG,"%lf %lf %lf",&(G+i)->x,&(G+i)->y,&(G+i)->z)!=EOF)
  {i++;}
 
 fclose(FG);
 
 return(G);
}
