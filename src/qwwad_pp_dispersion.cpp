/**
 * \file    qwwad_pp_dispersion.cpp
 * \brief   PseudoPotential EigenEnergies
 * \author  Paul Harrison  <p.harrison@shu.ac.uk>
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program gathers together the eigenenergies from the series of files
 *          Ek#.r and gathers them together in a format suitable for plotting the
 *          dispersion curves with xmgr.
 *
 *          Input files:
 *		Ek#.r		expansion coefficients of eigenvectors
 *
 *          Output files:
 *
 *		Ek.r		dispersion curves
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "struct.h"
#include "maths.h"

int main(int argc,char *argv[])
{
double	*Ek;
int	ik;		/* index over k-points				*/
int	in;		/* index over bands				*/
int	Nk;		/* number of k-points 				*/
int	nbands;		/* number of bands				*/
int	n_min;		/* first band					*/
int	n_max;		/* last band					*/
char	filename[9];	/* name of Ek#.r file				*/
FILE	*Fk;		/* pointer to k points file k.r			*/
FILE	*FEkr;		/* results file					*/
FILE	**FEk;		/* pointer to k points file k.r			*/
vector	k;		/* carrier momentum				*/

/* default values	*/

n_min=0;
n_max=3;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'n':
	   n_min=::atoi(argv[2])-1;         /* Note -1=>top VB=4, CB=5 */
           break;
  case 'm':
	   n_max=::atoi(argv[2])-1;         /* Note -1=>top VB=4, CB=5 */
           break;
  default :
           printf("Usage:  ppee [-n first band \033[1m1\033[0m][-m last band \033[1m4\033[0m]\n");
           exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

nbands=n_max-n_min+1;

Fk=fopen("k.r","r");			/* open k.r file for reading	*/

Nk=0;					/* count number of k points	*/
while(fscanf(Fk,"%*f %*f %*f")!=EOF)
  (Nk)++;
rewind(Fk);

Ek=(double *)calloc(Nk,(nbands+1)*sizeof(double));	/* allocate memory	*/
if (Ek==0) 						/* note column allocation*/
 {fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

FEk=(FILE **)calloc(Nk,sizeof(FILE));

for(ik=0;ik<Nk;ik++)
{
 int n_read = fscanf(Fk,"%lf %lf %lf",&k.x,&k.y,&k.z);
 if (n_read != 3)
 {
   fprintf(stderr, "Data missing in k.r");
   exit(EXIT_FAILURE);
 }

 *(Ek+ik*(nbands+1))=vmod(k);
 sprintf(filename,"Ek%i.r",ik);
 *(FEk+ik)=fopen(filename,"r");		/* open each file in turn	*/
 fseek(*(FEk+ik),n_min*11,0);		/* move file pointer on n_min lines */
 for(in=0;in<nbands;in++)
 {
 int n_read = fscanf(*(FEk+ik),"%lf",(Ek+ik*(nbands+1)+in+1));
 if (n_read != 1)
 {
   fprintf(stderr, "Data missing in k.r");
   exit(EXIT_FAILURE);
 }
 }
 fclose(*(FEk+ik));			/* close each individual file	*/

}

FEkr=fopen("Ek.r","w");		/* open file Ek.r for output	*/
for(ik=0;ik<Nk;ik++)
{
 fprintf(FEkr,"+%lf",*(Ek+ik*(nbands+1)));	/* output |k|	*/
 for(in=0;in<nbands;in++)
  fprintf(FEkr," %lf",*(Ek+ik*(nbands+1)+in+1));/* output each E	*/
 fprintf(FEkr,"\n");				/* end of line	*/
}
fclose(FEkr);


fclose(Fk);
free(Ek);

return EXIT_SUCCESS;
}/* end main */

