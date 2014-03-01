/*====================================================================
                     srrad Scattering Rate---RADiative
  ====================================================================*/

/* This program calculates the spontaneous radiative lifetime between
   two quantum well subbands.  The user is able to choose between the
   3D and 2D forms (see Smet, J. Appl. Phys. 79 (9305) (1996).

   Input files:       wf_e1.r
                      wf_e2.r    etc.
		      Ee.r             Note standard setup for
					electrons only!

   Output files:    

   Paul Harrison

                                                                */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_math.h>
#include "struct.h"
#include "maths.h"
#include "qclsim-constants.h"

int main(int argc, char *argv[])
{
double	E[9];		/* Energy levels			*/
double	gamma;		/* constant, see Smet			*/
double	lambda;		/* wavelength of radiation		*/
double	M();		/* Momentum matrix element <i|d/dz|f>	*/
double	mstar;		/* effective mass			*/
double	ri();		/* refractive index			*/
double	OS;		/* Oscillator Strength			*/
double	omega;		/* angular velocity of radiation	*/
double	rad;		/* Radiative lifetime, 
				Berger SST9, 1493 (1994)	*/
double Wmz;		/* cavity length			*/
int	N_1;		/* number of lines in data files	*/
int	N_2;		/* number of lines in data files	*/
int	state_i;	/* Initial state			*/
int	state_f;	/* Final state 				*/
void	read_E();	/* Reads energy levels			*/
char	filename_1[9];	/* character string for wavefunction 
                   		input file			*/
char	filename_2[9];	/* character string for wavefunction 
                 		input file			*/
char	p;		/* particle (e, h, or l)		*/
bool	D3_flag;	/* if set, use 3D form of lifetime	*/
data11	*read_data();
data11	*start_wf1;	/* start address of data		*/
data11	*start_wf2;	/* start address of data		*/


/* default values */

D3_flag=false;
state_i=2;
state_f=1;
mstar=0.067*me;
p='e';
Wmz=3e-6;		/* default from Smet			*/

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case '3':
	   D3_flag=true;
	   argv--;
	   argc++;
	   break;
  case 'f':
	   state_f=atoi(argv[2]);
	   break;
  case 'i':
	   state_i=atoi(argv[2]);
	   break;
  case 'm':
	   mstar=atof(argv[2])*me;
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  srrad [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
	   }
	   break;
  case 'W':
	   Wmz=atof(argv[2])*1e-6;
	   break;
  default :
           printf("Usage:  srrad [-i initial state \033[1m2\033[0m][-f final state \033[1m1\033[0m][-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("              [-m mass (\033[1m0.067\033[0mm0)][-3 3D value \033[1mfalse\033[0m]\n");
	   printf("              [-W cavity length (\033[1m3\033[0mmicrons)]\n");
	   exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

sprintf(filename_1,"wf_%c%i.r",p,state_i);
sprintf(filename_2,"wf_%c%i.r",p,state_f);

printf("state_i=%i state_f=%i\n",state_i,state_f);

read_E(E);

omega=(E[state_i-1]-E[state_f-1])/hBar;

lambda=2*pi*c/omega;

printf("frequency %le\n",c/lambda);

start_wf1=read_data(filename_1,&N_1);
start_wf2=read_data(filename_2,&N_2);

if(N_1!=N_2)
 {printf("Error: number of lines in %s and %s are not equal!\n",
         filename_1,filename_2);exit(0);
 }

OS=2*hBar/(mstar*omega)*gsl_pow_2(M(start_wf1,start_wf2,N_1));

printf("Oscillator strength %20.17le\n",OS);

if(D3_flag)gamma=ri(lambda)*gsl_pow_2(e)*gsl_pow_2(omega)/(6*pi*eps0*mstar*c*c*c);
else gamma=gsl_pow_2(e)*omega/(4*eps0*Wmz*c*c*mstar);

rad=1/(gamma*OS);

printf("Refractive index at %20.17le meV %20.17le\n",(E[state_i-1]-E[state_f-1])/(1e-3*e),ri(lambda));

printf("Radiative lifetime %20.17le s\n",rad);

free(start_wf1);
free(start_wf2);

return EXIT_SUCCESS;
}        /* end main */



double
M(pointer_wf1,pointer_wf2,N)

/* This function calculates the momentum matrix element <i|d/dz|f> */

data11 *pointer_wf1;       /* pointer to first wavefunction data */
data11 *pointer_wf2;       /* pointer to first wavefunction data */
int     N;                 /* number of wavefucntion points      */
{
 double delta_z;
 double matrix_element=0;
 int    i;                 /* index                              */

 delta_z=((pointer_wf1+1)->a)-(pointer_wf1->a);
 pointer_wf1++;pointer_wf2++;     /* Increment pointers to second z point */
 for(i=1;i<N-1;i++)
 {
  matrix_element+=(pointer_wf1->b)*
                    (((pointer_wf2+1)->b)-((pointer_wf2-1)->b))/(2*delta_z);
  pointer_wf1++;pointer_wf2++;
 }
 matrix_element*=delta_z;

 return(matrix_element);
}



data11
*read_data(filename,N)

/* This function reads the potential and masses into memory and returns
   the start address of this block of memory and the number of lines   */
 
char       filename[];
int        *N;              /* number of lines in data files     */

{
 FILE      *Fwf;            /* pointer to v.r                    */
 data11    *pointer_wf;     /* pointer to the data structure     */
 data11    *start_wf;       /* start address of data             */
 
 if((Fwf=fopen(filename,"r"))==0)
 {
  fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);exit(0);
 }

 *N=0;
 while(fscanf(Fwf,"%*e %*e")!=EOF)
  (*N)++;
 rewind(Fwf);

 start_wf=(data11 *)calloc(*N,sizeof(data11));
 if(start_wf==0)  
 {
  fprintf(stderr,"Cannot allocate memory!\n");exit(0);
 }

 pointer_wf=start_wf;

 while(fscanf(Fwf,"%lf %lf",&(pointer_wf->a),&(pointer_wf->b))!=EOF)
  pointer_wf++;

 fclose(Fwf);
 
 return(start_wf);
}



void
read_E(E)

/* This function reads the energy level file */

double	E[];
{
int	i=0;	/* index					*/
FILE	*FE;	/* file pointer to Ee.r				*/

FE=fopen("Ee.r","r");

while(fscanf(FE,"%*i %le",&E[i])!=EOF)
  {
   E[i]*=1e-3*e;                       /* eV->J */
   i++;
  }

fclose(FE);
}



double 
ri(lambda)

/* This function calculates the refractive index of the semiconductor
   using the expressions of Sellmeier and the data of Seraphin, see Adachi,
   `GaAs and related materials p423                          
                                                                */

double lambda;

{
 double A;     /*  Coefficients of empirical fit */
 double B; 
 double C;

 A=8.720;B=2.053;C=sqrt(0.358);

 /* Note the 1e-6 factor converts lambda into microns */

 return(sqrt(A+B*(gsl_pow_2(lambda/1e-6)/(gsl_pow_2(lambda/1e-6)-gsl_pow_2(C)))));

}
