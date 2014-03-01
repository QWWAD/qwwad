/*==================================================================
              efwf   Envelope Function Wave Functions
  ==================================================================*/

/* This program uses a shooting technique to calculate the
   uncorrelated one particle energies of any user supplied
   potential.  The potential is read from the file v.r

   Paul Harrison, July 1992                  

   The program has been updated to be run entirely from the command 
   line and stripped down to calculate the energies only.  It now also
   includes support for non-parabolicity.

   This program is an adapted version of efshoot to calculate the 
   wavefunctions which are read in from the user specified file E?.r

   Paul Harrison, December 1996                                   */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <malloc.h>
#include <gsl/gsl_math.h>
#include "ef-helpers.h"
#include "struct.h"
#include "qclsim-constants.h"
#include "maths.h"

int main(int argc,char *argv[])
{
double	read_delta_z();
double	wf();		/* calculates wavefunctions		*/
files	*read_data();	/* reads potential file into memory	*/

double	delta_z;	/* z separation of input potentials	*/
double	E;		/* electron and hole energies		*/
double	N;		/* normalization integral		*/
int	i_state;	/* electron/hole state index		*/
int	i;		/* index				*/
int	n;		/* length of potential file		*/
char	p;		/* particle				*/
char	filename[9];	/* output filename			*/
bool	T_flag;		/* Hamiltonian flag def.=1=>D(1/m)D	*/
bool	np_flag;	/* non-parabolicity flag		*/
FILE	*FE;		/* outputfile for el. energy states	*/
FILE	*Fwf;		/* outputfile for wavefunctions		*/
files	*data_start;	/* start address of potential		*/
data11	*data_m0Eg;	/* start address of m(0) and Eg		*/
data11	*data_zwf;	/* start address of z and wf		*/

/* default values */

np_flag=false;
T_flag=true;
p='e';

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
	   np_flag=true;
	   argv--;
	   argc++;
	   break;
  case 'k':
	   T_flag=false;
	   argv--;
	   argc++;
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  efshoot [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
	   }
	   break;
  default :
	   printf("Usage:  efwf [-a include non-parabolicity \033[1mfalse\033[0m]\n");
	   printf("             [-k use alternative KE operator (1/m)PP \033[1mP(1/m)P\033[0m]\n");
	   printf("             [-p particle (\033[1me\033[0m, h, or l)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}


data_start=read_data(&n);			/* reads potential file	*/
if(np_flag)data_m0Eg=read_Egdata(n,data_start);	/* reads bandgap data	*/	

delta_z=read_delta_z(data_start);

/* allocate memory for wavefunction */

data_zwf=(data11 *)calloc(n,sizeof(data11));
if(data_zwf==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* copy z data into wavefunction structure	*/

for(i=0;i<n;i++)
 (data_zwf+i)->a=(data_start+i)->z;

/* open appropriate energy file and read in energies sequentially	*/

sprintf(filename,"E%c.r",p);
FE=fopen(filename,"r");

while(fscanf(FE,"%i %le",&i_state,&E)!=EOF)
{
 E*=1e-3*e;	/* convert E from meV->J	*/
 N=wf(E,delta_z,data_start,data_m0Eg,data_zwf,n,np_flag,T_flag);

 sprintf(filename,"wf_%c%i.r",p,i_state);
 Fwf=fopen(filename,"w");

 for(i=0;i<n;i++)
  fprintf(Fwf,"%20.17le %20.17le\n",(data_zwf+i)->a,((data_zwf+i)->b)/sqrt(N));

 fclose(Fwf);

}/* end while */

fclose(FE);
free(data_start);
if(np_flag)free(data_m0Eg);
free(data_zwf);

return EXIT_SUCCESS;
} /* end main */




double
read_delta_z(fdata)

/* This function opens the external file v.r and calculates
   the separation along the z (growth) direction of the
   user supplied potentials                                        */

files *fdata;
{
 double z[2];           /* displacement along growth direction     */

 z[0]=fdata->z;
 fdata++;
 z[1]=fdata->z;
 return(z[1]-z[0]);
}



files
*read_data(n)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */

int	*n;

{
 FILE 	*Fv;            /* file pointer to potential file          */
 FILE 	*Fm;            /* file pointer to potential file          */
 files  *fdata;		/* temporary pointer to potential	   */
 files  *data_start;	/* start address of potential		   */

 if((Fm=fopen("m.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'm.r'!\n");exit(0);}

 if((Fv=fopen("v.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'v.r'!\n");exit(0);}

 *n=0;
 while(fscanf(Fv,"%*e %*e")!=EOF)
  (*n)++;
 rewind(Fv);

 data_start=(files *)calloc(*n,sizeof(files));
 if(data_start==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 fdata=data_start;

 while(fscanf(Fv,"%le %le",&(fdata->z),&(fdata->V))!=EOF)
 {
  int n_read = fscanf(Fm,"%*e %le",&(fdata->mstar));

  if (n_read == 2)
    fdata++;
 }

 fclose(Fm);
 fclose(Fv);

 return(data_start);

}

double
wf(E,delta_z,fdata,data_m0Eg,data_zwf,n,np_flag,T_flag)     

/* This function returns the value of the wavefunction (psi)
   at +infinity for a given value of the energy.  The solution
   to the energy occurs for psi(+infinity)=0.                      */

double E;
double delta_z;
files  *fdata;
data11 *data_m0Eg;
data11 *data_zwf;
int    n;
bool   np_flag;
bool   T_flag;
{
 double alpha;		     /* non-parabolicity parameter   */
 double N=0;		     /* normalization integral       */
 double psi[3];              /* wavefunction at z-delta_z,
                                z and z+delta_z              */
 int	i;		     /* index			     */

 /* account for non-parabolicity with mstar(E)	*/

 if(np_flag)
 {
  for(i=0;i<n;i++)
  {
   alpha=gsl_pow_2(1-((data_m0Eg+i)->a)/me)/((data_m0Eg+i)->b);
   (fdata+i)->mstar=((data_m0Eg+i)->a)*(1+alpha*(E-((fdata+i)->V)));
  }
 }

 /* boundary conditions */
 
 psi[0]=0.0;                 
 psi[1]=1.0;

 (data_zwf)->b=psi[0];
 (data_zwf+1)->b=psi[1];

 N+=gsl_pow_2(psi[0]);
 N+=gsl_pow_2(psi[1]);

 fdata++;                    /* ignore data corresponding to psi[0] */

 if (T_flag)
 {
  for(i=1;i<(n-1);i++)              /* last potential not used */
  {
   psi[2]=(
           (2*gsl_pow_2(delta_z/hBar)*(fdata->V-E)+
	    2/(fdata->mstar+(fdata+1)->mstar)+
	    2/(fdata->mstar+(fdata-1)->mstar))*psi[1]
           -2/(fdata->mstar+(fdata-1)->mstar)*psi[0]
          )
           *(fdata->mstar+(fdata+1)->mstar)/2;
   (data_zwf+i+1)->b=psi[2];
   N+=gsl_pow_2(psi[2]);
   psi[0]=psi[1];
   psi[1]=psi[2];
   fdata++;
  } 
 }
 else
 {
  for(i=1;i<(n-1);i++)              /* last potential not used */
  {
   psi[2]=(2*fdata->mstar*(fdata->V-E)*gsl_pow_2(delta_z/hBar)+2)*psi[1]
          -psi[0];
   (data_zwf+i+1)->b=psi[2];
   N+=gsl_pow_2(psi[2]);
   psi[0]=psi[1];
   psi[1]=psi[2];
   fdata++;
  }
 }

 return(N*delta_z);
}



double 
V_min(fdata,n)       

/* This function opens the external file v.r and finds     
   the minimum value for the potential energy, this value
   is used as the initial energy estimate.                         */

files *fdata;		/* pointer to potential			   */
int   n;		/* number of steps in potential		   */

{
 double min;            /* minimum value of potential energy       */
 int  i;                /* index                                   */
 
 min=1;

 for(i=0;i<n;i++)
 {
  if(fdata->V<min)
  {
   min=fdata->V;
  }
  fdata++;
 }
 return(min);
}

