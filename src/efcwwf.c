/*==================================================================
              efcwwf   Envelope Function Circular Wire Wave Functions
  ==================================================================*/

/* This program uses a shooting technique to calculate the
   uncorrelated one particle wavefunctions of any user supplied
   radial potential.  The potential is read from the file v.r
   The eigenenergies have been calculated previously with `efcwire'
   and are stored in the file E?.r

   Paul Harrison, December 1998                                   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <malloc.h>
#include "struct.h"
#include "const.h"
#include "maths.h"
#include "bools.h"
#include "io.h"


typedef
struct	{
 double	z;		/* z value of files    		  */
 double	V;		/* electron and hole values   	  */
 double	mstar;		/* electron and hole values   	  */
} files;

main(int argc,char *argv[])
{
double	read_delta_z();
double	wf();		/* calculates wavefunctions		*/
files	*read_data();	/* reads potential file into memory	*/
data11	*read_Egdata();	/* reads bandgap data			*/

double	delta_z;	/* z separation of input potentials	*/
double	E;		/* electron and hole energies		*/
double	N;		/* normalization integral		*/
int	i_state;	/* electron/hole state index		*/
int	i;		/* index				*/
int	n;		/* length of potential file		*/
char	p;		/* particle				*/
char	filename[9];	/* output filename			*/
boolean	np_flag;	/* non-parabolicity flag		*/
FILE	*FE;		/* outputfile for el. energy states	*/
FILE	*Fwf;		/* outputfile for wavefunctions		*/
files	*data_start;	/* start address of potential		*/
data11	*data_m0Eg;	/* start address of m(0) and Eg		*/
data11	*data_zwf;	/* start address of z and wf		*/

/* default values */

np_flag=false;
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
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  efcwwf [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
	   }
	   break;
  default :
	   printf("Usage:  efcwwf [-a include non-parabolicity \033[1mfalse\033[0m]\n");
	   printf("               [-p particle (\033[1me\033[0m, h, or l)]\n");
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
 E*=1e-3*e_0;	/* convert E from meV->J	*/
 N=wf(E,delta_z,data_start,data_m0Eg,data_zwf,n,np_flag,i_state);

 sprintf(filename,"wf_%c%i.r",p,i_state);
 Fwf=fopen(filename,"w");

 for(i=0;i<n;i++)
  fprintf(Fwf,"%20.17le %20.17le\n",(data_zwf+i)->a,((data_zwf+i)->b)/sqrt(N));

 fclose(Fwf);

}/* end while */

fclose(FE);
free(data_start);
free(data_m0Eg);
free(data_zwf);


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
 while(fscanf(Fv,"%*le %*le")!=EOF)
  (*n)++;
 rewind(Fv);

 data_start=(files *)calloc(*n,sizeof(files));
 if(data_start==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 fdata=data_start;

 while(fscanf(Fv,"%le %le",&(fdata->z),&(fdata->V))!=EOF)
 {
  fscanf(Fm,"%*le %le",&(fdata->mstar));
  fdata++;
 }

 fclose(Fm);
 fclose(Fv);

 return(data_start);

}



data11
*read_Egdata(n,data_start)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines         */

int     n;
files  *data_start;    /* start address of potential              */

{
 int	i;		/* index				*/
 data11	*data_m0Eg;	/* pointer to m0 and Eg data		*/
 FILE   *FEg;		/* file pointer to potential file	*/

 if((FEg=fopen("Eg.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'Eg.r'!\n");exit(0);}

 data_m0Eg=(data11 *)calloc(n,sizeof(data11));
 if(data_m0Eg==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 for(i=0;i<n;i++)
 {
  fscanf(FEg,"%*le %le",&(data_m0Eg+i)->b);
  ((data_m0Eg+i)->a)=(data_start+i)->mstar;
 }

 fclose(FEg);

 return(data_m0Eg);

}



double
wf(E,delta_z,fdata,data_m0Eg,data_zwf,n,np_flag,i_state)     

/* This function returns the value of the wavefunction (psi)
   at +infinity for a given value of the energy.  The solution
   to the energy occurs for psi(+infinity)=0.                      */

double E;
double delta_z;
files  *fdata;
data11 *data_m0Eg;
data11 *data_zwf;
int    n;
boolean np_flag;
int	i_state;
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
   alpha=sqr(1-((data_m0Eg+i)->a)/m0)/((data_m0Eg+i)->b);
   (fdata+i)->mstar=((data_m0Eg+i)->a)*(1+alpha*(E-((fdata+i)->V)));
  }
 }

 /* boundary conditions */
 
 psi[0]=1.0;
 psi[1]=1.0;

 (data_zwf)->b=psi[0];
 (data_zwf+1)->b=psi[1];

 N+=sqr(psi[0]);
 N+=sqr(psi[1]);

 fdata++;                    /* ignore data corresponding to psi[0] */

  for(i=1;i<(n-1);i++)              /* last potential not used */
  {
   psi[2]=(
           2*(fdata->z)*(2*(fdata->mstar)*sqr(delta_z/hbar)*(fdata->V-E)+2)*
	   psi[1]+
           (-2*(fdata->z)+delta_z)*psi[0]
          )
           /(2*(fdata->z)+delta_z);
   (data_zwf+i+1)->b=psi[2];
   N+=sqr(psi[2]);
   psi[0]=psi[1];
   psi[1]=psi[2];
   fdata++;
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

