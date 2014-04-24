/*==================================================================
              efcwire  Envelope Function Circular WIRE
  ==================================================================*/

/* This program uses a shooting technique to calculate the
   uncorrelated one particle energies of any user supplied
   radial potential.  The potential is read from the file v.r

   This program has been butchered from `efshoot', there may be 
   a little redundant code lying around...

   Paul Harrison, December 1998                                   */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <gsl/gsl_math.h>
#include "ef-helpers.h"
#include "struct.h"
#include "maths.h"
#include "qclsim-constants.h"

static double psi_at_inf(const double  E,
                         const double  delta_z,
                         files        *fdata,
                         const data11 *data_m0Eg,
                         const int     n,
                         const bool    np_flag);

int main(int argc,char *argv[])
{
double read_delta_z();
double V_min();
files  *read_data();	/* reads potential file into memory  */

double d_E;		/* infinitesmal energy               */
double delta_E;		/* small but finite energy           */
double delta_z;		/* z separation of input potentials  */
double dy;		/* derivative of function            */
double E;		/* electron and hole energies        */
double E_start;		/* initial energies (if not default) */
double x;		/* independent variable (energy)     */
double y;		/* function (psi at infinity)        */
double y1;		/* temporary y value                 */
double y2;		/* temporary y value                 */
int    i_state;		/* electron/hole state index         */
int    n;		/* length of potential file		 */
int    state;		/* electron and hole output states   */
char   p;		/* particle				 */
char   filename[9];	/* output filename			 */
bool   np_flag;         /* Hamiltonian flag def.=1=>D(1/m)D  */
FILE   *FE;		/* outputfile for el. energy states  */
files  *data_start;	/* start address of potential	 */
data11 *data_m0Eg=NULL;	/* start address of m(0) and Eg		*/

/* default values */

delta_E=1e-3*e;
E_start=0e-3*e;
np_flag=false;
p='e';
state=1;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
	   np_flag=true;
	   argv--;
	   argc++;
	   break;
  case 'd':
	   delta_E=atof(argv[2])*1e-3*e;
	   break;
  case 'e':
	   E_start=atof(argv[2])*1e-3*e;
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  efcwire [-p particle (e, h, or l)]\n");
                      exit(0);
	   }
	   break;
  case 's':
	   state=atoi(argv[2]);
	   break;
  default :
           printf("Usage:  efcwire [-a include non-parabolicity \033[1mfalse\033[0m]\n");
	   printf("                [-d energy step (\033[1m1\033[0mmeV)][-e initial energy (\033[1m0\033[0mmeV)]\n");
	   printf("                [-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("                [-s number of states \033[1m1\033[0m]\n");
	   exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

d_E=1e-8*e;

sprintf(filename,"E%c.r",p);
FE=fopen(filename,"w");

data_start=read_data(&n);			/* reads potential file	*/

if(np_flag)
    data_m0Eg=read_Egdata(n,data_start);	/* reads bandgap data	*/	

delta_z=read_delta_z(data_start);

if(fabs(E_start)<1e-3*1e-3*e) x=V_min(data_start,n);
 else x=E_start;                       /* initial energy estimate */

for(i_state=1;i_state<=state;i_state++)  
{
 /* increment energy-search for f(x)=0 */
 y2=psi_at_inf(x,delta_z,data_start,data_m0Eg,n,np_flag);

 do
 {
  y1=y2;
  x+=delta_E;
  y2=psi_at_inf(x,delta_z,data_start,data_m0Eg,n,np_flag);
 }while(y1*y2>0);

/* improve estimate using midpoint rule */

 x-=fabs(y2)/(fabs(y1)+fabs(y2))*delta_E;

/* implement Newton-Raphson method */

 do
 {
  y=psi_at_inf(x,delta_z,data_start,data_m0Eg,n,np_flag);
  dy=(psi_at_inf(x+d_E,delta_z,data_start,data_m0Eg,n,np_flag)-
      psi_at_inf(x-d_E,delta_z,data_start,data_m0Eg,n,np_flag))/
     (2.0*d_E);
  x-=y/dy;
 }while(fabs(y/dy)>1e-12*e);

 E=x;

 fprintf(FE,"%i %24.17le\n",i_state,E/(1e-3*e));

 x+=delta_E;    /* clears x from solution */

 }/* end i_state */

fclose(FE);
free(data_start);
if(np_flag)free(data_m0Eg);

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

/**
 * Finds the value of the wavefunction (\f$\psi\f$)
 * at \f$+\infty\f$ for a given energy. The solution
 * to the energy occurs for \f$\psi(+\infty)=0\f$.
 *
 * \param[in]      E         energy
 * \param[in]      delta_z
 * \param[in, out] fdata
 * \param[in]      data_m0Eg
 * \param[in]      n
 * \param[in]      np_flag   Whether to use non-parabolic effective mass
 *
 * \returns The value of the wavefunction at \f$\infty\f$
 */
static double psi_at_inf(const double  E,
                         const double  delta_z,
                         files        *fdata,
                         const data11 *data_m0Eg,
                         const int     n,
                         const bool    np_flag)
{
 double psi[3];              /* wavefunction at z-delta_z,
                                z and z+delta_z              */
 int	i;		     /* index			     */

 /* account for non-parabolicity with mstar(E)	*/

 if(np_flag)
 {
  for(i=0;i<n;i++)
  {
   /* Find nonparabolicity parameter using Eq. 3.77, QWWAD3 */
   const double alpha=gsl_pow_2(1-data_m0Eg[i].a/me)/data_m0Eg[i].b;

   /* Find effective mass at the desired energy using Eq. 3.76, QWWAD3 */
   fdata[i].mstar=data_m0Eg[i].a*(1.0+alpha*(E-fdata[i].V));
  }
 }

 /* boundary conditions: Eq. 8.55, QWWAD3 */
 psi[0]=1.0;
 psi[1]=1.0;
 psi[2]=0;

 fdata++;                    /* ignore data corresponding to psi[0] */

  for(i=1;i<(n-1);i++)              /* last potential not used */
  {
   /* Find wavefunction at next point, using Eq. 8.54, QWWAD3 */
   psi[2]=(
           2*(fdata->z)*(2*(fdata->mstar)*gsl_pow_2(delta_z/hBar)*(fdata->V-E)+2)*
	   psi[1]+
           (-2*(fdata->z)+delta_z)*psi[0]
          )
           /(2*(fdata->z)+delta_z);

   /* Shift previous values along ready for next iteration */
   psi[0]=psi[1];
   psi[1]=psi[2];
   fdata++;
  } 

 return psi[2];
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

