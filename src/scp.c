/*==================================================================
                    scp Self-Consistent Poisson 
  ==================================================================*/

/* This program calculates the potential due to a one-dimensional 
   distribution of dopants, in the file `d.r' and free carriers, 
   described by the wave functions `wf_pn.r'.  This potential is added 
   to the original band edge potential to produce a new combined potential
   suitable for continued solution with `efshoot'.

   Unnormalised subband populations in the file `Nu.r' are normalised 
   by integration the doping density over the z-domain.

	Input files:		
				wf_pn.r	p=particle n=state
				Nu.r	unnormalised subband populations

	Output files:		F.r	field strength, E(z)
				N.r	subband populations
				sigma.r	areal charge density
				V.r	Poisson potential


    Paul Harrison, March 1997                                   */

#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <malloc.h>
#include <gsl/gsl_math.h>
#include "struct.h"
#include "const.h"
#include "maths.h"

static double * read_E(char p, int *s);
static double * calc_field(const double  epsilon,
                           const double *sigma,
                           const double *z,
                           const size_t  n);

int main(int argc,char *argv[])
{
double	*calc_potn();	/* calculates the potential due to charge 	*/
double	*calc_sigma();	/* calculates the net areal charge density	*/
double	*read_N();	/* reads subband populations from file		*/
double	*read_Nda();	/* reads doping levels from file		*/
double 	*read_z();	/* reads z co-ordinates				*/
double	*read_wf();	/* reads wavefunctions into memory		*/
void	comb_potn();	/* combines `v.r' with charge potential		*/
void	write_field();	/* output field E(z) to file			*/
void	write_potn();	/* output potential to file			*/
void	write_sigma();	/* output charge density to file		*/

double	delta_z;	/* mesh length along growth (z-) axis		*/
double	epsilon;	/* low frequency dielectric constant		*/
double	*F;		/* electric field along z-axis			*/
double	*N;		/* pointer to subband populations		*/
double	*Nda;		/* dopant concentration in m^-3			*/
double	Ntotal;		/* total carrier concentration yielded by dopant*/
double	Nutotal;	/* total unnormalised carrier concentration	*/
double	q;		/* charge on carrier +/-e_0			*/
double	*sigma;		/* net areal charge density in m^-2		*/
double	*V;		/* potential due to charge distribution		*/
double	**wfs;		/* pointer to start addresses of wave functions	*/
double	*z;		/* pointer to z co-ordinates			*/
int	i;		/* general index				*/
int	is;		/* index over states				*/
int	n;		/* length of wavefunctions file			*/
int	s;		/* number of states as specified by Ep.r	*/
char	p;		/* particle					*/
FILE	*FN;		/* pointer to output file `N.r'			*/

/* default values */

epsilon=13.18*epsilon_0;/* low frequency dielectric constant for GaAs	*/
p='e';			/* electron			*/
q=-e_0;

/* default values for numerical calculations	*/

/* calculate step lengths	*/


while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'e':
	   epsilon=atof(argv[2])*epsilon_0;
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': q=-e_0;break;
	    case 'h': q=+e_0;break;
	    case 'l': q=+e_0;break;
	    default:  printf("Usage:  scp [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
	   }
	   break;
  default :
	   printf("Usage:  scp [-e permittivity (\033[1m13.18\033[0mepsilon_0]\n");
	   printf("            [-p particle (\033[1me\033[0m, h, or l)]\n");
	   exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

double *E=read_E(p, &s); /* read subband energies. Actually, we're only interested in number of states */
(void)E; /* Suppress compiler warning about unused variable */

N=read_N(s);		/* read in unnormalised subband populations	*/

z=read_z(&n,p);		/* read in z co-ordinates	*/
delta_z=*(z+1)-*(z);	/* calculate z step length, assume constant	*/

Nda=read_Nda(n);	/* read in dopant concentration	*/
 
/* Allocate memory for structure of pointers to wave function pointers */

wfs=(double **)calloc(s,sizeof(double *));
if(wfs==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* Now read in wave functions	*/

for(is=0;is<s;is++)
{
 *(wfs+is)=read_wf(is+1,n,p);
} 

/* Calculate total carrier concentration per unit area, yielded by dopant */

Ntotal=0;		/* Initialise sum	*/
for(i=0;i<n;i++)
{
 Ntotal+=*(Nda+i)*delta_z;
}
 
/* Calculate sum of unnormalised subband populations	*/

Nutotal=0;for(is=0;is<s;is++)Nutotal+=*(N+is);

/* And hence normalize populations, writing also to `N.r'	*/

FN=fopen("N.r","w");
for(is=0;is<s;is++)
{
 *(N+is)*=Ntotal/Nutotal;
 fprintf(FN,"%i %le\n",is+1,*(N+is)/(1e+4*1e+10)); /* m^-2->10^10cm^-2 */
}
fclose(FN);

/* Calculate `net' areal charge density	and output to file*/

sigma=calc_sigma(delta_z,N,Nda,wfs,n,s);
write_sigma(sigma,z,n);

/* Calculate electric field and output to file	*/

F=calc_field(epsilon,sigma,z,n);
write_field(F,z,n);

/* Calculate potential due to charge distribution and output to file	*/

V=calc_potn(delta_z,F,q,n);
write_potn(V,z,n);

/* Combine band edge potential with potential due to charge distribution */

comb_potn(V,z,n);

free(F);free(N);free(Nda);free(sigma);free(V);free(wfs);free(z);

return EXIT_SUCCESS;
} /* end main */


/**
 * Calculates the electric field along the z axis
 *
 * \param[in] epsilon permittivity
 * \param[in] sigma   areal charge density [number/m^2]
 * \param[in] z       spatial locations
 * \param[in] n       number of spatial samples
 *
 * \returns Electric field [V/m]
 */
static double * calc_field(const double  epsilon,
                           const double *sigma,
                           const double *z,
                           const size_t  n)
{
 double	      *F;      /* electric field as a function of z-	*/
 unsigned int  iz;     /* index over z co-ordinates	*/
 unsigned int  izdash; /* index over z' co-ordinates	*/

 /* Allocate memory for wave function */

 F=(double *)calloc(n,sizeof(double));
 if(F==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 /* Find total field using Eq. 3.83, QWWAD3 */
 for(iz=0;iz<n;iz++)
 {
  F[iz]=0; /* initialise F at each z co-ordinate */

  for(izdash=0;izdash<n;izdash++)
  {
   /* Note sigma is a number density per unit area, needs to be converted
      to Couloumb per unit area						*/
   F[iz] += sigma[izdash]*e_0*(float)GSL_SIGN(z[iz]-z[izdash])/epsilon;
  }
 }

return(F);
}



double
*calc_potn(delta_z,F,q,n)

/* This function calculates the potential (energy actually)	*/

double	delta_z;
double	*F;
double	q;
int	n;

{
 double	*V;		/* electric field as a function of z-	*/
 int	iz;		/* index over z co-ordinates	*/

 /* Allocate memory for wave function */

 V=(double *)calloc(n,sizeof(double));
 if(V==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 /* Calculate the potential, defining the first point as zero	*/

 *V=0;
 for(iz=1;iz<n;iz++)
 {
  *(V+iz)=*(V+iz-1)-q*(*(F+iz))*delta_z;
 }

return(V);
}



double
*calc_sigma(delta_z,N,Nda,wfs,n,s)

/* This function calculates `net' areal charge density	*/

double	delta_z;
double	*N;
double	*Nda;
double	**wfs;
int	n;
int	s;

{
 double	*sigma;
 int	i;		/* index over z co-ordinates	*/
 int	is;		/* index over states		*/

 /* Allocate memory for wave function */

 sigma=(double *)calloc(n,sizeof(double));
 if(sigma==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 for(i=0;i<n;i++)
 {
  *(sigma+i)=0;		/* initialise sigma at each z co-ordinate	*/
  for(is=0;is<s;is++)	/* sum over all subbands			*/
  {
   *(sigma+i)+=*(N+is)*gsl_pow_2(*(*(wfs+is)+i))*delta_z;
  }
  /* n-type dopants give -ve *(N+is) representing electrons, hence 
     addition of +ve ionised donors requires -*(Nda+i), note Nda is still a
     volume density, the delta_z converts it to an areal density	*/

  *(sigma+i)-=*(Nda+i)*delta_z;	
 }

return(sigma);
}



void
comb_potn(V,z,n)

/* This function combines the band edge potential with the potential
   due to the charge distribution				*/

double	*V;
double	*z;
int	n;

{
 double	v;		/* band edge potential			*/
 int	i;		/* index over z co-ordinates		*/
 int	nv;		/* number of lines in potential file	*/
 FILE	*Fv;		/* pointer to output file `v.r'	*/
 FILE	*Fv1;		/* pointer to input file `v1.r'	*/

 /* If the file `v1.r' exists, open it, if not create it first	*/

 if((Fv1=fopen("v1.r","r"))==0)
 {
  if((system("cp v.r v1.r"))==0) Fv1=fopen("v1.r","r");
  else 
   {printf("Error: Cannot open file 'v1.r' or find file 'v.r'!\n");exit(0);}
 }

 Fv=fopen("v.r","w");

 nv=0;
 while(fscanf(Fv1,"%*f %*f")!=EOF)
  nv++;
 rewind(Fv1);

 if(nv!=n)
 {printf("Error: number of lines in `v.r' not equal to those in `wf_??.r!\n");
  exit(0);}

 for(i=0;i<n;i++)
 {
  int n_read = fscanf(Fv1,"%*f %lf",&v);
  if(n_read != 2)
    error(EXIT_FAILURE, 0, "Data missing in v.r");
  fprintf(Fv,"%20.17le %20.17le\n",*(z+i),*(V+i)+v);
 }

 fclose(Fv);
 fclose(Fv1);

}

/**
 * Reads the subband minima into memory and returns
 * the number of defined states, and the start address of the data
 */
static double * read_E(char p, int *s)
{
 double	*E;
 int	i=0;		/* index over the energies			*/
 char	filename[9];	/* filename string				*/
 FILE 	*FE;		/* file pointer to energy data 			*/

 sprintf(filename,"E%c.r",p);
 if((FE=fopen(filename,"r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);
   exit(0);
 }

 *s=0;
 while(fscanf(FE,"%*i %*e")!=EOF)
  (*s)++;
 rewind(FE);

 E=(double *)calloc(*s,sizeof(double));
 if (E==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FE,"%*i %le",E+i)!=EOF)
 {
  *(E+i)*=1e-3*e_0;		/*convert meV->J		*/
  i++;
 }

 fclose(FE);

 return(E);

}

double
*read_N(s)

/* This function reads the unormalised subband populations into memory 
   and returns the start address of this block of memory	   */

int	s;

{
 double	*N;
 int	i=0;		/* index over the energies			*/
 FILE 	*FN;		/* file pointer to subband populations data	*/

 if((FN=fopen("Nu.r","r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file 'Nu.r'!\n");
   exit(0);
 }

 N=(double *)calloc(s,sizeof(double));
 if (N==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FN,"%*i %le",N+i)!=EOF)
 {
  *(N+i)*=1e+10*1e+4;	/*convert from units of 10^10cm^-2->m^-2 */
  i++;
 }

 fclose(FN);

 /* Now check number of defined populations equals number of subband minima */

 if(s!=i){
  fprintf(stderr,"Error: Number of defined populations in file 'Nu.r' is \n");
  fprintf(stderr,"       not equal to number of subband minima in Ep.r'\n");
  exit(0);
 }

 return(N);
}


double
*read_Nda(n)

/* This function reads the volume density of donors or acceptors 
   into memory and returns the start address 	   */

int	n;

{
 int	i;
 FILE 	*Fd;		/* file pointer to dopants file	*/
 double	*Nda;		/* pointer to dopant structure	*/

 /* Allocate memory for wave function */

 Nda=(double *)calloc(n,sizeof(double));
 if(Nda==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 /* Read in wave function	*/

  if((Fd=fopen("d.r","r"))==0)
   {fprintf(stderr,"Error: Cannot open input file '%s'!\n","d.r");exit(0);}

  i=0;			
  while(fscanf(Fd,"%*e %le",Nda+i)!=EOF)
   i++;

 fclose(Fd);			

 return(Nda);

}



double
*read_wf(is,n,p)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */

int	is;		/* particular state index	*/
int	n;
char	p;

{
 char	filename[9];	/* input filename			*/
 int	i;
 FILE 	*Fwf;		/* file pointer to wavefunction file	*/
 double	*wf;		/* pointer to wave function structure	*/

 /* Allocate memory for wave function */

 wf=(double *)calloc(n,sizeof(double));
 if(wf==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 /* Read in wave function	*/

  sprintf(filename,"wf_%c%i.r",p,is);	/* Open each file	*/
  if((Fwf=fopen(filename,"r"))==0)
   {fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);exit(0);}

  i=0;					
  while(fscanf(Fwf,"%*e %le",wf+i)!=EOF)
   i++;

 fclose(Fwf);				

 return(wf);

}



double
*read_z(n,p)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */

int	*n;
char	p;

{
 double	*z;		/* pointer to z co-ordinates		*/
 int	i;		/* index				*/
 char	filename[9];	/* input filename			*/
 FILE 	*Fwf;		/* file pointer to wavefunction file	*/

 /* Open first file and count number of lines */

 sprintf(filename,"wf_%c1.r",p);
 if((Fwf=fopen(filename,"r"))==0)
  {fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);exit(0);}
 
 *n=0;	
 while(fscanf(Fwf,"%*e %*e")!=EOF)
  (*n)++;
 rewind(Fwf);

 /* Allocate memory for z co-ordinates	*/

 z=(double *)calloc(*n,sizeof(double));
 if(z==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 /* Read in data 	*/

 i=0;						/* Read in each file	*/
 while(fscanf(Fwf,"%le %*e",z+i)!=EOF)
  i++;

 fclose(Fwf);					/* Close each file	*/


 return(z);

}



void
write_field(F,z,n)

/* This function writes the electric field data to a file	*/

double	*F;
double	*z;
int	n;

{
 int	i;		/* index over z co-ordinates		*/
 FILE	*FF;		/* pointer to output file `sigma.r'	*/

 FF=fopen("F.r","w");

 for(i=0;i<n;i++)
  fprintf(FF,"%20.17le %20.17le\n",*(z+i),*(F+i));

 fclose(FF);

}



void
write_potn(V,z,n)

/* This function writes the potential data to a file	*/

double	*V;
double	*z;
int	n;

{
 int	i;		/* index over z co-ordinates		*/
 FILE	*FV;		/* pointer to output file `V.r'	*/

 FV=fopen("V.r","w");

 for(i=0;i<n;i++)
  fprintf(FV,"%20.17le %20.17le\n",*(z+i),*(V+i));

 fclose(FV);

}



void
write_sigma(sigma,z,n)

/* This function writes the `net' areal charge density to a file	*/

double	*sigma;
double	*z;
int	n;

{
 int	i;		/* index over z co-ordinates		*/
 FILE	*Fsigma;	/* pointer to output file `sigma.r'	*/

 Fsigma=fopen("sigma.r","w");

 for(i=0;i<n;i++)
  fprintf(Fsigma,"%20.17le %20.17le\n",*(z+i),*(sigma+i));

 fclose(Fsigma);

}
