/*==================================================================
                  gde General Diffusion Equation
  ==================================================================*/

/* This program produces the general solution to the diffusion
   equation

           dn=d (D dn)
           -- --(  --)
           dt dx(  dx)

   for n=n(x,t) and D=D(x,t,n).

   Input files:
     x.r           initial (t=0) concentration profile versus z  

   Output files:
     X.r           final (diffused) concentration profile 

   Paul Harrison, July 1994
	
   Modifications June 1998                                       */

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <malloc.h>
#include <math.h>
#include "struct.h"
#include "const.h"
#include "maths.h"
#include "bools.h"
#include "io.h"

#include "dox.c"

typedef
struct	{
 double	x;		    /* concentration of diffusant 	  */
 double	z;		    /* z value of files    		  */
} data;
 
typedef
struct	{
 double	D;		    /* diffusion coefficient 	          */
 double	z;		    /* z value of files    		  */
} diff;

main(int argc,char *argv[])
{

data   *read_data();        /* reads data from external files    */
diff   *create_D();         /* creates diffusion coefficient     */
double read_delta_z();      /* deduces interval along z-axis     */
void   calculate_D();       /* calculates D for D=D(x)           */
void   diffuse();           /* diffuse for time delta_t          */
void   read_D();

double D0;                  /* the constant value of D           */
double delta_t;             /* time interval between iterations  */
double delta_z;             /* interval along growth (z-) axis   */
double t;                   /* time                              */
double t_final;             /* final time of diffusion           */
int    form_of_D;           /* dependencies of D                 */
int    i;
int    n;		    /* length of potential file		 */
char   c;                   /* general character                 */
FILE   *FX;                 /* file pointer to X.r               */
data   *data_start;    	    /* start address of data             */
data   *ddata;              /* pointer to data                   */
diff   *D_start;            /* start address of data             */

/* default values */

D0=1e-20;	/* 1 Angstrom^2/s	*/
delta_t=0.01;	
form_of_D=0;
t_final=1.0;


while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
           form_of_D=atoi(argv[2]);
	   switch(form_of_D)
	   {
	    case 1: break;
	    case 2: break;
	    case 3: break;
	    default: printf("Usage: gde [-a form of diffusion coefficient (\033[1m0\033[0m, 1, 2 or 3)\n");
		     printf("       \033[1m0\033[0m=>D=D0 a constant, specified by the -D option\n");
		     printf("       1=>D=D(z) only, read in from file D.r\n");
		     printf("       2=>D=D(D,x,z,t), general dependency, completely defined in dox.c\n");
		     printf("       3=>D=D(z,t), initial z-dependence from D.r, t-dependence defined in dox.c\n");
	             exit(0);
	   }
           break;
  case 'd':
	   delta_t=atof(argv[2]);
	   break;
  case 'D':
	   D0=atof(argv[2])*1e-20;	/* convert from A^2/s-> m^2/s	*/
	   break;
  case 't':
	   t_final=atof(argv[2]);
	   break;
  default:
	   printf("Usage: gde [-a form of diffusion coefficient (\033[1m0\033[0m, 1, 2 or 3)\n");
	   printf("           [-D the constant coefficient, if a=0 selected (\033[1m1\033[0mA^2/s)\n");
	   printf("           [-d time interval dt (\033[1m0.01\033[0ms)\n");
	   printf("           [-t final time (\033[1m1\033[0ms)\n");
	   exit(0);
	  
 }
 argv++;
 argv++;
 argc--;
 argc--;
}


data_start=read_data(&n);               /* reads all external data files */

delta_z=read_delta_z(data_start);       /* calculates interval delta_z   */

D_start=create_D(n);                    /* creates structure for D       */

FX=fopen("X.r","w");                    /* opens output file for writing */

switch(form_of_D)
{
 case 0:
	for(i=0;i<n;i++)(D_start+i)->D=D0;	/* set constant D	*/
	for(t=delta_t;t<=t_final;t+=delta_t)
        {
         diffuse(data_start,D_start,FX,delta_t,delta_z,n);
        }break;

 case 1:
        read_D(D_start);			/* read D from file	*/
        for(t=delta_t;t<=t_final;t+=delta_t)
        {
         diffuse(data_start,D_start,FX,delta_t,delta_z,n);
        }break;
 case 2:
        for(t=delta_t;t<=t_final;t+=delta_t)
        {
         calculate_D(data_start,D_start,t,n);	/* calculate D using `dox.c'	*/
         diffuse(data_start,D_start,FX,delta_t,delta_z,n);
        }break;
 case 3:
        read_D(D_start);			/* read D at t=0 from file	*/
        for(t=delta_t;t<=t_final;t+=delta_t)
        {
         calculate_D(data_start,D_start,t,n);	/* calculate subsequent D */
         diffuse(data_start,D_start,FX,delta_t,delta_z,n);
        }break;
}

ddata=data_start;
for(i=0;i<n;i++)
{
 fprintf(FX,"%20.17le %le %le\n",(ddata->z),(ddata->x),0.0);
 ddata++;
}

fclose(FX);

free(data_start);
free(D_start);

} /* end main */





void
calculate_D(ddata,dD,t,n)

/* This function recalculates the diffusion coefficient for all
   points along the z-axis when D is a function of the
   concentration.                                                */

data *ddata;    /* pointer to data, initialised to start address */
diff *dD;       /* pointer to diffusion structure, initialised
                   to start address                              */
double t;       /* time                                          */
int   n;
{
 extern double D_of_x();    /* concentration dependence of D=D(x)*/
 int i_z;                   /* index to all points along z-axis  */

 for(i_z=0;i_z<n;i_z++)     /* note limits cover all data points */
 {
  dD->D=D_of_x(dD->D,ddata->x,ddata->z,t);
  ddata++;
  dD++;
 }

}



void
diffuse(data_start,dD,FX,delta_t,delta_z,n)

/* This function projects the diffusant profile a short
   time interval delta_t into the future                         */

data   *data_start;         /* start address of data             */
diff   *dD;                 /* pointer to diffusion coefficient  */
FILE   *FX;                 /* file pointer to ABC.r             */
double delta_t;
double delta_z;
int    n;
{
 data   *ddata;              /* pointer to data structure         */
 data   *dnewdata;           /* pointer to data structure         */
 data   *newdata_start;      /* start address of data             */
 int    i;
 int    i_z;                 /* index through z values            */

 newdata_start=(data *)calloc(n,sizeof(data));
 if (newdata_start==0)  
 {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }
 
 ddata=data_start;
 dnewdata=newdata_start;

 ddata++;                      /* increment pointers to            */
 dnewdata++;                   /* z (x in theory!), x (n in        */ 
 dD++;                         /* theory!) and new x               */

 for(i_z=1;i_z<n-1;i_z++)     /* note limits                       */
 {
  dnewdata->x=delta_t*
              (
               (((dD+1)->D)-((dD-1)->D))*
               (((ddata+1)->x)-((ddata-1)->x))/sqr(2*delta_z)
              +(dD->D)*
               (((ddata+1)->x)-2*(ddata->x)+((ddata-1)->x))/sqr(delta_z)
              )
              +(ddata->x);
  dnewdata++;
  ddata++;
  dD++;
 }

/* Impose `closed-system boundary conditions */

 (newdata_start->x)=((newdata_start+1)->x);         
 ((newdata_start+n-1)->x)=((newdata_start+n-2)->x); 


 ddata=data_start;
 dnewdata=newdata_start;
 for(i=0;i<n;i++)
 {
  (ddata->x)=(dnewdata->x);
  ddata++;
  dnewdata++;
 }
 
free(newdata_start);


}



double 
read_delta_z(dp)

/* This function calculates the separation along the z (growth) 
   direction of the user supplied concentrations                   */

data *dp;               /* pointer to data                         */
{
 double z[2];           /* displacement along growth direction     */

 z[0]=dp->z;
 dp++;
 z[1]=dp->z;
 return(z[1]-z[0]);
}



diff
*create_D(n)

/* This function creates the memory space for the diffusion D
   versus z structure.                                             */

int	n;
{
 diff   *D_start;	/* start address of concentration          */

 D_start=(diff *)calloc(n,sizeof(diff));
 if (D_start==0)  
 {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 return(D_start);
}



void
read_D(dD)

/* This function reads the diffusion coefficient as a function 
   of z from the external file D.r                                 */

diff   *dD;		/* pointer to concentration                */
{
 FILE 	*FD;            /* file pointer to concentration file      */

 if((FD=fopen("D.r","r"))==0)
 {
  fprintf(stderr,"Error: Cannot open input file 'D.r'!\n");
  exit(0);
 }

 while(fscanf(FD,"%le %le", &(dD->z), &(dD->D))!=EOF)
  dD++;

 fclose(FD);

}



data
*read_data(n)

/* This function reads the initial (t=0) diffusant concentration
   profile into memory and returns the start
   address of this block of memory and the number of lines	   */

int	*n;

{
 FILE 	*Fx;            /* file pointer to concentration file      */
 data   *dt;		/* temporary pointer to concentration      */
 data   *data_start;	/* start address of concentration          */

 if((Fx=fopen("x.r","r"))==0)
 {
  fprintf(stderr,"Error: Cannot open input file 'x.r'!\n");
  exit(0);
 }
 *n=0;
 while(fscanf(Fx,"%*le %*le %*le")!=EOF)
  (*n)++;
 rewind(Fx);

 data_start=(data *)calloc(*n,sizeof(data));
 if (data_start==0)  
 {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }
 dt=data_start;

 /* Ignore the second alloy component data y(z) in input file `x.r'	*/

 while(fscanf(Fx,"%le %le %*le", &(dt->z), &(dt->x))!=EOF)
  dt++;

 fclose(Fx);

 return(data_start);

}

