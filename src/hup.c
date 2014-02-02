/*====================================================================
                  hup Heisenberg's Uncertainty Principle
  ====================================================================*/

/* This program calculates the uncertainty in the position and the
   momentum of any user supplied wave function.  The data is written 
   to the standard output.

   Input files:       wf_ps.r		p=particle (e, h or l)
					s=state

   Output files:    

   Paul Harrison, April 1998						*/
                                                               

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "const.h"

main(int argc, char *argv[])
{
double	ev_z;		/* Expectation Value of z			*/
double	ev_zsqr;	/* Expectation Value of sqr(z)			*/
double	ev_p;		/* Expectation Value of p			*/
double	ev_psqr;	/* Expectation Value of sqr(p)			*/
double	Delta_p;	/* Uncertainty in p 				*/
double	Delta_z;	/* Uncertainty in z 				*/
double	delta_z;	/* z step length				*/
int	N;		/* number of lines in data files		*/
int	i;		/* index					*/
int	state;		/* principal quantum number			*/
char	filename[9];	/* character string for wavefunction input file	*/
char	p;		/* particle (e, h, or l)			*/
data11	*read_data();
data11	*psi;		/* start address of wavefunction		*/


/* default values */

p='e';
state=1;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': break;
            case 'h': break;
            case 'l': break;
            default:  printf("Usage:  hup [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
           }
           break;
  case 's':
           state=atoi(argv[2]);
           break;
  default :
           printf("Usage:  hup [-p particle (\033[1me\033[0m, h, or l)][-s state \033[1m1\033[0m]\n");
           exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

sprintf(filename,"wf_%c%i.r",p,state);

psi=read_data(filename,&N);

ev_z=0;ev_zsqr=0;ev_p=0;ev_psqr=0;

/* Momentum algorithm requires omission of end points.
   Although step length calculated on each pass, algorithms do
   indeed assume uniform mesh					*/

for(i=1;i<N-1;i++)	
{
 delta_z=((psi+i+1)->a)-((psi+i)->a);
 ev_z+=((psi+i)->b)*((psi+i)->a)*((psi+i)->b)*delta_z;
 ev_zsqr+=((psi+i)->b)*sqr((psi+i)->a)*((psi+i)->b)*delta_z;

/* Note leave <p> in units of hbar and <p^2> in hbar^2 */

 ev_p-=((psi+i)->b)*
       (((psi+i+1)->b)-((psi+i-1)->b))/(2*delta_z)
       *delta_z;
 ev_psqr-=((psi+i)->b)*
          (((psi+i+1)->b)+((psi+i-1)->b)-2*((psi+i)->b))/sqr(delta_z)
          *delta_z;
}

Delta_z=sqrt(ev_zsqr-sqr(ev_z));
Delta_p=sqrt(ev_psqr-sqr(ev_p));

printf("<z>\t\t\t\t\t%20.17le\n",ev_z);
printf("<z^2>\t\t\t\t\t%20.17le\n",ev_zsqr);
printf("Delta_z=sqrt(<z^2>-<z>^2)\t\t%20.17le\n",Delta_z);

printf("\n");

printf("<p>/hbar\t\t\t\t%20.17le i\n",ev_p);
printf("<p^2>/sqr(hbar)\t\t\t\t%20.17le\n",ev_psqr);
printf("Delta_p/hbar=sqrt(<p^2>-<p>^2)/hbar\t%20.17le\n",Delta_p);

printf("\n");

printf("Delta_z*Delta_p\t\t\t\t%20.17le hbar\n",Delta_z*Delta_p);

free(psi);

}        /* end main */




data11
*read_data(filename,N)

/* This function reads the wave function into memory and returns
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
 while(fscanf(Fwf,"%*le %*le")!=EOF)
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


