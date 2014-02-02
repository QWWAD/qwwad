/*====================================================================
                        ovl  OVerLap integral
  ====================================================================*/

/* This program simply calculates the overlap integral of two user
   supplied wavefunctions, the filenames of which are read in as
   command line arguments.

   Input files:       wf_e1.r
                      wf_e2.r    etc.

   Output files:    

                                                                */

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "const.h"

main(int argc, char *argv[])
{
double    overlap();       /* calculates overlap integral       */
double    overlapmod();    /* calculates overlap integral       */
double    O;               /* overlap integral                  */
double    Omod;            /* overlap of the mod integral       */
int       N_1;             /* number of lines in data files     */
int       N_2;             /* number of lines in data files     */
char      filename_1[9];   /* character string for wavefunction 
                              input file                        */
char      filename_2[9];   /* character string for wavefunction 
                              input file                        */
data11    *read_data();
data11    *start_wf1;      /* start address of data             */
data11    *start_wf2;      /* start address of data             */


if(argc!=3)
 {printf("Usage: ovl [first wavefunction filename][second wavefunction filename]\n");exit(0);
 }

(void)strcpy(filename_1,argv[1]);
(void)strcpy(filename_2,argv[2]);

start_wf1=read_data(filename_1,&N_1);
start_wf2=read_data(filename_2,&N_2);

if(N_1!=N_2)
 {printf("Error: number of lines in %s and %s are not equal!\n",
         filename_1,filename_2);exit(0);
 }

O=overlap(start_wf1,start_wf2,N_1);
Omod=overlapmod(start_wf1,start_wf2,N_1);

/* Write result to standard output in a form suitable for processing
  with awk, for example							*/

printf("<Y1|Y2> %20.17le <|Y1|||Y2|> %20.17le\n",O,Omod);

free(start_wf1);
free(start_wf2);

}        /* end main */



double
overlap(pointer_wf1,pointer_wf2,N)

data11 *pointer_wf1;       /* pointer to first wavefunction data */
data11 *pointer_wf2;       /* pointer to first wavefunction data */
int     N;                 /* number of wavefucntion points      */
{
 double delta_z;
 double overlap=0;
 int    i;                 /* index                              */

 delta_z=((pointer_wf1+1)->a)-(pointer_wf1->a);
 for(i=0;i<N;i++)
 {
  overlap+=(pointer_wf1->b)*(pointer_wf2->b);
  pointer_wf1++;pointer_wf2++;
 }
 overlap*=delta_z;

 return(overlap);
}



double
overlapmod(pointer_wf1,pointer_wf2,N)

data11 *pointer_wf1;       /* pointer to first wavefunction data */
data11 *pointer_wf2;       /* pointer to first wavefunction data */
int     N;                 /* number of wavefucntion points      */
{
 double delta_z;
 double overlap=0;
 int    i;                 /* index                              */

 delta_z=((pointer_wf1+1)->a)-(pointer_wf1->a);
 for(i=0;i<N;i++)
 {
  overlap+=fabs(pointer_wf1->b)*fabs(pointer_wf2->b);
  pointer_wf1++;pointer_wf2++;
 }
 overlap*=delta_z;

 return(overlap);
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


