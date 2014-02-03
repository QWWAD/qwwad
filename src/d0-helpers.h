/**
 * \file  d0-helpers.h
 * \brief Functions that are common to all donor calculations
 */
#ifndef D0_HELPERS
#define D0_HELPERS

#include "struct.h"

data11 * read_v      (size_t       *n);
double   read_delta_z(data11       *Vp);
double   V_min       (data11       *Vp,
                      const size_t  n);

/**
 * Reads the potential into memory and returns the start
 * address of this block of memory and the number of lines
 *
 * \param[out] n The number of lines in the v.r file
 *
 * \returns The potential profile
 */
data11 * read_v(size_t *n)
{
 FILE   *fp;            /* file pointer to potential file          */
 data11 *Vp;           /* temporary pointer to potential          */
 data11 *Vstart;       /* start address of potential              */

 if((fp=fopen("v.r","r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file 'v.r'!\n");
   exit(0);
 }

 *n=0;
 while(fscanf(fp,"%*e %*e")!=EOF)
  (*n)++;
 rewind(fp);

 Vstart = (data11 *)calloc(*n,sizeof(data11));
 if (Vstart==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }
 Vp = Vstart;

 while(fscanf(fp,"%le %le", &(Vp->a), &(Vp->b))!=EOF)
  Vp++;

 fclose(fp);

 return Vstart;
}

/**
 * Calculates the separation along the z (growth)
 * direction of the user supplied potentials
 *
 * \param Vp Potential profile
 *
 * \returns Spatial separation
 */
double read_delta_z(data11 *Vp)
{
 double z[2];           /* displacement along growth direction     */

 z[0] = Vp->a;
 Vp++;
 z[1] = Vp->a;
 return(z[1]-z[0]);
}

/**
 * This function opens the external file v.r and finds     
 * the minimum value for the potential energy, this value
 * is used as the initial energy estimate.
 *
 * \param[in] Vp potential profile
 * \param[in] n  number of steps in potential profile
 *
 * \returns The minimum potential
 */
double V_min(data11       *Vp,
             const size_t  n)
{
 double       min = 1; /* minimum value of potential energy       */
 unsigned int i;       /* index                                   */

 for(i=0; i<n; i++)
 {
  if(Vp->b<min)
  {
   min=Vp->b;
  }
  Vp++;
 }

 return min;
}
#endif
