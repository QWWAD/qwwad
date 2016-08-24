/**
 * \file  ef-helpers.h
 * \brief Functions that are widely used in Envelope Function calculations
 */
#ifndef EF_HELPERS_H
#define EF_HELPERS_H

#include "struct.h"

typedef
struct	{
 double	z;		/* z value of files    		  */
 double	V;		/* electron and hole values   	  */
 double	mstar;		/* electron and hole values   	  */
} files;

data11 * read_Egdata(const size_t  n,
                     files        *data_start);

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */
files * read_data(int *n)
{
 FILE 	*Fv;            /* file pointer to potential file          */
 FILE 	*Fm;            /* file pointer to potential file          */
 files  *fdata;	/* start address of potential		   */
 int i;

 if((Fm=fopen("m.r","r"))==0)
 {
     fprintf(stderr,"Error: Cannot open input file 'm.r'!\n");
     exit(0);
 }

 if((Fv=fopen("v.r","r"))==0)
 {
     fprintf(stderr,"Error: Cannot open input file 'v.r'!\n");
     exit(0);
 }

 *n=0;
 while(fscanf(Fv,"%*e %*e")!=EOF)
  (*n)++;
 rewind(Fv);

 fdata=(files *)calloc(*n,sizeof(files));

 if(fdata==0){
     fprintf(stderr,"Cannot allocate memory!\n");
     exit(0);
 }

 for(i=0; i<*n; ++i)
 {
     int n_read = fscanf(Fv,"%le %le",&(fdata[i].z),&(fdata[i].V));
     n_read = fscanf(Fm,"%*e %le",&(fdata[i].mstar));
 }

 fclose(Fm);
 fclose(Fv);

 return(fdata);
}

/**
 * Reads the potential into memory and returns the start address of this
 * block of memory
 *
 * \param[in]  n          Number of lines to read
 * \param[out] data_start Start address of potential
 */
data11 * read_Egdata(const size_t  n,
                     files        *data_start)
{
 unsigned int  i;         /* index				*/
 data11	      *data_m0Eg; /* pointer to m0 and Eg data		*/
 FILE         *FEg;	  /* file pointer to potential file	*/

 if((FEg=fopen("Eg.r","r"))==0)
 {
   fprintf(stderr, "Cannot open input file 'Eg.r'!");
   exit(EXIT_FAILURE);
 }

 data_m0Eg=(data11 *)calloc(n, sizeof(data11));
 if(data_m0Eg==0)
 {
   fprintf(stderr, "Cannot allocate memory!");
   exit(EXIT_FAILURE);
 }

 for(i=0;i<n;i++)
 {
  int n_read = fscanf(FEg,"%*e %le",&(data_m0Eg+i)->b);
  if(n_read != 2)
  {
   fprintf(stderr, "Data missing on line %d of %s", i, "Eg.r");
   exit(EXIT_FAILURE);
  }

  data_m0Eg[i].a = data_start[i].mstar;
 }

 fclose(FEg);

 return data_m0Eg;
}

/**
 * \brief Calculates separation between spatial points
 *
 * \param[in] files User-supplied 2-column data
 *
 * \returns Spatial separation
 *
 * \details Reads the separation along the z (growth) direction from a
 *          user-supplied data table
 */
double read_delta_z(files *fdata)
{
 double z[2];           /* displacement along growth direction     */

 z[0]=fdata->z;
 fdata++;
 z[1]=fdata->z;
 return(z[1]-z[0]);
}

/**
 * \brief Finds minimum value for potential energy
 *
 * \param[in] fdata pointer to potential
 * \param[in] n     number of steps in potential
 *
 * \returns minimum value of potential
 *
 * \details This function opens the external file v.r and finds     
 *          the minimum value for the potential energy, this value
 *          is used as the initial energy estimate.
 */
double V_min(files *fdata, int n)       
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
#endif /* EF_HELPERS_H */
