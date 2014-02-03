/**
 * \file  ef-helpers.h
 * \brief Functions that are widely used in Envelope Function calculations
 */
#ifndef EF_HELPERS_H
#define EF_HELPERS_H

#include <error.h>
#include "struct.h"

typedef
struct	{
 double	z;		/* z value of files    		  */
 double	V;		/* electron and hole values   	  */
 double	mstar;		/* electron and hole values   	  */
} files;

data11 * read_Egdata(const size_t  n,
                     files        *data_start);

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
   error(EXIT_FAILURE, 0, "Cannot open input file 'Eg.r'!");

 data_m0Eg=(data11 *)calloc(n, sizeof(data11));
 if(data_m0Eg==0)
   error(EXIT_FAILURE, 0, "Cannot allocate memory!");

 for(i=0;i<n;i++)
 {
  int n_read = fscanf(FEg,"%*e %le",&(data_m0Eg+i)->b);
  if(n_read != 2)
    error(EXIT_FAILURE, 0, "Data missing on line %d of %s", i, "Eg.r");

  data_m0Eg[i].a = data_start[i].mstar;
 }

 fclose(FEg);

 return data_m0Eg;
}
#endif /* EF_HELPERS_H */
