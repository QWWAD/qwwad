/*==========================================================================
                                   struct.h
  ==========================================================================*/
#ifndef STRUCT_H
#define STRUCT_H

#ifdef __cplusplus
# include <complex>
#else
# include <complex.h>
#endif

typedef struct
{
 double x;
 double y;
 double z;
}vector;

typedef struct
{
 double a;
 double b;
}data11;
#endif /* STRUCT_H */
