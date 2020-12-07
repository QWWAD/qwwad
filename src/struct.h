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

struct vector
{
 double x;
 double y;
 double z;
};

struct data11
{
 double a;
 double b;
};
#endif /* STRUCT_H */
