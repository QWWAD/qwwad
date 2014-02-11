/*==========================================================================
                                   struct.h
  ==========================================================================*/
#ifndef STRUCT_H
#define STRUCT_H

#include <complex.h>

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

typedef struct
{
 double a;
 double b[2];
}data12;

typedef struct
{
#if __cplusplus
  std::complex<double> M[2][2];
#else
  complex double M[2][2];
#endif
}cmat2x2;
#endif /* STRUCT_H */
