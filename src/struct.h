/*==========================================================================
                                   struct.h
  ==========================================================================*/

typedef struct
{
 double x;
 double y;
 double z;
}vector;

typedef struct
{
 double re;
 double im;
}complex;

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
 complex M[2][2];
}cmat2x2;
