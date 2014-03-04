/*************************************/
/*  include file for math functions  */
/*************************************/

#ifndef MATHS_H
#define MATHS_H
#include "struct.h"
double sec   (const double x);
double cosec (const double x);
double coth  (const double x);
vector vadd  (const vector A,
              const vector B);
vector vsub  (const vector A,
              const vector B);
vector vmult (const vector A,
              const double c);
double vmod  (const vector A);
vector vvprod(const vector A,
              const vector B);
double vsprod(const vector A,
              const vector B);
double Theta (const double x);

/**
 * The secant of a number
 *
 * \param[in] x The number for which to find the secant (radians)
 *
 * \return The secant
 */
double sec(const double x)
{
 return 1.0/cos(x);
}

/**
 * The cosecant of a number
 *
 * \param[in] x The number for which to find the cosecant (radians)
 *
 * \return The cosecant
 */
double cosec(const double x)
{
 return 1/sin(x);
}

/**
 * cotangens hyperbolicus
 */
double coth(const double x)
{
 return 1/tanh(x);
}

/**
 * add two vectors
 */
vector vadd(const vector A, const vector B)
{
 vector	C;
 
 C.x=A.x+B.x;
 C.y=A.y+B.y;
 C.z=A.z+B.z;

 return C;
}

/**
 * subtracts on vector from another
 */
vector vsub(const vector A, const vector B)
{
 vector	C;
 
 C.x=A.x-B.x;
 C.y=A.y-B.y;
 C.z=A.z-B.z;

 return C;
}

/**
 * multiply a vector by a constant
 */
vector vmult(const vector A, const double c)
{ 
 vector	B;

 B.x=c*A.x;
 B.y=c*A.y;
 B.z=c*A.z;

 return B;
}

/**
 * returns modulus of a vector
 */
double vmod(const vector A)
{
 return sqrt((A.x)*(A.x)+(A.y)*(A.y)+(A.z)*(A.z));
}

/**
 * returns vector product of two vectors
 */
vector vvprod(const vector A, const vector B)
{
 vector	C;

 C.x=A.y*B.z-A.z*B.y;
 C.y=-(A.x*B.z-A.z*B.x);
 C.z=A.x*B.y-A.y*B.x;

 return C;
}
 
/**
 * returns scalar product of two vectors
 */
double vsprod(const vector A, const vector B)
{
 return A.x*B.x+A.y*B.y+A.z*B.z;
}

/**
 * Heaviside Step function
 */
double Theta(const double x)
{
 return (x<0)?0:1;
}
#endif
