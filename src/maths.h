/*************************************/
/*  include file for math functions  */
/*************************************/

double sqr(const double x);
double cub();
double sec();
double cosec();
double cot();
int sign();
double coth();
double modulus();
int Nint();
vector	vadd();		/* adds two vectors together	*/
vector	vsub();		/* subtracts one vector from another	*/
vector	vmult();	/* multiplies a vector by a constant	*/
double	vmod();		/* modulus of vector		*/
vector	vvprod();	/* vector product		*/
double	vsprod();	/* vector product		*/
double Theta();		/* Heaviside unit step function	*/


/**
 * square of a number
 *
 * \param[in] x The number to be squared
 *
 * \return The square of the number
 */
double sqr(const double x)
{
  return x*x;
}

/**************** cube function *********************/
double cub(x)
double x;
{
return(x*x*x);
}

/********************** secans **********************/
double
sec(x)
double x;
{
 return(1.0/cos(x));
}


/********************** cosecans ********************/
double
cosec(x)
double x;
{
 return(1/sin(x));
}



/********************** cotangent ********************/
double
cot(x)
double x;
{
 return(1/tan(x));
}


/*********************** sign ***********************/
int
sign(x)
double x;
{
 int i;

 if (x<0) i=-1;
 else i=1;
 return(i);
}


/*************** cotangens hyperbolicus *************/
double
coth(x)
double x;
{
 return(1/tanh(x));
}




/**************** Takes the modulus of a number ***************/
double modulus(x)

double x;
{
if(x<0)x*=-1;
return(x);
}

/****************** add two vectors ***********************************/

vector vadd(A,B)
vector	A;
vector	B;
{
 vector	C;
 
 C.x=A.x+B.x;
 C.y=A.y+B.y;
 C.z=A.z+B.z;

 return(C);
}

/****************** subtracts on vector from another ***********************/

vector vsub(A,B)
vector	A;
vector	B;
{
 vector	C;
 
 C.x=A.x-B.x;
 C.y=A.y-B.y;
 C.z=A.z-B.z;

 return(C);
}


/****************** multiply a vector bya constant ********************/

vector vmult(A,c)
vector	A;
double	c;
{ 
 vector	B;

 B.x=c*A.x;
 B.y=c*A.y;
 B.z=c*A.z;

 return(B);
}

/****************** returns modulus of a vector ***********************/

double vmod(A)
vector A;
{
 return(sqrt((A.x)*(A.x)+(A.y)*(A.y)+(A.z)*(A.z)));
}


/****************** returns vector product of two vectors *************/

vector vvprod(A,B)
vector	A,B;
{
 vector	C;

 C.x=A.y*B.z-A.z*B.y;
 C.y=-(A.x*B.z-A.z*B.x);
 C.z=A.x*B.y-A.y*B.x;

 return(C);
}
 
/****************** returns scalar product of two vectors *************/

double vsprod(A,B)
vector	A,B;
{
 return(A.x*B.x+A.y*B.y+A.z*B.z);
}

/**************** Heaviside Step function *****************************/
double Theta(x)

double x;
{
 return((x<0)?0:1);
}

