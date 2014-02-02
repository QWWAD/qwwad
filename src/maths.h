/*************************************/
/*  include file for math functions  */
/*************************************/

double sqr();
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
double	cmod();		/* modulus of complex number	*/
complex	cadd();		/* complex addition		*/
complex	csub();		/* complex subtraction		*/
complex	cmult();	/* complex mulitplication	*/
complex	cdiv();		/* complex division		*/
complex	cconj();	/* complex conjugate		*/
complex	cexp();		/* complex exponential		*/
double Theta();		/* Heaviside unit step function	*/


/**************** square of function ****************/
double sqr(x)
double x;
{
return(x*x);
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

/**************** returns nearest integer to a number **************/
int Nint(x)
double  x;
{
double  r;
double  y;
int	i;

y=modf(x,&r);

if (y>=0.5)
  i=(int)ceil(x);
if (y>=0 && y<0.5)
  i=(int)floor(x);
if (y<0 && y>-0.5)
  i=(int)ceil(x);
if (y<=-0.5)
  i=(int)floor(x);

return(i);
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
/****************** returns modulus of a complex number ****************/

double cmod(c)
complex c;
{
 return(sqrt((c.re)*(c.re)+(c.im)*(c.im)));
}


/****************** multiply two complex numbers ***********************/

complex cmult(c1,c2)
complex c1,c2;
{
 complex	c;
 c.re=c1.re*c2.re-c1.im*c2.im;
 c.im=c1.re*c2.im+c1.im*c2.re;
 return(c);
}


/****************** divide two complex numbers ***********************/

complex cdiv(c1,c2)
complex c1,c2;
{
 complex	c;
 c.re=(c1.re*c2.re+c1.im*c2.im)/(c2.re*c2.re+c2.im*c2.im);
 c.im=(c2.re*c1.im-c1.re*c2.im)/(c2.re*c2.re+c2.im*c2.im);
 return(c);
}



/****************** complex conjugate ***********************/

complex cconj(c)
complex c;
{
 complex	cstar;
 cstar.re=c.re;
 cstar.im=-1*c.im;
 return(cstar);
}


/****************** complex exponential ***********************/

complex cexp(d)
double	d;
{
 complex	c;
 c.re=cos(d);
 c.im=sin(d);
 return(c);
}



/****************** add two complex numbers ***********************/

complex cadd(c1,c2)
complex c1,c2;
{
 complex	c;
 c.re=c1.re+c2.re;
 c.im=c1.im+c2.im;
 return(c);
}


/****************** subtract two complex numbers ***********************/

complex csub(c1,c2)
complex c1,c2;
{
 complex	c;
 c.re=c1.re-c2.re;
 c.im=c1.im-c2.im;
 return(c);
}

/**************** Heaviside Step function *****************************/
double Theta(x)

double x;
{
 return((x<0)?0:1);
}

