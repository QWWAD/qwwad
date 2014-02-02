/* Need this second maths include file for linking with the NAG libraries
 * only.  NAG use their propriety type `Complex', rather than  `complex'
 * 									*/

Complex	Cadd();		/* complex addition		*/
Complex	Csub();		/* complex subtraction		*/
Complex	Cmult();	/* Complex mulitplication	*/
Complex	Cconj();	/* complex conjugate		*/
Complex	Cexp();		/* complex exponential		*/
double	Cmod();		/* modulus of complex number	*/


/****************** multiply two complex numbers ***********************/

Complex Cmult(c1,c2)
Complex c1,c2;
{
 Complex	c;
 c.re=c1.re*c2.re-c1.im*c2.im;
 c.im=c1.re*c2.im+c1.im*c2.re;
 return(c);
}



/****************** complex conjugate ***********************/


Complex Cconj(c)
Complex c;
{
 Complex	cstar;
 cstar.re=c.re;
 cstar.im=-1*c.im;
 return(cstar);
}



/****************** add two complex numbers ***********************/


Complex Cadd(c1,c2)
Complex c1,c2;
{
 Complex	c;
 c.re=c1.re+c2.re;
 c.im=c1.im+c2.im;
 return(c);
}


/****************** subtract two complex numbers ***********************/


Complex Csub(c1,c2)
Complex c1,c2;
{
 Complex	c;
 c.re=c1.re-c2.re;
 c.im=c1.im-c2.im;
 return(c);
}


/****************** complex exponential ***********************/

Complex Cexp(d)
double	d;
{
 Complex	c;
 c.re=cos(d);
 c.im=sin(d);
 return(c);
}



/****************** returns modulus of a complex number ****************/

double Cmod(c)
Complex c;
{
 return(sqrt((c.re)*(c.re)+(c.im)*(c.im)));
}

