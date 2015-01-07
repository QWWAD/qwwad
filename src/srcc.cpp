/*==================================================================
              srcc  Scattering Rate Carrier-Carrier
  ==================================================================*/

/* This program calculates the carrier-carrier scattering rate for
   both intra- and intersubband events.  The required rates are provided
   by the user in the file `rr.r'.  The other necessary inputs are listed
   below.

	Input files:	rr.r	contains required rates
			wf_xy.r	x=particle y=state
			N.r	subband populations
			Ex.r	x=particle, energies
			Ef.r	subband Fermi energies

	Output files:	ccABCD.r	cc rate versus Ei for each mechanism


    Paul Harrison, March 1997
    Modifications, February 1999					*/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <gsl/gsl_math.h>
#include "struct.h"
#include "qclsim-constants.h"
#include "qclsim-subband.h"

using namespace Leeds;
using namespace constants;

typedef
struct	{
 double	a;		/* z value of files			*/
 double	b[4];		/* wave function			*/
} data14;

static void output_ff(const double       W,
                      const std::vector<Subband> &subbands,
                      const unsigned int i,
                      const unsigned int j,
                      const unsigned int f,
                      const unsigned int g);

data11 * ff_table(const double                 Deltak0sqr,
                  const Subband               &isb,
                  const Subband               &jsb,
                  const Subband               &fsb,
                  const Subband               &gsb,
                  const std::valarray<double> &V,
                  const size_t                 nq);

data11 * PI_table(data11        *Aijfg,
                  const Subband &isb,
                  const double   T,
                  const size_t   nq,
                  const bool     S_flag);

double lookup_ff(data11       *Aijfg,
                 const double  q_perp,
                 const size_t  nq);

double lookup_PI(data11       *PIii,
                 const double  q_perp,
                 const size_t  nq);

double A(const std::valarray<double> &z,
         const double                 q_perp,
         const data14                *wf);

double PI(const Subband &isb,
          const size_t   q_perp,
          const double   T);

int main(int argc,char *argv[])
{
double	alpha;		/* angle between ki and kj			*/
double	Deltak0sqr;	/* twice the change in KE, see Smet (55)	*/
double	dalpha;		/* step length for alpha integration		*/
double	dtheta;		/* step length for theta integration		*/
double	dki;		/* step length for loop over ki 		*/
double	dkj;		/* step length for kj integration		*/
double	epsilon;	/* low frequency dielectric constant		*/
double	ki;		/* carrier momentum				*/
double	kij;		/* (vector)kj-(vector)(ki)			*/
double	kimax;		/* maximum value of ki				*/
double	kj;		/* carrier momentum				*/
double	kjmax;		/* maximum value of kj				*/
double	m;		/* carrier effective mass			*/
double	P;		/* probability factor, see Smet			*/
double	q_perp;		/* in-plane momentum, |ki-kf|			*/
double	q_perpsqr4;	/* 4*q_perp*q_perp, checks vailidity of q_perp	*/
double	T;		/* temperature					*/
double	theta;		/* angle between kij and kfg			*/
double	W;		/* arbitrary well width, soley for output	*/
double	Wbar;		/* FD weighted mean of Wijfg			*/
double	Wijfg;		/* the carrier-carrier scattering rate		*/
int	ialpha;		/* index over alpha				*/
int	iki;		/* index over ki				*/
int	ikj;		/* index over kj				*/
int	itheta;		/* index over theta				*/
int	nalpha;		/* number of strips in alpha integration	*/
int	nki;		/* number of ki calculations			*/
int	nkj;		/* number of strips in |kj| integration		*/
int	nq;		/* number of q_perp values for lookup table	*/
int	ntheta;		/* number of strips in theta integration	*/
char	filename[9];	/* character string for output filename		*/
char	p;		/* particle					*/
bool	ff_flag;	/* form factor flag, output to file if true	*/
bool	S_flag;		/* screening flag, include screening if true	*/
data11	*Aijfg;		/* form factor over all 4 wave functions	*/
data11	*PIii;		/* pointer to screening function versus q_perp	*/
FILE	*Fcc;		/* pointer to output file			*/
FILE	*FccABCD;	/* pointer to weighted mean output file		*/

/* default values */

epsilon=13.18*eps0;/* low frequency dielectric constant for GaAs	*/
ff_flag=false;		/* don't output formfactors	*/
m=0.067*me;		/* GaAs electron value		*/
p='e';			/* electron			*/
T=300;			/* temperature			*/
W=250e-10;		/* a well width, same as Smet	*/
S_flag=true;		/* include screening by default	*/

/* default values for numerical calculations	*/

nalpha=100;
ntheta=100;
nki=100;
nkj=100;
nq=100;

/* calculate step lengths	*/

dalpha=2*pi/((float)nalpha);
dtheta=2*pi/((float)ntheta);

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
	   ff_flag=true;
	   argv--;
           argc++;
           break;
  case 'e':
	   epsilon=atof(argv[2])*eps0;
	   break;
  case 'm':
	   m=atof(argv[2])*me;
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  srcc [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
	   }
	   break;
  case 'S':
	   S_flag=false;
	   argv--;
           argc++;
           break;
  case 'T':
	   T=atof(argv[2]);
	   break;
  case 'w':
	   W=atof(argv[2])*1e-10;
	   break;
  default :
	   printf("Usage:  srcc [-a generate form factors \033[1mfalse\033[0m][-e permittivity (\033[1m13.18\033[0mepsilon_0)]\n");
	   printf("             [-m mass (\033[1m0.067m0\033[0m)][-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("             [-S screening \033[1mtrue\033[0m]\n");
	   printf("             [-T temperature (\033[1m300\033[0mK)][-w a well width (\033[1m250\033[0mA)]\n");
	   exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

std::ostringstream E_filename; // Energy filename string
E_filename << "E" << p << ".r";
std::ostringstream wf_prefix;  // Wavefunction filename prefix
wf_prefix << "wf_" << p;

// Read data for all subbands from file
std::vector<Subband> subbands = Subband::read_from_file(E_filename.str(),
                                                        wf_prefix.str(),
                                                        ".r",
                                                        m);

// Read and set carrier distributions within each subband
std::valarray<double>       Ef;      // Fermi energies [J]
std::valarray<double>       N;       // Subband populations [m^{-2}]
std::valarray<unsigned int> indices; // Subband indices (garbage)
read_table_xy("Ef.r", indices, Ef);
Ef *= e/1000.0; // Rescale to J
read_table_xy("N.r", indices, N);	// read populations
N *= 1e+10*1e+4; // convert from units of 10^10cm^-2->m^-2

for(unsigned int isb = 0; isb < subbands.size(); ++isb)
    subbands[isb].set_distribution(Ef[isb], N[isb]);

// Read potential profile
std::valarray<double> z;
std::valarray<double> V;
read_table_xy("v.r", z, V);

// Read list of wanted transitions
std::valarray<unsigned int> i_indices;
std::valarray<unsigned int> j_indices;
std::valarray<unsigned int> f_indices;
std::valarray<unsigned int> g_indices;

read_table_xyzu("rr.r", i_indices, j_indices, f_indices, g_indices);

FccABCD=fopen("ccABCD.r","w");	/* open file for output of weighted means */

// Loop over all desired transitions
for(unsigned int itx = 0; itx < i_indices.size(); ++itx)
{
    // State indices for this transition (NB., these are indexed from 1)
    unsigned int i = i_indices[itx];
    unsigned int j = j_indices[itx];
    unsigned int f = f_indices[itx];
    unsigned int g = g_indices[itx];

    // Convenience labels for each subband (NB., these are indexed from 0)
    const Subband isb = subbands[i-1];
    const Subband jsb = subbands[j-1];
    const Subband fsb = subbands[f-1];
    const Subband gsb = subbands[g-1];

    // Subband minima
    const double Ei = isb.get_E();
    const double Ej = jsb.get_E();
    const double Ef = fsb.get_E();
    const double Eg = gsb.get_E();

    // Output form-factors if desired
    if(ff_flag)
        output_ff(W,subbands,i,j,f,g);

    /* Generate filename for particular mechanism and open file	*/
    sprintf(filename,"cc%i%i%i%i.r",i,j,f,g);
    Fcc=fopen(filename,"w");

    // Calculate Delta k0^2 [QWWAD3, Eq. 10.228]
    if(i+j == f+g) Deltak0sqr=0;
    else
        Deltak0sqr=4*m*(Ei + Ej - Ef - Eg)/(hBar*hBar);	

    Aijfg = ff_table(Deltak0sqr, isb, jsb, fsb, gsb, V,nq);
    PIii  = PI_table(Aijfg,isb,T,nq,S_flag);

    /* calculate maximum value of ki & kj and hence kj step length	*/
    kimax=sqrt(2*m*(V.max()-Ei))/hBar;	/* sqr(hBar*kimax)/2m=Vmax-Ei	*/
    dki=kimax/((float)nki);

    kjmax=sqrt(2*m*(V.max()-Ej))/hBar;	/* sqr(hBar*kjmax)/2m=Vmax-Ej	*/
    dkj=kjmax/((float)nkj);

 Wbar=0;			/* initialise integral sum */

 for(iki=0;iki<nki;iki++)	/* calculate c-c rate for all ki	*/
 {
 ki=dki*(float)iki;
 Wijfg=0;			/* Initialize for integration	*/
 for(ikj=0;ikj<nkj;ikj++)	/* integrate over |kj|		*/
 {
  kj=dkj*(float)ikj;

  // Find Fermi-Dirac occupation at kj
  P=jsb.f_FD_k(kj,T);

  for(ialpha=0;ialpha<nalpha;ialpha++)	/* Integral over alpha	*/
  {
   alpha=dalpha*(float)ialpha;
 
   kij=sqrt(ki*ki+kj*kj-2*ki*kj*cos(alpha));	/* calculate kij	*/

   for(itheta=0;itheta<ntheta;itheta++)		/* Integral over theta	*/
   {
    theta=dtheta*(float)itheta;	/* move origin to avoid 0*/
  
    /* calculate argument of sqrt function=4*q_perp*q_perp, see (8.179),
       to check for imaginary q_perp, if argument is positive, q_perp is
       real and hence calculate scattering rate, otherwise ignore and move
       onto next q_perp							*/

    q_perpsqr4=2*kij*kij+Deltak0sqr-2*kij*sqrt(kij*kij+Deltak0sqr)*cos(theta);

    /* recall areal element in plane polars is kj.dkj.dalpha */

    if(q_perpsqr4>=0) 
    {
     q_perp=sqrt(q_perpsqr4)/2;
     Wijfg+=gsl_pow_2(lookup_ff(Aijfg,q_perp,nq)/
                (q_perp+2*pi*e*e*lookup_PI(PIii,q_perp,nq)
                 *lookup_ff(Aijfg,q_perp,nq)/(4*pi*epsilon)
                )
               )*P*kj;
    }

    /* Note screening term has been absorbed into the above equation
       to avoid pole at q_perp=0				*/

   } /* end theta */

  } /* end alpha */

 } /* end kj   */
 
 Wijfg*=dtheta*dalpha*dkj;	/* multiply by all step lengths	*/

 Wijfg*=gsl_pow_2(e*e/(hBar*4*pi*epsilon))*m/(pi*hBar);

 /* output scattering rate versus carrier energy=subband minima+in-plane
    kinetic energy						*/
 fprintf(Fcc,"%20.17le %20.17le\n",(Ei+gsl_pow_2(hBar*ki)/(2*m))/
                                   (1e-3*e),Wijfg);

 /* calculate Fermi-Dirac weighted mean of scattering rates over the 
    initial carrier states, note that the integral step length 
    dE=2*sqr(hBar)*ki*dki/(2m)					*/
 Wbar+=Wijfg*ki*isb.f_FD_k(ki, T);
 } /* end ki	*/

 Wbar*=dki/(pi*isb.get_pop());

 fprintf(FccABCD,"%i %i %i %i %20.17le\n", i,j,f,g,Wbar);

 fclose(Fcc);	/* close output file for this mechanism	*/

 free(Aijfg);
 free(PIii);

} /* end while over states */

fclose(FccABCD);	/* close weighted mean output file	*/

return EXIT_SUCCESS;
} /* end main */

/* This function calculates the overlap integral over all four carrier
   states		*/
double A(const double   q_perp,
         const Subband &isb,
         const Subband &jsb,
         const Subband &fsb,
         const Subband &gsb)
{
 double	A=0;	// integral over z and hence form factor
 double	B;	/* integral over z'	*/
 const std::valarray<double> z = isb.z_array();
 const size_t nz = z.size();
 const double dz = z[1] - z[0];

 // Convenience labels for wave-functions in each subband
 const std::valarray<double> psi_i = isb.psi_array();
 const std::valarray<double> psi_j = jsb.psi_array();
 const std::valarray<double> psi_f = fsb.psi_array();
 const std::valarray<double> psi_g = gsb.psi_array();

 // Integral of i(=0) and f(=2) over z
 for(unsigned int iz=0;iz<nz;iz++)
 {
  B=0;

  // Integral of |j> and |g> over z'
  for(unsigned int izd=0;izd<nz;izd++)
      B += (psi_j[izd]*psi_g[izd])*exp(-q_perp*fabs(z[iz] - z[izd]));

  B *= dz;

  A += psi_i[iz] * psi_f[iz] * B;
 }

 A *= dz;

 return(A);
}



/**
 * \brief returns the screening factor, referred to by Smet as e_sc
 */
double PI(const Subband &isb,
          const double   q_perp,
          const double   T)
{
    const double m    = isb.get_md_0();    // Effective mass at band-edge [kg]
    const double ki_F = isb.get_k_fermi(); // Fermi wave-vector [1/m]

    // Equation 43 of Smet, QWWAD3, 10.236
    //
    // TODO: This is incorrect... need to bring it INSIDE the integral, because we should
    //       be using the wave-vector k(Ek), NOT the Fermi wave-vector!!!
    //
    //       See QWWAD4 for corrected calculation
    double	P = m/(pi*hBar*hBar);

    if(q_perp>2*ki_F)
        P -= m/(pi*hBar*hBar)*sqrt(1-gsl_pow_2(2*ki_F/q_perp));

    // Now perform the integration, equation 44 of Smet [QWWAD3, 10.238]
    double mu=isb.get_E();
    const double dmu=1e-3*e;
    double integral=0;

    double	dI; // Intervals in I

    do
    {
        dI=1/(4*kB*T*gsl_pow_2(cosh((isb.get_Ef() - mu)/(2*kB*T))));
        integral+=dI*dmu;
        mu+=dmu;
    }while(dI>integral/100); // continue until integral converges

    P*=integral;

    return P;
}

/* This function creates the polarizability table */
data11 * PI_table(data11        *Aijfg,
                  const Subband &isb,
                  const double   T,
                  const size_t   nq,
                  const bool     S_flag)
{
 data11	*PIii = new data11[nq];

 for(unsigned int iq=0;iq<nq;iq++)
 {
     const double q_perp = Aijfg[iq].a;
     PIii[iq].a = Aijfg[iq].a; // Use same scattering-vectors as matrix element table

     // Allow screening to be turned off	*/
     if(S_flag)
         PIii[iq].b = PI(isb, q_perp, T);
     else PIii[iq].b = 0;
 }

 return PIii;
}

/* This function creates a form factor look up table	*/
data11 * ff_table(const double  Deltak0sqr,
                  const Subband               &isb,
                  const Subband               &jsb,
                  const Subband               &fsb,
                  const Subband               &gsb,
                  const std::valarray<double> &V,
                  const size_t  nq)
{
 double	dq;	/* interval in q_perp			*/
 double	q_perp;	/* In-plane wave vector			*/
 double	q_perp_max;	/* maximum in-plane wave vector	*/
 data11	*Aijfg;

 const double vmax=V.max(); // Maximum potential [J]
 const double kimax=isb.k(vmax - isb.get_E()); // Max value of ki [1/m]
 const double kjmax=jsb.k(vmax - jsb.get_E()); // Max value of kj [1/m]

 Aijfg=(data11 *)calloc(nq,sizeof(data11));
  if (Aijfg==0)  {
   fprintf(stderr,"Cannot allocate memory!\n");
   exit(0);
  }

 q_perp_max=sqrt(2*gsl_pow_2(kimax+kjmax)+Deltak0sqr+2*(kimax+kjmax)*
                 sqrt(gsl_pow_2(kimax+kjmax)+Deltak0sqr))/2;

 dq=q_perp_max/((float)(nq-1));
 for(unsigned int iq=0;iq<nq;iq++)
 {
  q_perp=dq*(float)iq;
  Aijfg[iq].a=q_perp;
  Aijfg[iq].b=A(q_perp, isb, jsb, fsb, gsb);
 }

 return(Aijfg);
}

/**
 *  \brief looks up the form factor Aijfg in the table generated by ff_table
 */
double lookup_ff(data11       *Aijfg,
                 const double  q_perp,
                 const size_t  nq)
{
 double	A;	/* The looked up form factor	*/
 int	iq=0;	/* index over q_perp in table	*/

 if(q_perp>((Aijfg+nq-1)->a))
 {printf("q_perp>maximum allowed q!\n");exit(0);}

 while(q_perp>=((Aijfg+iq)->a))
 {
  iq++;		/* retreive the ff above q_perp	*/
 }

 /* Linearly interpolate the form factor from between the values
    directly above and below q_perp				*/

 A=((Aijfg+iq-1)->b)+(((Aijfg+iq)->b)-((Aijfg+iq-1)->b))*
                   (q_perp-((Aijfg+iq-1)->a))/
                   (((Aijfg+iq)->a)-((Aijfg+iq-1)->a));

 return(A);

}

/**
 * \brief looks up the screening function PIii in the table generated by PI_table
 */
double lookup_PI(data11       *PIii,
                 const double  q_perp,
                 const size_t  nq)
{
 double	P;	/* The looked up screening function	*/
 int	iq=0;	/* index over q_perp in table		*/

 if(q_perp>((PIii+nq-1)->a))
 {printf("q_perp>maximum allowed q!\n");exit(0);}

 while(q_perp>=((PIii+iq)->a))
 {
  iq++;		/* retreive the ff above q_perp	*/
 }

 /* Linearly interpolate the form factor from between the values
    directly above and below q_perp				*/

 P=((PIii+iq-1)->b)+(((PIii+iq)->b)-((PIii+iq-1)->b))*
                   (q_perp-((PIii+iq-1)->a))/
                   (((PIii+iq)->a)-((PIii+iq-1)->a));

 return(P);

}

/* This function outputs the formfactors into files	*/
static void output_ff(const double        W, // Arbitrary well width to generate q
                      const std::vector<Subband> &subbands,
                      const unsigned int  i,
                      const unsigned int  j,
                      const unsigned int  f,
                      const unsigned int  g)
{
 char	filename[9];	/* output filename				*/
 FILE	*FA;		/* output file for form factors versus q_perp	*/

 /* First generate filename and then open file	*/

 sprintf(filename,"A%i%i%i%i.r", i, j, f, g);	
 if((FA=fopen(filename,"w"))==0)
  {fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);exit(0);}

 // Convenience labels for each subband (NB., these are indexed from 0)
 const Subband isb = subbands[i-1];
 const Subband jsb = subbands[j-1];
 const Subband fsb = subbands[f-1];
 const Subband gsb = subbands[g-1];

 for(unsigned int iq=0;iq<100;iq++)
 {
  const double q_perp=6*iq/(100*W); // In-plane scattering vector
  const double Aijfg=A(q_perp,isb,jsb,fsb,gsb);
  fprintf(FA,"%le %le\n",q_perp*W,gsl_pow_2(Aijfg));
 }

 fclose(FA);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
