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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "qclsim-constants.h"
#include "qclsim-subband.h"

using namespace Leeds;
using namespace constants;

static void output_ff(const double       W,
                      const std::vector<Subband> &subbands,
                      const unsigned int i,
                      const unsigned int j,
                      const unsigned int f,
                      const unsigned int g);

gsl_spline * FF_table(const double                 Deltak0sqr,
                      const double                 epsilon,
                      const Subband               &isb,
                      const Subband               &jsb,
                      const Subband               &fsb,
                      const Subband               &gsb,
                      const std::valarray<double> &V,
                      const double                 T,
                      const size_t                 nq,
                      const bool                   S_flag);

double PI(const Subband &isb,
          const size_t   q_perp,
          const double   T);

int main(int argc,char *argv[])
{
    double	q_perpsqr4;	/* 4*q_perp*q_perp, checks vailidity of q_perp	*/
    double	W;		/* arbitrary well width, soley for output	*/
    double	Wbar;		/* FD weighted mean of Wijfg			*/
    double	Wijfg;		/* the carrier-carrier scattering rate		*/
    char	filename[9];	/* character string for output filename		*/
    char	p;		/* particle					*/
    bool	ff_flag;	/* form factor flag, output to file if true	*/
    bool	S_flag;		/* screening flag, include screening if true	*/
    FILE	*Fcc;		/* pointer to output file			*/
    FILE	*FccABCD;	/* pointer to weighted mean output file		*/

    /* default values */
    double epsilon=13.18*eps0; // low frequency dielectric constant for GaAs
    ff_flag=false;		/* don't output formfactors	*/
    double m=0.067*me;	       // effective mass [kg]
    p='e';			/* electron			*/
    double T=300;	       // temperature [K]
    W=250e-10;		/* a well width, same as Smet	*/
    S_flag=true;		/* include screening by default	*/

    /* default values for numerical calculations	*/
    size_t nalpha=101; // number of strips in alpha integration
    size_t ntheta=101; // number of strips in theta integration
    size_t nki=101; // number of ki calculations
    size_t nkj=101; // number of strips in |kj| integration
    size_t nq=101;  // number of q_perp values for lookup table

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
                              exit(EXIT_FAILURE);
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
                exit(EXIT_SUCCESS);
        }
        argv++;
        argv++;
        argc--;
        argc--;
    }

    /* calculate step lengths	*/
    const double dalpha=2*pi/((float)nalpha - 1); // step length for alpha integration
    const double dtheta=2*pi/((float)ntheta - 1); // step length for theta integration

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

    if(V.size() != subbands[0].z_array().size())
    {
        std::cerr << "Potential and wavefunction arrays are different sizes: " << V.size() << " and " << subbands[0].z_array().size() << " respectively." << std::endl;
        exit(EXIT_FAILURE);
    }

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
        //   twice the change in KE, see Smet (55)
        double Deltak0sqr = 0;
        if(i+j != f+g)
            Deltak0sqr=4*m*(Ei + Ej - Ef - Eg)/(hBar*hBar);	

        gsl_spline *FF = FF_table(Deltak0sqr, epsilon, isb, jsb, fsb, gsb, V, T,nq,S_flag); // Form factor table
        gsl_interp_accel *acc = gsl_interp_accel_alloc(); // Creates accelerator for interpolation of FF

        /* calculate maximum value of ki & kj and hence kj step length	*/
        const double kimax=sqrt(2*m*(V.max()-Ei))/hBar;
        const double dki=kimax/((float)nki - 1); // step length for loop over ki

        const double kjmax=sqrt(2*m*(V.max()-Ej))/hBar;
        const double dkj=kjmax/((float)nkj - 1); // step length for kj integration

        Wbar=0;			/* initialise integral sum */

        // calculate c-c rate for all ki
        for(unsigned int iki=0;iki<nki;iki++)
        {
            const double ki=dki*(float)iki; // carrier momentum
            Wijfg=0;			/* Initialize for integration	*/

            // integrate over |kj|
            std::valarray<double> Wijfg_integrand_kj(nkj);

            for(unsigned int ikj=0;ikj<nkj;ikj++)
            {
                const double kj=dkj*(float)ikj; // carrier momentum

                // Find Fermi-Dirac occupation at kj
                const double P=jsb.f_FD_k(kj,T);

                // Integral over alpha
                std::valarray<double> Wijfg_integrand_alpha(nalpha);

                for(unsigned int ialpha=0;ialpha<nalpha;ialpha++)
                {
                    const double alpha=dalpha*(float)ialpha; // angle between ki and kj

                    // Compute (vector)kj-(vector)(ki) [QWWAD3, 10.221]
                    const double kij=sqrt(ki*ki+kj*kj-2*ki*kj*cos(alpha));

                    // Now perform innermost integral (over theta)
                    std::valarray<double> Wijfg_integrand_theta(ntheta);

                    for(unsigned int itheta=0;itheta<ntheta;itheta++)
                    {
                        const double theta=dtheta*(float)itheta; // angle between kij and kfg

                        /* calculate argument of sqrt function=4*q_perp*q_perp, see (8.179),
                           to check for imaginary q_perp, if argument is positive, q_perp is
                           real and hence calculate scattering rate, otherwise ignore and move
                           onto next q_perp							*/
                        q_perpsqr4=2*kij*kij+Deltak0sqr-2*kij*sqrt(kij*kij+Deltak0sqr)*cos(theta);

                        if(q_perpsqr4>=0) 
                        {
                            const double q_perp=sqrt(q_perpsqr4)/2; // in-plane momentum, |ki-kf|

                            // Find the form-factor at this wave-vector by looking it up in the
                            // spline we created earlier
                            Wijfg_integrand_theta[itheta] = gsl_spline_eval(FF, q_perp, acc);
                        }
                    } /* end theta */

                    Wijfg_integrand_alpha[ialpha] = integral(Wijfg_integrand_theta, dtheta);
                } /* end alpha */

                Wijfg_integrand_kj[ikj] = integral(Wijfg_integrand_alpha, dalpha) * P * kj;
            } /* end kj   */

            Wijfg = integral(Wijfg_integrand_kj,dkj);

            // Multiply be pre-factor [QWWAD3, 10.233]
            Wijfg *= m*e*e*e*e / (4*pi*hBar*hBar*hBar*(4*4*pi*pi*epsilon*epsilon));

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

        gsl_spline_free(FF);
        gsl_interp_accel_free(acc);
} /* end while over states */

fclose(FccABCD);	/* close weighted mean output file	*/

return EXIT_SUCCESS;
} /* end main */

/** Tabulate the matrix element defined as 
 *    C_if⁺(q,z') = ∫_{z'}^∞ dz ψ_i(z) ψ_f(z)/exp(qz)]
 *  for a given wavevector, with respect to position
 */
std::valarray<double> find_Cif_p(const std::valarray<double>& psi_if, 
                                 const std::valarray<double>& exp_qz,
                                 const std::valarray<double>& z)
{
    const size_t nz = z.size();
    std::valarray<double> Cif_p(nz);
    const double dz=z[1]-z[0];

    Cif_p[nz-1] = psi_if[nz-1] / exp_qz[nz-1] * dz;

    for(int iz = nz-2; iz >=0; iz--)
        Cif_p[iz] = Cif_p[iz+1] + psi_if[iz] / exp_qz[iz] * dz;

    return Cif_p;
}

/** 
 * \brief Tabulate scattering matrix element component.
 *
 * \details defined as:
 *    C_if⁻(q,z') = ∫_{-∞}^{z'} dz ψ_i(z) ψ_f(z) exp(qz)
 *  for a given wavevector, with respect to position
 *
 * Note that the upper limit has to be the point just BEFORE each z'
 * value so that we don't double count
 */
std::valarray<double> find_Cif_m(const std::valarray<double>& psi_if, 
                                 const std::valarray<double>& exp_qz,
                                 const std::valarray<double>& z)
{
    const size_t nz = z.size();
    std::valarray<double> Cif_m(nz);
    const double dz = z[1]-z[0];

    // Seed the first value as zero
    Cif_m[0] = 0;

    // Now, perform a block integration by summing on top of the previous
    // value in the array
    for(unsigned int iz = 1; iz < nz; iz++)
        Cif_m[iz] = Cif_m[iz-1] + psi_if[iz-1] * exp_qz[iz-1] * dz;

    return Cif_m;
}

/** 
 * \brief Create an array of exp(qz) with respect to position
 *
 * \param q[in]       Scattering vector [1/m]
 * \param exp_qz[out] Output array (should be initialised before calling)
 * \param z[in]       Spatial positions [m]
 *
 * \todo  This is also useful for e-e scattering. Push into library
 */
std::valarray<double> find_exp_qz(const double q, const std::valarray<double>& z)
{
    //const double Lp = z.max() - z.min();

    // Use the midpoint of the z array as the origin, so as to minimise the
    // magnitude of the exponential terms
    return exp(q * (z - z[0]));
}

/** Find the matrix element Iif at a given dopant location z'.
 *
 * The matrix element is defined as
 *  I_if(q,z') = ∫dz ψ_i(z) ψ_f(z) exp(-q|z-z'|),
 * where z is the electron location.  The numerical solution
 * can however be speeded up by replacing the modulus function
 * with the sum of two integrals.  We can say that
 *
 *  I_if(q,z') = C_if⁻(q,z')/exp(qz') + C_if⁺(q,z') exp(qz')',
 *
 * Therefore, we have separated the z' dependence from the
 * z dependence of the matrix element.
 */
double Iif(const unsigned int iz0,
           const std::valarray<double>& Cif_p,
           const std::valarray<double>& Cif_m, 
           const std::valarray<double>& exp_qz)
{
    return Cif_m[iz0]/exp_qz[iz0] + Cif_p[iz0]*exp_qz[iz0];
}

/* This function calculates the overlap integral over all four carrier
   states		*/
double A(const double   q_perp,
         const Subband &isb,
         const Subband &jsb,
         const Subband &fsb,
         const Subband &gsb)
{
 const std::valarray<double> z = isb.z_array();
 const size_t nz = z.size();
 const double dz = z[1] - z[0];

 // Convenience labels for wave-functions in each subband
 const std::valarray<double> psi_i = isb.psi_array();
 const std::valarray<double> psi_j = jsb.psi_array();
 const std::valarray<double> psi_f = fsb.psi_array();
 const std::valarray<double> psi_g = gsb.psi_array();

 // Products of wavefunctions can be computed in advance
 const std::valarray<double> psi_if = psi_i * psi_f;
 const std::valarray<double> psi_jg = psi_j * psi_g;

 const std::valarray<double> expTerm = find_exp_qz(q_perp, z);
 const std::valarray<double> Cjg_plus  = find_Cif_p(psi_jg, expTerm, z);
 const std::valarray<double> Cjg_minus = find_Cif_m(psi_jg, expTerm, z);

 std::valarray<double> Aijfg_integrand(nz);

 // Integral of i(=0) and f(=2) over z
 for(unsigned int iz=0;iz<nz;iz++)
 {
     const double Ijg = Iif(iz, Cjg_plus, Cjg_minus, expTerm);
     Aijfg_integrand[iz] = psi_if[iz] * Ijg;
 }

 const double Aijfg = integral(Aijfg_integrand, dz);

 return Aijfg;
}

/**
 * \brief returns the screening factor, referred to by Smet as e_sc
 */
double PI(const Subband &isb,
          const double   q_perp,
          const double   T)
{
    const double m    = isb.get_md_0();    // Effective mass at band-edge [kg]

    // Now perform the integration, equation 44 of Smet [QWWAD3, 10.238]
    double mu=isb.get_E();
    const double dmu=1e-3*e;
    double integral=0;

    double dI; // Intervals in I

    do
    {
        const double ki = isb.k(mu);

        // Find low-temperature polarizability *at this wave-vector*
        //
        // Equation 43 of Smet, QWWAD3, 10.236
        double P0 = m/(pi*hBar*hBar);

        if(q_perp>2*ki)
            P0 -= m/(pi*hBar*hBar)*sqrt(1-gsl_pow_2(2*ki/q_perp));

        dI = P0/(4*kB*T*gsl_pow_2(cosh((isb.get_Ef() - mu)/(2*kB*T))));
        integral += dI*dmu;
        mu+=dmu;
    }while(dI>integral/100); // continue until integral converges

    return integral;
}

/**
 *  \brief Compute the form factor [Aijfg/(esc q)]^2
 */
gsl_spline * FF_table(const double                 Deltak0sqr,
                      const double                 epsilon,
                      const Subband               &isb,
                      const Subband               &jsb,
                      const Subband               &fsb,
                      const Subband               &gsb,
                      const std::valarray<double> &V,
                      const double                 T,
                      const size_t                 nq,
                      const bool                   S_flag)
{
    const double vmax=V.max(); // Maximum potential [J]
    const double kimax=isb.k(vmax - isb.get_E()); // Max value of ki [1/m]
    const double kjmax=jsb.k(vmax - jsb.get_E()); // Max value of kj [1/m]

    // maximum in-plane wave vector
    const double q_perp_max=sqrt(2*gsl_pow_2(kimax+kjmax)+Deltak0sqr+2*(kimax+kjmax)*
                 sqrt(gsl_pow_2(kimax+kjmax)+Deltak0sqr))/2;

    const double dq=q_perp_max/((float)(nq-1));	// interval in q_perp

    std::valarray<double> q_perp(nq);
    std::valarray<double> FF(nq);

    for(unsigned int iq=0;iq<nq;iq++)
    {
        q_perp[iq] = iq*dq;

        // Scattering matrix element (all 4 states)
        const double _Aijfg = A(q_perp[iq], isb, jsb, fsb, gsb);

        double _PI    = 0.0; // Polarizability
        double _Aiiii = 0.0; // Matrix element for lowest subband

        // Allow screening to be turned off
        if(S_flag)
        {
            _PI    = PI(isb, q_perp[iq], T);
            _Aiiii = A(q_perp[iq], isb, isb, isb, isb);
        }

        // Screening permittivity * wave vector
        // Note that the pole at q_perp=0 is avoided as long as screening is included
        const double esc_q = q_perp[iq] + 2*pi*e*e/(4*pi*epsilon) * _PI * _Aiiii;
        FF[iq] = _Aijfg*_Aijfg / (esc_q * esc_q);
    }

    // Fix singularity by "clipping" the top off it:
    if(!S_flag)
        FF[0] = FF[1];

    // Pack the table of FF vs q into a cubic spline
    gsl_spline *q_FF = gsl_spline_alloc(gsl_interp_cspline, nq);
    gsl_spline_init(q_FF, &(q_perp[0]), &(FF[0]), nq);

    return q_FF;
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
 {
     fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);
     exit(EXIT_FAILURE);
 }

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
