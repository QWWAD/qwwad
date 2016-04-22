/*=========================================================
         effv Envelope Function Magnetic Field to v
  =========================================================*/

/* This program adds the Zeeman splitting term to the 
   potential in the file `v.r'.

	Input files

		x.r
		v0.r
	or	v.r

	Output files

		v0.r
		v.r

   Paul Harrison, February 1998				 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include "struct.h"
#include "qwwad/constants.h"
#include "qwwad/maths-helpers.h"

using namespace QWWAD;
using namespace constants;

/**
 * \brief Computes the effective temperature T0(x)
 *
 * \details  (J.A.Gaj, C.Bodin-Deshayes, P.Peyla, J.Cibert, G.Feuillet,
 *           Y.Merle d'Aubigne, R.Romestain and A.Wasiela, Proc. 21st
 *           Int. Conf. Phys. Semiconductors 1936 (1992))
 */
static double Teff(const double x)
{
    const double F = 40.7;
    const double G = 4.6;

    return F*x/(1+G*x);
}

/**
 * \brief Compute the effective spin <Sz(x)>
 *
 * \details Data taken from 
 *
 *  J.A.Gaj, C.Bodin-Deshayes, P.Peyla, J.Cibert, G.Feuillet,
 *  Y.Merle d'Aubigne, R.Romestain and A.Wasiela, Proc. 21st
 *  Int. Conf. Phys. Semiconductors 1936 (1992)               
 *
 *  Note Gaj gives	Delta E(saturation)=6A+6B
 * 
 *  But Delta E(saturation)=x(N0alpha+N0beta)<Sz>
 *
 *  Hence <Sz> follows from Gaj's expression
 */
static double Seff(const double N0alpha,
                   const double N0beta,
                   const double x)
{
    const double A = 2488.0*e*1e-3;
    const double B = -57880.0*e*1e-3;
    const double C = 152.7;
    const double D = -20760.0*e*1e-3;
    const double E = 8.083;
    const double s = (A+B*gsl_pow_2(x)/(1+C*gsl_pow_2(x))+D*x/(1+E*x))/(N0alpha+N0beta);

    return s;
}

/**
 * \brief Brillouin function
 */
static double B_J(const double J,
                  const double MF,
                  const double T,
                  const double T0)
{
    const double y = 2*J*mu_b*MF/(kB*(T+T0));

    return sf_brillouin(J,y);
}

int main(int argc,char *argv[])
{
    /* Define global defaults */
    double MF      = 0.0; // Magnetic field
    char   p       = 'e'; // Particle (e, h, or l)
    char   s       = '+'; // Electron spin state, up or down, +/-
    double T       = 1.8;     // Sample temperature
    double N0alpha = 0.220*e; // Magnetic parameter
    double N0beta  = 0.880*e; // Magnetic parameter
    double J       = 2.5;     // magnetic ion spin

    while((argc>1)&&(argv[1][0]=='-'))
    {
        switch(argv[1][1])
        {
            case 'B':
                MF=atof(argv[2]);            /* read magnetic field in Tesla */
                break;
            case 'p':
                p=*argv[2];
                switch(p)
                {
                    case 'e': break;
                    case 'h': break;
                    case 'l': break;
                    default:  printf("Usage:  efmfv [-p particle (\033[1me\033[0m, h, or l)]\n");
                              exit(0);
                }
                break;
            case 's':
                s=*argv[2];
                break;
            case 'T':
                T=atof(argv[2]);
                break;
            default: 
                printf("Usage: efmfv [-B magnetic field (\033[1m0\033[0mTesla)][-p particle (\033[1me\033[0m, h, or l)]\n");
                printf("             [-s spin-state (\033[1m+\033[0m/-)][-T temperature (\033[1m1.8\033[0mK)]\n");
                exit(0);
        }
        argv++;
        argv++;
        argc--;
        argc--;
    }


    /* Open zero field potential file */
    FILE	*Fv;		/* pointer to input potential file	*/
    if((Fv=fopen("v0.r","r"))==0)
    {
        if((system("cp v.r v0.r"))==0) Fv=fopen("v0.r","r");
        else 
        {printf("Error: Cannot open file 'v0.r' or find file 'v.r'!\n");exit(0);}
    }

    /* Open structure file `x.r'	*/
    FILE	*Fx;		/* pointer to x.r 			*/
    if((Fx=fopen("x.r","r"))==NULL)
    {printf("Error: Cannot open file 'x.r'!\n");exit(0);}

    /* Count number of lines in potential file */
    int n=0;
    while(fscanf(Fv,"%*e %*e")!=EOF)
        n++;
    rewind(Fv);

    /* Allocate memory for potential and alloy concentration structures */

    data11 *V0=(data11 *)calloc(n,sizeof(data11));
    if(V0==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

    data11 *X0=(data11 *)calloc(n,sizeof(data11));
    if(X0==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

    /* Read data into structure */

    int i=0;
    while(fscanf(Fv,"%le %le",&((V0+i)->a),&((V0+i)->b))!=EOF)
    {
        /* Note line below assumes magnetic ion in `x' column of `s.r'	*/

        int n_read = fscanf(Fx,"%le %le %*e",&((X0+i)->a),&((X0+i)->b));
        if(n_read == 3)
            i++;
    }

    fclose(Fv);
    fclose(Fx);

    /* Add field in the form +/- 3A, +/- 3B or +/- B depending on spin */
    for(int i=0;i<n;i++)
    {
        const double x=(X0+i)->b;
        const double T0 = Teff(x); // Effective temperature
        const double Sz = Seff(N0alpha,N0beta,x);

        const double A=x*N0alpha*Sz*B_J(J,MF,T,T0)/6;
        const double B=x*N0beta*Sz*B_J(J,MF,T,T0)/6;

        /* calculate change in potential DeltaV due to field	*/
        double DeltaV=3*A; /* Default value */
        switch(p)
        {
            case 'e': 
                switch(s)
                {
                    case '+':DeltaV=3*A;break;
                    case '-':DeltaV=-3*A;break;
                }
                break;
            case 'h':
                switch(s)
                {
                    case '+':DeltaV=3*B;break;
                    case '-':DeltaV=-3*B;break;
                }
                break;
            case 'l': 
                switch(s)
                {
                    case '+':DeltaV=B;break;
                    case '-':DeltaV=B;break;
                }
                break;
        }

        ((V0+i)->b)+=DeltaV;
    }

    /* Write data to file 'v.r' */
    FILE	*FvF;		/* pointer to v.r with field		*/
    if((FvF=fopen("v.r","w"))==NULL)
        printf("Error: Cannot open file 'v.r' to output data!\n");

    for(i=0;i<n;i++)
        fprintf(FvF,"%20.17le %20.17le\n",(V0+i)->a,(V0+i)->b);

    /* Close files */

    fclose(FvF);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
