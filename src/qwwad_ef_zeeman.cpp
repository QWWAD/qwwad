/**
 * \file   qwwad_ef_zeeman.cpp
 * \brief  Adds Zeeman splitting term to an input potential
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavania@leeds.ac.uk>
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include "struct.h"
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
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
                    default:  printf("Usage:  qwwad_ef_zeeman [-p particle (\033[1me\033[0m, h, or l)]\n");
                              exit(EXIT_FAILURE);
                }
                break;
            case 's':
                s=*argv[2];
                break;
            case 'T':
                T=atof(argv[2]);
                break;
            default: 
                printf("Usage: qwwad_ef_zeeman [-B magnetic field (\033[1m0\033[0mTesla)][-p particle (\033[1me\033[0m, h, or l)]\n");
                printf("             [-s spin-state (\033[1m+\033[0m/-)][-T temperature (\033[1m1.8\033[0mK)]\n");
                exit(0);
        }
        argv++;
        argv++;
        argc--;
        argc--;
    }

    // Open baseline potential and alloy files
    std::valarray<double> z;  // Spatial locations [m]
    std::valarray<double> Vb; // Band-edge potential [J]
    read_table("v_b.r", z, Vb);

    const auto nz = z.size(); // Number of spatial points

    // TODO: Check that alloy file uses same coordinates
    std::valarray<double> z_tmp;
    std::valarray<double> x;     // Alloy fraction at each point
    read_table("x.r", z_tmp, x);

    std::valarray<double> V_zeeman(nz);

    /* Add field in the form +/- 3A, +/- 3B or +/- B depending on spin */
    for(unsigned int iz=0; iz<nz; iz++)
    {
        const double T0 = Teff(x[iz]); // Effective temperature
        const double Sz = Seff(N0alpha,N0beta,x[iz]);

        const double A=x[iz]*N0alpha*Sz*B_J(J,MF,T,T0)/6;
        const double B=x[iz]*N0beta*Sz*B_J(J,MF,T,T0)/6;

        // calculate change in potential due to magnetic field
        switch(p)
        {
            case 'e': 
                switch(s)
                {
                    case '+': V_zeeman[iz] =  3*A;break;
                    case '-': V_zeeman[iz] = -3*A;break;
                }
                break;
            case 'h':
                switch(s)
                {
                    case '+': V_zeeman[iz] =  3*B;break;
                    case '-': V_zeeman[iz] = -3*B;break;
                }
                break;
            case 'l': 
                switch(s)
                {
                    case '+': V_zeeman[iz] = B;break;
                    case '-': V_zeeman[iz] = B;break;
                }
                break;
        }
    }

    const std::valarray<double> V_total = Vb + V_zeeman;

    // Write data to file
    write_table("v.r", z, V_total);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
