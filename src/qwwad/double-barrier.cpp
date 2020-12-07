/**
 * \file   double-barrier.cpp
 * \brief  Calculate transmission coefficient for double barrier structure
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cmath>
#include <armadillo>
#include "constants.h"

using namespace arma;

namespace QWWAD {
using namespace constants;

auto get_transmission_coefficient(const double E,
                                    const double m_w,
                                    const double m_b,
                                    const double V,
                                    const double L1,
                                    const double L2,
                                    const double L3) -> double
{
    // Calculate interfaces
    const double I2=L1;
    const double I3=L1+L2;
    const double I4=L1+L2+L3;

    double Tx = 0; // Transmission coefficient

    // Find wave vector in well and decay constant in barrier
    const double k=sqrt(2*m_w*E)/hBar;
    const double K=sqrt(2*m_b*(V-E))/hBar;

    // Define transfer matrices
    cx_mat M1(2,2);
    M1(0,0) = 1;
    M1(0,1) = 1;
    M1(1,0) = cx_double(0.0,  k/m_w);
    M1(1,1) = cx_double(0.0, -k/m_w);

    cx_mat M2(2,2);
    M2(0,0) = 1;
    M2(0,1) = 1;
    M2(1,0) = +K/m_b;
    M2(1,1) = -K/m_b;

    cx_mat M3(2,2);
    M3(0,0) = exp(+K*I2);
    M3(0,1) = exp(-K*I2);
    M3(1,0) =  K*exp(+K*I2)/m_b;
    M3(1,1) = -K*exp(-K*I2)/m_b;

    cx_mat M4(2,2);
    M4(0,0) = cx_double(cos(k*I2),         +sin(k*I2));
    M4(0,1) = cx_double(cos(k*I2),         -sin(k*I2));
    M4(1,0) = cx_double(-k*sin(+k*I2)/m_w, +k*cos(k*I2)/m_w);
    M4(1,1) = cx_double(+k*sin(-k*I2)/m_w, -k*cos(k*I2)/m_w);

    cx_mat M5(2,2);
    M5(0,0) = cx_double(cos(k*I3),         +sin(k*I3));
    M5(0,1) = cx_double(cos(k*I3),         -sin(k*I3));
    M5(1,0) = cx_double(-k*sin(+k*I3)/m_w, +k*cos(k*I3)/m_w);
    M5(1,1) = cx_double(+k*sin(-k*I3)/m_w, -k*cos(k*I3)/m_w);

    cx_mat M6(2,2);
    M6(0,0) = exp(+K*I3);
    M6(0,1) = exp(-K*I3);
    M6(1,0) =  K*exp(+K*I3)/m_b;
    M6(1,1) = -K*exp(-K*I3)/m_b;

    cx_mat M7(2,2);
    M7(0,0) = exp(+K*I4);
    M7(0,1) = exp(-K*I4);
    M7(1,0) =  K*exp(+K*I4)/m_b;
    M7(1,1) = -K*exp(-K*I4)/m_b;

    cx_mat M8(2,2);
    M8(0,0) = cx_double(cos(k*I4),         +sin(k*I4));
    M8(0,1) = cx_double(cos(k*I4),         -sin(k*I4));
    M8(1,0) = cx_double(-k*sin(+k*I4)/m_w, +k*cos(k*I4)/m_w);
    M8(1,1) = cx_double(+k*sin(-k*I4)/m_w, -k*cos(k*I4)/m_w);

    // Little hack to stop nonsense output when E = 0
    if (gsl_fcmp(E, 0, 1e-9*e) == 1)
    {
        cx_mat M = inv(M1) * M2 * inv(M3) * M4 * inv(M5) * M6 * inv(M7) * M8;
        Tx = 1/(norm(M(1,1))); // Transmission coeff
    }

    return Tx;
}
}// namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
