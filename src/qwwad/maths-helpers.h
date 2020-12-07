/** 
 * \file   maths-helpers.h
 * \brief  Declarations for mathematical utility functions
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_MATHS_HELPERS_H
#define QWWAD_MATHS_HELPERS_H

#include <cstddef>
#include <stdexcept>
#include <sstream>

#include <gsl/gsl_math.h>

#include <armadillo>

namespace QWWAD
{
/**
 * \brief Integrate using Simpson's rule
 *
 * \param[in] y  Samples of the function to be integrated
 * \param[in] dx Spatial step between samples
 *
 * \details The number of samples must be odd, and >= 3
 *
 * \returns The integral
 */
template <class complex_type, class real_type>
auto simps(const arma::Col<complex_type>& y, const real_type dx) -> complex_type
{
    const size_t n = y.size();

    if(n < 3)
        throw std::runtime_error("Not enough points for Simpson's rule");

    if(GSL_IS_EVEN(n))
    {
        std::ostringstream oss;
        oss << "Simpson's rule needs odd number of points: " << n << " received.";
        throw std::length_error(oss.str());
    }

    complex_type ans=0;
    for(unsigned int i=0; i<n-2; i+=2)
        ans += y[i] + 4.0*y[i+1] + y[i+2];

    ans *= dx/3.0;
    return ans;
}

/**
 * \brief Integrate using the trapezium rule
 *
 * \param[in] y  Samples of the function to be integrated
 * \param[in] dx Spatial step between samples
 *
 * \details The number of samples must be >= 2
 *
 * \returns The integral
 */
template <class complex_type, class real_type>
auto trapz(const arma::Col<complex_type>& y, const real_type dx) -> complex_type
{
    const size_t n = y.size();

    if(n < 2)
        throw std::runtime_error("Need at least two points for trapezium rule");

    complex_type ans=0;

    for(unsigned int i=0; i<n-1; i++)
        ans += (y[i] + y[i+1])/2.0;

    ans *= dx;

    return ans;
}

/**
 * \brief Compute a numerical integral using a sensible solver
 *
 * \param[in] y  Samples of the function to be integrated
 * \param[in] dx Spatial step between samples
 *
 * \details If the number of samples is odd, and >=3, a Simpson's rule solver
 *          will be used. This gives higher precision than the fallback
 *          trapezium rule solver. The method is selected automatically, based
 *          on the number of samples in the function. This should only result
 *          in a small time overhead, so it should generally be fine to use
 *          this function instead of directly calling trapz or simps.
 *
 *          Complex functions of real variables are supported
 *
 * \todo Replace with arma::trapz when available in all distros (Armadillo >= 6.5)
 */
template <class complex_type, class real_type>
auto integral(const arma::Col<complex_type>& y, const real_type dx) -> complex_type
{
    const size_t n = y.size();

    if(n < 2)
        throw std::runtime_error("Need at least two points for numerical integration.");

    if(GSL_IS_ODD(n) && n >= 3)
        return simps(y, dx);
    else
        return trapz(y, dx);
}

auto lookup_y_from_x(const arma::vec &x_values,
                       const arma::vec &y_values,
                       const double     x0) -> double;

auto lin_interp(const double y0,
                  const double y1,
                  const double x,
                  const double b=0) -> double;

auto cot(const double x) -> double;

auto coth(const double x) -> double;

auto Theta(const double x) -> unsigned int;

auto sf_brillouin(const double J,
                    const double x) -> double;
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
