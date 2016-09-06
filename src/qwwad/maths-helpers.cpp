/**
 * \file   maths-helpers.cpp
 * \brief  Mathematical utility functions
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include "maths-helpers.h"

namespace QWWAD
{
/**
 * \brief      Interpolates y=f(x) between f(0) and f(1)
 *
 * \param[in]  y0 y-value for f(x=0)
 * \param[in]  y1 y-value for f(x=1)
 * \param[in]  x  x-value for which to find y=f(x)
 * \param[in]  b  An optional "bowing factor"
 *
 * \details  x must lie in the range 0 < x < 1.
 *           The interpolation is computed using
 *
 *           f(x) = (1-x)f(0) + xf(1) + x(1-x)b
 *
 * \return y = f(x)
 */
double lin_interp(const double y0,
                  const double y1,
                  const double x,
                  const double b){
    if(x < 0 or x > 1)
        throw std::domain_error("x value out of range");

    return y0*(1.0-x) + y1*x + b*x*(1.0-x);
}

/**
 * Looks up a y value in a table of the form y=f(x)
 *
 * \param[in] x_values  The x values, in ascending order
 * \param[in] y_values  The y values, in ascending order
 * \param[in] x0        The desired x value for which to find y
 * 
 * \returns The y-value that corresponds to x0
 *
 * \details Linear interpolation is used to find the value more accurately
 *
 * \todo This is a very inefficient way of doing it.  Use splines!
 */
double lookup_y_from_x(const arma::vec &x_values,
                       const arma::vec &y_values,
                       const double x0)
{
    if (x0 > x_values.max() or x0 < x_values.min())
    {
        std::ostringstream oss;
        oss << "Desired x value: " << x0 << " is out of range (" << x_values.min() << ", " << x_values.max() << ").";
        throw std::domain_error(oss.str());
    }

    unsigned int ix=0;
    while (x0 >= x_values[ix])
        ix++;

    if (ix==0)
        return y_values[0];
    else
        return y_values[ix-1] + (y_values[ix] - y_values[ix-1]) * (x0 - x_values[ix-1])/(x_values[ix] - x_values[ix-1]);
}

/**
 * \brief The cotangent of a number
 *
 * \param[in] x The number for which to find the cotangent (radians)
 *
 * \return The cotangent
 */
double cot(const double x)
{
    return 1.0/tan(x);
}

/**
 * \brief The hyperbolic cotangent of a number
 *
 * \param[in] x The number for which to find the hyperbolic cotangent
 *
 * \return The hyperbolic cotangent
 */
double coth(const double x)
{
    return 1.0/tanh(x);
}

/**
 * \brief The Heaviside step function
 *
 * \param[in] x The argument of the step function
 *
 * \return 1 if \f$x > 0\r$, 0 otherwise
 */
unsigned int Theta(const double x)
{
    if(x > 0) return 1;
    else      return 0;
}

/**
 * \brief The Brillouin function y = B_J(x)
 *
 * \param[in] J The order of the Brillouin function
 * \param[in] x The argument of the function
 *
 * \return The value of the function, y
 */
double sf_brillouin(const double J,
                    const double x)
{
    const auto tJ   = 2.0*J;
    const auto tJp1 = tJ + 1.0;
    const auto tJp1_tJ = tJp1/tJ;
    const auto y = tJp1_tJ * coth(tJp1_tJ * x) - 1.0/tJ*coth(x/tJ);
    return y;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
