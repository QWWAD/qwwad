/** 
 * \file   qclsim-maths.h
 * \brief  Declarations for mathematical utility functions
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date   2012-03-14
 */

#ifndef QCLSIM_MATHS_H
#define QCLSIM_MATHS_H

#include <cstddef>
#include <valarray>

namespace Leeds
{

double simps(const std::valarray<double>& y, const double dx);

double trapz(const std::valarray<double>& y, const double dx);

double lookup_y_from_x(const std::valarray<double>& x_values,
        const std::valarray<double>& y_values,
        const double x0);

double lin_interp(const double y0,
                  const double y1,
                  const double x,
                  const double b=0);

double cot(const double x);

/**
 * A numerical solver for Laplace transforms
 *
 * \details This uses the Stehfest algorithm to compute transforms as described
 *          at: http://www.codeproject.com/Articles/25189/Numerical-Laplace-Transforms-and-Inverse-Transform
 *
 * \todo Replace this with a call to an external library if/when one becomes
 *       available.  This isn't widely used in QCLsim, so we probably don't
 *       want to maintain it ourselves!
 */
class Laplace
{
private:
    /**
     * Number of summation terms to use in algorithm.
     *
     * This must be an even number.
     *
     * TODO: Make this configurable?
     */
    static const size_t N = 20;

    /**
     * Stehfest algorithm coefficients
     *
     * More terms gives a more accurate solution, but is slower
     */
    std::valarray<double> V;
    
public:
    Laplace();

    double inverse_transform(double (*F)(const double), const double t) const;
};

} // namespace Leeds
#endif

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
