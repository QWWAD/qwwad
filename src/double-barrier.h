/**
 * \file   double-barrier.h
 * \brief  Calculate transmission coefficient for double barrier structure
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef DOUBLE_BARRIER_H
#define DOUBLE_BARRIER_H
namespace Leeds {
double get_transmission_coefficient(const double E,
                                    const double m_w,
                                    const double m_b,
                                    const double V,
                                    const double L1,
                                    const double L2,
                                    const double L3);
};

#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
