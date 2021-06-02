/**
 * \file   double-barrier.h
 * \brief  Calculate transmission coefficient for double barrier structure
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_DOUBLE_BARRIER_H
#define QWWAD_DOUBLE_BARRIER_H
namespace QWWAD {
auto get_transmission_coefficient(double E,
                                  double m_w,
                                  double m_b,
                                  double V,
                                  double L1,
                                  double L2,
                                  double L3) -> double;
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
