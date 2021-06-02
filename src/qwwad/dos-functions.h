#ifndef DOS_FUNCTIONS_H
#define DOS_FUNCTIONS_H
#include <valarray>

namespace QWWAD
{
auto calculate_dos_3D(double mass,
                      double energy,
                      double V       = 0,
                      double alpha   = 0) noexcept -> double;

auto calculate_dos_2D(double                       mass,
                      double                       E_carrier,
                      const std::valarray<double> &E_subbands,
                      double                       V = 0,
                      double                       alpha = 0) noexcept -> double;

auto calculate_dos_1D(double                       mass,
                      double                       E_carrier,
                      const std::valarray<double> &E_subbands) noexcept -> double;
}
#endif // DOS_FUNCTIONS_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
