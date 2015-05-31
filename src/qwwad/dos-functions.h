#ifndef DOS_FUNCTIONS_H
#define DOS_FUNCTIONS_H
#include <valarray>

namespace QWWAD
{
double calculate_dos_3D(const double mass,
                        const double energy,
                        const double V     = 0,
                        const double alpha = 0);

double calculate_dos_2D(const double                 mass,
                        const double                 E_carrier,
                        const std::valarray<double> &E_subbands,
                        const double                 V = 0,
                        const double                 alpha = 0);

double calculate_dos_1D(const double                 mass,
                        const double                 E_carrier,
                        const std::valarray<double> &E_subbands);
}
#endif // DOS_FUNCTIONS_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
