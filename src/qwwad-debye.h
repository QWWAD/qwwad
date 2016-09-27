/**
 * \file   qwwad-debye.h
 *
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  A Debye model of specific heat capacity
 */

#ifndef QWWAD_DEBYE
#define QWWAD_DEBYE

#include <cstddef>

namespace QWWAD
{
class DebyeModel
{
public:
    DebyeModel(const double T_D,
               const double M,
               const size_t natoms);

    double get_internal_energy(const double T) const;
    double get_cp(const double T) const;
    double get_cp_approx(const double T) const;
    double get_cp_low_T(const double T) const;
    double get_cp_high_T() const;

private:
    double _T_D;    ///< Debye temperature [K]
    double _M;      ///< Molar mass [kg/mol]
    size_t _natoms; ///< Number of atoms per molecular unit

    static double find_U(double T, void *params);
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
