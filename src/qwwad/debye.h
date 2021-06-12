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
    DebyeModel(double T_D,
               double M,
               size_t natoms);

    [[nodiscard]] auto get_internal_energy(double T) const -> double;
    [[nodiscard]] auto get_cp(double T) const -> double;
    [[nodiscard]] auto get_cp_approx(double T) const -> double;
    [[nodiscard]] auto get_cp_low_T(double T) const -> double;
    [[nodiscard]] auto get_cp_high_T() const -> double;

private:
    double _T_D;    ///< Debye temperature [K]
    double _M;      ///< Molar mass [kg/mol]
    size_t _natoms; ///< Number of atoms per molecular unit

    static auto find_U(double T, void *params) -> double;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
