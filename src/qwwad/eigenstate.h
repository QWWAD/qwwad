/**
 * \file   eigenstate.h
 * \brief  A pure 1D eigenstate of a Hamiltonian
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_EIGENSTATE
#define QWWAD_EIGENSTATE

#include <string>
#include <valarray>
#include <vector>

namespace QWWAD {

/**
 * A pure 1D eigenstate of a Hamiltonian
 */
class Eigenstate {
private:
    double _E; ///< The energy of the state [J]

    std::valarray<double> _z;   ///< Spatial sampling positions [m]
    std::valarray<double> _psi; ///< Wave function [m^{-0.5}]

    double get_total_probability() const;
    void normalise();

public:
    Eigenstate(decltype(_E)   E,
               decltype(_z)   z,
               decltype(_psi) psi);

    inline double get_energy() const {return _E;}
    inline double get_wavefunction_at_index(const unsigned int iz) const {return _psi[iz];}
    inline decltype(_psi) get_wavefunction_samples() const {return _psi;}
    inline decltype(_psi) get_PD() const {return _psi*_psi;}
    inline decltype(_z)   get_position_samples() const {return _z;}

    static double psi_squared_max(const std::vector<Eigenstate> &EVP);

    static std::vector<Eigenstate> read_from_file(const std::string &Eigenval_name,
                                                  const std::string &Eigenvect_prefix,
                                                  const std::string &Eigenvect_ext,
                                                  const double       eigenvalue_scale    = 1.0,
                                                  const bool         ignore_first_column = false);

    static void write_to_file(const std::string             &Eigenval_name,
                              const std::string             &Eigenvect_prefix,
                              const std::string             &Eigenvect_ext,
                              const std::vector<Eigenstate> &states,
                              const bool                     with_num=false);

    // TODO: Should probably be part of an Operator class
    double get_expectation_position() const;

    // TODO: Should probably be part of an Operator class
    static double get_position_matrix_element(const Eigenstate &i,
                                              const Eigenstate &j);
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
