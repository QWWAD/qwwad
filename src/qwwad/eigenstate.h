/**
 * \file   eigenstate.h
 * \brief  A pure 1D eigenstate of a Hamiltonian
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_EIGENSTATE
#define QWWAD_EIGENSTATE

#include <string>
#include <armadillo>

namespace QWWAD {

/**
 * A pure 1D eigenstate of a Hamiltonian
 */
class Eigenstate {
private:
    double E_; ///< The energy of the state [J]

    arma::vec _z;      ///< Spatial sampling positions [m]
    arma::cx_vec _psi; ///< Wave function [m^{-0.5}]

    [[nodiscard]] auto get_total_probability() const -> double;
    void normalise();

public:
    explicit Eigenstate(decltype(E_)   E,
                        decltype(_z)   z,
                        decltype(_psi) psi);

    [[nodiscard]] inline auto get_energy()                                     const {return E_;}
    [[nodiscard]] inline auto get_wavefunction_at_index(const unsigned int iz) const {return _psi[iz];}
    [[nodiscard]] inline auto get_wavefunction_samples()                       const {return _psi;}
    [[nodiscard]] inline auto get_position_samples()                           const {return _z;}
    [[nodiscard]] inline auto get_PD()                                         const -> arma::vec {return square(abs(_psi));}

    static auto psi_squared_max(const std::vector<Eigenstate> &states) -> double;

    static auto read_from_file(const std::string &Eigenval_name,
                               const std::string &Eigenvect_prefix,
                               const std::string &Eigenvect_ext,
                               double             eigenvalue_scale    = 1.0,
                               bool               ignore_first_column = false) -> std::vector<Eigenstate>;

    static void write_to_file(const std::string             &Eigenval_name,
                              const std::string             &Eigenvect_prefix,
                              const std::string             &Eigenvect_ext,
                              const std::vector<Eigenstate> &states,
                              bool                           with_num=false);

    // TODO: Should probably be part of an Operator class
    [[nodiscard]] auto get_expectation_position() const -> double;

    // TODO: Should probably be part of an Operator class
    static auto get_position_matrix_element(const Eigenstate &i,
                                            const Eigenstate &j) -> double;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
