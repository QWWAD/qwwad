/**
 * \file   schroedinger-solver.h
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_H
#define QWWAD_SCHROEDINGER_SOLVER_H

#include "eigenstate.h"

namespace QWWAD
{
/**
 * Abstract base class for any Schroedinger-equation solver
 *
 * \todo Cache the solutions and only calculate on first request
 */
class SchroedingerSolver
{
private:
    // Options for specifying cut-off energy
    double E_min_ = 0.0;       ///< Lower cut-off energy for solutions [J]
    double E_max_ = 0.0;       ///< Upper cut-off energy for solutions [J]
    bool   E_min_set_ = false; ///< True if lower cut-off energy has been set
    bool   E_max_set_ = false; ///< True if upper cut-off energy has been set

    size_t _nst_max = 0; ///< Maximum number of states to find (0 means all states)

    ///< Set of solutions to the Schroedinger equation
    std::vector<Eigenstate> _solutions;

    arma::vec    _z; ///< Spatial points [m]
    arma::vec    V_; ///< Confining potential [J]

protected:
    [[nodiscard]] auto get_E_min_set() const -> bool {return E_min_set_;}
    [[nodiscard]] auto get_E_max_set() const -> bool {return E_max_set_;}

    inline void set_nst_max(const size_t nst_max) {_nst_max = nst_max;}
    [[nodiscard]] auto get_nst_max() const -> decltype(_nst_max) {return _nst_max;}

    [[nodiscard]] auto get_E_search_min() const -> double;
    [[nodiscard]] auto get_E_search_max() const -> double;

    [[nodiscard]] auto energy_above_range(const double E) const -> bool;
    [[nodiscard]] auto energy_below_range(const double E) const -> bool;

    /**
     * \brief Calculate all eigenstates
     *
     * \details This needs to be implemented by any derived classes
     *
     * \return A set of all the eigenstates in the system
     */
    virtual auto calculate() -> decltype(_solutions) = 0;

    /**
     * \brief Set the potential at each point
     *
     * \details Subclasses should always call this in their constructor
     *
     * \todo Make sub-classes define a vfunc to populate the V & z
     *       vectors
     *
     * \param[in] V The potential at each point [J]
     */
    inline void set_V(const decltype(V_) &V) {V_ = V;}

    /**
     * \brief Set the position at each point
     *
     * \details Subclasses should always call this in their constructor
     *
     * \todo Make sub-classes define a vfunc to populate the V & z
     *       vectors
     *
     * \param[in] z The position at each point [z]
     */
    inline void set_z(const decltype(_z) &z) {_z = z;}

    void refresh_solutions();

public:
    auto get_solutions(const bool convert_to_meV=false) -> std::vector<Eigenstate>;

    [[nodiscard]] auto get_z() const -> decltype(_z);
    [[nodiscard]] auto get_V() const -> decltype(V_);

    virtual auto get_name() -> std::string = 0;

    void set_E_min(const double E_min);
    void set_E_max(const double E_max);

    /**
     * \brief Turn off filtering of solutions by energy
     */
    inline void unset_E_cutoff()
    {
        E_min_set_ = false;
        E_max_set_ = false;
    }
};
} // namespace QWWAD

#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
