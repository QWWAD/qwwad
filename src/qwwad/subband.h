/**
 * \file   qclsim-subband.h
 * \brief  A subband in a 2D system
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date   2013-01-10
 */

#ifndef QCLSIM_SUBBAND_H
#define QCLSIM_SUBBAND_H

#if HAVE_CONFIG_H
# include "config.h"
#endif //HAVE_CONFIG_H

#include <string>
#include "linear-algebra.h"

namespace QWWAD
{
class Subband
{
private:
    State  _ground_state;      ///< State at bottom of subband

    /// TODO: The spatial profile should really be part of the State class
    std::valarray<double> _z; ///< Spatial profile of subband [m]

    double _m;                ///< Effective mass at subband minimum (for dispersion) [kg]
    double _alpha;            ///< In-plane nonparabolicity parameter [1/J]
    double _V;                ///< Conduction band edge [J]

    // Carrier distribution parameters
    bool   _dist_known;       ///< True if the carrier distribution is set
    double _Ef;               ///< Quasi-Fermi energy [J]
    double _Te;               ///< Temperature of carrier distribution [K]
    double _N;                ///< TOTAL sheet-density of carriers [m^{-2}]

public:
    Subband(State                 ground_state,
            double                m,
            std::valarray<double> z);
        
    Subband(State                 ground_state,
            double                m,
            std::valarray<double> z,
            double                alpha,
            double                V);

    void set_distribution_from_Ef_Te(const double Ef,
                                     const double Te);

    inline State get_ground() const {return _ground_state;}

    /**
     * \brief find the total length of the spatial region [m]
     *
     * \todo This should be part of State class
     */
    inline std::valarray<double>       z_array()    const {return _z;}
    inline double                      get_dz()     const {return _z[1]-_z[0];}
    inline double                      get_length() const {return _z[_z.size()-1]-_z[0];}

    /** Find expectation position for the ground state [m] */
    inline double                      get_z_av_0() const {return z_av(_ground_state, _z);}

    inline double                      get_Ef()     const {return _Ef;}

    /**
     * \brief Find energy of subband edge
     */
    inline double                      get_E_min()  const {return _ground_state.get_E();}

    double                             get_total_population()    const;

    inline std::valarray<double>       psi_array()  const {return _ground_state.psi_array();}
    inline double                      get_condband_edge() const {return _V;}

    double                             get_k_fermi() const;

    static std::vector<Subband> read_from_file(const std::string &energy_input_path,
                                               const std::string &wf_input_prefix,
                                               const std::string &wf_input_ext,
                                               const std::string &m_filename);

    static std::vector<Subband> read_from_file(const std::string &energy_input_path,
                                               const std::string &wf_input_prefix,
                                               const std::string &wf_input_ext,
                                               const double       m);

    static std::vector<Subband> read_from_file(const std::string &energy_input_path,
                                               const std::string &wf_input_prefix,
                                               const std::string &wf_input_ext,
                                               const std::string &m_filename,
                                               const std::string &alpha_filename,
                                               const std::string &potential_filename);

    static std::vector<Subband> read_from_file(const std::string &energy_input_path,
                                               const std::string &wf_input_prefix,
                                               const std::string &wf_input_ext,
                                               const double       m,
                                               const double       alpha,
                                               const double       V);

    double get_Ek_at_k(const double k) const;
    double get_k_at_Ek(const double Ek) const;
    double get_k_max(const double Te) const;

    /// Return total energy of carrier at a given wave-vector
    inline double get_E_total_at_k(const double k) const {return get_E_min() + get_Ek_at_k(k);}

    decltype(_m) get_effective_mass    (const double E = 0.0) const;
    decltype(_m) get_effective_mass_dos(const double E = 0.0) const;

    /// Return nonparabolicity parameter
    inline double get_alpha() const {return _alpha;}

    double get_density_of_states(const double E = 0.0) const;

    double get_occupation_at_E_total(const double E) const;
       
    double get_occupation_at_k(const double k) const;

    double get_population_at_k(const double k) const;
};
} // namespace
#endif // QCLSIM_SUBBAND_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
