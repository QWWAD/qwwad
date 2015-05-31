/**
 * \file   subband.cpp
 * \brief  A subband in a 2D system
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include "subband.h"
#include "file-io.h"
#include "maths-helpers.h"
#include "constants.h"
#include "fermi.h"

namespace QWWAD
{
using namespace constants;

/**
 * \brief Create a parabolic subband
 *
 * \param[in] ground_state The quantum state at the edge of the subband
 * \param[in] m            The effective mass for the subband [kg]
 */
Subband::Subband(const Eigenstate &ground_state,
                 const double      m) :
    _ground_state(ground_state),
    _m(m),
    _alpha(0.0),
    _V(0.0),
    _dist_known(false),
    _Ef(ground_state.get_energy()),
    _Te(0.0),
    _N(0.0)
{}

/**
 * \brief Create a non-parabolic subband
 *
 * \param[in] ground_state The quantum state at the edge of the subband
 * \param[in] m            The effective mass for the subband [kg]
 * \param[in] alpha        The nonparabolicity parameter [1/J]
 * \param[in] V            The edge of the bulk band that contains this subband [J]
 */
Subband::Subband(const Eigenstate &ground_state,
                 const double      m,
                 const double      alpha,
                 const double      V) :
    _ground_state(ground_state),
    _m(m),
    _alpha(alpha),
    _V(V),
    _dist_known(false),
    _Ef(ground_state.get_energy()),
    _Te(0.0),
    _N(0.0)
{}

/**
 * \brief Sets the carrier distribution function in the subband
 *
 * \param[in] Ef Quasi-Fermi energy [J]
 * \param[in] Te Carrier temperature [K]
 *
 * \details A Fermi-Dirac distribution is assumed
 */
void Subband::set_distribution_from_Ef_Te(const double Ef,
                                          const double Te)
{
    if(Te <= 0.0)
        throw "Carrier temperature must be positive";

    _dist_known = true;
    _Ef         = Ef;
    _Te         = Te;
}

/**
 * \brief Find Fermi wave-vector
 *
 * \details Uses QWWAD3, Eq. 10.237
 *          Note that the degeneracy is set as 1
 *
 * \returns Fermi wave-vector [1/m]
 */
double Subband::get_k_fermi() const
{
    if(!_dist_known)
        throw std::runtime_error("Distribution has not been set");

    const auto N = get_total_population();

    return sqrt(2.0*pi*N);
}

/**
 * Reads a set of subbands from data files (not including nonparabolic dispersion)
 *
 * \param[in] energy_input_path     Path for energy and wavefunction files
 * \param[in] wf_input_prefix       Prefix for wavefunction filenames
 * \param[in] wf_input_ext          Extension for wavefunction filenames
 * \param[in] m_d_filename          Name of data file containing density-of-states
 *                                  effective mass.
 */
std::vector<Subband> Subband::read_from_file(const std::string& energy_input_path,
                                             const std::string& wf_input_prefix,
                                             const std::string& wf_input_ext,
                                             const std::string& m_d_filename)
{
    // Read ground state data
    const auto ground_state = Eigenstate::read_from_file(energy_input_path,
                                                         wf_input_prefix,
                                                         wf_input_ext,
                                                         1000.0/e,
                                                         true);
    
    const size_t nst = ground_state.size();

    if(nst == 0)
        throw std::runtime_error("No states found in file");

    // Read effective mass table
    std::valarray<double> z; // m
    std::valarray<double> m_d_z; // kg
    read_table(m_d_filename.c_str(), z, m_d_z);

    // Copy subband data to vector
    std::vector<Subband> subbands;

    for(unsigned int ist = 0; ist < nst; ist++)
    {
        // Find the expectation-value of in-plane effective mass using QWWAD 4, 12.22
        // TODO: Note that the value used for transitions between a PAIR of subbands
        //       should use the inverse-mass matrix element; not this expectation value
        const auto z  = ground_state[ist].get_position_samples();
        const auto dz = z[1] - z[0];

        const std::valarray<double> mass_integrand = 1.0 / m_d_z * ground_state[ist].get_PD();
        const auto mass = 1.0 / integral(mass_integrand, dz);

        subbands.push_back(Subband(ground_state[ist], mass));
    }

    return subbands;
}

/**
 * Reads a set of subbands from data files (not including nonparabolic dispersion)
 *
 * \param[in] energy_input_path     Path for energy and wavefunction files
 * \param[in] wf_input_prefix       Prefix for wavefunction filenames
 * \param[in] wf_input_ext          Extension for wavefunction filenames
 * \param[in] m_d                   Density-of-states effective mass [kg]
 */
std::vector<Subband> Subband::read_from_file(const std::string& energy_input_path,
                                             const std::string& wf_input_prefix,
                                             const std::string& wf_input_ext,
                                             const double       m_d)
{
    const auto ground_state = Eigenstate::read_from_file(energy_input_path,
                                                         wf_input_prefix,
                                                         wf_input_ext,
                                                         1000.0/e,
                                                         true);

    const size_t nst = ground_state.size();

    if(nst == 0)
        throw std::runtime_error("No states found in file");

    // Copy subband data to vector
    std::vector<Subband> subbands;

    for (unsigned int ist = 0; ist < nst; ist++)
        subbands.push_back(Subband(ground_state[ist], m_d));

    return subbands;
}

/**
 * Reads a set of subbands from data files (including nonparabolic dispersion)
 *
 * \param[in] energy_input_path     Path for energy and wavefunction files
 * \param[in] wf_input_prefix       Prefix for wavefunction filenames
 * \param[in] wf_input_ext          Extension for wavefunction filenames
 * \param[in] m_d_filename          Name of data file containing density-of-states effective mass.
 * \param[in] alphad_filename       Name of data file containing dispersion nonparabolicity profile
 * \param[in] potential_filename    Name of data file containing band edge potential.
 */
std::vector<Subband> Subband::read_from_file(const std::string& energy_input_path,
                                             const std::string& wf_input_prefix,
                                             const std::string& wf_input_ext,
                                             const std::string& m_filename,
                                             const std::string& alpha_filename,
                                             const std::string& potential_filename)
{
    // Read ground state data
    const auto ground_state = Eigenstate::read_from_file(energy_input_path,
                                                         wf_input_prefix,
                                                         wf_input_ext,
                                                         1000.0/e,
                                                         true);

    const size_t nst = ground_state.size();

    if(nst == 0)
        throw std::runtime_error("No states found in file");

    // Read effective mass table
    std::valarray<double> z; // m
    std::valarray<double> m; // kg
    read_table(m_filename.c_str(), z, m);

    // Read non-parabolicity parameter
    std::valarray<double> alpha; // [1/J]
    read_table(alpha_filename.c_str(), z, alpha);

    // Read band-edge potential
    std::valarray<double> V; // [J]
    read_table(potential_filename.c_str(), z, V);

    // Copy subband data to vector
    std::vector<Subband> subbands;

    for(unsigned int ist = 0; ist < nst; ist++)
    {
        // Find the expectation-value of in-plane effective mass using QWWAD 4, 12.22
        // TODO: Note that the value used for transitions between a PAIR of subbands
        //       should use the inverse-mass matrix element; not this expectation value
        const auto z  = ground_state[ist].get_position_samples();
        const auto dz = z[1] - z[0];

        const std::valarray<double> mass_integrand = 1.0 / m * ground_state[ist].get_PD();
        const auto mass = 1.0 / integral(mass_integrand, dz);

        // Get "expectation values" for potential and non-parabolicity too
        // TODO: Check whether these make sense!
        const std::valarray<double> V_integrand = V * ground_state[ist].get_PD();
        const auto V_exp = integral(V_integrand, dz);

        const std::valarray<double> alpha_integrand = alpha * ground_state[ist].get_PD();
        const auto alpha_exp = integral(alpha_integrand, dz);

        Subband sb(ground_state[ist],
                   mass,
                   alpha_exp,
                   V_exp);

        subbands.push_back(sb);
    }
    
    return subbands;
}

/**
 * Reads a set of subbands from data files (including nonparabolic dispersion)
 *
 * \param[in] energy_input_path     Path for energy and wavefunction files
 * \param[in] wf_input_prefix       Prefix for wavefunction filenames
 * \param[in] wf_input_ext          Extension for wavefunction filenames
 * \param[in] m_d                   Density-of-states effective mass [kg]
 * \param[in] alphad                Dispersion nonparabolicity profile [1/J]
 * \param[in] V                     Band edge potential [J]
 */
std::vector<Subband> Subband::read_from_file(const std::string& energy_input_path,
                                             const std::string& wf_input_prefix,
                                             const std::string& wf_input_ext,
                                             const double       m_d,
                                             const double       alphad,
                                             const double       V)
{
    // Read ground-state data
    const auto ground_state = Eigenstate::read_from_file(energy_input_path,
                                                         wf_input_prefix,
                                                         wf_input_ext,
                                                         1000.0/e,
                                                         true);

    const size_t nst = ground_state.size();

    if(nst == 0)
        throw std::runtime_error("No states found in file");

    // Copy subband data to vector
    std::vector<Subband> subbands; // Output structure
    for (unsigned int ist = 0; ist < nst; ist++)
        subbands.push_back(Subband(ground_state[ist], m_d, alphad, V));

    return subbands;
}

/**
 * \begin Find the kinetic energy at a given wavevector
 *
 * \param[in] k In-plane wave vector [1/m]
 *
 * \return Kinetic energy [J]
 */
double Subband::get_Ek_at_k(double k) const
{
    double Ek;

    // Check if subband is initialised as being nonparabolic
    if(_alpha == 0.0)
        Ek = hBar*hBar*k*k/(2.0*_m);
    else
    {
        const auto En = get_E_min();
        const auto b       = 1.0 + _alpha*(En - _V);
        const auto four_ac = 4.0*_alpha*(-hBar*hBar*k*k)/(2.0*_m);

        // Check solveable
        if(four_ac > b*b)
        {
            std::ostringstream oss;
            oss << "No real energy solution exists at wavevector k = " << k*1.0e-9 << " nm^{-1}.";
            throw std::domain_error(oss.str());
        }

        const auto root = sqrt(b*b - four_ac);
        if(root >= b)
            Ek = (-b + root)/(2.0*_alpha);
        else
        {
            std::ostringstream oss;
            oss << "Negative energy found at wavevector k = " << k*1.0e-9 << " nm^{-1}.";
            throw std::domain_error(oss.str());
        }
    }

    return Ek;
}

/**
 * Return the wavevector at some energy above subband minima
 *
 * \param[in] Ek Kinetic energy [J]
 * 
 * \return Wavevector [m^{-1}]
 */
double Subband::get_k_at_Ek(const double Ek) const
{
    if(Ek < 0.0)
    {
        std::ostringstream oss;
        oss << "Cannot find wavevector at negative kinetic energy, Ek = " << Ek/e*1000 << " meV.";
        throw std::domain_error(oss.str());
    }

    // Find the energy-dependent effective mass, including nonparabolicity
    const auto E_total = Ek + get_E_min();
    const auto m = get_effective_mass(E_total);

    const auto k = sqrt(Ek*2.0*m)/hBar;
    
    return k;
} 

/**
 * \brief Return 2D density of states for subband
 *
 * \param[in] E Absolute energy of state [J]
 */
double Subband::get_density_of_states(const double E) const
{
    double rho = 0.0;

    // Density of states is zero unless we're above the
    // subband edge
    if(E > get_E_min())
    {
        const auto m = get_effective_mass_dos(E);

        // QWWAD 4, Eq. 2.63
        rho = m/(pi*hBar*hBar);
    }

    return rho;
}

/**
 * \brief Find the occupation probability of a state at a given energy
 *
 * \param[in] E_total The energy of the state [J]
 *
 * \details Note that E is the TOTAL energy of the state, specified on the same
 *          absolute scale as the Fermi energy and the subband minimum
 */
double Subband::get_occupation_at_E_total(const double E) const
{
    if(!_dist_known)
        throw std::runtime_error("Distribution has not been set");

    return f_FD(_Ef, E, _Te);
}

/**
 * \brief Find the occupation probability of a state at a given wave vector
 *
 * \param[in] k  The in-plane wave-vector of the state [1/m]
 */
double Subband::get_occupation_at_k(const double k) const
{
    if(!_dist_known)
        throw std::runtime_error("Distribution has not been set");

    const auto E = get_E_total_at_k(k);
    return get_occupation_at_E_total(E);
}

/**
 * \brief Get the total population of the subband
 *
 * \returns Population [m^{-2}]
 */
double Subband::get_total_population() const
{
    if(!_dist_known)
        throw std::runtime_error("Distribution has not been set");

    const auto E = get_E_min();
    const auto N = find_pop(E, _Ef, _m, _Te, _alpha, _V);
    return N;
}

/**
 * \brief find the wave-vector below which 99% of population lies
 *
 * \todo This is a pretty crude approximation that assumes Maxwell-Boltzmann
 *       statistics (and adds a correction for subbands with low populations).
 *       Would be better to integrate over the band to find the correct solution.
 */
double Subband::get_k_max(const double Te) const
{
    if(!_dist_known)
        throw std::runtime_error("Distribution has not been set");

    double Ek_max = 5.0*kB*Te;

    if(get_E_min() < get_Ef())
        Ek_max += get_Ef();

    return get_k_at_Ek(Ek_max);
}

/**
 * \brief Return effective mass at a given energy [kg]
 *
 * \param[in] E The absolute energy of the state [J]
 *
 * \details This is the mass that should be used for finding quantised states
 *          in a system. Note that the energy is specified on the same absolute
 *          scale as the subband minimum.
 */
double Subband::get_effective_mass(const double E) const
{
    const auto m = _m*(1.0 + _alpha*(E-_V));

    return m;
}

/**
 * \brief Return the density-of-states effective mass [kg]
 *
 * \param[in] E The absolute energy of the state [J]
 *
 * \details This is the mass that should be used for finding density of states
 *          in a system. Note that the energy is specified on the same absolute
 *          scale as the subband minimum.
 */
double Subband::get_effective_mass_dos(const double E) const
{
    // QWWAD 4, Eq. 2.65
    const auto m = _m*(1.0 + 2.0*_alpha*(E-_V));

    return m;
}

/**
 * \brief Find the population of the subband at wavevector k
 *
 * \param[in] k Wave vector [1/m]
 */
double Subband::get_population_at_k(const double k) const
{
    const auto E    = get_E_total_at_k(k);
    const auto rho  = get_density_of_states(E);
    const auto f_FD = get_occupation_at_k(k);

    return rho*f_FD;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
