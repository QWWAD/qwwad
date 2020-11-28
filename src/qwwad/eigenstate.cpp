#include "eigenstate.h"
#include <sstream>
#include <utility>

#include "maths-helpers.h"
#include "file-io.h"

namespace QWWAD {
Eigenstate::Eigenstate(decltype(_E)   E,
                       decltype(_z)   z,
                       decltype(_psi) psi) :
    _E(E),
    _z(std::move(z)),
    _psi(std::move(psi))
{
    normalise();
}

/**
 * \brief Find the total probability of the state over all space
 */
double Eigenstate::get_total_probability() const
{
    const auto PD = get_PD();
    const auto dz = _z[1] - _z[0];
    const auto probability = integral(PD, dz);

    return probability;
}

/**
 * \brief Normalise the wave function amplitude
 *
 * \details This normalises the amplitude so that the total
 *          probability = 1
 */
void Eigenstate::normalise()
{
    const auto P = get_total_probability();

    // Normalisation factor
    const auto A = sqrt(P);

    _psi *= 1.0/A;
}

/** 
 * \brief Read a set of eigenstates from file.
 *
 * \param[in]   Eigenval_name       The name of the file which holds the eigenvalues in a single column
 * \param[in]   Eigenvect_prefix    Prefix of the files holding the eigenvectors
 * \param[in]   Eigenvect_ext       Extension of the files holding the eigenvectors
 * \param[in]   eigenvalue_scale    Value by which all eigenvalues will be divided upon
 *                                  read
 * \param[in]   ignore_first_column True if first column of eigenvalue file should be
 *                                  ignored
 * 
 * \returns  A vector containing the eigenstates
 *
 * \details Reads in eigenstates from files into a vector
 */
std::vector<Eigenstate>
Eigenstate::read_from_file(const std::string &Eigenval_name,
                           const std::string &Eigenvect_prefix,
                           const std::string &Eigenvect_ext,
                           const double       eigenvalue_scale,
                           const bool         ignore_first_column)
{
    std::vector<Eigenstate> states;

    // Read eigenvalues into tempory memory
    arma::vec E_temp;

    if(ignore_first_column)
    {
        arma::vec indices;
        read_table(Eigenval_name.c_str(), indices, E_temp);
    }
    else
        read_table(Eigenval_name.c_str(), E_temp);

    E_temp /= eigenvalue_scale;

    // Set number of states
    const size_t nst = E_temp.size();
    if(nst==0)
    {
        std::ostringstream oss;
        oss << Eigenval_name << " appears to be empty. Is this the correct eigenvalue input file?.";
        throw std::runtime_error(oss.str());
    }

    // Read first eigenvector into tempory memory to get size of vectors
    arma::vec z_temp;
    arma::vec psi_temp;
    std::string Eigenvect_name = Eigenvect_prefix + "1" + Eigenvect_ext;
    read_table(Eigenvect_name.c_str(), z_temp, psi_temp);

    if(z_temp.size() == 0)
    {
        std::ostringstream oss;
        oss << "No data found in " << Eigenvect_name << ". Is this the correct eigenvector input file?";
        throw std::runtime_error(oss.str());
    }

    // Resize permanent store of eigen-solutions to correct size and copy in first eigen-solution
    const auto psi_size = z_temp.size();
    states.emplace_back(E_temp[0], z_temp, psi_temp);

    // Read in remaining eigenvectors and copy into permanent store
    for(unsigned int ist=1; ist<nst; ist++){
        std::stringstream Eigenvect_name_sstream;
        Eigenvect_name_sstream << Eigenvect_prefix << ist+1 << Eigenvect_ext;
        Eigenvect_name = Eigenvect_name_sstream.str();
        read_table(Eigenvect_name.c_str(), z_temp, psi_temp, psi_size);
        states.emplace_back(E_temp[ist], z_temp, psi_temp);
    }

    return states;
}
        
/** 
 * \brief Write a set of eigenstates to file
 *
 * \param[in]  Eigenval_name     The name of the file to which the eigenvalues will be written
 * \param[in]  Eigenvect_prefix  Prefix of the files holding the eigenvectors
 * \param[in]  Eigenvect_ext     Extension of the files holding the eigenvectors
 * \param[in]  states            Set of eigenstates
 * \param[in]  with_num          First column of eigenvalue file should contain state index
 *
 * \details The eigenvectors are saved in separate files with the form <prefix>i<extension>
 *          where i is the index of the state (starting from 1).
 */
void Eigenstate::write_to_file(const std::string             &Eigenval_name,
                               const std::string             &Eigenvect_prefix,
                               const std::string             &Eigenvect_ext,
                               const std::vector<Eigenstate> &states,
                               const bool                     with_num)
{
    // Output eigenvalues
    arma::vec E_temp(states.size());
    for(unsigned int ist=0; ist < states.size(); ist++)
        E_temp[ist] = states[ist].get_energy();

    write_table(Eigenval_name.c_str(), E_temp, with_num, 17);

    // Output eigenvectors
    for(unsigned int ist=0; ist < states.size(); ist++)
    {
        std::stringstream Eigenvect_name_sstream;
        Eigenvect_name_sstream << Eigenvect_prefix << ist+1 << Eigenvect_ext;
        std::string Eigenvect_name = Eigenvect_name_sstream.str();
        const auto z   = states[ist].get_position_samples();
        const auto psi = states[ist].get_wavefunction_samples();
        write_table(Eigenvect_name.c_str(), z, psi, false, 17);
    }
}

/**
 * \brief Find the expectation position for a given state
 *
 * \param[in] i State
 * \param[in] z Spatial coordinates [m]
 *
 * \return Expectation position [m]
 */
double Eigenstate::get_expectation_position() const
{
    const auto dz = _z[1] - _z[0];
    const decltype(_psi) dz_av = _psi * _psi * _z;

    return integral(dz_av, dz);
}

/** 
 * \brief Find dipole matrix element between a pair of eigenvectors
 *
 * \param[in] i Initial state
 * \param[in] j Final state
 *
 * \return Dipole matrix element [m]
 *
 * \todo FIXME: When |i> == |j>, this should just return the expectation
 *       position.  At the moment, we have no way of figuring this out
 *       and the "pivoting" code below generates the wrong value.  It's
 *       not too important however, because the matrix element for an
 *       intrasubband transition is never really used!
 */
double mij(const Eigenstate &i, const Eigenstate &j)
{
    // FIXME: Currently it is assumed that both states use same spatial grid
    const auto z = i.get_position_samples();
    const double dz = z[1] - z[0];

    /* Because we have a nonparabolic effective mass, the Schroedinger solutions
     * are NOT part of an orthonormal set. As such, we need to do something to
     * ensure spatial invariance when we calculate dipole matrix element.  We
     * therefore define a pivot point halfway between the expectation positions
     * of the electron in each state.
     */
    const auto z_exp_i = i.get_expectation_position();
    const auto z_exp_j = j.get_expectation_position();
    const double z0 = 0.5 * (z_exp_i + z_exp_j);

    const auto psi_i = i.get_wavefunction_samples();
    const auto psi_j = j.get_wavefunction_samples();

    const arma::vec dmij = psi_i * (z - z0) * psi_j;

    return integral(dmij, dz);
}

/**
 * \brief Find the largest probability density at any point in a set of eigenstates
 */
double Eigenstate::psi_squared_max(const std::vector<Eigenstate> &states)
{
    double PDmax = 0.0;

    // Loop through all states
    for(auto st : states)
    {
        // If this state has highest probability so far, store its value
        const auto PD = st.get_PD();
        PDmax = GSL_MAX_DBL(PDmax, PD.max());
    }

    return PDmax;
}

} //namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
