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
#include "qclsim-linalg.h"

namespace Leeds {

class Subband
{
    public:
        Subband(State                 ground_state,
                double                m_d,
                std::valarray<double> z);
        
        Subband(State                 ground_state,
                double                m_d,
                std::valarray<double> z,
                double                alphad,
                double                condband_edge);

        void set_distribution(const double Ef,
                              const double population);

        inline State                       get_ground() const {return ground_state;}
        inline std::valarray<double>       z_array()    const {return _z;}
        inline double                      get_dz()     const {return _z[1]-_z[0];}

        /** Find the total length of the spatial region [m] */
        inline double                      get_length() const {return _z[_z.size()-1]-_z[0];}

        /** Find expectation position for the ground state [m] */
        inline double                      get_z_av_0() const {return z_av(ground_state, _z);}

        inline double                      get_Ef()     const {return Ef;}
        inline double                      get_E()      const {return ground_state.get_E();}
        inline double                      get_pop()    const {return population;}
        inline std::valarray<double>       psi_array()  const {return ground_state.psi_array();}
        inline double                      get_condband_edge() const {return condband_edge;}

        double                             get_k_fermi() const;

        static std::vector<Subband> read_from_file(const std::string &energy_input_path,
                                                   const std::string &wf_input_prefix,
                                                   const std::string &wf_input_ext,
                                                   const std::string &m_d_filename);

        static std::vector<Subband> read_from_file(const std::string &energy_input_path,
                                                   const std::string &wf_input_prefix,
                                                   const std::string &wf_input_ext,
                                                   const double       m_d);

        static std::vector<Subband> read_from_file(const std::string &energy_input_path,
                                                   const std::string &wf_input_prefix,
                                                   const std::string &wf_input_ext,
                                                   const std::string &m_d_filename,
                                                   const std::string &alphad_filename,
                                                   const std::string &potential_filename);

        static std::vector<Subband> read_from_file(const std::string &energy_input_path,
                                                   const std::string &wf_input_prefix,
                                                   const std::string &wf_input_ext,
                                                   const double       m_d,
                                                   const double       alphad,
                                                   const double       V);

        double Ek(double k) const;
        double k(double Ek) const;
        
        /// Return d.o.s mass at bottom of subband 
        inline double get_md_0() const {return md_0;}
       
        /// Return nonparabolicity parameter
        inline double get_alphad() const {return alphad;}
        
        /// Get in plane effective mass at some absolute energy.
        inline double get_m_d(double E) const {return md_0*(1.0 + alphad*(E - condband_edge));}
        
        double rho(const double E) const;
        
        double f_FD(const double E,
                    const double Te) const;
        
        /// Find the population of the subband at wavevector k
        inline double N(double k, double Te) {return rho(get_E()+Ek(k))*f_FD(get_E()+Ek(k),Te);}

    private:
        State  ground_state;        ///< State at bottom of subband
        std::valarray<double> _z;   ///< Spatial profile of subband
        double md_0;                ///< Density of states effective mass
        double alphad;              ///< In-plane nonparabolicity parameter
        double condband_edge;       ///< Energy of the conduction band located below 'most populated' region of subband.

        bool   distribution_known;  ///< True if the carrier distribution is set
        double Ef;                  ///< Quasi-Fermi energy [J]
        double population;          ///< Sheet-density of carriers [m^{-2}]
};

} // namespace Leeds

#endif // QCLSIM_SUBBAND_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
