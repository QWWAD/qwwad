/**
 * \file    heterostructure.cpp
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date    2012-08-03
 *
 * \brief   A description of layers in a heterostructure
 */

#include "heterostructure.h"

#include <gsl/gsl_sf_erf.h>
#include <error.h>
#include <fstream>

#ifdef DEBUG
# include <iostream>
#endif

#include "qclsim-fileio.h"

/**
 * Gauss error function
 */
#define erf gsl_sf_erf

namespace Leeds {
/**
 * \brief Create a heterostructure using specifications of each layer parameter
 *
 * \param[in] x_layer   Alloy fraction of each layer
 * \param[in] W_layer   Thickness of each layer [m]
 * \param[in] n3D_layer Volume doping in each layer [m^{-3}]
 * \param[in] dz_max    The maximum permissible spatial separation between points in the
 *                      structure
 * \param[in] n_periods Number of periods of the structure to generate
 * \param[in] L_diff    Alloy diffusion length [m]
 *
 * \details The Heterostructure is created using an automatically-determined spatial
 *          resolution.  You specify an upper limit on the separation between points
 *          by using the dz_max and the class works out the a suitable value.  We do
 *          things this way so that we can guarantee that each period in the structure
 *          has the same, fixed number of points.
 */
Heterostructure::Heterostructure(const std::valarray<double>& x_layer,
                                 const std::valarray<double>& W_layer,
                                 const std::valarray<double>& n3D_layer,
                                 const double                 dz_max,
                                 const size_t                 n_periods,
                                 const double                 L_diff) :
    _x_layer(x_layer),
    _W_layer(W_layer),
    _n3D_layer(n3D_layer),
    _n_periods(n_periods),
    _L_diff(L_diff),
    _nz_1per(ceil( ceil(_W_layer.sum()*n_periods/dz_max) / n_periods )),
    _layer_top_index(_x_layer.size() * n_periods),
    _z(_nz_1per*n_periods),
    _x_nominal(_nz_1per*n_periods),
    _x_diffuse(_nz_1per*n_periods),
    _n3D(_nz_1per*n_periods)
{
    const double dz = _W_layer.sum()/_nz_1per;

    // TODO: Check that no layer is thinner than dz and throw an error if it is

    for (unsigned int iL = 0; iL < _x_layer.size() * _n_periods; ++iL)
    {
        const double height = get_height_at_top_of_layer(iL);

        // Divide by the spatial separation to find the total height of the
        // stack up to the top of the desired layer
        _layer_top_index[iL] = (unsigned int)round(height / dz);
    }

    // Calculate spatial points
    for (unsigned int iz = 0; iz < _nz_1per * n_periods; ++iz)
    {
        _z[iz] = iz*dz; // Calculate the spatial location [m]

        const unsigned int iL = get_layer_from_height(_z[iz]);
        _n3D[iz]       = get_n3D_in_layer(iL);
        _x_nominal[iz] = get_x_in_layer_nominal(iL);
        _x_diffuse[iz] = calculate_x_annealed_at_point(iz);
    }
}

/**
 * Create a heterostructure using data from an input file containing data for each layer
 *
 * \param[in] layer_filename Name of input file
 * \param[in] unit           The unit of measurement for the layer thicknesses
 * \param[in] dz_max         The maximum permissible spatial separation between points in the
 *                           structure
 * \param[in] n_periods      Number of periods to generate
 *
 * \return A new heterostructure object for the system.  Remember to delete it after use!
 */
Heterostructure* Heterostructure::read_from_file(const std::string &layer_filename,
                                                 const Unit         thickness_unit,
                                                 const double       dz_max,
                                                 const size_t       n_periods,
                                                 const double       L_diff)
{
    std::valarray<double> x_layer;   // Alloy fraction for each layer
    std::valarray<double> W_layer;   // Thickness of each layer
    std::valarray<double> n3D_layer; // Doping density of each layer

    read_table_xyz(layer_filename.c_str(), W_layer, x_layer, n3D_layer);

    const size_t nL = W_layer.size();

    // Check data integrity
    for(unsigned int iL = 0; iL < nL; iL++)
    {
        check_positive(&W_layer[iL]);
        check_c_interval_0_1(&x_layer[iL]);
        check_not_negative(&n3D_layer[iL]);
    }

    // Scale layer widths according to unit
    switch(thickness_unit)
    {
        case UNIT_NM:
            W_layer *= 1.0e-9;
            break;
        case UNIT_ANGSTROM:
            W_layer *= 1.0e-10;
    }

    n3D_layer *= 1000000.0; // Scale doping to m^{-3}

    // Pack input data into a heterostructure object
    return new Heterostructure(x_layer, W_layer, n3D_layer, dz_max, n_periods, L_diff);
}

/**
 * \brief Return the doping concentration in a given layer
 *
 * \param[in] iL The index of the layer
 * 
 * \return The nominal doping density [m\f$^{-3}\f$]
 *
 * \todo Implement dopant diffusion
 */
double Heterostructure::get_n3D_in_layer(const unsigned int iL) const
{
    if(iL > _n3D_layer.size() * _n_periods)
        throw std::domain_error("Tried to access the doping concentration in a layer outside the heterostructure.");

    return _n3D_layer[iL%_n3D_layer.size()];
}

/** Get the doping concentration at a given point in the structure */
double Heterostructure::get_n3D_at_point(const unsigned int iz) const
{
    if(iz > _n3D.size())
        throw std::domain_error("Tried to access the doping concentration at a point outside the heterostructure.");

    return _n3D[iz];
}

/**
 * \brief Return the nominal alloy fraction in a given layer
 *
 * \param[in] iL The index of the layer
 *
 * \return The nominal alloy fraction
 */
double Heterostructure::get_x_in_layer_nominal(const unsigned int iL) const
{
    if(iL > _x_layer.size() * _n_periods)
        throw std::domain_error("Tried to access the alloy fraction in a layer outside the heterostructure.");

    return _x_layer[iL%_x_layer.size()];
}

/**
 * \brief     Finds the index of the point at the top of a layer
 *
 * \param[in] iL        The layer number
 * 
 * \returns The index of the last point in the specified layer.
 *
 * \details Normally, if the layer number exceeds the number of layers in the
 *          structure, the function will throw an error.  However, if the
 *          "periodic" flag is set as true, the function will return the index
 *          of the point, assuming that the structure is infinitely long and
 *          periodic.
 */
unsigned int Heterostructure::get_layer_top_index(const unsigned int iL) const
{
    if(iL >= _W_layer.size()*_n_periods)
        throw std::domain_error("Cannot find layer top: layer index is outside the heterostructure");

    return _layer_top_index[iL];
}

/**
 * Get the height of the structure up to the top of a given layer
 *
 * \param[in] iL The index of the layer
 *
 * \return The height of the structure at the top of that layer
 *
 * \details If the heterostructure is periodic, then the specified layer index
 *          can be greater than the number in the heterostructure.
 */
double Heterostructure::get_height_at_top_of_layer(const unsigned int iL) const
{
    if(iL >= _W_layer.size() * _n_periods)
        throw std::domain_error("Tried to find top of a layer that is outside the heterostructure.");

    // Find the height of the highest **incomplete** period
    double height = _W_layer[std::slice(0, iL%_W_layer.size() + 1, 1)].sum();

    // Add on the height of all the **complete** periods below this layer
    if(_n_periods > 1)
        height += _W_layer.sum() * floor(static_cast<double>(iL)/_W_layer.size());

    return height;
}

/**
 * \brief Find the alloy fraction at a given position, when annealing is included
 *
 * \param[in] z     The spatial position at which to find alloy fraction [m]
 *
 * \return The alloy fraction in the annealed structure at the given point
 *
 * \details See [E. H. Li et. al., IEEE J. Quantum Electron. 32, 1399 (1996)]
 *          for details.  I have replaced the infinitely thick barrier
 *          criterion with periodic boundaries.  This might mess up if very
 *          short structures are used (thinner than the diffusion length) as I
 *          only included adjacent periods in the calculation.
 */
double Heterostructure::calculate_x_annealed_at_point(const unsigned int iz) const
{
    const double Lp  = _W_layer.sum();  // Length of period [m]
    const double z   = fmod(_z[iz], Lp); // Position within period [m]
    const double dz  = _z[1] - _z[0];   // Spatial step [m]
    double _x = 0.0; // Ge fraction

    if (_L_diff > 0)
    {
        // Find contribution from each layer
        for(unsigned int iL = 0; iL < _W_layer.size(); iL++){
            // Find top of layer
            const double zi = dz*get_layer_top_index(iL);

            // Find bottom of layer
            double zi_1 = 0;
            if(iL != 0)
                zi_1 = dz*get_layer_top_index(iL-1);

            _x += 0.5 * _x_layer[iL] * (erf((z-zi_1)/_L_diff) - erf((z-zi)/_L_diff));

            // Add contribution from layers in previous period
            _x += 0.5 * _x_layer[iL] * (erf((z-(zi_1-Lp))/_L_diff)
                    - erf((z-(zi-Lp))/_L_diff));

            // Add contributions from layers in next period
            _x += 0.5 * _x_layer[iL] * (erf((z-(zi_1+Lp))/_L_diff)
                    - erf((z-(zi+Lp))/_L_diff));
        }
    }
    else
        _x = _x_nominal[iz];

    return _x;
}

/**
 * Determine whether a point is in a layer
 *
 * \param[in] z  A spatial point [m]
 * \param[in] iL Index of the layer
 *
 * \return True if the point is in the layer
 *
 * \details The point is considered to be within the layer
 *          if z_top[iL-1] <= z < z_top[iL]
 */
bool Heterostructure::point_is_in_layer(const double       z,
                                        const unsigned int iL) const
{
    if(iL >= get_n_layers_total() or z > get_total_length())
    {
        std::ostringstream oss;
        oss << "Tried to check whether the layer with index " << iL << " contains the point at " << z*1e9 << " nm."
            " However, the structure only contains " << get_n_layers_total() << " layers and has a total height of "
            << get_total_length()*1e9 << " nm.";

        throw std::domain_error(oss.str());
    }

    double       top_of_previous_layer = 0;
    const double top_of_this_layer     = get_height_at_top_of_layer(iL);

    if(iL > 0)
        top_of_previous_layer = get_height_at_top_of_layer(iL-1);

#ifdef DEBUG
    std::cout << z*1e9 << " nm.  Layer = " << iL << "." << std::endl
              << " * Top of previous layer = " << top_of_previous_layer*1e9 << " nm." << std::endl
              << " * Top of this layer     = " << top_of_this_layer*1e9     << " nm." << std::endl;
#endif

    const bool result = (z >= top_of_previous_layer and z < top_of_this_layer);

#ifdef DEBUG
    std::cout << "Point " << (result ? "is" : "is not") << " in this layer" << std::endl;
#endif

    return result;
}

/**
 * Find the index of the layer that contains a given point
 *
 * \param[in] z height [m]
 *
 * \returns The index of the layer containing point z
 *
 * \details A point lies within layer i if z_(i-1) <= z < z_i. 
 */
unsigned int Heterostructure::get_layer_from_height(const double z) const
{
    if (z > get_period_length()*_n_periods)
    {
        std::ostringstream oss;
        oss << "Tried to find layer index at a height of " << z*1e9 << " nm, but total structure length is "
            << get_period_length()*_n_periods << " nm.";
        throw std::domain_error(oss.str());
    }

    unsigned int iL = 0; // index of layer

    while(!point_is_in_layer(z,iL))
        iL++;

    return iL;
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
