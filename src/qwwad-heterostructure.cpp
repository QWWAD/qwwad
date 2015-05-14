/**
 * \file    qwwad-heterostructure.cpp
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief   A description of layers in a heterostructure
 */

#include "qwwad-heterostructure.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <fstream>

#ifdef DEBUG
# include <iostream>
#endif

#include "qwwad/file-io.h"

/**
 * Gauss error function
 */
#define erf gsl_sf_erf

namespace QWWAD
{
/**
 * \brief Create a heterostructure using specifications of each layer parameter
 *
 * \param[in] x_layer   Alloy fraction of each layer
 * \param[in] W_layer   Thickness of each layer [m]
 * \param[in] n3D_layer Volume doping in each layer [m^{-3}]
 * \param[in] nz_1per   The number of points to be used in mapping out a single
 *                      period of the structure
 * \param[in] n_periods Number of periods of the structure to generate
 *
 * \details The number of points in the entire structure is given by
 *          N = (nz_{\text{1per}} - 1) \times n_{\text{period}} + 1
 */
Heterostructure::Heterostructure(const alloy_vector          &x_layer,
                                 const std::valarray<double> &W_layer,
                                 const std::valarray<double> &n3D_layer,
                                 const size_t                 nz_1per,
                                 const size_t                 n_periods) :
    _n_alloy(x_layer.at(0).size()),
    _x_layer(x_layer),
    _W_layer(W_layer),
    _n3D_layer(n3D_layer),
    _n_periods(n_periods),
    _nz_1per(nz_1per),
    _layer_top_index(_x_layer.size() * n_periods),
    _z((_nz_1per - 1)*n_periods + 1),
    _x(_z.size(), std::valarray<double>(_n_alloy)),
    _n3D(_z.size())
{
    const double Lp = _W_layer.sum();    // Length of one period [m]
    const size_t n_cell_1per = _nz_1per - 1; // Number of spatial intervals in a period
    _dz = Lp/n_cell_1per; // Separation between points [m]

    // TODO: Check that no layer is thinner than dz and throw an error if it is

    // Calculate spatial points
    unsigned int iL_cache = 0;
    unsigned int idx = 0;
    const size_t nz = _z.size(); // Total points in structure

    for (unsigned int iz = 0; iz < nz; ++iz)
    {
        _z[iz] = iz*_dz; // Calculate the spatial location [m]

        const double zp = fmod(_z[iz], Lp); // Position within period [m]
        unsigned int iL = 0; // Layer index

        // If this is the last point in the structure, insist that we stay
        // inside the last layer
        if (iz < nz-1)
            iL = get_layer_from_height(zp); // Layer index
        else
        {
            iL = iL_cache;
            _layer_top_index[_x_layer.size()-1] = _z.size();
        }

        // Only recalculate if we're in a new period
        if(iL != iL_cache)
        {
            // TODO: Add some kind of bounds checking
            _layer_top_index[idx++] = iz;
            iL_cache = iL;
        }

        _n3D[iz]       = get_n3D_in_layer(iL);

        for(unsigned int ialloy = 0; ialloy < _x_layer[0].size(); ++ialloy)
            _x.at(iz)[ialloy] = get_x_in_layer(iL, ialloy);
    }
}

void Heterostructure::read_layers_from_file(const std::string     &filename,
                                            alloy_vector          &x_layer,
                                            std::valarray<double> &W_layer,
                                            std::valarray<double> &n3D_layer,
                                            const Unit             thickness_unit)
{
    std::valarray<double> coltemp; // A temporary vector containing the first column of the input file
    read_table(filename.c_str(), coltemp);
    const size_t nL = coltemp.size();

    std::valarray<double> rowtemp; // A temporary vector containing the first row of the input file

    std::ifstream stream(filename.c_str());

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << filename;
        throw std::runtime_error(oss.str());
    }

    size_t n_alloy = 0; // Number of alloy components in material system

    for(unsigned int i = 0; i < nL; ++i)
    {
        int result = 0;

        if(i == 0)
        {
            read_line_array_u(rowtemp, stream);
            // Find number of alloy components (noting that first element on line is the
            // thickness and the last is the doping density)
            n_alloy = rowtemp.size() - 2;

            x_layer.resize(nL, std::valarray<double>(n_alloy));
            W_layer.resize(nL);
            n3D_layer.resize(nL);
        }
        else
        {
            if(!stream)
                throw std::runtime_error("Could not read stream");

            result = read_line_array(rowtemp, n_alloy + 2, stream); 
        }

        if(result == 0)
        {
            W_layer[i] = rowtemp[0];

            for(unsigned int j = 0; j < n_alloy; ++j)
                x_layer.at(i)[j] = rowtemp[1+j];

            n3D_layer[i] = rowtemp[n_alloy+1];
        }
        else
            throw std::runtime_error("Could not read stream");
    }

    stream.close();

    // Check data integrity
    for(unsigned int iL = 0; iL < nL; iL++)
    {
        check_positive(&W_layer[iL]);

        for(unsigned int ialloy = 0; ialloy < n_alloy; ++ialloy)
            check_c_interval_0_1(&x_layer.at(iL)[ialloy]);

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
}

/**
 * Create a heterostructure using data from an input file, with the number of spatial points
 * per period calculated automatically
 *
 * \param[in] layer_filename Name of input file
 * \param[in] thickness_unit The unit of measurement for the layer thicknesses
 * \param[in] n_periods      Number of periods to generate
 * \param[in] dz_max         The maximum allowable separation between spatial points [m]
 *
 * \return A new heterostructure object for the system.  Remember to delete it after use!
 */
Heterostructure* Heterostructure::create_from_file_auto_nz(const std::string &layer_filename,
                                                           const Unit         thickness_unit,
                                                           const size_t       n_periods,
                                                           const double       dz_max)
{
    alloy_vector          x_layer;   // Alloy fraction for each layer
    std::valarray<double> W_layer;   // Thickness of each layer
    std::valarray<double> n3D_layer; // Doping density of each layer

    read_layers_from_file(layer_filename, x_layer, W_layer, n3D_layer, thickness_unit);

    const double period_length = W_layer.sum();
    const size_t nz_1per = ceil(period_length/dz_max) + 1;

    // Pack input data into a heterostructure object
    return new Heterostructure(x_layer, W_layer, n3D_layer, nz_1per, n_periods);
}

/**
 * Create a heterostructure using data from an input file containing data for each layer
 *
 * \param[in] layer_filename Name of input file
 * \param[in] thickness_unit The unit of measurement for the layer thicknesses
 * \param[in] nz_1per        The number of points to be used in mapping out a single
 *                           period of the structure
 * \param[in] n_periods      Number of periods to generate
 *
 * \return A new heterostructure object for the system.  Remember to delete it after use!
 */
Heterostructure* Heterostructure::create_from_file(const std::string &layer_filename,
                                                   const Unit         thickness_unit,
                                                   const size_t       nz_1per,
                                                   const size_t       n_periods)
{
    alloy_vector          x_layer;   // Alloy fraction for each layer
    std::valarray<double> W_layer;   // Thickness of each layer
    std::valarray<double> n3D_layer; // Doping density of each layer

    read_layers_from_file(layer_filename, x_layer, W_layer, n3D_layer, thickness_unit);

    // Pack input data into a heterostructure object
    return new Heterostructure(x_layer, W_layer, n3D_layer, nz_1per, n_periods);
}

/**
 * \brief Return the doping concentration in a given layer
 *
 * \param[in] iL The index of the layer
 * 
 * \return The doping density [m\f$^{-3}\f$]
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
 * \brief Return the alloy fraction in a given layer
 *
 * \param[in] iL The index of the layer
 *
 * \return The alloy fraction
 */
double Heterostructure::get_x_in_layer(const unsigned int iL,
                                               const unsigned int ialloy) const
{
    if(iL >= _x_layer.size() * _n_periods)
        throw std::domain_error("Tried to access the alloy fraction in a layer outside the heterostructure.");

    return _x_layer.at(iL%_x_layer.size())[ialloy];
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
    if(iL >= _W_layer.size()*_n_periods + 1)
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
//    if(iL >= _W_layer.size() * _n_periods + 10)
//        throw std::domain_error("Tried to find top of a layer that is outside the heterostructure.");

    // Find the height of the highest **incomplete** period
    double height = _W_layer[std::slice(0, iL%_W_layer.size() + 1, 1)].sum();

    // Add on the height of all the **complete** periods below this layer
    if(_n_periods > 1)
        height += _W_layer.sum() * floor(static_cast<double>(iL)/_W_layer.size());

    return height;
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
    /*
    if(iL >= get_n_layers_total() + 1 or gsl_fcmp(z, get_total_length()*2, get_dz()/10) == 1)
    {
        std::ostringstream oss;
        oss << "Tried to check whether the layer with index " << iL << " contains the point at " << z*1e9 << " nm."
            " However, the structure only contains " << get_n_layers_total() << " layers and has a total height of "
            << get_total_length()*1e9 << " nm.";

        throw std::domain_error(oss.str());
    }*/

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
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
