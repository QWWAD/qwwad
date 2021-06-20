/**
 * \file    mesh.cpp
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief   A description of layers in a heterostructure
 */

#include "mesh.h"

#include <fstream>
#include <utility>


#ifdef DEBUG
# include <iostream>
#endif

#include "data-checker.h"
#include "file-io.h"

namespace QWWAD
{
/**
 * \brief Create a Mesh using specifications of each layer parameter
 *
 * \param[in] x_layer    Alloy fractions in each layer
 * \param[in] W_layer    Thickness of each layer [m]
 * \param[in] n3D_layer  Volume doping in each layer [m^{-3}]
 * \param[in] ncell_1per The number of cells to be used in the mesh for each period
 * \param[in] n_periods  Number of periods of the structure to generate
 */
Mesh::Mesh(const decltype(_x_layer)    &x_layer,
           decltype(_W_layer)           W_layer,
           decltype(_n3D_layer)         n3D_layer,
           const decltype(_ncell_1per)  ncell_1per,
           const decltype(_n_periods)   n_periods) :
    _n_alloy(x_layer.at(0).size()),
    _x_layer(x_layer),
    _W_layer(std::move(W_layer)),
    _n3D_layer(std::move(n3D_layer)),
    _n_periods(n_periods),
    _ncell_1per(ncell_1per),
    _layer_top_index(_x_layer.size() * n_periods),
    _z(_ncell_1per*_n_periods),
    _x(_z.size(), std::valarray<double>(_n_alloy)),
    _n3D(_z.size()),
    _Lp(sum(_W_layer)),
    _dz(_Lp/_ncell_1per)
{
    const auto n_layer_1per = _W_layer.size(); // Number of layers in one period
    const auto n_layer      = n_layer_1per * _n_periods; // Total number of layers in system

    // Check that no layer is thinner than dz and throw an error if it is
    if (_W_layer.min() < _dz)
    {
        std::ostringstream oss;
        oss << "Layer with " << _W_layer.min()*1e10 << " angstrom width is thinner than spatial separation " << _dz << " angstrom." << std::endl;
        throw std::runtime_error(oss.str());
    }

    // Find the index at the top of each layer
    for(unsigned int iL = 0; iL < n_layer; ++iL)
    {
        _layer_top_index[iL] = get_layer_top_index(iL);

        unsigned int _previous_layer_top_index = 0;

        if (iL > 0) {
            _previous_layer_top_index = _layer_top_index[iL-1];
        }

        // Now fill in the properties for each cell within the layer
        for (unsigned int icell = _previous_layer_top_index;
                          icell < _layer_top_index[iL];
                          ++icell)
        {
            _z[icell]   = (icell+0.5) * _dz; // Set location to middle of cell
            _n3D[icell] = get_n3D_in_layer(iL);

            // Copy all the alloy fractions for this layer
            for(unsigned int ialloy = 0; ialloy < _n_alloy; ++ialloy) {
                _x.at(icell)[ialloy] = _x_layer.at(iL%n_layer_1per)[ialloy];
            }
        }
    }
}


void Mesh::read_layers_from_file(const std::string &filename,
                                 alloy_vector      &x_layer,
                                 arma::vec         &W_layer,
                                 arma::vec         &n3D_layer)
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
            if(!stream) {
                throw std::runtime_error("Could not read stream");
            }

            result = read_line_array(rowtemp, n_alloy + 2, stream); 
        }

        if(result == 0)
        {
            W_layer[i] = rowtemp[0];

            for(unsigned int j = 0; j < n_alloy; ++j) {
                x_layer.at(i)[j] = rowtemp[1+j];
            }

            n3D_layer(i) = rowtemp[n_alloy+1];
        } else {
            throw std::runtime_error("Could not read stream");
        }
    }

    stream.close();

    // Check data integrity
    DataChecker::check_positive(W_layer);
    DataChecker::check_not_negative(n3D_layer);

    // Check data integrity
    for(unsigned int iL = 0; iL < nL; iL++)
    {
        for(unsigned int ialloy = 0; ialloy < n_alloy; ++ialloy) {
            check_c_interval_0_1(&x_layer.at(iL)[ialloy]);
        }
    }

    // Scale layer widths angstrom -> metre
    W_layer *= 1.0e-10;

    n3D_layer *= 1000000.0; // Scale doping to m^{-3}
}

/**
 * Create a Mesh using data from an input file, with the number of spatial points
 * per period calculated automatically
 *
 * \param[in] layer_filename Name of input file
 * \param[in] n_periods      Number of periods to generate
 * \param[in] dz_max         The maximum allowable width of each cell [m]
 *
 * \return A new Mesh object for the system.  Remember to delete it after use!
 */
auto Mesh::create_from_file_auto_nz(const std::string &layer_filename,
                                     const size_t       n_periods,
                                     const double       dz_max) -> Mesh*
{
    alloy_vector x_layer;   // Alloy fraction for each layer
    arma::vec    W_layer;   // Thickness of each layer
    arma::vec    n3D_layer; // Doping density of each layer

    read_layers_from_file(layer_filename, x_layer, W_layer, n3D_layer);

    const double period_length = sum(W_layer);

    // Round up to get the required number of cells per period
    const size_t nz_1per = ceil(period_length/dz_max);

    // Pack input data into a Mesh object
    return new Mesh(x_layer, W_layer, n3D_layer, nz_1per, n_periods);
}

/**
 * Create a Mesh using data from an input file containing data for each layer
 *
 * \param[in] layer_filename Name of input file
 * \param[in] ncell_1per     The number of cells to be used in each period
 * \param[in] n_periods      Number of periods to generate
 *
 * \return A new Mesh object for the system.  Remember to delete it after use!
 */
auto Mesh::create_from_file(const std::string &layer_filename,
                             const size_t       ncell_1per,
                             const size_t       n_periods) -> Mesh*
{
    alloy_vector x_layer;   // Alloy fraction for each layer
    arma::vec    W_layer;   // Thickness of each layer
    arma::vec    n3D_layer; // Doping density of each layer

    read_layers_from_file(layer_filename, x_layer, W_layer, n3D_layer);

    // Pack input data into a Mesh object
    return new Mesh(x_layer, W_layer, n3D_layer, ncell_1per, n_periods);
}

/**
 * \brief Return the doping concentration in a given layer
 *
 * \param[in] iL The index of the layer
 * 
 * \return The doping density [m\f$^{-3}\f$]
 */
auto Mesh::get_n3D_in_layer(const unsigned int iL) const -> double
{
    if(iL > _n3D_layer.size() * _n_periods) {
        throw std::domain_error("Tried to access the doping concentration in a layer outside the heterostructure.");
    }

    return _n3D_layer[iL%_n3D_layer.size()];
}

/** Get the doping concentration at a given point in the structure */
auto Mesh::get_n3D_at_point(const unsigned int iz) const -> double
{
    if(iz > _n3D.size()) {
        throw std::domain_error("Tried to access the doping concentration at a point outside the heterostructure.");
    }

    return _n3D[iz];
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
auto Mesh::get_layer_top_index(const unsigned int iL) const -> unsigned int
{
    // First figure out how many complete periods precede this layer
    const auto previous_periods = iL / _W_layer.size(); // Integer division

    // ...and hence how many cells in those underlying periods
    const auto previous_period_cells = _ncell_1per * previous_periods;

    // Now work within this (incomplete) period
    const auto iL_per = iL % _W_layer.size(); // Index of layer WITHIN period
    const auto z_at_top  = get_height_at_top_of_layer(iL_per);
    unsigned int iz_at_top = round(z_at_top / _dz); // Round to nearest layer

    // Fix possible rounding error at top of period
    if(iz_at_top > _ncell_1per) {
        iz_at_top = _ncell_1per;
    }

    return iz_at_top + previous_period_cells;
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
auto Mesh::get_height_at_top_of_layer(const unsigned int iL) const -> double
{
    const arma::vec layer_tops = cumsum(_W_layer);

    // Find the height of the highest **incomplete** period
    double height = layer_tops(iL%_W_layer.size());

    // Add on the height of all the **complete** periods below this layer
    if(_n_periods > 1)
    {
        const auto n_previous_periods = iL / _W_layer.size(); // Integer division
        height += sum(_W_layer) * n_previous_periods;
    }

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
auto Mesh::point_is_in_layer(const double       z,
                             const unsigned int iL) const -> bool
{
    double       top_of_previous_layer = 0;
    const double top_of_this_layer     = get_height_at_top_of_layer(iL);

    if(iL > 0) {
        top_of_previous_layer = get_height_at_top_of_layer(iL-1);
    }

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
auto Mesh::get_layer_from_height(const double z) const -> unsigned int
{
    if (z > get_period_length()*_n_periods)
    {
        std::ostringstream oss;
        oss << "Tried to find layer index at a height of " << z*1e9 << " nm, but total structure length is "
            << get_period_length()*_n_periods << " nm.";
        throw std::domain_error(oss.str());
    }

    unsigned int iL = 0; // index of layer

    while(!point_is_in_layer(z,iL)) {
        iL++;
    }

    return iL;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
