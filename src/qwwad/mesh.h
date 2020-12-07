/**
 * \file   mesh.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for heterostructure array generation
 */

#ifndef MESH_H
#define MESH_H

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <armadillo>
#include <string>
#include <valarray>
#include <vector>

namespace QWWAD
{
/// Convenience wrapper for a list of vector components in each layer
using alloy_vector = std::vector<std::valarray<double>>;

/**
 * \brief A stack of layers making up a quantum heterostructure
 */
class Mesh
{
private:
    // Parameters for each individual layer of the structure
    size_t       _n_alloy;    ///< Number of alloy components
    alloy_vector _x_layer;    ///< Alloy fractions in each layer
    arma::vec    _W_layer;    ///< Width of each layer [m]
    arma::vec    _n3D_layer;  ///< Donor density in each layer [m^{-3}]

    static void read_layers_from_file(const std::string    &filename,
                                      decltype(_x_layer)   &x_layer,
                                      decltype(_W_layer)   &W_layer,
                                      decltype(_n3D_layer) &n3D_layer);


    size_t                _n_periods;  ///< Number of periods in the structure
    size_t                _ncell_1per; ///< Number of cells in each period of the mesh

    std::valarray<unsigned int> _layer_top_index; ///< Index of the last cell in each layer

    // Parameters for each point in the entire, expanded structure
    std::valarray<double> _z;   ///< Spatial position at the middle of each cell [m]
    alloy_vector          _x;   ///< Alloy fractions at the middle of each cell
    std::valarray<double> _n3D; ///< Volume doping at the middle of each cell [m^{-3}]
    double                _Lp;  ///< Length of one period [m]
    double                _dz;  ///< Width of each cell [m]

public:
    Mesh(const decltype(_x_layer)   &x_layer,
         decltype(_W_layer)          W_layer,
         decltype(_n3D_layer)        n3D_layer,
         const decltype(_ncell_1per) ncell_1per,
         const decltype(_n_periods)  n_periods = 1);

    static auto create_from_file_auto_nz(const std::string &layer_filename,
                                         const size_t       n_periods,
                                         const double       dz_max = 1e-10) -> Mesh*;

    static auto create_from_file(const std::string &layer_filename,
                                 const size_t       nz_1per,
                                 const size_t       n_periods) -> Mesh*;

    /** Return the number of cells in one period of the mesh */
    [[nodiscard]] inline auto get_ncell_1per() const {return _ncell_1per;}

    /** Return the total number of sampling points in the entire structure */
    [[nodiscard]] inline auto get_ncell() const {return _z.size();}

    [[nodiscard]] inline auto get_z()                const {return _z;}
    [[nodiscard]] inline auto get_z(unsigned int iz) const {return _z[iz];}
    [[nodiscard]] inline auto get_dz()               const {return _dz;}

    /** Return the number of alloy components in the structure */
    [[nodiscard]] inline auto get_n_alloy()      const {return _n_alloy;}
    [[nodiscard]] inline auto get_layer_widths() const {return _W_layer;}
    [[nodiscard]] inline auto get_x_array()      const {return _x;}

    [[nodiscard]] auto get_n3D_in_layer(const unsigned int iL) const -> double;
    [[nodiscard]] auto get_n3D_at_point(const unsigned int iz) const -> double;

    /**
     * Return the entire array of doping at each point
     */
    [[nodiscard]] inline auto get_n3D_array() const {return _n3D;}

    /**
     * \brief Return the number of layers in a single period
     *
     * \returns Number of layers in the structure
     */
    [[nodiscard]] inline auto get_n_layers_per_period() const {return _W_layer.size();}

    /// Return the number of layers in the entire structure
    [[nodiscard]] inline auto get_n_layers_total() const {return _W_layer.size()*_n_periods;}

    [[nodiscard]] auto get_layer_from_height(const double z) const -> unsigned int;

    [[nodiscard]] auto point_is_in_layer(const double       z,
                                         const unsigned int iL) const -> bool;

    [[nodiscard]] auto get_height_at_top_of_layer(const unsigned int iL) const -> double;

    [[nodiscard]] auto get_layer_top_index(const unsigned int iL) const -> unsigned int;
    [[nodiscard]] inline auto get_layer_top_indices() const {return _layer_top_index;}

    /// Return the length of a single period of the structure
    [[nodiscard]] inline auto get_period_length() const {return sum(_W_layer);}

    /// Return the entire length of the structure
    [[nodiscard]] inline auto get_total_length() const {return sum(_W_layer)*_n_periods;}
};
} // namespace
#endif // MESH_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
