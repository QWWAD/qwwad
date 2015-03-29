#ifndef QCLSIM_MATERIAL_SPECIFICATION_H
#define QCLSIM_MATERIAL_SPECIFICATION_H

namespace Glib {
class ustring;
}

class Material;
class MaterialLibrary;

/**
 * Conduction band valleys
 */
enum Valley
{
    GAMMA,   ///< The \f$\Gamma\f$ valley
    DELTA_2, ///< The 2 \f$\Delta\f$ valleys in the (001) direction
    DELTA_4, ///< The 4 \f$\Delta\f$ valleys in the (001) plane
    L_1,     ///< The \f$L\f$ valley in the (111) direction
    L_3      ///< The other 3 \f$L\f$ valleys
};

/** Crystallographic orientations */
enum Orientation
{
    ORIENT_001, ///< The (001) growth direction
    ORIENT_111  ///< The (111) growth direction
};

/**
 * A complete specification of the material that should be passed to any of the
 * functions here
 */
class MaterialSpecification
{
public:
    MaterialSpecification();

    MaterialSpecification(const char          *mat_name,
                          const Valley        &valley,
                          const Orientation   &orientation,
                          const double         alloy,
                          const double         doping);

    MaterialSpecification(const Glib::ustring &mat_name,
                          const Valley        &valley,
                          const Orientation   &orientation,
                          const double         alloy,
                          const double         doping);

    MaterialSpecification(const MaterialSpecification &spec);

    ~MaterialSpecification();

    MaterialSpecification & operator=(const MaterialSpecification &rhs);

    double get_prop_val_0(const char    *prop_name) const;
    double get_prop_val_0(Glib::ustring &prop_name) const;

    double get_prop_val_x(const char    *prop_name) const;
    double get_prop_val_x(Glib::ustring &prop_name) const;

    const Glib::ustring & get_prop_text(const char    *prop_name) const;
    const Glib::ustring & get_prop_text(Glib::ustring &prop_name) const;

    MaterialLibrary *lib;         ///< The library used to obtain most of the data
    Material        *xml;         ///< Material data read from the external library
    Valley           valley;      ///< The conduction band valley to consider (if applicable)
    Orientation      orientation; ///< The crystal orientation (if applicable)
    double           alloy;       ///< The alloy fraction (if applicable)
    double           n3d;         ///< The doping density (if applicable)
};
#endif // QCLSIM_MATERIAL_SPECIFICATION_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
