#ifndef QCLSIM_MATERIAL_SPECIFICATION_H
#define QCLSIM_MATERIAL_SPECIFICATION_H

namespace Glib {
class ustring;
}

namespace QWWAD {
class Material;
class MaterialLibrary;

/**
 * A complete specification of the material that should be passed to any of the
 * functions here
 */
class MaterialSpecification
{
public:
    MaterialSpecification();

    MaterialSpecification(const char          *mat_name,
                          const double         alloy,
                          const double         doping);

    MaterialSpecification(const Glib::ustring &mat_name,
                          const double         alloy,
                          const double         doping);

    MaterialSpecification(const MaterialSpecification &spec);

    ~MaterialSpecification();

    MaterialSpecification & operator=(const MaterialSpecification &rhs);

    double get_prop_val(const char    *prop_name,
                        const double   x) const;
    double get_prop_val(Glib::ustring &prop_name,
                        const double   x) const;

    double get_prop_val_0(const char    *prop_name) const;
    double get_prop_val_0(Glib::ustring &prop_name) const;

    double get_prop_val_x(const char    *prop_name) const;
    double get_prop_val_x(Glib::ustring &prop_name) const;

    const Glib::ustring & get_prop_text(const char    *prop_name) const;
    const Glib::ustring & get_prop_text(Glib::ustring &prop_name) const;

    MaterialLibrary *lib;         ///< The library used to obtain most of the data
    Material        *xml;         ///< Material data read from the external library
    double           alloy;       ///< The alloy fraction (if applicable)
    double           n3d;         ///< The doping density (if applicable)
};
} // end namespace
#endif // QCLSIM_MATERIAL_SPECIFICATION_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
