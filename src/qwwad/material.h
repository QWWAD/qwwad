#ifndef QWWAD_MATERIAL
#define QWWAD_MATERIAL

#include <boost/ptr_container/ptr_map.hpp>
#include <libxml++/libxml++.h>

namespace xmlpp {
class Element;
}

namespace Glib {
class ustring;
}

namespace QWWAD {
class MaterialProperty;
class MaterialPropertyNumeric;

/** Wrapper for XML data for a material */
class Material {
public:
    Material(const Material *mat);
    Material(xmlpp::Element *elem);

    [[nodiscard]] const Glib::ustring & get_name() const;
    [[nodiscard]] const Glib::ustring & get_description() const;

    MaterialProperty const * get_property(const char          *property_name) const;
    [[nodiscard]] MaterialProperty const * get_property(const Glib::ustring &property_name) const;

    MaterialPropertyNumeric const * get_numeric_property(const char    *property_name) const;
    MaterialPropertyNumeric const * get_numeric_property(Glib::ustring &property_name) const;

    double get_property_value(const char   *property_name,
                              const double  x = 0) const;

    double get_property_value(Glib::ustring &property_name,
                              const double   x = 0) const;

private:
    /// Cached set of material properties
    boost::ptr_map<Glib::ustring, MaterialProperty> properties;

    Glib::ustring          name;           ///< The name of the material
    Glib::ustring          description;    ///< The description of the material
};
} // end namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
