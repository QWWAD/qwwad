#ifndef QCLSIM_MATERIAL_PROPERTY_H
#define QCLSIM_MATERIAL_PROPERTY_H

namespace xmlpp {
class Element;
}

namespace Glib {
class ustring;
}

class MaterialPropertyImpl;

/**
 * The way in which the material property is specified
 *
 * \todo This should probably be implemented as subclasses
 *       of MaterialProperty
 */
enum MaterialPropertyType {
    MATERIAL_PROPERTY_CONSTANT,
    MATERIAL_PROPERTY_INTERP,
    MATERIAL_PROPERTY_POLY
};

/** Wrapper for XML data describing a physical property of a material */
class MaterialProperty {
public:
    MaterialProperty(xmlpp::Element *elem);
    ~MaterialProperty();

    const Glib::ustring & get_name() const;
    const Glib::ustring & get_unit() const;
    const Glib::ustring & get_text() const;

    double get_val(const double x = 0) const;

    MaterialPropertyType get_type() const;

    double get_interp_y0()   const;
    double get_interp_y1()   const;
    double get_interp_b()    const;
    bool   interp_is_bowed() const;
private:
    MaterialPropertyImpl *priv; ///< Private members of the class
};
#endif // QCLSIM_MATERIAL_PROPERTY_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
