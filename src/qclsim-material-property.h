#ifndef QCLSIM_MATERIAL_PROPERTY_H
#define QCLSIM_MATERIAL_PROPERTY_H

#include <vector>
#include <glibmm/ustring.h>

namespace xmlpp {
class Element;
}

/**
 * The way in which the material property is specified
 *
 * \todo This should probably be implemented as subclasses
 *       of MaterialProperty
 */
enum MaterialPropertyType {
    MATERIAL_PROPERTY_CONSTANT,
    MATERIAL_PROPERTY_POLY
};

/** Wrapper for XML data describing a physical property of a material */
class MaterialProperty {
public:
    MaterialProperty(){}
    MaterialProperty(xmlpp::Element *elem);

    const Glib::ustring & get_name() const;
    const Glib::ustring & get_unit() const;
    const Glib::ustring & get_text() const;

    virtual double get_val(const double x = 0) const;

    MaterialPropertyType get_type() const;

protected:
    xmlpp::Element       *elem; ///< Underlying XML data
    MaterialPropertyType  type; ///< The way in which the property is specified

    double        _constant; ///< The value of a constant property
    Glib::ustring _text;     ///< Text representation of a constant property

    // Values for polynomial properties
    std::vector<int>    _poly_index; ///< The index of each term
    std::vector<double> _poly_coeff; ///< The coefficient of each term

    Glib::ustring _name; ///< The name of the property
    Glib::ustring _unit; ///< The unit of the property
};
#endif // QCLSIM_MATERIAL_PROPERTY_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
