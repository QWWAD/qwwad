#ifndef QCLSIM_MATERIAL_PROPERTY_H
#define QCLSIM_MATERIAL_PROPERTY_H

#include <vector>
#include <glibmm/ustring.h>

namespace xmlpp {
class Element;
}

/** Wrapper for XML data describing a physical property of a material */
class MaterialProperty {
public:
    MaterialProperty(){}
    MaterialProperty(xmlpp::Element *elem);

    const Glib::ustring & get_name() const;
    const Glib::ustring & get_unit() const;
    const Glib::ustring & get_text() const;

    virtual double get_val(const double x = 0) const;

protected:
    xmlpp::Element       *elem; ///< Underlying XML data

    double        _constant; ///< The value of a constant property
    Glib::ustring _text;     ///< Text representation of a constant property

    Glib::ustring _name; ///< The name of the property
    Glib::ustring _unit; ///< The unit of the property
};
#endif // QCLSIM_MATERIAL_PROPERTY_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
