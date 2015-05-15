#ifndef QCLSIM_MATERIAL_PROPERTY_H
#define QCLSIM_MATERIAL_PROPERTY_H

#include <glibmm/ustring.h>

namespace xmlpp {
class Element;
}

/** Wrapper for XML data describing a physical property of a material */
class MaterialProperty {
public:
    MaterialProperty(){}
    MaterialProperty(xmlpp::Element *elem);

    virtual const Glib::ustring & get_name() const;

protected:
    xmlpp::Element       *elem; ///< Underlying XML data

    Glib::ustring _name; ///< The name of the property
};
#endif // QCLSIM_MATERIAL_PROPERTY_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
