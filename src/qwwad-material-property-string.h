#ifndef QWWAD_MATERIAL_PROPERTY_STRING_H
#define QWWAD_MATERIAL_PROPERTY_STRING_H

#include "qclsim-material-property.h"

class MaterialPropertyString : public MaterialProperty {
public:
    MaterialPropertyString(){}
    MaterialPropertyString(xmlpp::Element *elem);

    const Glib::ustring & get_text() const;

private:
    Glib::ustring _text; ///< The text value of the property
};
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
