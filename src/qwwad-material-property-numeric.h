#ifndef QWWAD_MATERIAL_PROPERTY_NUMERIC_H
#define QWWAD_MATERIAL_PROPERTY_NUMERIC_H

#include "qclsim-material-property.h"

/**
 * A physical property of a material that can be described by a numerical value
 */
class MaterialPropertyNumeric : public MaterialProperty {
public:
    MaterialPropertyNumeric(){}
    MaterialPropertyNumeric(xmlpp::Element *elem);

    const Glib::ustring & get_unit() const;

    virtual double get_val(const double x = 0) const = 0;

protected:
    Glib::ustring _unit; ///< The unit of the property
};
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
