#ifndef QWWAD_MATERIAL_PROPERTY_CONSTANT_H
#define QWWAD_MATERIAL_PROPERTY_CONSTANT_H

#include "qwwad-material-property-numeric.h"

class MaterialPropertyConstant : public MaterialPropertyNumeric {
public:
    MaterialPropertyConstant(){}
    MaterialPropertyConstant(xmlpp::Element *elem);

    double get_val(const double x = 0) const;

private:
    double _constant; ///< The constant value of the property
};
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
