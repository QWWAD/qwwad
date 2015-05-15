#ifndef QWWAD_MATERIAL_PROPERTY_CONSTANT_H
#define QWWAD_MATERIAL_PROPERTY_CONSTANT_H

#include "qwwad-material-property-numeric.h"

class MaterialPropertyConstant : public MaterialPropertyNumeric {
private:
    double _constant; ///< The constant value of the property

public:
    MaterialPropertyConstant(){}
    MaterialPropertyConstant(xmlpp::Element *elem);
    MaterialPropertyConstant(decltype(_name)        name,
                             decltype(_description) description,
                             decltype(_reference)   reference,
                             decltype(_unit)        unit,
                             decltype(_constant)    value);

    decltype(_constant) get_val(const double x = 0) const;
};
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
