#ifndef QWWAD_MATERIAL_PROPERTY_POLY_H
#define QWWAD_MATERIAL_PROPERTY_POLY_H

#include "material-property-numeric.h"

namespace QWWAD {
class MaterialPropertyPoly : public MaterialPropertyNumeric {
private:
    std::map<int, double> _poly_coeffs; ///< The terms in the polynomial

public:
    MaterialPropertyPoly(xmlpp::Element *elem);
    MaterialPropertyPoly(decltype(_name)        name,
                         decltype(_description) description,
                         decltype(_reference)   reference,
                         decltype(_unit)        unit,
                         decltype(_poly_coeffs) poly_coeffs);

    virtual MaterialPropertyPoly * clone() const;

    double get_val(const double x = 0) const;
};
} // end namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
