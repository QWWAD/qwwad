#ifndef QWWAD_MATERIAL_PROPERTY_POLY_H
#define QWWAD_MATERIAL_PROPERTY_POLY_H

#include <map>
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

    [[nodiscard]] auto clone() const -> MaterialPropertyPoly * override;

    [[nodiscard]] auto get_val(const double x = 0) const -> double override;
};
} // end namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
