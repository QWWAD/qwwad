#ifndef QWWAD_MATERIAL_PROPERTY_POLY_H
#define QWWAD_MATERIAL_PROPERTY_POLY_H

#include "qwwad-material-property-numeric.h"

namespace QWWAD {
class MaterialPropertyPoly : public MaterialPropertyNumeric {
public:
    MaterialPropertyPoly(xmlpp::Element *elem);
    double get_val(const double x = 0) const;

private:
    // Values for polynomial properties
    std::map<int, double> _poly_coeffs; ///< The terms in the polynomial
};
} // end namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
