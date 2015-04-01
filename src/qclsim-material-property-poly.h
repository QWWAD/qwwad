#ifndef QCLSIM_MATERIAL_PROPERTY_POLY_H
#define QCLSIM_MATERIAL_PROPERTY_POLY_H

#include "qclsim-material-property.h"

class MaterialPropertyPoly : public MaterialProperty {
public:
    MaterialPropertyPoly(xmlpp::Element *elem);
    double get_val(const double x = 0) const;

private:
    // Values for polynomial properties
    std::vector<int>    _poly_index; ///< The index of each term
    std::vector<double> _poly_coeff; ///< The coefficient of each term
};
#endif // QCLSIM_MATERIAL_PROPERTY_POLY_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
