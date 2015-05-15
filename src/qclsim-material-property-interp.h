#ifndef QCLSIM_MATERIAL_PROPERTY_INTERP_H
#define QCLSIM_MATERIAL_PROPERTY_INTERP_H

#include "qwwad-material-property-numeric.h"

class MaterialPropertyInterp : public MaterialPropertyNumeric {
public:
    MaterialPropertyInterp(){}
    MaterialPropertyInterp(xmlpp::Element *elem);

    double get_interp_y0()   const;
    double get_interp_y1()   const;
    double get_interp_b()    const;
    bool   interp_is_bowed() const;
    double get_val(const double x = 0) const;

private:
    // Values for interpolated properties
    double _interp_y0; ///< The value when input variable is 0
    double _interp_y1; ///< The value when input variable is 1

    bool   _interp_is_bowed; ///< Whether a bowing factor is used
    double _interp_b;  ///< The optional bowing factor

    double _xmin; ///< The minimum value of x that can be used in value lookup
    double _xmax; ///< The maximum value of x that can be used in value lookup
};
#endif // QCLSIM_MATERIAL_PROPERTY_INTERP_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
