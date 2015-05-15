#ifndef QCLSIM_MATERIAL_PROPERTY_INTERP_H
#define QCLSIM_MATERIAL_PROPERTY_INTERP_H

#include "qwwad-material-property-numeric.h"

/**
 * A MaterialProperty that is interpolated between two values
 *
 * \details Values \f$y_0\f$ and \f$y_1\f$ are stored for the property, which correspond to \f$x=0\f$ and \f$x=1\f$
 *          respectively.
 *          An optional bowing parameter $b$ may be provided.
 *          The value of the property is given by:
 *          \f[
 *          y(x) = y_0(1-x) + y_1(x) + bx(1-x)
 *          \f]
 */
class MaterialPropertyInterp : public MaterialPropertyNumeric {
private:
    // Values for interpolated properties
    double _y0; ///< The value when input variable is 0
    double _y1; ///< The value when input variable is 1
    double _b;  ///< The optional bowing factor

    double _xmin; ///< The minimum value of x that can be used in value lookup
    double _xmax; ///< The maximum value of x that can be used in value lookup

public:
    MaterialPropertyInterp(){}
    MaterialPropertyInterp(xmlpp::Element *elem);
    MaterialPropertyInterp(decltype(_name)        name,
                           decltype(_description) description,
                           decltype(_reference)   reference,
                           decltype(_unit)        unit,
                           decltype(_y0)          y0,
                           decltype(_y1)          y1,
                           decltype(_b)           b = 0.0);

    void set_limits(const decltype(_xmin) xmin,
                    const decltype(_xmax) xmax);

    void get_limits(decltype(_xmin) &xmin,
                    decltype(_xmax) &xmax);

    inline decltype(_y0) get_interp_y0() const {return _y0;}
    inline decltype(_y1) get_interp_y1() const {return _y1;}
    inline decltype(_b)  get_interp_b()  const {return _b;}
    decltype(_y0) get_val(const double x = 0) const;
};
#endif // QCLSIM_MATERIAL_PROPERTY_INTERP_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
