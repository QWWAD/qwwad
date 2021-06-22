#ifndef QWWAD_MATERIAL_PROPERTY_INTERP_H
#define QWWAD_MATERIAL_PROPERTY_INTERP_H

#include "material-property-numeric.h"

namespace QWWAD {
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
    MaterialPropertyInterp() = default;
    MaterialPropertyInterp(xmlpp::Element *elem);
    MaterialPropertyInterp(const decltype(_name)        &name,
                           const decltype(_description) &description,
                           const decltype(_reference)   &reference,
                           decltype(_unit)               unit,
                           decltype(_y0)                 y0,
                           decltype(_y1)                 y1,
                           decltype(_b)                  b = 0.0);

    [[nodiscard]] auto clone() const -> MaterialPropertyInterp * override;

    void set_limits(decltype(_xmin) xmin,
                    decltype(_xmax) xmax);

    void get_limits(decltype(_xmin) &xmin,
                    decltype(_xmax) &xmax) const;

    [[nodiscard]] inline auto get_interp_y0() const {return _y0;}
    [[nodiscard]] inline auto get_interp_y1() const {return _y1;}
    [[nodiscard]] inline auto get_interp_b()  const {return _b;}

    [[nodiscard]] auto get_val(double x = 0) const -> decltype(_y0) override;
};
} // end namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
