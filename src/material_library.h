/**
 * \file     material_library.h 
 * \brief    Fundamental material data and interpolation functions
 *
 * \details  Note that this is only supposed to contain the kind of
 *           things that would appear in a table of parameters in a
 *           paper.  Derived properties like the strain tensor
 *           elements for a specific structure really belong somewhere
 *           else!
 *
 * \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef MATERIAL_LIBRARY_H
#define MATERIAL_LIBRARY_H

#include <boost/ptr_container/ptr_map.hpp>
#include <libxml++/libxml++.h>

class Material;
class MaterialProperty;

namespace Glib {
class ustring;
}

/** Library of material data */
class MaterialLibrary {
public:
    MaterialLibrary(const Glib::ustring &filename);
    ~MaterialLibrary();

    Material * get_material(const char          *mat_name);
    Material * get_material(const Glib::ustring &mat_name);

    MaterialProperty * get_property(Glib::ustring &mat_name,
                                    Glib::ustring &property_name);

    double get_val(Glib::ustring &mat_name,
                   Glib::ustring &property_name);

    const Glib::ustring & get_property_unit(Glib::ustring &mat_name,
                                            Glib::ustring &property_name);
private:
    boost::ptr_map<Glib::ustring, Material>  materials;
    xmlpp::Node::NodeList                    material_nodes;
    xmlpp::DomParser                        *parser;
    xmlpp::Document                         *doc;
    xmlpp::Element                          *root_element;
};
#endif //MATERIAL_LIBRARY_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
