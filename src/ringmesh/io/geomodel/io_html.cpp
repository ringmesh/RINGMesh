/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
 * Applications (ASGA) All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <geologyjs/main_export.h>

namespace {

class HTMLIOHandler final : public GeoModelIOHandler<3> {
  public:
    void load(const std::string &filename, GeoModel &geomodel) final {
        throw RINGMeshException(
            "I/O",
            "Geological model loading of a from HTML mesh not yet implemented");
        return false;
    }

    void save(const GeoModel &geomodel, const std::string &filename) final {
        GEOLOGYJS::JSWriter js(filename);
        js.build_js_gui_ = true;

        save_all_lines(geomodel, js);
        save_interfaces(geomodel, js);

        // Check validity and write
        std::string error_message;
        if (js.check_validity(error_message)) {
            js.write();
        } else {
            throw RINGMeshException("I/O", error_message);
        }
    }

  private:
    void save_all_lines(const GeoModel &geomodel,
                        GEOLOGYJS::JSWriter &js) const {
        std::vector<std::vector<double>> xyz;
        xyz.resize(geomodel.nb_lines());
        for (const auto &line : geomodel.lines()) {
            index_t line_id = line.index();
            xyz[line_id].reserve(3 * line.nb_vertices());
            for (index_t v_itr : range(line.nb_vertices())) {
                xyz[line_id].push_back(line.vertex(v_itr).x);
                xyz[line_id].push_back(line.vertex(v_itr).y);
                xyz[line_id].push_back(line.vertex(v_itr).z);
            }
        }
        js.add_lines("all_lines", xyz);
    }

    void save_interfaces(const GeoModel &geomodel,
                         GEOLOGYJS::JSWriter &js) const {
        for (index_t interface_itr : range(geomodel.nb_geological_entities(
                 Interface::type_name_static()))) {
            const GeoModelGeologicalEntity &cur_interface =
                geomodel.geological_entity(Interface::type_name_static(),
                                           interface_itr);
            if (!GeoModelGeologicalEntity::is_stratigraphic_limit(
                    cur_interface.geological_feature()) &&
                !GeoModelGeologicalEntity::is_fault(
                    cur_interface.geological_feature())) {
                continue;
            }

            index_t nb_vertices = 0;
            index_t nb_triangles = 0;
            for (index_t surf_itr : range(cur_interface.nb_children())) {
                const Surface &cur_surface =
                    geomodel.surface(cur_interface.child(surf_itr).index());
                nb_vertices += cur_surface.nb_vertices();
                nb_triangles += cur_surface.nb_mesh_elements();
            }

            std::vector<double> xyz;
            xyz.reserve(3 * nb_vertices);
            std::vector<index_t> indices;
            indices.reserve(3 * nb_triangles);

            index_t vertex_count = 0;
            for (index_t surf_itr : range(cur_interface.nb_children())) {
                const Surface &cur_surface =
                    geomodel.surface(cur_interface.child(surf_itr).index());

                for (index_t v_itr : range(cur_surface.nb_vertices())) {
                    xyz.push_back(cur_surface.vertex(v_itr).x);
                    xyz.push_back(cur_surface.vertex(v_itr).y);
                    xyz.push_back(cur_surface.vertex(v_itr).z);
                }

                for (index_t p_itr : range(cur_surface.nb_mesh_elements())) {
                    for (index_t v_itr : range(3)) {
                        indices.push_back(vertex_count +
                                          cur_surface.mesh_element_vertex_index(
                                              p_itr, v_itr));
                    }
                }

                vertex_count += cur_surface.nb_vertices();
            }
            js.add_surface(cur_interface.name(), xyz, indices);
        }
    }
};

} // namespace
