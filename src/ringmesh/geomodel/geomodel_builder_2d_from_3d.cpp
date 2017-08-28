/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/geomodel/geomodel_builder_2d_from_3d.h>

namespace {
    using namespace RINGMesh;

    static const MeshEntityType projectable_entity_types[3] = {
        Corner3D::type_name_static(), Line3D::type_name_static(),
        Surface3D::type_name_static() };

    void copy_geomodel_3d_topology(
        const GeoModel3D& geomodel_from,
        GeoModelBuilder2D& geomodel_to_builder )
    {
        for( const auto& entity_type : projectable_entity_types ) {
            geomodel_to_builder.topology.create_mesh_entities( entity_type,
                geomodel_from.nb_mesh_entities( entity_type ) );
        }
        //@todo boundary relations
    }

    std::vector< vec2 > compute_projected_vertices(
        const GeoModelMeshEntity3D& entity,
        Geometry::Plane projection_plane )
    {
        std::vector< vec2 > projected_vertices;
        projected_vertices.resize( entity.nb_vertices() );
        for( const auto v : range( entity.nb_vertices() ) ) {
            vec3 projection_3d_on_plane;
            std::tie( std::ignore, projection_3d_on_plane ) =
                Distance::point_to_plane( entity.vertex( v ), projection_plane );
        }
        return projected_vertices;
    }

    void project_geomodel_3d_mesh_entities(
        const GeoModel3D& geomodel_from,
        const Geometry::Plane& projection_plane,
        GeoModelBuilder2D& geomodel_to_builder )
    {
        for( const auto& entity_type : projectable_entity_types ) {
            for( const auto entity_id : range(
                geomodel_from.nb_mesh_entities( entity_type ) ) ) {
                std::vector< vec2 > projected_vertices = compute_projected_vertices(
                    geomodel_from.mesh_entity( entity_type, entity_id ),
                    projection_plane );
                geomodel_to_builder.geometry.set_mesh_entity_vertices( gmme_id {
                    entity_type, entity_id }, projected_vertices );
            }
        }
    }
}

namespace RINGMesh {

    void GeoModelBuilder2DProjection::build_geomodel()
    {
        copy_geomodel_3d_topology( geomodel3d_from_, *this );
        project_geomodel_3d_mesh_entities( geomodel3d_from_, plane_, *this );
    }

}
