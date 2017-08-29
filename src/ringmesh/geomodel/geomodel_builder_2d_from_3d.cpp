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

}

namespace RINGMesh {

    void GeoModelBuilder2DProjection::build_geomodel()
    {
        copy_geomodel_3d_topology();
        project_geomodel_3d_mesh_entities();
    }

    void GeoModelBuilder2DProjection::copy_geomodel_3d_topology()
    {
        for( const auto& entity_type : projectable_entity_types ) {
            topology.create_mesh_entities( entity_type,
                geomodel3d_from_.nb_mesh_entities( entity_type ) );
        }
        //@todo boundary relations
    }

    void GeoModelBuilder2DProjection::project_geomodel_3d_mesh_entities()
    {
        for( const auto& entity_type : projectable_entity_types ) {
            for( const auto entity_id : range(
                geomodel3d_from_.nb_mesh_entities( entity_type ) ) ) {
                std::vector< vec2 > projected_vertices = compute_projected_vertices(
                    geomodel3d_from_.mesh_entity( entity_type, entity_id ) );
                geometry.set_mesh_entity_vertices(
                    gmme_id { entity_type, entity_id }, projected_vertices, false );
            }
        }
    }

    std::vector< vec2 > GeoModelBuilder2DProjection::compute_projected_vertices(
        const GeoModelMeshEntity3D& entity )
    {
        std::vector< vec2 > projected_vertices;
        projected_vertices.resize( entity.nb_vertices() );
        for( const auto v : range( entity.nb_vertices() ) ) {
            projected_vertices.push_back( get_2d_coord( entity.vertex( v ) ) );

        }
        return projected_vertices;
    }

}
