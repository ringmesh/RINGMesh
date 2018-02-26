/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <array>
#include <iomanip>
#include <iostream>

#include <geogram/basic/command_line.h>
#include <geogram/basic/progress.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_entity.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>


#include <ringmesh/tetrahedralize/tetra_gen.h>

/*!
 * @file Set of high level API functions
 */

namespace RINGMesh
{
    template < index_t DIMENSION >
    void geomodel_tools_api copy_geomodel(
        const GeoModel< DIMENSION >& from, GeoModel< DIMENSION >& to )
    {
        GeoModelBuilder< DIMENSION > to_builder( to );
        to_builder.topology.copy_topology( from );
        to_builder.geometry.copy_meshes( from );
        to_builder.geology.copy_geology( from );
    }

    template < index_t DIMENSION >
    void translate( GeoModel< DIMENSION >& geomodel,
        const vecn< DIMENSION >& translation_vector )
    {
        GeoModelBuilder< DIMENSION > builder( geomodel );
        for( auto v : range( geomodel.mesh.vertices.nb() ) )
        {
            // Coordinates are not directly modified to
            // update the matching vertices in geomodel entities
            const auto& p = geomodel.mesh.vertices.vertex( v );
            builder.geometry.set_mesh_entity_vertex( v, p + translation_vector );
        }
    }

    void rotate( GeoModel3D& geomodel,
        const vec3& origin,
        const vec3& axis,
        double angle,
        bool degrees )
    {
        if( length( axis ) < geomodel.epsilon() )
        {
            Logger::err( "GeoModel",
                "Rotation around an epsilon length axis is impossible" );
            return;
        }
        GEO::Matrix< 4, double > rot_mat{ rotation_matrix_about_arbitrary_axis(
            origin, axis, angle, degrees ) };

        GeoModelBuilder3D builder( geomodel );
        for( auto v : range( geomodel.mesh.vertices.nb() ) )
        {
            const vec3& p = geomodel.mesh.vertices.vertex( v );
            std::array< double, 4 > old{ { p[0], p[1], p[2], 1. } };
            std::array< double, 4 > new_p{ { 0, 0, 0, 1. } };
            GEO::mult( rot_mat, old.data(), new_p.data() );
            ringmesh_assert( std::fabs( new_p[3] - 1. ) < global_epsilon );

            builder.geometry.set_mesh_entity_vertex( v, vec3{ new_p.data() } );
        }
    }

    void tetrahedralize(
        GeoModel3D& geomodel, index_t region_id, bool add_steiner_points )
    {
        std::vector< std::vector< vec3 > > internal_vertices(
            geomodel.nb_regions() );
        tetrahedralize(
            geomodel, region_id, add_steiner_points, internal_vertices );
    }

    void tetrahedralize( GeoModel3D& geomodel,
        index_t region_id,
        bool add_steiner_points,
        const std::vector< std::vector< vec3 > >& internal_vertices )
    {
        const std::string method{ GEO::CmdLine::get_arg( "algo:tet" ) };
        if( region_id == NO_ID )
        {
            Logger::out( "Info", "Using ", method );
            GEO::ProgressTask progress( "Compute", geomodel.nb_regions() );
            for( auto i : range( geomodel.nb_regions() ) )
            {
                tetrahedralize(
                    geomodel, i, add_steiner_points, internal_vertices );
                progress.next();
            }
        }
        else
        {
            std::unique_ptr< TetraGen > tetragen{ TetraGen::create(
                geomodel, region_id, method ) };
            tetragen->set_boundaries(
                geomodel.region( region_id ), geomodel.wells() );
            tetragen->set_internal_points( internal_vertices[region_id] );
            bool status{ Logger::instance()->is_quiet() };
            Logger::instance()->set_quiet( true );
            tetragen->tetrahedralize( add_steiner_points );
            Logger::instance()->set_quiet( status );
        }

        // The GeoModelMesh should be updated, just erase everything
        // and it will be re-computed during its next access.
        geomodel.mesh.vertices.clear();
    }

    template void geomodel_tools_api copy_geomodel( const GeoModel2D&, GeoModel2D& );
    template void geomodel_tools_api translate( GeoModel2D&, const vec2& );

    template void geomodel_tools_api copy_geomodel( const GeoModel3D&, GeoModel3D& );
    template void geomodel_tools_api translate( GeoModel3D&, const vec3& );

} // namespace RINGMesh
