/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geometry.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/tetra_gen.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/mesh/mesh_geometry.h>

namespace RINGMesh {

    double model_element_size( const GeoModelElement& E )
    {
        double result = 0. ;

        if( E.nb_children() ) {
            /// If this element has children sum up their sizes
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                result += model_element_size( E.child( i ) ) ;
            }
            return result ;
        } else {

            /// Else it is a base element and its size is computed

            switch( E.gme_id().type ) {

                // If this is a region
                case GeoModelElement::REGION: {
                    const Region& R = dynamic_cast< const Region& >( E ) ;
                    // Compute the volume if this is a region
                    for( index_t i = 0; i < R.nb_boundaries(); i++ ) {
                        const Surface& surface =
                            dynamic_cast< const Surface& >( R.boundary( i ) ) ;

                        for( index_t t = 0; t < surface.nb_cells(); t++ ) {
                            const vec3& p0 = surface.vertex( t, 0 ) ;
                            for( index_t v = 1;
                                v + 1 < surface.nb_vertices_in_facet( t ); ++v ) {
                                double cur_volume = ( dot( p0,
                                    cross( surface.vertex( t, v ),
                                        surface.vertex( t, v + 1 ) ) ) )
                                    / static_cast< double >( 6 ) ;
                                R.side( i ) ? result -= cur_volume : result +=
                                                  cur_volume ;
                            }
                        }
                    }
                    return fabs( result ) ;
                }
                case GeoModelElement::SURFACE: {
                    const Surface& S = dynamic_cast< const Surface& >( E ) ;
                    const GEO::Mesh& mesh = S.mesh() ;
                    for( index_t i = 0; i < S.nb_cells(); i++ ) {
                        result += GEO::Geom::mesh_facet_area( mesh, i ) ;
                    }
                    return result ;
                }
                case GeoModelElement::LINE: {
                    const Line& L = dynamic_cast< const Line& >( E ) ;
                    for( index_t i = 1; i < L.nb_vertices(); ++i ) {
                        result += GEO::Geom::distance( L.vertex( i ),
                            L.vertex( i - 1 ) ) ;
                    }
                    return result ;
                }
                case GeoModelElement::CORNER: {
                    return 0 ;
                }
            }
            ringmesh_assert_not_reached;
            return result ;
        }
    }

    double model_element_cell_size( const GeoModelElement& E, index_t c )
    {
        double result = 0. ;

        switch( E.gme_id().type ) {
            case GeoModelElement::REGION: {
                const Region& R = dynamic_cast< const Region& >( E ) ;
                const GEO::Mesh& mesh = R.mesh() ;
                return RINGMesh::mesh_cell_volume( mesh, c ) ;
            }
            case GeoModelElement::SURFACE: {
                const Surface& S = dynamic_cast< const Surface& >( E ) ;
                const GEO::Mesh& mesh = S.mesh() ;
                return GEO::Geom::mesh_facet_area( mesh, c ) ;
            }
            case GeoModelElement::LINE: {
                const Line& L = dynamic_cast< const Line& >( E ) ;
                const GEO::Mesh& mesh = L.mesh() ;
                index_t v0 = mesh.edges.vertex( c, 0 ) ;
                index_t v1 = mesh.edges.vertex( c, 1 ) ;
                return GEO::Geom::distance( L.vertex( v0 ), L.vertex( v1 ) ) ;
            }
        }
        ringmesh_assert_not_reached ;
        return result ;
    }

    vec3 model_element_center( const GeoModelElement& E )
    {
        vec3 result( 0., 0., 0. ) ;
        index_t nb_vertices = 0 ;

        if( E.nb_children() ) {
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                const GeoModelMeshElement& F =
                    dynamic_cast< const GeoModelMeshElement& >( E.child( i ) ) ;
                nb_vertices += F.nb_vertices() ;
                result += model_element_center( F ) * F.nb_vertices() ;
            }
            return result / static_cast< double >( nb_vertices ) ;
        } else {
            const GeoModelMeshElement& F =
                dynamic_cast< const GeoModelMeshElement& >( E ) ;
            for( index_t v = 0; v < F.nb_vertices(); v++ ) {
                result += F.vertex( v ) ;
            }
            return result / static_cast< double >( F.nb_vertices() ) ;
        }
    }

    vec3 model_element_cell_center( const GeoModelElement& E, index_t c )
    {
        vec3 result( 0., 0., 0. ) ;
        index_t nb_vertices = 0 ;

        switch( E.gme_id().type ) {
            case GeoModelElement::REGION: {
                const Region& R = dynamic_cast< const Region& >( E ) ;
                const GEO::Mesh& mesh = R.mesh() ;
                return RINGMesh::mesh_cell_center( mesh, c ) ;
            }
            case GeoModelElement::SURFACE: {
                const Surface& S = dynamic_cast< const Surface& >( E ) ;
                const GEO::Mesh& mesh = S.mesh() ;
                return GEO::Geom::mesh_facet_center( mesh, c ) ;
            }
            case GeoModelElement::LINE: {
                const Line& L = dynamic_cast< const Line& >( E ) ;
                const GEO::Mesh& mesh = L.mesh() ;
                index_t v0 = mesh.edges.vertex( c, 0 ) ;
                index_t v1 = mesh.edges.vertex( c, 1 ) ;
                return 0.5 * ( L.vertex( v0 ), L.vertex( v1 ) ) ;
            }
            case GeoModelElement::CORNER: {
                const Corner& C = dynamic_cast< const Corner& >( E ) ;
                return C.vertex() ;
            }
        }

        ringmesh_assert_not_reached ;
        return result ;
    }

    void translate( GeoModel& M, const vec3& translation_vector )
    {
        // Note: if the translation is null, do nothing.
        if( translation_vector == vec3( 0, 0, 0 ) ) {
            return ;
        }

        for( index_t v = 0; v < M.mesh.vertices.nb(); ++v ) {
            vec3 p = M.mesh.vertices.vertex( v ) ;
            for( index_t i = 0; i < 3; i++ ) {
                p[i] += translation_vector[i] ;
            }
            M.mesh.vertices.update_point( v, p ) ;
        }
    }

    void rotate(
        GeoModel& M,
        const vec3& origin,
        const vec3& axis,
        float64 theta,
        bool degrees )
    {
        // Note: Rotation is impossible about an axis with null length.
        ringmesh_debug_assert( axis != vec3() ) ;
        if( theta == 0. ) {
            return ;
        }

        GEO::Matrix< float64, 4 > rot_mat ;
        rotation_matrix_about_arbitrary_axis( origin, axis, theta, degrees,
            rot_mat ) ;

        for( index_t v = 0; v < M.mesh.vertices.nb(); ++v ) {
            const vec3& p = M.mesh.vertices.vertex( v ) ;

            float64 old[4] = { p[0], p[1], p[2], 1. } ;
            float64 new_p[4] = { 0, 0, 0, 1. } ;
            GEO::mult( rot_mat, old, new_p ) ;
            ringmesh_debug_assert( new_p[3] == 1. ) ;

            M.mesh.vertices.update_point( v, vec3( new_p[0], new_p[1], new_p[2] ) ) ;
        }
    }

    void tetrahedralize(
        GeoModel& M,
        const std::string& method,
        index_t region_id,
        bool add_steiner_points,
        std::vector< std::vector< vec3 > >& internal_vertices )
    {
        if( internal_vertices.empty() ) internal_vertices.resize( M.nb_regions() ) ;
        GEO::Logger::out( "Info" ) << "Using " << method << std::endl ;
        if( region_id == NO_ID ) {
            GEO::ProgressTask progress( "Compute", M.nb_regions() ) ;
            for( index_t i = 0; i < M.nb_regions(); i++ ) {
                TetraGen_var tetragen = TetraGen::create( M.region( i ).mesh(),
                    method ) ;
                tetragen->set_boundaries( M.region( i ), M.wells() ) ;
                tetragen->set_internal_points( internal_vertices[i] ) ;
                GEO::Logger::instance()->set_quiet( true ) ;
                tetragen->tetrahedralize( add_steiner_points ) ;
                GEO::Logger::instance()->set_quiet( false ) ;
                progress.next() ;
            }
        } else {
            TetraGen_var tetragen = TetraGen::create( M.region( region_id ).mesh(),
                method ) ;
            tetragen->set_boundaries( M.region( region_id ), M.wells() ) ;
            tetragen->set_internal_points( internal_vertices[region_id] ) ;
            GEO::Logger::instance()->set_quiet( true ) ;
            tetragen->tetrahedralize( add_steiner_points ) ;
            GEO::Logger::instance()->set_quiet( false ) ;
        }
    }

}
