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

#include <ringmesh/geomodel/geomodel_builder_from_mesh.h>

#include <stack>

/*!
 * @file ringmesh/geomodel/geomodel_builder.cpp
 * @brief Implementation of the classes to build GeoModel from various meshes
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    /*! @details Adds separately each connected component of the mesh
     *          as a Surface of the geomodel under construction.
     *          All the facets of the input mesh are visited and added to a
     *          Surface of the GeoModel.
     *          Connected components of the mesh are determined with a
     *          propagation (or "coloriage" algorithm) using the adjacent_facet
     *          information provided on the input GEO::Mesh.
     *
     * @todo Old code - old building - to delimit connected components
     * vertices are duplicated in the input mesh
     *
     */
    void GeoModelBuilderSurfaceMesh::build_polygonal_surfaces_from_connected_components()
    {
        std::vector< index_t > global_vertex_id_to_id_in_cc( mesh_.vertices.nb(),
            NO_ID ) ;

        std::vector< bool > visited( mesh_.facets.nb(), false ) ;
        for( index_t i = 0; i < mesh_.facets.nb(); i++ ) {
            if( !visited[i] ) {
                std::vector< index_t > cc_corners ;
                std::vector< index_t > cc_facets_ptr ;
                std::vector< vec3 > cc_vertices ;

                /// @todo Review : This should not be necessary as each vertex should
                /// be in one and only one connected component. To test. [JP]
                std::fill( global_vertex_id_to_id_in_cc.begin(),
                    global_vertex_id_to_id_in_cc.end(), NO_ID ) ;

                // First facet begin at corner 0
                cc_facets_ptr.push_back( 0 ) ;

                // Propagate from facet #i 
                std::stack< index_t > S ;
                S.push( i ) ;
                while( !S.empty() ) {
                    index_t f = S.top() ;
                    S.pop() ;
                    visited[f] = true ;

                    for( index_t c = mesh_.facets.corners_begin( f );
                        c < mesh_.facets.corners_end( f ); ++c ) {
                        index_t v = mesh_.facet_corners.vertex( c ) ;
                        if( global_vertex_id_to_id_in_cc[v] == NO_ID ) {
                            index_t index =
                                static_cast< index_t >( cc_vertices.size() ) ;
                            global_vertex_id_to_id_in_cc[v] = index ;
                            cc_vertices.push_back( mesh_.vertices.point( v ) ) ;
                        }
                        cc_corners.push_back( global_vertex_id_to_id_in_cc[v] ) ;

                        index_t n = mesh_.facet_corners.adjacent_facet( c ) ;
                        if( n != NO_ID && !visited[n] ) {
                            visited[n] = true ;
                            S.push( n ) ;
                        }
                    }
                    index_t nb_cc_corners =
                        static_cast< index_t >( cc_corners.size() ) ;
                    cc_facets_ptr.push_back( nb_cc_corners ) ;
                }

                gmme_t surface_gme = topology.create_mesh_entity< Surface >() ;
                geometry.set_surface_geometry( surface_gme.index(), cc_vertices,
                    cc_corners, cc_facets_ptr ) ;
            }
        }
    }

} // namespace
