/*
*  Copyright (c) 2012-2014, Bruno Levy
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*
*  * Redistributions of source code must retain the above copyright notice,
*  this list of conditions and the following disclaimer.
*  * Redistributions in binary form must reproduce the above copyright notice,
*  this list of conditions and the following disclaimer in the documentation
*  and/or other materials provided with the distribution.
*  * Neither the name of the ALICE Project-Team nor the names of its
*  contributors may be used to endorse or promote products derived from this
*  software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
*  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
*  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
*  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
*  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
*  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
*  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
*  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
*  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
*  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: Bruno Levy
*
*     Bruno.Levy@inria.fr
*     http://www.loria.fr/~levy
*
*     ALICE Project
*     LORIA, INRIA Lorraine,
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX
*     FRANCE
*
*
*  MODIFIED BY 
* 
* Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
* All rights reserved.
*
*     http://www.ring-team.org
*
*     RING Project
*     Ecole Nationale Superieure de Geologie - GeoRessources
*     2 Rue du Doyen Marcel Roubault - TSA 70605
*     54518 VANDOEUVRE-LES-NANCY
*     FRANCE
*/
#ifndef __RINGMESH_GEOGRAM_MESH_REPAIR__
#define __RINGMESH_GEOGRAM_MESH_REPAIR__

#include <algorithm>

#include <geogram/basic/memory.h>
#include <geogram/mesh/mesh.h>

#include <ringmesh/common.h>


/* 
 * @file High level repair operations on GEO::Mesh modified from geogram/mesh/mesh_repair.cpp
 * @author Jeanne Pellerin
 */

namespace RINGMesh {

    /**
    * \brief Predicate to disable connection of facets in repair_connect_facet
    */
    class MeshEdgesOnBorder
    {
    public:
        typedef std::pair< index_t, index_t > Edge ;
       /*!
        * @param[in] mesh_edges_on_border Pair of mesh vertex indices on a boundary.
        *            In each pair first < second.
        */
        MeshEdgesOnBorder( const std::vector<Edge>& mesh_edges_on_border )
            : edges_on_border_( mesh_edges_on_border.begin(), mesh_edges_on_border.end() )
        {}

        bool operator()( index_t v0, index_t v1 ) const
        {
            Edge edge( std::min( v0, v1 ), std::max( v0, v1 ) ) ;
            return edges_on_border_.find( edge ) != edges_on_border_.end() ;
        }
        void debug( index_t v0, index_t v1 ) const
        {
            ringmesh_assert_not_reached ;
        }
    private:
        const std::set< Edge > edges_on_border_ ;
    } ;


    /*!
    * @brief Connects the facets in a TRIANGULATED GEO::Mesh.
    * @details Reconstructs the corners.adjacent_facet links.
    *          Orientation is not checked
    * @note Modified from geogram to take into account a predicate to disconnect facets
    *       along identified edges [JP]
    *       The predicate should implement
    *            bool operator() (index_t v1, index_t v2) const ;
    *            void debug(index_t v1, index_t v2) ;
    *
    * @note It is template and inline on purpose - Theorem from Bruno Levy.
    * 
    * @warning DIFFICULT NOT CLEAN CODE
    */
    template< typename P > inline
    void repair_connect_facets( GEO::Mesh& M, P is_border )
    {
        using GEO::index_t ;

        const index_t NO_FACET = index_t( -1 ) ;
        const index_t NO_CORNER = index_t( -1 ) ;
        const index_t NON_MANIFOLD = index_t( -2 ) ;

        // Reset all facet-facet adjacencies.
        for( index_t c = 0; c < M.facet_corners.nb(); ++c ) {
            M.facet_corners.set_adjacent_facet( c, NO_FACET ) ;
        }

        // For each vertex v, v2c[v] gives the index of a 
        // corner incident to vertex v.
        GEO::vector< index_t > v2c( M.vertices.nb(), NO_CORNER ) ;

        // For each corner c, next_c_around_v[c] is the 
        // linked list of all the corners incident to 
        // vertex v.
        GEO::vector< index_t > next_c_around_v( M.facet_corners.nb(), NO_CORNER ) ;

        // Compute v2c and next_c_around_v
        for( index_t c = 0; c < M.facet_corners.nb(); ++c ) {
            index_t v = M.facet_corners.vertex( c ) ;
            next_c_around_v[ c ] = v2c[ v ] ;
            v2c[ v ] = c ;
        }

        for( index_t f1 = 0; f1 < M.facets.nb(); ++f1 ) {
            for( index_t c1 = M.facets.corners_begin( f1 );
                 c1 < M.facets.corners_end( f1 ); ++c1 ) {

                if( M.facet_corners.adjacent_facet( c1 ) == NO_FACET ) {
                    index_t adj_corner = NO_CORNER ;
                    index_t v1 = M.facet_corners.vertex( c1 ) ;
                    index_t v2 = M.facet_corners.vertex(
                        M.facets.next_corner_around_facet( f1, c1 ) ) ;

                    index_t c2 = v2c[ v1 ] ;

                    // Lookup candidate adjacent edges from incident
                    // edges list.
                    while( c2 != NO_CORNER ) {
                        if( c2 != c1 ) {
                            index_t f2 = c2 / 3 ;
                            index_t c3 = M.facets.prev_corner_around_facet( f2, c2 ) ;
                            index_t v3 = M.facet_corners.vertex( c3 ) ;
                            // Check with standard orientation.
                            if( v3 == v2 ) {
                                if( !is_border( v1, v2 ) ) {
                                    if( adj_corner == NO_CORNER ) {
                                        adj_corner = c3 ;
                                    } else {
                                        // Unexpected non-manifold edge store it
                                        is_border.debug( v1, v2 ) ;
                                        adj_corner = NON_MANIFOLD ;
                                    }
                                }
                            } else {
                                // Check with the other ("wrong") orientation
                                c3 = M.facets.next_corner_around_facet( f2, c2 ) ;
                                v3 = M.facet_corners.vertex( c3 ) ;
                                if( v3 == v2 ) {
                                    if( !is_border( v1, v2 ) ) {
                                        if( adj_corner == NO_CORNER ) {
                                            adj_corner = c2 ;
                                        } else {
                                            // Unexpected non-manifold edge store it
                                            is_border.debug( v1, v2 ) ;
                                            adj_corner = NON_MANIFOLD ;
                                        }
                                    }
                                }
                            }
                        }
                        c2 = next_c_around_v[ c2 ] ;
                    }
                    if( adj_corner != NO_CORNER && adj_corner != NON_MANIFOLD ) {
                        M.facet_corners.set_adjacent_facet( adj_corner, f1 ) ;
                        index_t f2 = adj_corner / 3 ;
                        M.facet_corners.set_adjacent_facet( c1, f2 ) ;
                    }
                }
            }
        }
    }


    /*!
     * @brief Connects the facets of a Mesh except along edges flagged so.
     * @note inline otherwise linking trouble on Windows. No idea why [JP]
     */
    inline void connect_mesh_facets( GEO::Mesh& mesh,
        const std::vector<std::pair<index_t, index_t> >& border_edges )
    {
        MeshEdgesOnBorder edges_to_keep_disconnected( border_edges ) ;
        return repair_connect_facets< MeshEdgesOnBorder >( mesh, edges_to_keep_disconnected ) ;
    }

    inline void connect_mesh_facets_except_on_mesh_edges( GEO::Mesh& mesh )
    {
        index_t nb_edges = mesh.edges.nb() ;
        std::vector< std::pair<index_t, index_t > > edges( nb_edges ) ;
        for( index_t i = 0; i< nb_edges; ++i ) {
            index_t v0 = mesh.edges.vertex( i, 0 );
            index_t v1 = mesh.edges.vertex( i, 1 );
            if( v1 < v0 ) {
                std::swap( v0, v1 ) ;
            }
            edges[ i ] = std::pair<index_t, index_t>( v0, v1 ) ;
        }
        connect_mesh_facets( mesh, edges ) ;
    }
   

    /*!
    * @brief Returns true if there are colocated vertices in the Mesh
    * @details This is a wrapper around Geogram colocate functions
    */
    bool RINGMESH_API has_mesh_colocate_vertices( const GEO::Mesh& M, double tolerance ) ;

    /*!
    * @brief Merges the vertices of a mesh that are at the same geometric location
    * @note Copied from geogram/mes/mesh_repair.cpp. No choice since BL will not give access to it.
    */
    void RINGMESH_API repair_colocate_vertices( GEO::Mesh& M, double tolerance ) ;

}

#endif