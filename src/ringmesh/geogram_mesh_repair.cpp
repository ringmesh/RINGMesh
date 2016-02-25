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

#include <ringmesh/geogram_mesh_repair.h>

#include <algorithm>

#include <geogram/basic/memory.h>
#include <geogram/mesh/mesh.h>
#include <geogram/points/colocate.h>


/*!
* @file Implementation of high level repairing functions on GEO::Mesh 
* @note Code modified from geogram/mesh/mesh_repair.cpp
* @author Jeanne Pellerin
*/

namespace {
    using RINGMesh::index_t ;
    using RINGMesh::NO_ID ;

    struct Edge {
        Edge( index_t v0, index_t v1 ) :
            vertices_( std::min( v0, v1 ), std::max( v0, v1 ) )
        {}
        Edge( std::pair<index_t, index_t> edge ) :
            Edge( edge.first, edge.second )
        {}
        Edge() : vertices_( NO_ID, NO_ID )
        {} 
        bool operator<( const Edge& rhs ) const
        {
            if( v0() == rhs.v0() ) {
                return v1() < rhs.v1() ;
            } else {
                return v0() < rhs.v0() ;
            }
        }
        index_t v0() const
        {
            return vertices_.first ;
        }
        index_t v1() const
        {
            return vertices_.second ;
        }
       
    private:
        std::pair<index_t, index_t> vertices_ ;
    };

    /*!
    * @brief
    */
    class MeshEdgesOnBorder {
    public:
        /*!
        * @param[in] mesh_edges_on_border Pair of mesh vertex indices on a boundary.
        *            In each pair first < second.
        */
        MeshEdgesOnBorder( const std::vector< index_t >& mesh_edges_on_border )
        {            
            for( index_t i = 0; i+1 < mesh_edges_on_border.size(); i+=2 ) {
                index_t v0 = mesh_edges_on_border[ i ];
                index_t v1 = mesh_edges_on_border[ i+1 ];
                edges_on_border_.insert( Edge( v0, v1 ) );
            }
        }

        bool operator()( index_t v0, index_t v1 ) const
        {
            Edge edge( v0, v1 ) ;
            return edges_on_border_.find( edge ) != edges_on_border_.end() ;
        }
        void non_manifold( index_t v0, index_t v1 ) 
        {
            non_manifold_edges_.insert( Edge( v0, v1 ) ) ;
        }
        bool non_manifold_edges( std::vector< index_t >& edges ) const
        {
            edges.clear() ;
            if( !non_manifold_edges_.empty() ) {
                edges.resize( 2*non_manifold_edges_.size() ) ;
                index_t count = 0 ;
                for( std::set<Edge>::const_iterator itr( non_manifold_edges_.begin() );
                     itr!= non_manifold_edges_.end(); ++itr ) {
                    edges[ count ] = itr->v0() ;
                    edges[ count+1 ] = itr->v1();
                    count += 2 ;
                }
                return true ;
            } else {
                return false ;
            }
        }
    private:
        std::set< Edge > edges_on_border_ ;
        std::set< Edge > non_manifold_edges_ ;
    } ;
   
    class DoNothing {
    public:
        bool operator()( index_t v0, index_t v1 ) const
        {
            return false ;
        }
        void non_manifold( index_t v0, index_t v1 ) const
        {}
    } ;

    /*!
    * @brief Connects the facets in a TRIANGULATED GEO::Mesh.
    * @details Reconstructs the corners.adjacent_facet links.
    *          Orientation is not checked
    * @note Modified from geogram to take into account a predicate to disconnect facets
    *       along identified edges [JP]
    *       The predicate should implement
    *            bool operator() (index_t v1, index_t v2) const ;
    *            void non_manifold(index_t v1, index_t v2) ;
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
                                        is_border.non_manifold( v1, v2 ) ;
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
                                            is_border.non_manifold( v1, v2 ) ;
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

    void get_mesh_edges( const GEO::Mesh& mesh, std::vector< index_t >& edges )
    {
        index_t nb_edges = mesh.edges.nb() ;
        edges.resize( 2*nb_edges ) ;
        for( index_t i = 0; i < nb_edges; ++i ) {
            edges[2*i]= mesh.edges.vertex( i, 0 );
            edges[2*i+1]= mesh.edges.vertex( i, 1 );            
        }
    }
} // anonymous namespace

namespace RINGMesh {
    
    void connect_mesh_facets( GEO::Mesh& mesh,
        const std::vector< index_t >& border_edges )
    {
        MeshEdgesOnBorder edges_to_keep_disconnected( border_edges ) ;
        return repair_connect_facets< MeshEdgesOnBorder >( mesh, edges_to_keep_disconnected ) ;
    }

    void connect_mesh_facets_except_on_mesh_edges( GEO::Mesh& mesh )
    {
        std::vector<index_t> edges ;
        get_mesh_edges( mesh, edges ) ;
        connect_mesh_facets( mesh, edges ) ;
    }

    void connect_mesh_facets_except_on_mesh_edges( GEO::Mesh& mesh,
        std::vector<index_t>& non_manifold_edges )
    {
        std::vector< index_t > edges ;
        get_mesh_edges( mesh, edges ) ;
        
        MeshEdgesOnBorder border_edges( edges ) ;
        repair_connect_facets( mesh, border_edges ) ;
        border_edges.non_manifold_edges( non_manifold_edges ) ;        
    }                                                   

    void connect_facets( GEO::Mesh& mesh )
    {
        DoNothing predicate ;
        return repair_connect_facets( mesh, predicate ) ;
    }


    /*******************************************************************************/

    index_t detect_mesh_colocated_vertices(
        const GEO::Mesh& M, double tolerance, GEO::vector< index_t >& old2new )
    {
        index_t nb_unique_vertices = 0 ;
        if( tolerance == 0.0 ) {
            nb_unique_vertices = GEO::Geom::colocate_by_lexico_sort(
                M.vertices.point_ptr( 0 ), 3, M.vertices.nb(), old2new,
                M.vertices.dimension() ) ;
        } else {
            nb_unique_vertices = GEO::Geom::colocate(
                M.vertices.point_ptr( 0 ), 3, M.vertices.nb(), old2new,
                tolerance, M.vertices.dimension() ) ;
        }
        index_t nb_colocated_vertices( M.vertices.nb() - nb_unique_vertices ) ;
        return nb_colocated_vertices ;
    }

    bool has_mesh_colocate_vertices( const GEO::Mesh& M, double tolerance )
    {
        GEO::vector< index_t > old2new ;
        index_t nb_colocated_vertices = detect_mesh_colocated_vertices( M, tolerance, old2new ) ;
        if( nb_colocated_vertices == 0 ) {
            return false ;
        } else {
            return true ;
        }
    }

    void update_mesh_edges_vertices( GEO::Mesh& M, const GEO::vector<index_t>& old2new )
    {
        for( index_t e = 0; e < M.edges.nb(); ++e ) {
            M.edges.set_vertex( e, 0, old2new[ M.edges.vertex( e, 0 ) ] );
            M.edges.set_vertex( e, 1, old2new[ M.edges.vertex( e, 1 ) ] );
        }
    }

    void update_mesh_facets_vertices( GEO::Mesh& M, const GEO::vector<index_t>& old2new )
    {
        for( index_t c = 0; c < M.facet_corners.nb(); ++c ) {
            M.facet_corners.set_vertex( c, old2new[ M.facet_corners.vertex( c ) ] );
        }
    }

    void update_mesh_cells_vertices( GEO::Mesh& M, const GEO::vector<index_t>& old2new )
    {
        for( index_t ce = 0; ce < M.cells.nb(); ++ce ) {
            for( index_t c = M.cells.corners_begin( ce );
                 c<M.cells.corners_end( ce ); ++c
                 ) {
                M.cell_corners.set_vertex( c, old2new[ M.cell_corners.vertex( c ) ] );
            }
        }
    }

    void delete_colocated_vertices( GEO::Mesh& M, GEO::vector< index_t >& old2new )
    {
        for( index_t i = 0; i < old2new.size(); i++ ) {
            if( old2new[ i ] == i ) {
                old2new[ i ] = 0;
            } else {
                old2new[ i ] = 1;
            }
        }
        M.vertices.delete_elements( old2new );
    }


    void repair_colocate_vertices( GEO::Mesh& M, double colocate_epsilon )
    {
        GEO::vector<index_t> old2new;
        index_t nb_colocated_vertices = detect_mesh_colocated_vertices( M, colocate_epsilon, old2new ) ;
        if( nb_colocated_vertices == 0 ) {
            return ;
        }

        GEO::Logger::out( "GeoModel" ) << "Removing "
            << nb_colocated_vertices
            << " duplicated vertices" << std::endl;

        update_mesh_edges_vertices( M, old2new ) ;
        update_mesh_facets_vertices( M, old2new ) ;
        update_mesh_cells_vertices( M, old2new ) ;

        delete_colocated_vertices( M, old2new ) ;
    }
}