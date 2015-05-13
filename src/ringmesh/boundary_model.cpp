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
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*! \author Jeanne Pellerin and Arnaud Botella */

#include <ringmesh/boundary_model.h>
#include <ringmesh/boundary_model_builder.h>
#include <ringmesh/utils.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/points/colocate.h>
#include <geogram/mesh/triangle_intersection.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_intersection.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_io.h>


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <set>
#include <map>



namespace {
    using namespace GEO ;
    using namespace RINGMesh ;
    using GEO::index_t ;
    using GEO::vec3 ;

    /*---------------------------------------------------------------------------*/
    /* ----- Some code below is copied/modified from geogram\mesh\mesh_intersection.cpp ---*/
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
    */


    /*********************************************************************/

    /** 
    * \brief Computes the intersection between two triangular facets in
    *  a mesh
    * \param[in] M the mesh
    * \param[in] f1 index of the first facet
    * \param[in] f2 index of the second facet
    * \param[out] sym symbolic representation of the intersection (if any)
    * \return true if facets \p f1 and \p f2 have an intersection, false
    *  otherwise
    */
    bool triangles_intersect(
        const Mesh& M, index_t f1, index_t f2,
        vector<TriangleIsect>& sym
        )
    {
        geo_debug_assert( M.facets.nb_vertices( f1 ) == 3 );
        geo_debug_assert( M.facets.nb_vertices( f2 ) == 3 );
        index_t c1 = M.facets.corners_begin( f1 );
        const vec3& p1 = GEO::Geom::mesh_vertex( M, M.facet_corners.vertex( c1 ) );
        const vec3& p2 = GEO::Geom::mesh_vertex( M, M.facet_corners.vertex( c1 + 1 ) );
        const vec3& p3 = GEO::Geom::mesh_vertex( M, M.facet_corners.vertex( c1 + 2 ) );
        index_t c2 = M.facets.corners_begin( f2 );
        const vec3& q1 = GEO::Geom::mesh_vertex( M, M.facet_corners.vertex( c2 ) );
        const vec3& q2 = GEO::Geom::mesh_vertex( M, M.facet_corners.vertex( c2 + 1 ) );
        const vec3& q3 = GEO::Geom::mesh_vertex( M, M.facet_corners.vertex( c2 + 2 ) );
        return triangles_intersections( p1, p2, p3, q1, q2, q3, sym );
    }

    /*!
     * @brief Returns the Line identification if the points given define 
     *       an edge of one of this Line of the model
     */
    BME::bme_t is_edge_on_line(
        const BoundaryModel& model,
        index_t v0,
        index_t v1 )
    {
        const std::vector< BoundaryModelVertices::VertexInBME >&
            v0_bme = model.vertices.bme_vertices( v0 ) ;
        const std::vector< BoundaryModelVertices::VertexInBME >&
            v1_bme = model.vertices.bme_vertices( v1 ) ;

        // Get the local indices of the vertices in 
        // a common Line if any 
        BME::bme_t result ;
        index_t lv0 = NO_ID ;
        index_t lv1 = NO_ID ;
        for( index_t i = 0; i < v0_bme.size() ; ++i ) {
            if( v0_bme[ i ].bme_id.type == BME::LINE ) {
                for( index_t j = 0; j < v1_bme.size() ; ++j ) {
                    if(
                        v1_bme[ j ].bme_id.type == BME::LINE &&
                        v0_bme[ i ].bme_id.index == v1_bme[ j ].bme_id.index
                        ) {
                        if( lv0 == NO_ID ) {
                            lv0 = v0_bme[ i ].v_id ;
                            lv1 = v1_bme[ j ].v_id ;
                            result = v0_bme[ i ].bme_id ;
                        } else {
                            // The two points should be corners ...
                            // I am not completely sure (JP) - If they define an edge 
                            // model topology is not correct 
                            return BME::bme_t() ;
                        }
                    }
                }
            }
        }

        // There is an edge between the 2 points if their indices
        // in the Line of i and i+1
        if( abs( (int)lv0 - (int)lv1 ) == 1 ) {
            return result ;
        } else {
            return BME::bme_t() ;
        }
    }

    BME::bme_t is_edge_on_line(
        const BoundaryModel& model,
        const vec3& p0,
        const vec3& p1 )
    {
        // Get the ids in the model of these 2 points
        index_t v0 = model.vertices.vertex_index( p0 ) ;
        index_t v1 = model.vertices.vertex_index( p1 ) ;
        ringmesh_assert( v0 != NO_ID && v1 != NO_ID ) ;

        return is_edge_on_line( model, v0, v1 ) ;
    }

   


    /*! 
     * @brief Returns true if the facets of the mesh share an edge
     *       that is on one Line of the boundary model
     * @pre the mesh is triangulated
     */
    bool facets_share_line_edge(
        const Mesh& M,
        const BoundaryModel& BM,
        index_t f1,
        index_t f2 )
    {
        geo_debug_assert( M.facets.nb_vertices( f1 ) == 3 );
        geo_debug_assert( M.facets.nb_vertices( f2 ) == 3 );
        
        // I only want to test the edges that are on boundary 
        for( index_t i = 0; i < 3; ++i ) {
            if( M.facets.adjacent( f1, i ) == NO_ID ) {
                for( index_t j = 0; j < 3; ++j ) {
                    if( M.facets.adjacent( f2, j ) == NO_ID ) {
                        const vec3& p10 = M.vertices.point( M.facets.vertex( f1, i ) ) ;
                        const vec3& p11 = M.vertices.point( M.facets.vertex( f1, i==2 ? 0 : i+1 ) ) ;

                        const vec3& p20 = M.vertices.point( M.facets.vertex( f2, j ) ) ;
                        const vec3& p21 = M.vertices.point( M.facets.vertex( f2, j==2 ? 0 : j+1 ) ) ;

                        index_t v10 = BM.vertices.vertex_index( p10 ) ;
                        index_t v11 = BM.vertices.vertex_index( p11 ) ;
                        ringmesh_assert( v10 != NO_ID && v11 != NO_ID ) ;

                        index_t v20 = BM.vertices.vertex_index( p20 ) ;
                        index_t v21 = BM.vertices.vertex_index( p21 ) ;

                        if( v10 == v20 && v11 == v21 && is_edge_on_line( BM, p20, p21 ).is_defined() ) {
                            return true ;
                        }
                        if( v10 == v21 && v11 == v20 && is_edge_on_line( BM, p20, p21 ).is_defined() ) {
                            return true ;
                        }
                    }
                }
            }
        }

        return false ;
    }


    /**
    * \brief Tests whether two facets are adjacent
    * \details Two facets are adjacents if they share an edge
    *          In a Surface two facets are adjacent if they are stored as such
    *          in the Mesh, but they can also share an edge along the boundary of the
    *          Surface - checked with the global model indices
    *
    * \param[in] M the mesh
    * \param[in] f1 index of the first facet
    * \param[in] f2 index of the second facet
    * \return true if facets \p f1 and \p f2 share an edge, false
    *  otherwise
    */
    bool facets_are_adjacent(
        const Mesh& M, 
        index_t f1, index_t f2 )
    {
        if( f1 == f2 ) {
            return true;
        }                 
        for( index_t c = M.facets.corners_begin( f1 );
             c != M.facets.corners_end( f1 ); ++c ) {
            if( M.facet_corners.adjacent_facet( c ) == f2 ) {
                return true;
            }
        }        
        return false;
    }

    /**
    * \brief Action class for storing intersections when traversing
    *  a AABBTree.
    */
    class StoreIntersections {
    public:
        /**
        * \brief Constructs the StoreIntersections
        * \param[in] M the mesh
        * \param[out] has_isect the flag that indicates for each facet
        *  whether it has intersections
        */
        StoreIntersections(
            const Mesh& M, const BoundaryModel& BM, vector<index_t>& has_isect
            ) :
            M_(M), BM_(BM),
            has_intersection_( has_isect )
        {
            has_intersection_.assign( M.facets.nb(), 0 );
        }

        /**
        * \brief Determines the intersections between two facets
        * \details It is a callback for AABBTree traversal
        * \param[in] f1 index of the first facet
        * \param[in] f2 index of the second facet
        */
        void operator() ( index_t f1, index_t f2 )
        {
            if( f1 != f2 &&
                !facets_are_adjacent( M_, f1, f2 ) &&                
                !facets_share_line_edge( M_, BM_, f1, f2 ) && 
                triangles_intersect( M_, f1, f2, sym_ )
                ) {
                has_intersection_[ f1 ] = 1;
                has_intersection_[ f2 ] = 1;
            }
        }

    private:
        const Mesh& M_;
        const BoundaryModel& BM_ ;
        vector<index_t>& has_intersection_;
        vector<TriangleIsect> sym_;
    };


    /**
    * \brief Detect intersecting facets in a mesh TRIANGULATED !!
    * \param[in] M the mesh
    * \return number of intersecting facets
    */
    index_t detect_intersecting_facets( 
        const BoundaryModel& model, 
        Mesh& M )
    {
        geo_assert( M.vertices.dimension() >= 3 );

        vector<index_t> has_intersection;
        StoreIntersections action( M, model, has_intersection );
        MeshFacetsAABB AABB( M );
        AABB.compute_facet_bbox_intersections( action );

        return std::count( has_intersection.begin(), has_intersection.end(), 1 ) ;
    }

    /*********************************************************************/


    /*---------------------------------------------------------------------------*/
    /* ----- Code copied and modified from geogram\mesh\mesh_repair.cpp ---*/

    /*!
    * @brief Trigger an assertion if several vertices of a mesh at the same geometric location
    * @note Code modified from geogram/mesh/mesh_repair.cpp
    * @param[in] M the mesh
    * @param[in] colocate_epsilon tolerance
    */
    void assert_no_colocate_vertices( GEO::Mesh& M, double colocate_epsilon )
    {
        GEO::vector<index_t> old2new;

        index_t nb_new_vertices = 0;
        if( colocate_epsilon == 0.0 ) {
            nb_new_vertices = GEO::Geom::colocate_by_lexico_sort(
                M.vertices.point_ptr( 0 ), 3, M.vertices.nb(),
                old2new, M.vertices.dimension()
                );
        } else {
            nb_new_vertices = GEO::Geom::colocate(
                M.vertices.point_ptr( 0 ), 3, M.vertices.nb(),
                old2new, colocate_epsilon, M.vertices.dimension()
                );
        }
        if( nb_new_vertices != M.vertices.nb() ) {
            geo_assert_not_reached;
        }
    }

    /*!
    * @brief Merges the vertices of a mesh that are at the same geometric location
    * @note Code modified from geogram/mesh/mesh_repair.cpp
    * @param[in] M the mesh
    * @param[in] colocate_epsilon tolerance for merging vertices
    * @param[out] old2new mapping from previous M.vertices to new M.vertices
    */
    void repair_colocate_vertices(
        GEO::Mesh& M,
        double colocate_epsilon,
        GEO::vector<index_t>& old2new )
    {
        old2new.clear();

        index_t nb_new_vertices = 0;
        if( colocate_epsilon == 0.0 ) {
            nb_new_vertices = GEO::Geom::colocate_by_lexico_sort(
                M.vertices.point_ptr( 0 ), 3, M.vertices.nb(),
                old2new, M.vertices.dimension()
                );
        } else {
            nb_new_vertices = GEO::Geom::colocate(
                M.vertices.point_ptr( 0 ), 3, M.vertices.nb(),
                old2new, colocate_epsilon, M.vertices.dimension()
                );
        }
        if( nb_new_vertices == M.vertices.nb() ) {
            return;
        }
        for( index_t c = 0; c < M.facet_corners.nb(); c++ ) {
            M.facet_corners.set_vertex( c, old2new[ M.facet_corners.vertex( c ) ] );
        }

        // Some index magic to flag the point to delete and the right 
        // mapping between old and new vertices of the mesh
        GEO::vector< index_t > to_delete( old2new.size() );
        for( index_t i = 0; i < old2new.size(); i++ ) {
            if( old2new[ i ] == i ) {
                to_delete[ i ] = 0;
            } else {
                to_delete[ i ] = 1;
            }
        }
        M.vertices.delete_elements( to_delete, false );

        // The to_delete vector is used for mapping in the delete_elements function
        // We need it to get the correct mapping
        for( index_t i = 0; i < old2new.size(); i++ ) {
            if( to_delete[ i ] != NO_ID ) {
                old2new[ i ] = to_delete[ i ];
            } else {
                old2new[ i ] = to_delete[ old2new[ i ] ];
            }
        }
    }

    /**
    * \brief Connects the facets in a triangulated mesh.
    * \details Reconstructs the corners.adjacent_facet links. 
    *          Orientation not checked 
    *
    * \note Modified from geogram to take into account a predicate that impose to disconnect facets
    *       along identified edges - Jeanne
    *       The predicate implements bool operator() (index_t v1, index_t v2) const ;
    */
    template< typename P >
    void repair_connect_facets(
        Mesh& M, P is_border
        )
    {
        const index_t NO_FACET = index_t( -1 );
        const index_t NO_CORNER = index_t( -1 );
        const index_t NON_MANIFOLD = index_t( -2 );

        // Reset all facet-facet adjacencies.
        for( index_t c = 0; c<M.facet_corners.nb(); ++c ) {
            M.facet_corners.set_adjacent_facet( c, NO_FACET );
        }

        // For each vertex v, v2c[v] gives the index of a 
        // corner incident to vertex v.
        vector<index_t> v2c( M.vertices.nb(), NO_CORNER );

        // For each corner c, next_c_around_v[c] is the 
        // linked list of all the corners incident to 
        // vertex v.
        vector<index_t> next_c_around_v( M.facet_corners.nb(), NO_CORNER );

       
        // Compute v2c and next_c_around_v
        for( index_t c = 0; c<M.facet_corners.nb(); ++c ) {
            index_t v = M.facet_corners.vertex( c );
            next_c_around_v[ c ] = v2c[ v ];
            v2c[ v ] = c;
        }

        for( index_t f1 = 0; f1<M.facets.nb(); ++f1 ) {
            for(
                index_t c1 = M.facets.corners_begin( f1 );
                c1<M.facets.corners_end( f1 ); ++c1
                ) {

                if( M.facet_corners.adjacent_facet( c1 ) == NO_FACET ) {
                    index_t adj_corner = NO_CORNER;
                    index_t v1 = M.facet_corners.vertex( c1 );
                    index_t v2 = M.facet_corners.vertex(
                        M.facets.next_corner_around_facet( f1, c1 )
                        );

                    index_t c2 = v2c[ v1 ];

                    // Lookup candidate adjacent edges from incident
                    // edges list.
                    while( c2 != NO_CORNER ) {
                        if( c2 != c1 ) {
                            index_t f2 = c2/3 ;
                            index_t c3 = 
                                M.facets.prev_corner_around_facet( f2, c2 );
                            index_t v3 = M.facet_corners.vertex( c3 );
                            // Check with standard orientation.
                            if( v3 == v2 ) {
                                if( !is_border( M.vertices.point( v1 ), 
                                                M.vertices.point( v2 ) ) )
                                { // Jeanne
                                    if( adj_corner == NO_CORNER ) {
                                        adj_corner = c3;
                                    } else {
                                        // Non-manifold edge
                                        is_border.debug( M.vertices.point( v1 ),
                                                         M.vertices.point( v2 ) ) ;
                                        adj_corner = NON_MANIFOLD;
                                    }
                                }
                            } else {
                                // Check with the other ("wrong") orientation
                                c3 = M.facets.next_corner_around_facet( f2, c2 );
                                v3 = M.facet_corners.vertex( c3 );
                                if( v3 == v2 ) {
                                    if( !is_border( M.vertices.point( v1 ),
                                                    M.vertices.point( v2 ) ) 
                                       ) { // Jeanne
                                        if( adj_corner == NO_CORNER ) {
                                            adj_corner = c2;
                                        } else {
                                            // Non-manifold edge
                                            is_border.debug( M.vertices.point( v1 ), M.vertices.point( v2 ) ) ;
                                            adj_corner = NON_MANIFOLD;
                                        }
                                    }
                                }
                            }
                        }
                        c2 = next_c_around_v[ c2 ];
                    }
                    if(
                        adj_corner != NO_CORNER &&
                        adj_corner != NON_MANIFOLD
                        ) {
                        M.facet_corners.set_adjacent_facet( adj_corner, f1 );
                        index_t f2 = adj_corner/3 ;
                        M.facet_corners.set_adjacent_facet( c1, f2 );
                    }
                }
            }
        }
    }

   
    /**
     * \brief Predicate to be used by the function setting facet adjacencies in the GEO::Mesh
     *  to force disconnection of facets on a Line edge and detect unexpected non-manifold edges
     */
    class EdgeOnLine {
    public:
        EdgeOnLine( const BoundaryModel& model, Mesh& non_manifold ) :
            M_( model ), non_manifold_( non_manifold )
        {} ;
        bool operator()( const vec3& p0, const vec3& p1 ) const
        {
            return is_edge_on_line( M_, p0, p1 ).is_defined() ;
        }
        void debug( const vec3& p0, const vec3& p1 )
        {
            index_t v0 = non_manifold_.vertices.create_vertex( p0.data() ) ;
            index_t v1 = non_manifold_.vertices.create_vertex( p0.data() ) ;
            non_manifold_.edges.create_edge( v0, v1 ) ;
        }
    private:
        const BoundaryModel& M_ ;
        Mesh& non_manifold_ ;
    } ;


    /*----------------------------------------------------------------------------*/

    /*!
     * @brief Build a Mesh from the model non-duplicated vertices
     *        and its Surface facets.     
     */
    void mesh_from_boundary_model( const BoundaryModel& model, Mesh& M )
    {
        // Clear the Mesh keeping the attributes, otherwise we crash
        M.clear( true ) ;

        // Set the vertices 
        index_t nbv = model.nb_vertices() ;
        M.vertices.create_vertices( nbv ) ;
        for( index_t i = 0; i < nbv; ++i ) {
            M.vertices.point( i ) = model.vertices.unique_vertex( i ) ;
        }

        // Set the facets  
        index_t begin_S = 0 ;
        for( index_t s = 0; s < model.nb_surfaces(); ++s ) {
            begin_S = M.facets.nb() ;

            const Surface& S = model.surface( s ) ;
            for( index_t f = 0; f < S.nb_cells(); ++f ) {
                index_t nbv = S.nb_vertices_in_facet( f ) ;
                GEO::vector< index_t > ids( nbv ) ;

                for( index_t v = 0; v < nbv; ++v ) {
                    ids[ v ] = S.model_vertex_id( f, v ) ;
                }
                M.facets.create_polygon( ids ) ;
            }
            for( index_t f = 0; f < S.nb_cells(); ++f ) {
                index_t nbv = S.nb_vertices_in_facet( f ) ;
                for( index_t v = 0; v < nbv; ++v ) {
                    index_t adj = S.adjacent( f, v ) == NO_ID ? NO_ID : S.adjacent( f, v ) + begin_S ;
                    M.facets.set_adjacent( f, v, adj ) ;
                }
            }
        }
    }

    /*! 
     * @brief Get the BMME defining the boundaries of an element
     */
    void boundary_bmme( 
        const BME& E,
        std::vector< BME::bme_t >& borders, 
        bool with_inside_borders )
    {
        borders.clear() ;

        BME::TYPE T = E.bme_id().type ;
        if( T == BME::CORNER ) {
            return ;
        }
        if( BME::parent_allowed( T ) ) {
            // We are dealing with basic elements 
            for( index_t i = 0; i < E.nb_boundaries(); ++i ) {

                if( with_inside_borders ||
                    ( !with_inside_borders &&
                    !E.boundary( i ).is_inside_border( E ) )
                  ) {
                    borders.push_back( E.boundary_id( i ) ) ;
                }
            }
        } else {
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                const BME& C = E.child( i ) ;
                for( index_t j = 0; j < C.nb_boundaries(); ++j ) {
                    if( with_inside_borders ||
                        ( !with_inside_borders &&
                          !C.boundary( j ).is_inside_border( C ) )
                    ) {
                        borders.push_back( E.child( i ).boundary_id( j ) ) ;
                    }
                }
            }
            std::sort( borders.begin(), borders.end() ) ;
            index_t nb = std::unique( borders.begin(), borders.end() )-borders.begin() ;
            borders.resize( nb ) ;
        }
    }

     
    /*!
     * @brief Get the elements BME in_boundary
     * @details For BMME, get the contents of in_boudnary vector
     *          For high level elements determine in_boundary high level elements
     */
    void in_boundary_bme( const BME& E, std::vector< BME::bme_t >& in_boundary )
    {
        in_boundary.clear() ;

        BME::TYPE T = E.bme_id().type ;
        if( T == BME::REGION || T == BME::LAYER ) {
            return ;
        }
        if( BME::parent_allowed( T ) ) {
            // We are dealing with basic elements 
            for( index_t i = 0; i < E.nb_in_boundary(); ++i ) {
                in_boundary.push_back( E.in_boundary_id( i ) ) ;
            }
        } else {
            // We are dealing with high level elements
            // Need to go through the children to get information
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                for( index_t j = 0; j < E.child( i ).nb_in_boundary(); ++j ) {
                    in_boundary.push_back( E.child( i ).in_boundary( j ).parent_id() ) ;
                }
            }
            // Remove duplicates
            std::sort( in_boundary.begin(), in_boundary.end() ) ;
            index_t nb = std::unique( in_boundary.begin(), 
                                      in_boundary.end() ) - in_boundary.begin() ;
            in_boundary.resize( nb ) ;
        }

    }




    /*!
     * @brief Build a Mesh from the boundaries of the given element
     * @details Inside borders are ignored.
     */
    void mesh_from_element_boundaries( const BME& E, Mesh& M )
    {
        M.clear();

        BME::TYPE T = E.bme_id().type ;
        if( T == BME::CORNER ) {
            return ;
        }
        else {
            std::vector< BME::bme_t > borders ;
            boundary_bmme( E, borders, false ) ;
            if( borders.size() == 0 ) {
                return ;
            } else {
                if( T == BME::LINE || T == BME::CONTACT ) {
                    // There are only points to add
                    M.vertices.create_vertices( borders.size() ) ;
                    for( index_t i = 0; i < borders.size(); ++i ) {
                        M.vertices.point( i ) = 
                            E.model().element( borders[ i ] ).vertex() ;
                    }                            
                } else {
                    // Put an attribute on the ModelVertices to know its index
                    // in this Mesh
                    const BoundaryModel& model = E.model() ;
                    GEO::Attribute< index_t > old2new ; 
                    old2new.bind( model.vertices.attribute_manager(), "old2new" ) ;
                    old2new.fill( NO_ID ) ;

                    // Add the vertices 
                    for( index_t i = 0; i < borders.size(); ++i ) {
                        const BME& b = model.element( borders[ i ] ) ;                       
                        for( index_t v = 0; v < b.nb_vertices(); ++v ) {
                            index_t global_v = b.model_vertex_id( v ) ;
                            if( old2new[ global_v ] == NO_ID ) {
                                old2new[ global_v ] = M.vertices.create_vertex(
                                    model.vertices.unique_vertex( global_v ).data() ) ;
                            }
                        }
                    }

                    if( T == BME::SURFACE || T == BME::INTERFACE ) {
                        // Build edges
                        for( index_t i = 0; i < borders.size(); ++i ) {
                            ringmesh_debug_assert( borders[ i ].type == BME::LINE ) ;
                            const Line& L = model.line( borders[ i ].index ) ;
                            index_t off = M.edges.create_edges( L.mesh().edges.nb() ) ;
                            for( index_t e = 0; e < L.mesh().edges.nb(); ++e ) {
                                M.edges.set_vertex( off+e, 0, old2new[
                                    L.model_vertex_id( L.mesh().edges.vertex( e, 0 ) ) ] );
                                M.edges.set_vertex( off+e, 1, old2new[
                                    L.model_vertex_id( L.mesh().edges.vertex( e, 1 ) ) ] );
                            }
                        }

                    } else if( T == BME::REGION ) {
                        // Build facets
                        index_t off = 0 ;                        
                        for( index_t i = 0; i < borders.size(); ++i ) {
                            ringmesh_debug_assert( borders[ i ].type == BME::SURFACE ) ;
                            const Surface& S = model.surface( borders[ i ].index ) ;
                            index_t off = M.facets.nb() ;
                            for( index_t f = 0; f < S.nb_cells(); ++f ) {
                                index_t nbv = S.nb_vertices_in_facet( f ) ;
                                GEO::vector< index_t > ids( nbv ) ;
                                for( index_t v = 0; v < nbv; ++v ) {
                                    ids[ v ] = old2new[ S.model_vertex_id( f, v ) ] ;
                                }
                                M.facets.create_polygon( ids ) ;
                            }
                            for( index_t f = 0; f < S.nb_cells(); ++f ) {
                                index_t nbv = S.nb_vertices_in_facet( f ) ;
                                for( index_t v = 0; v < nbv; ++v ) {
                                    index_t adj = S.adjacent( f, v ) == NO_ID ? NO_ID : S.adjacent( f, v ) + off ;
                                    M.facets.set_adjacent( f, v, adj ) ;
                                }
                            }
                        }
                    }
                    old2new.unbind() ;
                }                
            }
        }
    }


    bool is_region_valid( const BoundaryModelElement& region )
    {
        bool valid = true ;
        if( region.bme_id().type != BME::REGION ) {
            valid = false ;
        }
        Mesh mesh ;
        mesh_from_element_boundaries( region, mesh ) ;
        GEO::mesh_repair( mesh ) ;
      
        if( GEO::mesh_nb_connected_components( mesh ) != 1 ) {
            valid = false ;
        }
        if( GEO::mesh_nb_borders( mesh ) != 0 ) {
            valid = false ;
        }       
#ifdef RINGMESH_DEBUG
        if( !valid ) {
            std::ostringstream file ;
            file << "D:\\Programming\\DataTest\\debug\\region_border_"
                << region.bme_id().index << ".mesh"  ;
            GEO::mesh_save( mesh, file.str() ) ;
        }
#endif
        return valid;
    }

    /*********************************************************************/

    /*! 
     * @brief Check if element @param is of the @param model is in the 
     *        in_boundary vector of element @param in. 
     */
    bool is_in_in_boundary(
        const BoundaryModel& model,
        BME::bme_t is, BME::bme_t in )
    {
        const BME& E = model.element( in ) ;
        for( index_t i = 0; i < E.nb_in_boundary(); ++i ) {
            if( E.in_boundary_id( i ) == is ) {
                return true ;
            }
        }
        return false ;
    }


    /*! 
     * @brief Check the geometrical-topological consistency of the model
     * @details Verification is based on the information stored by the unique
     *          vertices of the model
     */
    bool check_model_points_validity( const BoundaryModel& M )
    {
        // For all the vertices of the model 
        // We check that the elements in which they are are consistent 
        // to have a valid model
        std::vector< bool > valid( M.vertices.nb_unique_vertices(), true ) ;

        for( index_t i = 0 ; i < M.vertices.nb_unique_vertices(); ++i ) {
            
            index_t corner = NO_ID ;
            std::vector< index_t > lines ;
            std::vector< index_t > surfaces ;
            
            bool valid_vertex = true ; 

            const std::vector< BoundaryModelVertices::VertexInBME >&
                bmes = M.vertices.bme_vertices( i ) ;

            for( index_t j = 0; j < bmes.size(); ++j ) {
                BME::TYPE T = bmes[ j ].bme_id.type ;
                index_t id = bmes[ j ].bme_id.index ;

                switch( T ) {
                    case BME::SURFACE:
                        surfaces.push_back( id ) ;
                        break ;
                    case BME::LINE:
                        lines.push_back( id ) ;
                        break ;
                    case BME::CORNER :
                        if( corner != NO_ID ) {
                            valid_vertex = false ;
                        }
                        else {
                            corner = id ;
                        }
                        break ;                    
                    default :
                        valid_vertex = false ;
                }
            }

            if( valid_vertex ) {
                if( corner == NO_ID && lines.empty() ) {
                    // This is a point on one SURFACE and only one
                    if( surfaces.size() != 1 ) {
                        valid_vertex = false ;
                    }
                }

                else if( corner == NO_ID && !lines.empty() ) {
                    // This is a point on one LINE 
                    if( lines.size() != 1 ) {
                        valid_vertex = false ;
                    } else {
                        // This point must also be in at least one SURFACE
                        if( surfaces.empty() ) {
                            valid_vertex = false ;
                        }
                        // Check that one point is no more than twice in a SURFACE
                        for( index_t k = 0; k < surfaces.size(); ++k ) {
                            index_t nb = std::count( surfaces.begin(), surfaces.end(), surfaces[ k ] ) ;
                            if( nb > 2 ) {
                                valid_vertex = false ;
                                break ;
                            }
                        }
                        // Check that all the surfaces are in in_boundary of all
                        // the lines 
                        for( index_t k = 0; k < surfaces.size(); ++k ) {
                            for( index_t l = 0; l < lines.size(); ++l ) {
                                if( !is_in_in_boundary( M,
                                    BME::bme_t( BME::SURFACE, surfaces[ k ] ),
                                    BME::bme_t( BME::LINE, lines[ l ] ) )
                                    ) {
                                    valid_vertex = false ;
                                }
                            }
                        }
                    }
                }
                else if( corner != NO_ID ) {
                    // This is one point at a CORNER
                    // It must be in at least one LINE
                    if( lines.empty() ) {
                        valid_vertex = false ;
                    }
                    else {
                        if( lines.size() < 2 ) {
                            valid_vertex = false ;
                        }
                        // Check that a point is no more than twice in a LINE
                        for( index_t k = 0; k < lines.size(); ++k ) {
                            index_t nb = std::count( lines.begin(), lines.end(), lines[ k ] ) ;
                            if( nb == 2 ) {
                                // The line must be closed
                                if( !M.line( lines[ k ] ).is_closed() ) {
                                    valid_vertex = false ;
                                }
                            }
                            if( nb > 2 ) {
                                valid_vertex = false ;
                                break ;
                            }
                        }   
                        // Check that all the lines are in in_boundary of this corner
                        for( index_t k = 0; k < lines.size(); ++k ) {
                            if( !is_in_in_boundary( M,
                                BME::bme_t( BME::LINE, lines[ k ] ),
                                BME::bme_t( BME::CORNER, corner ) )
                                ) {
                                valid_vertex = false ;
                            }
                        }

                    }
                    // It must also be in a least one surface ? perhaps 2
                    if( surfaces.empty() ) {
                        valid_vertex ;                    
                    }
                }
            }
            valid[ i ] = valid_vertex ;
        }
        return std::count( valid.begin(), valid.end(), false ) == 0 ;
    }

    void save_edges( 
        const std::string& file,
        const BoundaryModel& M, 
        const std::vector< index_t >& e )
    {
        std::ofstream out( file ) ;
        if( out.is_open() ) {
            out.precision( 16 ) ;
            for( index_t i = 0 ; i < e.size(); ++i ) {
                out << "v " <<  M.vertices.unique_vertex( e[i] ) << std::endl ;
            }
            for( index_t i = 0 ; i+1 < e.size() ; i+=2 ) {
                out << "s "<< i+1 << " "<< i+2 << std::endl ;
            }
            out.close() ;
        }
    }

    /*! 
     * @brief Check boundary of a surface 
     * @details All the edges on the boundary of a surface must be in a Line
     *          of the associated model
     *          The Line boundaries must form a closed manifold line.
     */
    bool surface_boundary_valid( const Surface& S )
    {
        std::vector< index_t > invalid_corners ;
        for( index_t f = 0; f < S.nb_cells(); ++f ) {
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); ++v ) {
                if( S.adjacent( f, v ) == NO_ID &&
                    !is_edge_on_line(
                        S.model(),
                        S.model_vertex_id( f, v ),
                        S.model_vertex_id( f, S.next_in_facet( f, v ) ) ).is_defined()
                 ) {
                    invalid_corners.push_back( S.model_vertex_id( f, v ) ) ;
                    invalid_corners.push_back( S.model_vertex_id( f, S.next_in_facet( f, v ) ) ) ;
                }
            }
        }
#ifdef RINGMESH_DEBUG
        if( !invalid_corners.empty() ) {
            std::ostringstream file ;
            file << "D:\\Programming\\DataTest\\debug\\invalid_border_edge_S"
                << S.bme_id().index << ".lin"  ;
            save_edges( file.str(), S.model(), invalid_corners ) ;
        }
#endif

        
        return invalid_corners.empty() ;
    }

} // anonymous namespace 



namespace RINGMesh {

    BoundaryModelVertices::~BoundaryModelVertices()
    {
        delete ann_ ;
    }

    void BoundaryModelVertices::initialize_unique_vertices()
    {
        index_t nb_corners = bm_.nb_corners();
        index_t nb_lines = bm_.nb_lines();
        index_t nb_surfaces = bm_.nb_surfaces();

        // Total number of vertices in the Corners - Lines and Surfaces of the model
        index_t nb = bm_.nb_corners();

        for( index_t l = 0; l < nb_lines; l++ ) {
            nb += bm_.line( l ).nb_vertices();
        }
        for( index_t s = 0; s < nb_surfaces; s++ ) {
            nb += bm_.surface( s ).nb_vertices();
        }

        // Get out if the BME has no vertex
        if( nb == 0 ) {
            return;
        }

        std::vector< vec3 > all_vertices( nb );
        index_t index = 0;
        for( index_t c = 0; c < nb_corners; c++ ) {
            all_vertices[ index++ ] = bm_.corner( c ).vertex();
        }
        for( index_t l = 0; l < nb_lines; l++ ) {
            const Line& line = bm_.line( l );
            for( index_t v = 0; v < line.nb_vertices(); v++ ) {
                all_vertices[ index++ ] = line.vertex( v );
            }
        }
        for( index_t s = 0; s < nb_surfaces; s++ ) {
            const Surface& surface = bm_.surface( s );
            for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
                all_vertices[ index++ ] = surface.vertex( v );
            }
        }

        unique_vertices_.vertices.create_vertices( all_vertices.size() );
        unique_vertices_.vertices.assign_points( all_vertices[ 0 ].data(), 3, all_vertices.size() );

        GEO::vector< index_t > old2new;
        repair_colocate_vertices( unique_vertices_, epsilon, old2new );

        // We do the same loop as above
        index = 0;
        for( index_t c = 0; c < nb_corners; c++ ) {
            Corner& C = const_cast<Corner&>( bm_.corner( c ) );
            C.set_model_vertex_id( old2new[ index++ ] );
            // I am crazy paranoid (JP)
            ringmesh_debug_assert( 
                length2( C.vertex() - unique_vertex(
                old2new[ index - 1 ] ) ) < epsilon_sq );
        }
        for( index_t l = 0; l < nb_lines; l++ ) {
            Line& L = const_cast<Line&>( bm_.line( l ) );
            for( index_t v = 0; v < L.nb_vertices(); v++ ) {
                L.set_model_vertex_id( v, old2new[ index++ ] );
                ringmesh_debug_assert( 
                    length2( L.vertex( v ) - unique_vertex( 
                    old2new[ index - 1 ] ) ) < epsilon_sq );
            }
        }
        for( index_t s = 0; s < nb_surfaces; s++ ) {
            Surface& S = const_cast<Surface&>( bm_.surface( s ) );
            for( index_t v = 0; v < S.nb_vertices(); v++ ) {
                S.set_model_vertex_id( v, old2new[ index++ ] );
                ringmesh_debug_assert( 
                    length2( S.vertex( v ) - unique_vertex(
                    old2new[ index - 1 ] ) ) < epsilon_sq );
            }
        }

        ann_ = new ColocaterANN( unique_vertices_, ColocaterANN::VERTICES );

#ifdef RINGMESH_DEBUG
        // Paranoia (JP)
        assert_no_colocate_vertices(unique_vertices_, epsilon);
#endif
    }


    void BoundaryModelVertices::initialize_reverse()
    {
        if (unique_vertices_.vertices.nb() == 0) {
            initialize_unique_vertices();
        }
        if (!unique2bme_.is_bound()) {
            unique2bme_.bind(attribute_manager(), "unique2bme");
        }

        for (index_t c = 0; c < bm_.nb_corners(); c++) {
            unique2bme_[c].push_back(VertexInBME(BME::bme_t(BME::CORNER, c), 0));
        }
        for (index_t l = 0; l < bm_.nb_lines(); l++) {
            for (index_t v = 0; v < bm_.line(l).nb_vertices(); v++) {
                VertexInBME cur(BME::bme_t(BME::LINE, l), v);
                unique2bme_[unique_vertex_id(cur)].push_back(cur);
            }
        }
        for (index_t s = 0; s < bm_.nb_surfaces(); s++) {
            for (index_t v = 0; v < bm_.surface(s).nb_vertices(); v++) {
                VertexInBME cur(BME::bme_t(BME::SURFACE, s), v);
                unique2bme_[unique_vertex_id(cur)].push_back(cur);
            }
        }
    }

    void BoundaryModelVertices::update_point(index_t v, const vec3& point) const
    {
        ringmesh_assert(v < nb_unique_vertices());
        const std::vector< VertexInBME >& bme_v = bme_vertices(v);
        for (index_t i = 0; i < bme_v.size(); i++) {
            const VertexInBME& info = bme_v[i];
            const_cast<BME&>(bm_.element(
                BME::bme_t(info.bme_id))).set_vertex(
                info.v_id, point, false);
        }
    }


    const std::vector< BoundaryModelVertices::VertexInBME >&
        BoundaryModelVertices::bme_vertices(index_t v) const
    {
        ringmesh_assert(v < nb_unique_vertices());
        if (!unique2bme_.is_bound()) {
            const_cast<BoundaryModelVertices*>(this)->initialize_reverse();
        }
        return unique2bme_[v];
    }

    index_t BoundaryModelVertices::add_unique_vertex(const vec3& point)
    {
        return unique_vertices_.vertices.create_vertex(point.data());
    }

    void BoundaryModelVertices::add_unique_to_bme(
        index_t unique_id,
        BME::bme_t type,
        index_t v_id)
    {
        /// The attribute unique2bme is bound if not already ? Good idea or not ? not sure ....
        if (!unique2bme_.is_bound()) {
            unique2bme_.bind(attribute_manager(), "unique2bme");
        }
        ringmesh_assert(unique_id < nb_unique_vertices());
        unique2bme_[unique_id].push_back(VertexInBME(type, v_id));
    }

    /*!
     * @brief Returns the index of the given vertex in the model
     * \todo Implement the function - Add a KdTree for geometrical request on model vertices
     *
     * @param[in] p input point coordinates
     * @return NO_ID
     */
    index_t BoundaryModelVertices::vertex_index(const vec3& p) const
    {
        if (unique_vertices_.vertices.nb() == 0) {
            const_cast<BoundaryModelVertices*>(this)->initialize_unique_vertices();
        }
        std::vector< index_t > result;
        if (ann_->get_colocated(p, result)) return result[0];
        return NO_ID;
    }

    index_t BoundaryModelVertices::nb_unique_vertices() const
    {
        if (unique_vertices_.vertices.nb() == 0) {
            const_cast<BoundaryModelVertices*>(this)->initialize_unique_vertices();
        }
        return unique_vertices_.vertices.nb();
    }


    index_t BoundaryModelVertices::unique_vertex_id(
        BME::bme_t t,
        index_t v) const
    {
        if (unique_vertices_.vertices.nb() == 0) {
            const_cast<BoundaryModelVertices*>(this)->initialize_unique_vertices();
        }
        ringmesh_assert(v < bm_.element(t).nb_vertices());
        return bm_.element(t).model_vertex_id(v);
    }


    index_t BoundaryModelVertices::unique_vertex_id(
        const VertexInBME& v) const
    {
        return unique_vertex_id(v.bme_id, v.v_id);
    }


    const vec3& BoundaryModelVertices::unique_vertex(index_t v) const
    {
        if (unique_vertices_.vertices.nb() == 0) {
            const_cast<BoundaryModelVertices*>(this)->initialize_unique_vertices();
        }
        ringmesh_assert(v < nb_unique_vertices());
        return unique_vertices_.vertices.point(v);
    }

    void BoundaryModelVertices::clear()
    {
        /// \todo Unbind all attributes !!!! otherwise we'll get a crash
        // For the moment 
        if (unique2bme_.is_bound()) unique2bme_.unbind();

        unique_vertices_.clear(true, true);

        // Clear the information for the Corner - Line - Surface
        for (index_t c = 0; c < bm_.nb_corners(); c++) {
            Corner& C = const_cast<Corner&>(bm_.corner(c));
            C.set_model_vertex_id(NO_ID);
        }
        for (index_t l = 0; l < bm_.nb_lines(); l++) {
            Line& L = const_cast<Line&>(bm_.line(l));
            for (index_t v = 0; v < L.nb_vertices(); v++) {
                L.set_model_vertex_id(v, NO_ID);
            }
        }
        for (index_t s = 0; s < bm_.nb_surfaces(); s++) {
            Surface& S = const_cast<Surface&>(bm_.surface(s));
            for (index_t v = 0; v < S.nb_vertices(); v++) {
                S.set_model_vertex_id(v, NO_ID);
            }
        }
    }

    /*******************************************************************************/


    BoundaryModel::~BoundaryModel()
    {
        for (index_t i = 0; i < corners_.size(); i++) {
            if (corners_[i]) delete corners_[i];
        }
        for (index_t i = 0; i < lines_.size(); i++) {
            if (lines_[i]) delete lines_[i];
        }
        for (index_t i = 0; i < surfaces_.size(); i++) {
            if (surfaces_[i]) delete surfaces_[i];
        }
        for (index_t i = 0; i < regions_.size(); i++) {
            if (regions_[i]) delete regions_[i];
        }
        for (index_t i = 0; i < contacts_.size(); i++) {
            if (contacts_[i]) delete contacts_[i];
        }
        for (index_t i = 0; i < interfaces_.size(); i++) {
            if (interfaces_[i]) delete interfaces_[i];
        }
        for (index_t i = 0; i < layers_.size(); i++) {
            if (layers_[i]) delete layers_[i];
        }
    }

    /*!
     * @brief Total number of facets in the model Surface s
     */
    index_t BoundaryModel::nb_facets() const
    {
        index_t result = 0;
        for (index_t i = 0; i < nb_surfaces(); ++i) {
            result += surface(i).nb_cells();
        }
        return result;
    }

    /*!
     * Copies a BoundaryModel in another one
     * @param[in] from BoundaryModel to copy
     */
    void BoundaryModel::copy(const BoundaryModel& from)
    {
        BoundaryModelBuilder builder(*this);
        builder.copy_macro_topology(from);
        builder.copy_meshes(from);
    }

    /*!
     * @brief Returns the index of the region neighboring the surface.
     * @param[in] surface_part_id Index of the Surface
     * @param[in] side Side of the Surface
     * @return The region index or NO_ID if none found.
     */
    index_t BoundaryModel::find_region(
        index_t surface_part_id,
        bool side) const
    {
        ringmesh_debug_assert(surface_part_id < nb_surfaces());
        BME::bme_t cur_surface(BME::SURFACE,
            surface_part_id);
        for (index_t r = 0; r < nb_regions(); r++) {
            const BME& cur_region = region(r);
            for (index_t s = 0; s < cur_region.nb_boundaries(); s++) {
                if (cur_region.side(s) == side
                    && cur_region.boundary_id(s) == cur_surface)
                {
                    return r;
                }
            }
        }
        return BME::NO_ID;
    }


    /*!
     * @brief Modify the model so that it is compatible with a Gocad Model3d
     *   and can be saved in .ml format
     *
     * @return True if this was a success, False if modifications could not be done.
     */
    bool BoundaryModel::check_model3d_compatibility()
    {
        BoundaryModelBuilder builder(*this);

        /// 1. Check that the Interfaces exist
        if (nb_interfaces() == 0 && nb_surfaces() > 0) {
            /// If not create one Interface per Surface
            for (index_t i = 0; i < surfaces_.size(); ++i) {
                // Set name, type, links with other elements
                std::ostringstream name;
                name << "surface_" << i;
                BME::bme_t id = builder.create_interface(name.str());
                builder.add_child(id, BME::bme_t(BME::SURFACE, i));
            }

            // Set links from surfaces_ toward interfaces_
            for (index_t i = 0; i < interfaces_.size(); ++i) {
                builder.set_parent(
                    one_interface(i).child_id(0),
                    BME::bme_t(BME::INTERFACE, i));
            }

            // Is it really useful to have contacts, let's hope not... I am not doing it
        }

        /// 2. Check that the Universe region exists
        /// \todo Write some code to create the universe (cf. line 805 to 834 de s2_b_model.cpp)
        if (universe_.name() != "Universe") {
            GEO::Logger::err("")
                <<
                "The region universe is not defined for the model. IMPLEMENTATION TO DO"
                << std::endl;
            return false;
        }

        /// 3. Check that each region has a name and valid surfaces
        for (index_t i = 0; i < regions_.size(); ++i) {
            const BME& region = this->region(i);

            if (region.name() == "") {
                std::ostringstream name;
                name << "region_" << i;
                builder.set_element_name(
                    BME::bme_t(BME::REGION, i),
                    name.str());
            }
            if (region.nb_boundaries() == 0) {
                GEO::Logger::err("") << "The region " << region.name()
                    << " has no Surfaces on its boundary" <<
                    std::endl;
                return false;
            }
        }

        /// 4. Check that all the surfaces_ of the model are triangulated
        /// \todo Implement a triangulation function in SurfaceMutator
        for (index_t s = 0; s < nb_surfaces(); s++) {
            if (!surface(s).is_triangulated()) {
                GEO::Logger::err("") << "Surface " << s <<
                    " is not triangulated" << std::endl;
                return false;
            }
        }
        return true;
    }

    /*!
     * @brief Check the validity of all individual elements 
     * @details Check that the elements belong to this model, 
     *          call the check validity for each element
     *          For regions, check that their boundary is a one connected component
     *          manifold closed surface.
     *
     */
    bool BoundaryModel::check_elements_validity() const
    {
        std::vector< bool > valid( nb_elements( BME::ALL_TYPES ), true ) ;
        for( index_t i = 0; i < nb_elements( BME::ALL_TYPES ); ++i ) {
            const BME& E = element( BME::bme_t( BME::ALL_TYPES, i ) ) ;
            // Verify that E points actually to this BoundaryModel
            if( &E.model() != this ) {
                valid[i] = false ;
                ringmesh_debug_assert( false ) ; 
                break ;
            }
            valid[ i ] = E.is_valid() ;
            
            if( valid[i] && E.bme_id().type == BME::REGION ) {
                // Check validity of region definition
                valid[ i ] = is_region_valid( E ) ;
            }
        }        
        index_t nb_invalid = std::count( valid.begin(), valid.end(), false) ;
        std::cout << nb_invalid << " elements are invalid " << std::endl ;
        return nb_invalid == 0 ;
    }

    /*! 
     * @brief Check geological validity -
     * @details Only a fault can have a free border and 
     *          an interface can on the boundary of maximum two layers      
     *          See Building and Editing a Sealed Geological Model,
     *          Caumon et al. 2004
     */          
    bool BoundaryModel::check_geology_validity() const
    {
        for( index_t i = 0; i < nb_lines(); ++i ) {
            if( line( i ).nb_in_boundary() == 1 ) {
                const BME& S = line( i ).in_boundary( 0 ) ;
                if( S.has_parent() && S.parent().has_geological_feature() &&
                    S.parent().geological_feature() != BME::FAULT ) {
                    return false ;
                }
            }
        }

        for( index_t i = 0; i < nb_interfaces(); ++i ) {
            std::vector< BME::bme_t > regions ;
            in_boundary_bme( one_interface( i ), regions ) ;
            if( regions.size() == 0 ) {
                return false ;
            }
            if( !one_interface(i).geological_feature() == BME::STRATI &&
                regions.size() > 2 
              ) {
                return false ;
            }            
        }
    }

    /*!
     * @brief Check consistency of the geometry of the model elements 
     *        with the stored connectivity information
     * @details Finite extension - Universe region exists - has no hole -
     *          its boundary is be a closed manifold surface - one connected component.
     *          
     *          No intersection between two different elements except along
     *          shared boundaries - that must be actual boundaries  
     *          Performed on Corners, Lines and Surfaces.
     */
    bool BoundaryModel::check_model_validity() const
    {
        if( !check_elements_validity() ) {
            return false ;
        }
        if( nb_interfaces() > 0 && !check_geology_validity() ) {
            return false ;
        }

        // Check that the model has a finite extension 
        // The boundary of the universe region is a one connected component 
        // manifold closed surface 
        if( !is_region_valid( universe() ) ) {
            return false ;
        }        
        /// @todo check that facet orientation is consistent ? useful really ?

        // Check geometrical-connectivity consistency
        if( !check_model_points_validity( *this ) ) {
            return false ;
        }

        // No edge of a Surface can be on the boundary of this Surface without
        // being in a Line
        for( index_t i = 0; i < nb_surfaces(); ++i ) {
            if( !surface_boundary_valid( surface( i ) ) ){
                return false ;
            }
        }
      

        // Check surface-surface intersections 
        // Build a global triangulated mesh corresponding
        // to this model
        GEO::Mesh model_mesh ;
        mesh_from_boundary_model( *this, model_mesh ) ;
        GEO::mesh_repair( model_mesh, MESH_REPAIR_TRIANGULATE ) ;

        // Check that non-manifold edges are all Lines
        GEO::Mesh non_manifold_edges ;
        EdgeOnLine P( *this, non_manifold_edges ) ;
        repair_connect_facets( model_mesh, P ) ;

        if( non_manifold_edges.vertices.nb() > 0 ) {
            return false ;
        }
        
        // Check there is no intersections except along Line boundaries
        if( detect_intersecting_facets( *this, model_mesh ) > 0 ) {
            return false ;
        }

        return true ;
    }




    /*!
     * @brief Write a region information in a stream
     * @details Used by function to save the Model in a .ml file
     *
     * @param[in] count Region index in the file
     * @param[in] region The region to save
     * @param[in,out] out The file output stream
     */
    void save_region(
        index_t count,
        const BoundaryModelElement& region,
        std::ostream& out )
    {
        out << "REGION " << count << "  " << region.name() << " " << std::endl ;
        index_t it = 0 ;

        for( index_t i = 0; i < region.nb_boundaries(); ++i ) {
            out << "  " ;
            if( region.side( i ) ) {
                out << "+" ;
            } else {
                out << "-" ;
            }
            out << region.boundary( i ).bme_id().index + 1 ;
            it++ ;
            if( it == 5 ) {
                out << std::endl ;
                it = 0 ;
            }
        }
        out << "  0" << std::endl ;
    }


    /*!
     * @brief Write information for on layer in a stream
     * @details Used by function to save the Model in a .ml file
     *
     * @param[in] count Index of the layer in the file
     * @param[in] offset Offset of region indices in the file
     * @param[in] layer The layer to write
     * @param[in,out] out The output file stream
     */
    void save_layer(
        index_t count,
        index_t offset,
        const BoundaryModelElement& layer,
        std::ostream& out )
    {
        out << "LAYER " << layer.name() << " " << std::endl ;
        index_t it = 0 ;

        for( index_t i = 0; i < layer.nb_children(); ++i ) {
            out << "  " << layer.child_id( i ).index + offset + 1 ;
            it++ ;
            if( it == 5 ) {
                out << std::endl ;
                it = 0 ;
            }
        }
        out << "  0" << std::endl ;
    }


    /*!
     * @brief Write basic header for Gocad coordinate system.
     * @param[in,out] out Output .ml file stream
     */
    void save_coordinate_system( std::ostream& out )
    {
        out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl << "NAME Default"
            << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
            << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl << "ZPOSITIVE Elevation"
            << std::endl << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl ;
    }


    /*!
     * @brief Save the model in a .ml file if it compatible
     *
     * @param[in,out] out Output file stream
     * @return false if the model is not compatible with a Gocad model
     */
    bool BoundaryModel::save_gocad_model3d( std::ostream& out )
    {
        out.precision( 16 ) ;
        if( !check_model3d_compatibility() ) {
            GEO::Logger::err( "" ) << "The BoundaryModel " << name_
                                   << " cannot be saved in .ml format " << std::endl ;
            return false ;
        }

        // Print Gocad Model3d headers
        out << "GOCAD Model3d 1" << std::endl << "HEADER {" << std::endl << "name:"
            << name() << std::endl << "}" << std::endl ;

        save_coordinate_system( out ) ;

        // Print the TSurf = Interface information
        for( index_t i = 0; i < nb_interfaces(); ++i ) {
            out << "TSURF " << one_interface( i ).name() << std::endl ;
        }

        index_t count = 1 ;

        // Print the TFace = Surface information
        for( index_t i = 0; i < nb_surfaces(); ++i ) {
            const Surface& s = surface( i ) ;
            out << "TFACE " << count << "  " ;
            out << BME::geol_name( s.geological_feature() ) ;
            out << " " << s.parent().name() << std::endl ;

            // Print the key facet points, whuich are simply the first three
            // vertices of the first facet
            out << "  " << s.vertex( 0, 0 ) << std::endl ;
            out << "  " << s.vertex( 0, 1 ) << std::endl ;
            out << "  " << s.vertex( 0, 2 ) << std::endl ;

            ++count ;
        }

        index_t offset_layer = count ;

        // Print universe, region, and layer information
        save_region( count, universe_, out ) ;
        ++count ;

        for( index_t i = 0; i < nb_regions(); ++i ) {
            save_region( count, region( i ), out ) ;
            ++count ;
        }

        for( index_t i = 0; i < nb_layers(); ++i ) {
            save_layer( count, offset_layer, layer( i ), out ) ;
            ++count ;
        }
        out << "END" << std::endl ;

        // Save the geometry of the Surfaces (TFace), Interface (TSurf) by Interface
        for( index_t i = 0; i < nb_interfaces(); ++i ) {
            const BME& tsurf = one_interface( i ) ;

            // Header
            out << "GOCAD TSurf 1" << std::endl << "HEADER {" << std::endl <<
            "name:"
                << tsurf.name() << std::endl << "name_in_model_list:" << tsurf.name()
                << std::endl << "}" << std::endl ;
            save_coordinate_system( out ) ;

            out << "GEOLOGICAL_FEATURE " << tsurf.name() << std::endl
                << "GEOLOGICAL_TYPE " ;
            out << BME::geol_name( tsurf.geological_feature() ) ;
            out << std::endl ;

            out << "PROPERTY_CLASS_HEADER Z {" << std::endl << "is_z:on" <<
            std::endl
                << "}" << std::endl ;

            // Save surfaces_ geometry
            index_t vertex_count = 1 ;
            index_t offset = vertex_count ;

            std::vector< index_t > bstones ;
            std::vector< index_t > next_vertex ;
            std::set< index_t > set_end_corners ;

            bstones.reserve( tsurf.nb_boundaries() ) ;
            next_vertex.reserve( tsurf.nb_boundaries() ) ;

            for( index_t j = 0; j < tsurf.nb_children(); ++j ) {
                offset = vertex_count ;

                const Surface& sp = dynamic_cast< const Surface& >( tsurf.child( j ) ) ;

                out << "TFACE" << std::endl ;
                for( index_t k = 0; k < sp.nb_vertices(); ++k ) {
                    out << "VRTX " << vertex_count << " " << sp.vertex( k ) <<
                    std::endl ;
                    vertex_count++ ;
                }

                for( index_t k = 0; k < sp.nb_cells(); ++k ) {
                    if( sp.nb_vertices_in_facet( k ) != 3 ) {
                        GEO::Logger::err( "I/O" ) << "Model is not triangulated"
                            << std::endl ;
                        return false ;
                    }
                    out << "TRGL " << sp.surf_vertex_id( k, 0 ) + offset << " "
                        << sp.surf_vertex_id( k, 1 ) + offset << " "
                        << sp.surf_vertex_id( k, 2 ) + offset << std::endl ;
                }

                // Gather information on Corners (BStones) and Lines (getting the next point on the line)
                for( index_t k = 0; k < sp.nb_boundaries(); ++k ) {
                    const Line& cp = dynamic_cast< const Line& >( sp.boundary( k ) ) ;

                    vec3 c = cp.vertex( 0 ) ;
                    vec3 next = cp.vertex( 1 ) ;

                    // To be sure that we have all corners we need to ensure
                    // that all corners at the end of lines are saved too
                    std::vector< index_t > result ;
                    sp.tools.ann().get_colocated( cp.vertex( cp.nb_vertices() - 1 ), result ) ;
                    ringmesh_debug_assert( !result.empty() ) ;
                    set_end_corners.insert( result[0] + offset ) ;

                    result.clear() ;
                    sp.tools.ann().get_colocated( c, result ) ;
                    ringmesh_debug_assert( !result.empty() ) ;
                    index_t c_id = result[0] ;
                    result.clear() ;
                    sp.tools.ann().get_colocated( next, result ) ;
                    ringmesh_debug_assert( !result.empty() ) ;
                    index_t next_id = result[0] ;

                    ringmesh_assert( c_id != NO_ID && next_id != NO_ID ) ;

                    bstones.push_back( c_id + offset ) ;
                    next_vertex.push_back( next_id + offset ) ;
                }
            }

            // Print Corners and Lines
            std::vector< index_t > end_corners(
                set_end_corners.begin(), set_end_corners.end() ) ;
            std::vector< bool > end_corner_to_print( end_corners.size(), true ) ;

            for( index_t j = 0; j < bstones.size(); ++j ) {
                out << "BSTONE " << bstones[ j ] << std::endl ;

                // Determine the corners at the end of the lines that are not saved
                for( index_t k = 0; k < end_corners.size(); k++ ) {
                    if( bstones[ j ] == end_corners[ k ] ) {
                        end_corner_to_print[ k ] = false ;
                        break ;
                    }
                }
            }

            // Print the corners that were at the beginning of none of the contacts
            // in this Interface
            for( index_t j = 0; j < end_corners.size(); j++ ) {
                if( end_corner_to_print[ j ] ) {
                    out << "BSTONE " << end_corners[ j ] << std::endl ;
                }
            }

            // Print the the information to build the lines :
            // index of the vertex at the corner and index of the second vertex on the line
            for( index_t j = 0; j < bstones.size(); ++j ) {
                out << "BORDER " << vertex_count << " " << bstones[ j ] << " "
                    << next_vertex[ j ] << std::endl ;
                vertex_count++ ;
            }
            out << "END" << std::endl ;
        }
        return true ;
    }


    /*! To save the attributes in a Graphite readable file, we need to write the correct
     * keyword for the attribute type - We restrict ourselves to the 3 types
     * int          "integer"
     * double       "real"
     * float        "real"
     * bool         "boolean"
     */
    inline std::string alias_name( const std::string& in )
    {
        if( in == "int" ) {return "integer" ;} else if( in == "index" ) {
            return "integer" ;
        } else if( in == "double" ) {
            return "real" ;
        } else if( in ==
                   "float" )
        {
            return "real" ;
        } else if( in ==
                   "bool" )
        {
            return "boolean" ;
        }
        ringmesh_assert_not_reached ;
        return "" ;
    }


    /*!
     * @brief DEBUG function - Save the surfaces of the model with their facet attributes into an .eobj file.
     * @details WARNING We assume that all Surface have the same attributes - if not this function will most
     *  certainly crash.
     *
     * @param[in] file_name Name of the file
     *
     * \todo Make this function const
     *
     */
    void BoundaryModel::save_as_eobj_file( const std::string& file_name ) const
    {
        std::ofstream out ;
        out.open( file_name.c_str() ) ;
        if( out.bad() ) {
            std::cout << "Error when opening the file: " << file_name.c_str() <<
            std::endl ;
            return ;
        }
        out.precision( 16 ) ;
        std::vector< index_t > offset( nb_surfaces(), 0 ) ;
        index_t cur_offset = 0 ;

        // Write vertices once for each surface
        for( index_t s = 0; s < nb_surfaces(); s++ ) {
            const Surface& S = surface( s ) ;
            offset[ s ] = cur_offset ;
            for( index_t p = 0; p < S.nb_vertices(); p++ ) {
                const vec3& V = S.vertex( p ) ;
                out << "v"
                    << " " << V.x
                    << " " << V.y
                    << " " << V.z
                    << std::endl ;
            }
            cur_offset += S.nb_vertices() ;
        }

        // Write the facets for a each surface
        for( index_t s = 0; s < nb_surfaces(); s++ ) {
            const Surface& S = surface( s ) ;
            for( index_t f = 0; f < S.nb_cells(); f++ ) {
                out << "f" << " " ;
                for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                    out << offset[ s ] + S.surf_vertex_id( f, v ) + 1 << " " ;
                }
                out << std::endl ;
            }
        }

        // Write facet attributes
        {
            // Get the attributes that can be saved on the first Surface
//            std::vector< SerializedAttribute< BME::FACET > > facet_attribs ;
//            get_serializable_attributes( surface(
//                    0 ).facet_attribute_manager(), facet_attribs ) ;
//
//            for( index_t i = 0; i < facet_attribs.size(); i++ ) {
//                // Output global information on the attribute
//                out << "# attribute " << facet_attribs[ i ].name() << " facet "
//                    << alias_name( facet_attribs[ i ].type_name() )
//                    << std::endl ;
//            }
//            if( facet_attribs.size() > 0 ) {
//                // Global counter for all the facets of all surfaces
//                index_t count = 0 ;
//                for( index_t s = 0; s < nb_surfaces(); s++ ) {
//                    const Surface& S = surface( s ) ;
//
//                    std::vector< SerializedAttribute< BME::FACET > > cur_attribs ;
//                    get_serializable_attributes(
//                        S.facet_attribute_manager(), cur_attribs ) ;
//
//                    ringmesh_assert( cur_attribs.size() == facet_attribs.size() ) ;
//                    for( index_t i = 0; i < facet_attribs.size(); ++i ) {
//                        ringmesh_assert(
//                            facet_attribs[ i ].type_name() ==
//                            cur_attribs[ i ].type_name() &&
//                            facet_attribs[ i ].name() == cur_attribs[ i ].name() ) ;
//                    }
//
//                    // Output attributes values
//                    for( index_t f = 0; f < S.nb_cells(); f++ ) {
//                        out << "# attrs f " << count + 1 << " " ;
//                        serialize_write_attributes( out, f, cur_attribs ) ;
//                        out << std::endl ;
//                        count++ ;
//                    }
//                }
//            }
        }
        // Write the attribute <index_t> on facets called ""chart" - because I want to have it
        {
            out << "# attribute " << "chart" << " facet "
                << "integer"
                << std::endl ;
            
            // Global counter for all the facets of all surfaces
            index_t count = 0 ;
            for( index_t s = 0; s < nb_surfaces(); s++ ) {
                const Surface& S = surface( s ) ;

                GEO::Attribute< index_t > A( S.cell_attribute_manager(), "chart") ;

                for( index_t f = 0; f < S.nb_cells(); f++ ) {
                    out << "# attrs f " << count + 1 << " " << A[f] 
                        << std::endl ;
                    count++ ;
                }
            }
        }
    }


    /*!
     * @brief Debug: Save a Surface of the model in the file OBJ format is used
     */
    void BoundaryModel::save_surface_as_obj_file(
        index_t s,
        const std::string& file_name ) const
    {
        std::ofstream out ;
        out.open( file_name.c_str() ) ;
        if( out.bad() ) {
            std::cout << "Error when opening the file: " << file_name.c_str() <<
            std::endl ;
            return ;
        }
        out.precision( 16 ) ;
        const Surface& S = surface( s ) ;
        for( index_t p = 0; p < S.nb_vertices(); p++ ) {
            const vec3& V = S.vertex( p ) ;
            out << "v"
                << " " << V.x
                << " " << V.y
                << " " << V.z
                << std::endl ;
        }
        for( index_t f = 0; f < S.nb_cells(); f++ ) {
            out << "f" << " " ;
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                out << S.surf_vertex_id( f, v ) + 1 << " " ;
            }
            out << std::endl ;
        }
    }


    /*!
     * @brief Write in the out stream things to save for CONTACT, INTERFACE and LAYERS
     */
    void save_high_level_bme(
        std::ofstream& out,
        const BoundaryModelElement& E )
    {
        /// First line:  TYPE - ID - NAME - GEOL
        out << BoundaryModelElement::type_name( E.bme_id().type ) << " "
            << E.bme_id().index << " " ;
        if( E.has_name() ) { out << E.name() << " " ;} else { out << "no_name " ;}
        out <<  BoundaryModelElement::geol_name( E.geological_feature() )
            << std::endl ;

        /// Second line:  IDS of children
        for( index_t j = 0; j < E.nb_children(); ++j ) {
            out << " " << E.child_id( j ).index ;
        }
        out << std::endl ;
    }


    /*!
     * @brief Save the BoundaryModel into a dedicated format bm
     */
    void BoundaryModel::save_bm_file( const std::string& file_name ) const
    {
        std::ofstream out ;
        out.open( file_name.c_str() ) ;
        if( out.bad() ) {
            std::cout << "Error when opening the file: " << file_name.c_str() <<
            std::endl ;
            return ;
        }
        out.precision( 16 ) ;

        out << "RINGMESH BOUNDARY MODEL" << std::endl ;
        out << "NAME " << name() << std::endl ;

        // Number of the different type of elements
        for( index_t i = BME::CORNER; i < BME::NO_TYPE; i++ ) {
            BME::TYPE type = (BME::TYPE) i ;
            out <<  "NB_" << BME::type_name( type ) << " " << nb_elements( type ) <<
            std::endl ;
        }

        // Write high-level elements
        for( index_t i = BME::CONTACT; i < BME::NO_TYPE; i++ ) {
            BME::TYPE type = (BME::TYPE) i ;
            index_t nb = nb_elements( type ) ;
            for( index_t j = 0; j < nb; ++j ) {
                save_high_level_bme( out, element( BME::bme_t( type, j ) ) ) ;
            }
        }

        // Regions
        for( index_t i = 0; i < nb_regions(); ++i ) {
            const BME& E = region( i ) ;

            // Save ID - NAME -
            out << BME::type_name( BME::REGION ) << " " << E.bme_id().index << " " ;
            if( E.has_name() ) {out << E.name() ;} else {out << "no_name" ;}
            out << std::endl ;

            // Second line Signed ids of boundary surfaces
            for( index_t j = 0; j < E.nb_boundaries(); ++j ) {
                if( E.side( j ) ) {out << "+" ;} else {out << "-" ;}
                out << E.boundary_id( j ).index << " " ;
            }
            out << std::endl ;
        }

        // Universe
        out << "UNIVERSE " << std::endl ;
        for( index_t j = 0; j < universe().nb_boundaries(); ++j ) {
            if( universe().side( j ) ) {out << "+" ;} else {out << "-" ;}
            out << universe().boundary_id( j ).index << " " ;
        }
        out << std::endl ;

//        // Vertices and attributes on vertices
//        out << "MODEL_VERTICES" << " " << nb_vertices() << std::endl ;
//        out << "MODEL_VERTEX_ATTRIBUTES " ;
//        std::vector< SerializedAttribute< VERTEX > > vertex_attribs ;
//        get_serializable_attributes( &vertex_attribute_manager_, vertex_attribs, out ) ;
//        for( index_t i = 0; i < nb_vertices(); ++i ) {
//            out << vertex( i )  << "  " ;
//            serialize_write_attributes( out, i, vertex_attribs ) ;
//            out << std::endl ;
//        }

        // Corners
        for( index_t i = 0; i < nb_corners(); ++i ) {
            out << BME::type_name( BME::CORNER ) << " "
                << corner( i ).bme_id().index << " " << corner( i ).vertex() <<
            std::endl ;
        }

        // Lines
        for( index_t i = 0; i < nb_lines(); ++i ) {
            const Line& L = line( i ) ;
            out << BME::type_name( BME::LINE ) << " " << L.bme_id().index << std::endl ;
            out << "LINE_VERTICES " << L.nb_vertices() << std::endl ;
            for( index_t j = 0; j < L.nb_vertices(); ++j ) {
                out << L.vertex( j ) << std::endl ;
            }
//            out << "LINE_VERTEX_ATTRIBUTES " ;
//            std::vector< SerializedAttribute< BME::VERTEX > > line_vertex_attribs ;
//            get_serializable_attributes(
//                L.vertex_attribute_manager(), line_vertex_attribs, out ) ;
//            for( index_t j = 0; j < L.nb_vertices(); ++j ) {
//                out << j << "  " ;
//                serialize_write_attributes( out, j, line_vertex_attribs ) ;
//                out << std::endl ;
//            }
//            out << "LINE_SEGMENT_ATTRIBUTES " ;
//            std::vector< SerializedAttribute< BME::FACET > > line_segments_attribs ;
//            get_serializable_attributes(
//                L.facet_attribute_manager(), line_segments_attribs, out ) ;
//            if( line_segments_attribs.size() > 0 ) {
//                for( index_t j = 0; j < L.nb_cells(); ++j ) {
//                    out << j << "  " ;
//                    serialize_write_attributes( out, j, line_segments_attribs ) ;
//                    out << std::endl ;
//                }
//            }
            out << "IN_BOUNDARY " ;
            for( index_t j = 0; j < L.nb_in_boundary(); ++j ) {
                out << L.in_boundary_id( j ).index << " " ;
            }
            out << std::endl ;
        }

        // Surfaces
        for( index_t i = 0; i < nb_surfaces(); ++i ) {
            const Surface& S = surface( i ) ;
            out << BME::type_name( BME::SURFACE ) << " " << S.bme_id().index << std::endl ;
            out << "SURFACE_VERTICES " << S.nb_vertices() << std::endl ;
            for( index_t j = 0; j < S.nb_vertices(); ++j ) {
                out << S.vertex( j ) << std::endl ;
            }
//            out << "SURFACE_VERTEX_ATTRIBUTES " ;
//            std::vector< SerializedAttribute< BME::VERTEX > > surface_vertex_attribs ;
//            get_serializable_attributes(
//                S.vertex_attribute_manager(), surface_vertex_attribs, out ) ;
//            for( index_t j = 0; j < S.nb_vertices(); ++j ) {
//                out << j << "  " ;
//                serialize_write_attributes( out, j, surface_vertex_attribs ) ;
//                out << std::endl ;
//            }

            out << "SURFACE_CORNERS " << S.nb_corners() << std::endl ;
            out << "SURFACE_FACETS " << S.nb_cells() << std::endl ;
//            out << "SURFACE_FACET_ATTRIBUTES " ;
//            std::vector< SerializedAttribute< BME::FACET > > surface_facet_attribs ;
//            get_serializable_attributes(
//                S.facet_attribute_manager(), surface_facet_attribs, out ) ;
//
            for( index_t j = 0; j < S.nb_cells(); ++j ) {
                out << S.nb_vertices_in_facet( j ) << " " ;
                for( index_t v = 0; v < S.nb_vertices_in_facet( j ); ++v ) {
                    out << S.surf_vertex_id( j, v ) << " " ;
                }
//                serialize_write_attributes( out, j, surface_facet_attribs ) ;
                out << std::endl ;
            }
        }
    }


    /*!
     * \brief Save the model in smesh format
     * \details No attributes and no boundary marker are transferred
     */
    void BoundaryModel::save_smesh_file( const std::string& file_name ) const 
    {
        std::ofstream out;
        out.open(file_name.c_str());
        if (out.bad()) {
            std::cout << "Error when opening the file: " << file_name.c_str() <<
                std::endl;
            return;
        }
        out.precision(16);

        /// 1. Write the unique vertices
        out << "# Node list" << std::endl;
        out << "# node count, 3 dim, no attribute, no boundary marker" << std::endl;
        out << vertices.nb_unique_vertices() << " 3 0 0" << std::endl;
        out << "# node index, node coordinates " << std::endl;
        for (index_t p = 0; p < vertices.nb_unique_vertices(); p++){
            const vec3& V = vertices.unique_vertex(p);
            out << p << " "
                << " " << V.x
                << " " << V.y
                << " " << V.z
                << std::endl;
        }

        /// 2. Write the triangles 
        out << "# Part 2 - facet list" << std::endl;
        out << "# facet count, no boundary marker" << std::endl;
        out << nb_facets() << "  0 " << std::endl;

        for (index_t i = 0; i < nb_surfaces(); ++i) {
            const Surface& S = surface(i);
            for (index_t f = 0; f < S.nb_cells(); f++){
                out << S.nb_vertices_in_facet(f) << " ";
                for (index_t v = 0; v < S.nb_vertices_in_facet(f); v++){
                    out << S.model_vertex_id( f, v ) << " ";
                }
                out << std::endl;
            }
        }

        // Do not forget the stupid zeros at the end of the file 
        out << std::endl << "0" << std::endl << "0" << std::endl;
    }

    signed_index_t BoundaryModel::find_interface( const std::string& name) const {
        for(index_t i = 0 ; i < nb_interfaces() ; i++ ) {
            if( one_interface(i).name() == name ) {
                return i ;
            }
        }
        GEO::Logger::err("") << "Surface name did not match with an actual\
                                interface name of the Boundary Model. Abort.. " 
                                << std::endl ;
        ringmesh_assert_not_reached ;
        return -1 ;
    }

    signed_index_t BoundaryModel::find_region( const std::string& name) const {
        for(index_t r = 0 ; r < nb_regions() ; r++ ) {
            if( region(r).name() == name ) {
                return r ;
            }
        }
        GEO::Logger::err("") << "Region name did not match with an actual\
                                region name of the Boundary Model. Abort.. " 
                                << std::endl ;
        ringmesh_assert_not_reached ;
        return -1 ;
    }


} // namespace
