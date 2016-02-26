/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/geo_model_validity.h>

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>

#include <geogram/basic/algorithm.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/string.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_intersection.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/triangle_intersection.h>

#include <geogram/points/colocate.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/geometry.h>
#include <ringmesh/geogram_mesh_repair.h>

/*!
 * @file ringmesh/geo_model_validity.cpp
 * @brief Implementation of functions to check the validity of GeoModels
 * @author Jeanne Pellerin
 */

namespace {

    using namespace GEO ;
    using namespace RINGMesh ;
    using GEO::index_t ;
    using GEO::vec3 ;

    typedef GeoModelMeshElement GMME ;

    /*---------------------------------------------------------------------------*/
    /*----- Some pieces of the code below are copied or modified from -----------*/
    /*----- geogram\mesh\mesh_intersection.cpp-----------------------------------*/
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

    /** \note Copied from geogram
     * \brief Computes the intersection between two triangular facets in
     *  a mesh
     * \param[in] M the mesh
     * \param[in] f1 index of the first facet
     * \param[in] f2 index of the second facet
     * \param[warn] sym symbolic representation of the intersection (if any)
     * \return true if facets \p f1 and \p f2 have an intersection, false
     *  otherwise
     */
    bool triangles_intersect(
        const Mesh& M,
        index_t f1,
        index_t f2,
        vector< TriangleIsect >& sym )
    {
        ringmesh_debug_assert( M.facets.nb_vertices( f1 ) == 3 ) ;
        ringmesh_debug_assert( M.facets.nb_vertices( f2 ) == 3 ) ;
        index_t c1 = M.facets.corners_begin( f1 ) ;
        const vec3& p1 = GEO::Geom::mesh_vertex( M, M.facet_corners.vertex( c1 ) ) ;
        const vec3& p2 = GEO::Geom::mesh_vertex( M,
            M.facet_corners.vertex( c1 + 1 ) ) ;
        const vec3& p3 = GEO::Geom::mesh_vertex( M,
            M.facet_corners.vertex( c1 + 2 ) ) ;

        index_t c2 = M.facets.corners_begin( f2 ) ;
        const vec3& q1 = GEO::Geom::mesh_vertex( M, M.facet_corners.vertex( c2 ) ) ;
        const vec3& q2 = GEO::Geom::mesh_vertex( M,
            M.facet_corners.vertex( c2 + 1 ) ) ;
        const vec3& q3 = GEO::Geom::mesh_vertex( M,
            M.facet_corners.vertex( c2 + 2 ) ) ;
        return triangles_intersections( p1, p2, p3, q1, q2, q3, sym ) ;
    }

    /*!
     * @brief Returns the Line identification if the given points define
     *       an edge of one of the Line of the model
     * @param model The GeoModel to consider
     * @param v0 Index of the first point in the model
     * @param v1 Index of the second point in the model
     */
    GME::gme_t is_edge_on_line( const GeoModel& model, index_t v0, index_t v1 )
    {
        const std::vector< GMEVertex >& v0_bme = model.mesh.vertices.gme_vertices(
            v0 ) ;
        const std::vector< GMEVertex >& v1_bme = model.mesh.vertices.gme_vertices(
            v1 ) ;

        // Get the local indices of the vertices in 
        // a common Line if any 
        GME::gme_t result ;
        index_t lv0 = NO_ID ;
        index_t lv1 = NO_ID ;

        // No sorting to optimize since 
        // v0_bme and v1_bme are very small sets ( < 10 elements ) [JP]
        for( index_t i = 0; i < v0_bme.size(); ++i ) {
            if( v0_bme[i].gme_id.type == GME::LINE ) {
                for( index_t j = 0; j < v1_bme.size(); ++j ) {
                    if( v1_bme[j].gme_id.type == GME::LINE
                        && v0_bme[i].gme_id.index == v1_bme[j].gme_id.index ) {
                        if( lv0 == NO_ID ) {
                            lv0 = v0_bme[i].v_id ;
                            lv1 = v1_bme[j].v_id ;
                            result = v0_bme[i].gme_id ;
                        } else {
                            if( !model.line( result.index ).is_closed() ) {
                                // Most certainly there is a problem (JP)
                                return GME::gme_t() ;
                            }

                        }
                    }
                }
            }
        }
        if( !result.is_defined() ) {
            // The two points are not on the same Line
            return GME::gme_t() ;
        } else {
            // Determine if the points define an edge
            if( lv0 > lv1 ) {
                std::swap( lv0, lv1 ) ;
            }
            // Casts are here to avoid a compiler warning [JP]
            int delta_i = static_cast< int >( lv1 ) - static_cast< int >( lv0 ) ;

            if( delta_i == 1 ) {
                // There is an edge if their indices in the Line are i and i+1
                return result ;
            } else if( model.line( result.index ).is_closed()
                && delta_i == model.line( result.index ).nb_vertices() - 2 ) {
                // If the Line is closed we can also have 0; n-2 or n-1; 1
                return result ;
            } else {
                // The two points are on the same line but
                // do not define an edge
                return GME::gme_t() ;
            }
        }
    }

    /*!
     * @brief Returns the Line identification if the given points define
     *       an edge of one of the Line of the model
     */
    GME::gme_t is_edge_on_line(
        const GeoModel& model,
        const vec3& p0,
        const vec3& p1 )
    {
        // Get the ids in the model of these 2 points
        index_t v0 = model.mesh.vertices.index( p0 ) ;
        index_t v1 = model.mesh.vertices.index( p1 ) ;
        ringmesh_debug_assert( v0 != NO_ID && v1 != NO_ID ) ;

        return is_edge_on_line( model, v0, v1 ) ;
    }

    /*!
     * @brief Returns true if the facets @param f1 and @param f2
     *        of the mesh @param M share an edge
     *        that is on one Line of the boundary model @param BM
     * @pre The mesh M is triangulated
     *
     */
    bool facets_share_line_edge(
        const Mesh& M,
        const GeoModel& BM,
        index_t f1,
        index_t f2 )
    {
        ringmesh_debug_assert( M.facets.nb_vertices( f1 ) == 3 ) ;
        ringmesh_debug_assert( M.facets.nb_vertices( f2 ) == 3 ) ;

        // I only want to test the edges that are on boundary 
        for( index_t i = 0; i < 3; ++i ) {
            if( M.facets.adjacent( f1, i ) == NO_ID ) {
                for( index_t j = 0; j < 3; ++j ) {
                    if( M.facets.adjacent( f2, j ) == NO_ID ) {
                        const vec3& p10 = M.vertices.point(
                            M.facets.vertex( f1, i ) ) ;
                        const vec3& p11 = M.vertices.point(
                            M.facets.vertex( f1, i == 2 ? 0 : i + 1 ) ) ;

                        const vec3& p20 = M.vertices.point(
                            M.facets.vertex( f2, j ) ) ;
                        const vec3& p21 = M.vertices.point(
                            M.facets.vertex( f2, j == 2 ? 0 : j + 1 ) ) ;

                        index_t v10 = BM.mesh.vertices.index( p10 ) ;
                        index_t v11 = BM.mesh.vertices.index( p11 ) ;
                        ringmesh_debug_assert( v10 != NO_ID && v11 != NO_ID ) ;

                        index_t v20 = BM.mesh.vertices.index( p20 ) ;
                        index_t v21 = BM.mesh.vertices.index( p21 ) ;

                        if( v10 == v20 && v11 == v21
                            && is_edge_on_line( BM, p20, p21 ).is_defined() ) {
                            return true ;
                        }
                        if( v10 == v21 && v11 == v20
                            && is_edge_on_line( BM, p20, p21 ).is_defined() ) {
                            return true ;
                        }
                    }
                }
            }
        }

        return false ;
    }

    /** \note Copied from geogram
     * \brief Tests whether two facets are adjacent
     * \details Two facets are adjacents if they share an edge
     *
     * \param[in] M the mesh
     * \param[in] f1 index of the first facet
     * \param[in] f2 index of the second facet
     * \return true if facets \p f1 and \p f2 share an edge, false
     *  otherwise
     */
    bool facets_are_adjacent( const Mesh& M, index_t f1, index_t f2 )
    {
        if( f1 == f2 ) {
            return true ;
        }
        for( index_t c = M.facets.corners_begin( f1 );
            c != M.facets.corners_end( f1 ); ++c ) {
            if( M.facet_corners.adjacent_facet( c ) == f2 ) {
                return true ;
            }
        }
        return false ;
    }

    /** \note Modified from geogram
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
            const Mesh& M,
            const GeoModel& BM,
            vector< index_t >& has_isect )
            : M_( M ), BM_( BM ), has_intersection_( has_isect )
        {
            has_intersection_.assign( M.facets.nb(), 0 ) ;
        }

        /**
         * \brief Determines the intersections between two facets
         * \details It is a callback for AABBTree traversal
         * \param[in] f1 index of the first facet
         * \param[in] f2 index of the second facet
         */
        void operator()( index_t f1, index_t f2 )
        {
            if( f1 != f2 && !facets_are_adjacent( M_, f1, f2 )
                && !facets_share_line_edge( M_, BM_, f1, f2 )
                && triangles_intersect( M_, f1, f2, sym_ ) ) {
                has_intersection_[f1] = 1 ;
                has_intersection_[f2] = 1 ;
            }
        }

    private:
        const Mesh& M_ ;
        const GeoModel& BM_ ;
        vector< index_t >& has_intersection_ ;
        vector< TriangleIsect > sym_ ;
    } ;

    /** \note Copied from geogram
     * \brief Detect intersecting facets in a TRIANGULATED mesh
     * \param[in] M the mesh
     * \return number of intersecting facets
     */
    index_t detect_intersecting_facets( const GeoModel& model, Mesh& M )
    {
        geo_assert( M.vertices.dimension() >= 3 ) ;

        vector< index_t > has_intersection ;
        StoreIntersections action( M, model, has_intersection ) ;
        MeshFacetsAABB AABB( M ) ;
        AABB.compute_facet_bbox_intersections( action ) ;

        index_t nb_intersections = static_cast< index_t >( std::count(
            has_intersection.begin(), has_intersection.end(), 1 ) ) ;

        if( nb_intersections > 0 ) {
            GEO::Mesh mesh ;
            for( index_t f = 0; f < has_intersection.size(); f++ ) {
                if( !has_intersection[f] ) continue ;
                GEO::vector< index_t > vertices ;
                vertices.reserve( 3 ) ;
                for( index_t v = 0; v < M.facets.nb_vertices( f ); v++ ) {
                    index_t id = mesh.vertices.create_vertex(
                        M.vertices.point_ptr( M.facets.vertex( f, v ) ) ) ;
                    vertices.push_back( id ) ;
                }
                mesh.facets.create_polygon( vertices ) ;
            }
            std::ostringstream file ;
            file << validity_errors_directory << "/intersected_facets.mesh" ;
            GEO::mesh_save( mesh, file.str() ) ;
            GEO::Logger::out( "I/O" ) << std::endl ;
        }
        return nb_intersections ;
    }

    /***************************************************************************/

    /*---------------------------------------------------------------------------*/
    /*----- Some pieces of the code below are copied or modified from -----------*/
    /*----- geogram\mesh\mesh_repair.cpp-----------------------------------------*/

    /*!
     * @brief Trigger an assertion if several vertices of a mesh at the same geometric location

     */
    void assert_no_colocate_vertices( const GEO::Mesh& M, double colocate_epsilon )
    {
        if( has_mesh_colocate_vertices( M, colocate_epsilon ) ) {
            geo_assert_not_reached;
        }
    }

    /*!
     * @brief Get the colocated vertices of a mesh, i.e. which have the same geometric location
     * @note Code modified from geogram/mesh/mesh_repair.cpp
     * @param[in] M the mesh
     * @param[in] colocate_epsilon tolerance for merging vertices
     * @param[out] old2new if old2new[i] == i, point is to keep; otherwise
     *             old2new[i] = j, j is the index of the matching point kept
     * @returns true if there are colocated vertices
     *
     * @todo replace by a call to  GEO::mesh_detect_colocated_vertices since it
     *       is now in the API
     *
     * @pre The mesh has no facet, cell or edges.
     */
    bool colocate_vertices(
        GEO::Mesh& M,
        double colocate_epsilon,
        GEO::vector< index_t >& old2new )
    {
        old2new.clear() ;

        if( M.edges.nb() > 0 || M.facets.nb() > 0 || M.cells.nb() > 0 ) {
            // This function is not sufficient to update the complete mesh.
            ringmesh_debug_assert( false ) ;
        }

        index_t nb_new_vertices = 0 ;
        if( colocate_epsilon == 0.0 ) {
            nb_new_vertices = GEO::Geom::colocate_by_lexico_sort(
                M.vertices.point_ptr( 0 ), 3, M.vertices.nb(), old2new,
                M.vertices.dimension() ) ;
        } else {
            nb_new_vertices = GEO::Geom::colocate( M.vertices.point_ptr( 0 ), 3,
                M.vertices.nb(), old2new, colocate_epsilon,
                M.vertices.dimension() ) ;
        }
        return nb_new_vertices != M.vertices.nb() ;
    }

 
    /*----------------------------------------------------------------------------*/

    /*!
     * @brief Get the BMME defining the boundaries of an element
     */
    void boundary_bmme(
        const GME& E,
        std::vector< GME::gme_t >& borders,
        bool with_inside_borders )
    {
        borders.clear() ;

        GME::TYPE T = E.type() ;
        if( T == GME::CORNER ) {
            return ;
        }
        if( GME::parent_allowed( T ) ) {
            // We are dealing with basic elements 
            for( index_t i = 0; i < E.nb_boundaries(); ++i ) {
                if( with_inside_borders
                    || ( !with_inside_borders
                        && !E.boundary( i ).is_inside_border( E ) ) ) {
                    borders.push_back( E.boundary_gme( i ) ) ;
                }
            }
        } else {
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                const GME& C = E.child( i ) ;
                for( index_t j = 0; j < C.nb_boundaries(); ++j ) {
                    if( with_inside_borders
                        || ( !with_inside_borders
                            && !C.boundary( j ).is_inside_border( C ) ) ) {
                        borders.push_back( E.child( i ).boundary_gme( j ) ) ;
                    }
                }
            }
            GEO::sort_unique( borders ) ;
        }
    }

    /*!
     * @brief Get the elements in the boundary of which @param E is
     * @details For BMME, get the contents of the in_boundary vector
     *          For high level elements, determine in_boundary high level elements
     */
    void in_boundary_bme( const GME& E, std::vector< GME::gme_t >& in_boundary )
    {
        in_boundary.clear() ;

        GME::TYPE T = E.type() ;
        if( T == GME::REGION || T == GME::LAYER ) {
            return ;
        }
        if( GME::parent_allowed( T ) ) {
            // We are dealing with basic elements 
            for( index_t i = 0; i < E.nb_in_boundary(); ++i ) {
                in_boundary.push_back( E.in_boundary_gme( i ) ) ;
            }
        } else {
            // We are dealing with high level elements
            // Need to go through the children to get information
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                for( index_t j = 0; j < E.child( i ).nb_in_boundary(); ++j ) {
                    in_boundary.push_back(
                        E.child( i ).in_boundary( j ).parent_id() ) ;
                }
            }
            // Remove duplicates
            GEO::sort_unique( in_boundary ) ;
        }
    }

    /*!
     * @brief Build a Mesh from the boundaries of the given element
     * @details Inside borders are ignored. Adjacencies are not set.
     * Client should call mesh repair functions afterwards.
     */
    void mesh_from_element_boundaries( const GME& E, Mesh& M )
    {
        M.clear() ;

        GME::TYPE T = E.type() ;
        if( T == GME::CORNER ) {
            return ;
        } else {
            std::vector< GME::gme_t > borders ;
            boundary_bmme( E, borders, false ) ;
            if( borders.size() == 0 ) {
                return ;
            } else {
                if( T == GME::LINE || T == GME::CONTACT ) {
                    // There are only points to add
                    M.vertices.create_vertices( borders.size() ) ;
                    for( index_t i = 0; i < borders.size(); ++i ) {
                        M.vertices.point( i ) =
                            E.model().mesh_element( borders[i] ).vertex() ;
                    }
                } else {
                    // Put an attribute on the ModelVertices to know its index
                    // in this Mesh
                    const GeoModel& model = E.model() ;
                    GEO::Attribute< index_t > old2new ;
                    old2new.bind( model.mesh.vertex_attribute_manager(),
                        "old2new" ) ;
                    old2new.fill( NO_ID ) ;

                    // Add the vertices 
                    for( index_t i = 0; i < borders.size(); ++i ) {
                        const GMME& b = model.mesh_element( borders[i] ) ;
                        for( index_t v = 0; v < b.nb_vertices(); ++v ) {
                            index_t global_v = b.model_vertex_id( v ) ;
                            if( old2new[global_v] == NO_ID ) {
                                old2new[global_v] = M.vertices.create_vertex(
                                    model.mesh.vertices.vertex( global_v ).data() ) ;
                            }
                        }
                    }

                    if( T == GME::SURFACE || T == GME::INTERFACE ) {
                        // Build edges
                        for( index_t i = 0; i < borders.size(); ++i ) {
                            ringmesh_debug_assert( borders[i].type == GME::LINE ) ;
                            const Line& L = model.line( borders[i].index ) ;
                            index_t off = M.edges.create_edges(
                                L.mesh().edges.nb() ) ;
                            for( index_t e = 0; e < L.mesh().edges.nb(); ++e ) {
                                M.edges.set_vertex( off + e, 0,
                                    old2new[L.model_vertex_id(
                                        L.mesh().edges.vertex( e, 0 ) )] ) ;
                                M.edges.set_vertex( off + e, 1,
                                    old2new[L.model_vertex_id(
                                        L.mesh().edges.vertex( e, 1 ) )] ) ;
                            }
                        }

                    } else if( T == GME::REGION ) {
                        // Build facets              
                        for( index_t i = 0; i < borders.size(); ++i ) {
                            ringmesh_debug_assert( borders[i].type == GME::SURFACE ) ;
                            const Surface& S = model.surface( borders[i].index ) ;
                            for( index_t f = 0; f < S.nb_cells(); ++f ) {
                                index_t nbv = S.nb_vertices_in_facet( f ) ;
                                GEO::vector< index_t > ids( nbv ) ;
                                for( index_t v = 0; v < nbv; ++v ) {
                                    ids[v] = old2new[S.model_vertex_id( f, v )] ;
                                }
                                M.facets.create_polygon( ids ) ;
                            }
                        }
                    }
                    old2new.unbind() ;
                }
            }
        }
    }

    /*
     * @brief Checks that boundary surfaces of @param region define
     *        a one connected component closed manifold surface
     * @details Builds a GEO::Mesh from the surface meshes, repairs it and analyses it.
     * @todo Put this function in Region class
     */
    bool is_region_valid( const GeoModelElement& region )
    {
        if( region.type() != GME::REGION ) {
            GEO::Logger::err( "GeoModel" ) << " Incorrect element type "
                << GME::type_name( region.type() ) << " for " << region.gme_id()
                << std::endl ;
            return false ;
        }
        if( region.nb_boundaries() == 0 ) {
            GEO::Logger::warn( "GeoModel" ) << region.gme_id()
                << " has no boundary Surface" << std::endl ;
            return false ;
        } else {
            Mesh mesh ;
            GEO::Logger::instance()->set_quiet( true ) ;
            mesh_from_element_boundaries( region, mesh ) ;
            GEO::mesh_repair( mesh ) ;
            GEO::Logger::instance()->set_quiet( false ) ;

            bool valid = true ;
            index_t nb_cc = GEO::mesh_nb_connected_components( mesh ) ;
            index_t nb_b = GEO::mesh_nb_borders( mesh ) ;
            if( nb_cc != 1 ) {
                GEO::Logger::warn( "GeoModel" ) << " Surface boundary of "
                    << region.gme_id() << " has " << nb_cc
                    << " connected components " << std::endl ;
                valid = false ;
            }
            if( nb_b != 0 ) {
                GEO::Logger::warn( "GeoModel" ) << " Surface boundary of "
                    << region.gme_id() << " has " << nb_b
                    << " border connected components " << std::endl ;
                valid = false ;
            }
            if( !valid ) {
                std::ostringstream file ;
                file << validity_errors_directory << "/boundary_surface_region_"
                    << region.index() << ".mesh" ;
                GEO::mesh_save( mesh, file.str() ) ;
                return false ;
            } else {
                return true ;
            }
        }
    }

    /*!
     * @brief Check if element @param is of the @param model is in the
     *        in_boundary vector of element @param in.
     */
    bool is_in_in_boundary( const GeoModel& model, GME::gme_t is, GME::gme_t in )
    {
        const GME& E = model.element( in ) ;
        for( index_t i = 0; i < E.nb_in_boundary(); ++i ) {
            if( E.in_boundary_gme( i ) == is ) {
                return true ;
            }
        }
        return false ;
    }

    void save_invalid_points(
        const std::string& file,
        const GeoModel& M,
        const std::vector< bool >& valid )
    {
        std::ofstream out ;
        out.open( file.c_str() ) ;
        if( out.bad() ) {
            GEO::Logger::err( "File" ) << "Failed to open file: " << file
                << std::endl ;
        } else {
            out.precision( 16 ) ;
            for( index_t i = 0; i < valid.size(); ++i ) {
                if( !valid[i] ) {
                    const vec3& V = M.mesh.vertices.vertex( i ) ;
                    out << "v" << " " << V.x << " " << V.y << " " << V.z
                    // Not in format file, but convenient
                        << " model index " << i << std::endl ;
                }
            }
        }
    }

    /*!
     * @brief Check the geometrical-topological consistency of the model
     * @details Verification is based on the information stored by the unique
     *          vertices of the model which validity must be checked beforehand
     * @todo Check that the model vertices are consistent with the model_vertex_ids
     *       stored at by the GMME
     * @todo Implementation for regions
     * @todo Split in smaller functions
     */
    bool check_model_points_validity( const GeoModel& M )
    {
        // For all the vertices of the model 
        // We check that the elements in which they are are consistent 
        // to have a valid B-Rep model
        std::vector< bool > valid( M.mesh.vertices.nb(), true ) ;
        for( index_t i = 0; i < M.mesh.vertices.nb(); ++i ) {
            bool valid_vertex = true ;

            // Get the mesh elements in which this vertex is            
            index_t corner = NO_ID ;
            std::vector< index_t > lines ;
            std::vector< index_t > surfaces ;
            std::vector< index_t > regions ;

            const std::vector< GMEVertex >& bmes = M.mesh.vertices.gme_vertices(
                i ) ;

            for( index_t j = 0; j < bmes.size(); ++j ) {
                GME::TYPE T = bmes[j].gme_id.type ;
                index_t id = bmes[j].gme_id.index ;
                switch( T ) {
                    case GME::REGION:
                        regions.push_back( id ) ;
                        break ;
                    case GME::SURFACE:
                        surfaces.push_back( id ) ;
                        break ;
                    case GME::LINE:
                        lines.push_back( id ) ;
                        break ;
                    case GME::CORNER:
                        if( corner != NO_ID ) {
                            GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                                << " is in at least 2 Corners" << std::endl ;
                            valid_vertex = false ;
                        } else {
                            corner = id ;
                        }
                        break ;
                    default:
                        GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                            << " is in no Element of the Model" << std::endl ;
                        valid_vertex = false ;
                        break ;
                }
            }

            if( valid_vertex ) {
                if( surfaces.empty() ) {
                    if( regions.size() != 1 ) {
                        GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                            << " is in " << regions.size() << " Regions: " ;
                        for( index_t j = 0; j < surfaces.size(); ++j ) {
                            GEO::Logger::warn( "GeoModel" ) << regions[j] << " ; " ;
                        }
                        GEO::Logger::warn( "GeoModel" ) << std::endl ;
                        valid_vertex = false ;
                    } /// @todo Implement the other conditions for Region point validity
                } else if( corner == NO_ID && lines.empty() ) {
                    // This is a point on one SURFACE and only one
                    if( surfaces.size() != 1 ) {
                        GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                            << " is in " << surfaces.size() << " Surfaces: " ;
                        for( index_t j = 0; j < surfaces.size(); ++j ) {
                            GEO::Logger::warn( "GeoModel" ) << surfaces[j] << " ; " ;
                        }
                        GEO::Logger::warn( "GeoModel" ) << std::endl ;
                        valid_vertex = false ;
                    }
                } else if( corner == NO_ID && !lines.empty() ) {
                    // This is a point on one LINE 
                    if( lines.size() != 1 ) {
                        GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                            << " is in " << lines.size() << " Lines " ;
                        for( index_t j = 0; j < lines.size(); ++j ) {
                            GEO::Logger::warn( "GeoModel" ) << lines[j] << " ; " ;
                        }
                        GEO::Logger::warn( "GeoModel" ) << std::endl ;
                        valid_vertex = false ;
                    } else {
                        // This point must also be in at least one SURFACE
                        if( surfaces.empty() ) {
                            GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                                << " is in a Line but in no Surface " << std::endl ;
                            valid_vertex = false ;
                        }
                        // Check that one point is no more than twice in a SURFACE
                        for( index_t k = 0; k < surfaces.size(); ++k ) {
                            index_t nb = static_cast< index_t >( std::count(
                                surfaces.begin(), surfaces.end(), surfaces[k] ) ) ;
                            if( nb > 2 ) {
                                GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                                    << " is " << nb << " times in Surface "
                                    << M.surface( surfaces[k] ).gme_id()
                                    << std::endl ;
                                valid_vertex = false ;
                            } else if( nb == 2 ) {
                                // If a point is twice in a SURFACE, it must be
                                // on an internal boundary Line.
                                bool internal_boundary = false ;
                                for( index_t l = 0; l < lines.size(); ++l ) {
                                    if( M.line( lines[l] ).is_inside_border(
                                        M.surface( surfaces[k] ) ) ) {
                                        internal_boundary = true ;
                                        break ;
                                    }
                                }
                                if( !internal_boundary ) {
                                    GEO::Logger::warn( "GeoModel" ) << " Vertex "
                                        << i << " appears " << nb
                                        << " times in Surface "
                                        << M.surface( surfaces[k] ).gme_id()
                                        << std::endl ;
                                    valid_vertex = false ;
                                }
                            }
                        }
                        // Check that all the surfaces are in in_boundary of all
                        // the lines 
                        for( index_t k = 0; k < surfaces.size(); ++k ) {
                            for( index_t l = 0; l < lines.size(); ++l ) {
                                GME::gme_t s_id( GME::SURFACE, surfaces[k] ) ;
                                GME::gme_t l_id( GME::LINE, lines[l] ) ;
                                if( !is_in_in_boundary( M, s_id, l_id ) ) {
                                    GEO::Logger::warn( "GeoModel" )
                                        << " Inconsistent Line-Surface connectivity "
                                        << " Vertex " << i << " shows that "
                                        << M.element( s_id ).gme_id()
                                        << " must be in the boundary of "
                                        << M.element( l_id ).gme_id() << std::endl ;
                                    valid_vertex = false ;
                                }
                            }
                        }
                    }
                } else if( corner != NO_ID ) {
                    // This is one point at a CORNER
                    // It must be in at least one LINE
                    if( lines.empty() ) {
                        GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                            << " is at a Corner but in no Line " << std::endl ;
                        valid_vertex = false ;
                    } else {
                        if( lines.size() < 2 ) {
                            GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                                << " is in at a Corner but in one Line only: "
                                << lines[0] << std::endl ;
                            valid_vertex = false ;
                        }
                        // Check that a point is no more than twice in a LINE
                        for( index_t k = 0; k < lines.size(); ++k ) {
                            index_t nb = static_cast< index_t >( std::count(
                                lines.begin(), lines.end(), lines[k] ) ) ;
                            if( nb == 2 ) {
                                // The line must be closed
                                if( !M.line( lines[k] ).is_closed() ) {
                                    GEO::Logger::warn( "GeoModel" ) << " Vertex "
                                        << i << " is twice in Line " << lines[k]
                                        << std::endl ;
                                    valid_vertex = false ;
                                }
                            }
                            if( nb > 2 ) {
                                GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                                    << " appears " << nb << " times in Line "
                                    << lines[k] << std::endl ;
                                valid_vertex = false ;
                                break ;
                            }
                        }
                        // Check that all the lines are in in_boundary of this corner
                        for( index_t k = 0; k < lines.size(); ++k ) {
                            GME::gme_t l_id( GME::LINE, lines[k] ) ;
                            GME::gme_t c_id( GME::CORNER, corner ) ;
                            if( !is_in_in_boundary( M, l_id, c_id ) ) {
                                GEO::Logger::warn( "GeoModel" )
                                    << " Inconsistent Line-Corner connectivity "
                                    << " vertex " << i << " shows that "
                                    << M.element( l_id ).gme_id()
                                    << " must be in the boundary of "
                                    << M.element( c_id ).gme_id() << std::endl ;
                                valid_vertex = false ;
                            }
                        }
                    }
                    // It must also be in a least one surface ? perhaps 2
                    if( surfaces.empty() ) {
                        GEO::Logger::warn( "GeoModel" ) << " Vertex " << i
                            << " is at a Corner but in no Surface " << std::endl ;
                        valid_vertex = false ;
                    }
                }
            }
            valid[i] = valid_vertex ;
        }
        index_t nb_invalid = static_cast< index_t >( std::count( valid.begin(),
            valid.end(), false ) ) ;

        if( nb_invalid > 0 ) {
            std::ostringstream file ;
            file << validity_errors_directory << "/invalid_global_vertices.pts" ;
            save_invalid_points( file.str(), M, valid ) ;

            GEO::Logger::warn( "GeoModel" ) << nb_invalid << " invalid vertices "
                << std::endl << "Saved in file: " << file.str() << std::endl ;
            return false ;
        } else {
            return true ;
        }
    }

    void save_edges(
        const std::string& file,
        const GeoModel& M,
        const std::vector< index_t >& e )
    {
        std::ofstream out( file.c_str() ) ;
        if( out.bad() ) {
            GEO::Logger::err( "File" ) << "Failed to open file: " << file
                << std::endl ;
        } else {
            out.precision( 16 ) ;
            GEO::Logger::out( "I/O" ) << "Saving " << file <<"...." << std::endl ;
            
            for( index_t i = 0; i < e.size(); ++i ) {
                out << "v " << M.mesh.vertices.vertex( e[i] ) << std::endl ;
            }
            for( index_t i = 0; i + 1 < e.size(); i += 2 ) {
                out << "s " << i + 1 << " " << i + 2 << std::endl ;
            }
            out.close() ;

            GEO::Logger::out( "I/O" ) << std::endl ;
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
                if( S.adjacent( f, v ) == NO_ID
                    && !is_edge_on_line( S.model(), S.model_vertex_id( f, v ),
                        S.model_vertex_id( f, S.next_in_facet( f, v ) ) ).is_defined() ) {
                    invalid_corners.push_back( S.model_vertex_id( f, v ) ) ;
                    invalid_corners.push_back(
                        S.model_vertex_id( f, S.next_in_facet( f, v ) ) ) ;
                }
            }
        }
        if( !invalid_corners.empty() ) {
            std::ostringstream file ;
            file << validity_errors_directory << "/invalid_boundary_surface_"
                << S.index() << ".lin" ;
            save_edges( file.str(), S.model(), invalid_corners ) ;

            GEO::Logger::warn( "GeoModel" ) << " Invalid surface boundary: "
                << invalid_corners.size() / 2 << " boundary edges of " << S.gme_id()
                << "  are in no line of the model " << std::endl
                << " Saved in file: " << file.str() << std::endl ;
            return false ;
        } else {
            return true ;
        }
    }
} // anonymous namespace 

namespace RINGMesh {
    /*!
    * @brief Save in a .lin file the
    */
    void debug_save_non_manifold_edges( const GeoModel& geomodel,
                                        const std::vector< index_t>& edge_vertices )
    {
        std::string extension( "lin" ) ;
        std::string what( "non_manifold_edges" ) ;
        std::string file_name( validity_errors_directory +"/" + what + "." + extension ) ;

        save_edges( file_name, geomodel, edge_vertices ) ;
    }

    bool are_geomodel_elements_valid( const GeoModel& GM )
    {
        std::vector< bool > valid( GM.nb_elements( GME::ALL_TYPES ), true ) ;
        for( index_t e = 0; e < GM.nb_elements( GME::ALL_TYPES ); ++e ) {
            const GME& E = GM.element( GME::gme_t( GME::ALL_TYPES, e ) ) ;
            // Verify that E points actually to this GeoModel
            if( &E.model() != &GM ) {
                GEO::Logger::err( "GeoModel" ) << "The model stored for "
                    << GME::type_name( E.type() ) << " " << E.index()
                    << " is not correct " << std::endl ;
                valid[ e ] = false ;
                // This is a major problem
                ringmesh_debug_assert( false ) ;
                break ;
            }
            valid[ e ] = E.is_valid() ;

            if( valid[ e ] && E.type() == GME::REGION ) {
                // Check validity of region definition
                valid[ e ] = is_region_valid( E ) ;
            }
        }
        index_t nb_invalid = static_cast< index_t >(std::count( valid.begin(),
            valid.end(), false )) ;
        if( nb_invalid != 0 ) {
            GEO::Logger::warn( "GeoModel" ) << nb_invalid
                << " individual elements of the model are invalid " << std::endl ;
        }
        return nb_invalid == 0 ;
    }

    bool is_geomodel_geology_valid( const GeoModel& GM )
    {
        bool valid = true ;
        for( index_t l = 0; l < GM.nb_lines(); ++l ) {
            if( GM.line( l ).nb_in_boundary() == 1 ) {
                const GME& S = GM.line( l ).in_boundary( 0 ) ;
                if( S.has_parent()
                    && !GME::is_fault( S.parent().geological_feature() ) ) {
                    GEO::Logger::warn( "GeoModel" ) << " Invalid free border: "
                        << GM.line( l ).gme_id() << " is in the boundary of Surface "
                        << S.gme_id() << " that is not a FAULT " << std::endl
                        << std::endl ;
                    valid = false ;
                }
            }
        }

        for( index_t i = 0; i < GM.nb_interfaces(); ++i ) {
            std::vector< GME::gme_t > layers ;
            in_boundary_bme( GM.one_interface( i ), layers ) ;
            if( layers.size() == 0 ) {
                GEO::Logger::warn( "GeoModel" ) << " Invalid interface: "
                    << GM.one_interface( i ).gme_id()
                    << " is in the boundary of no Layer " << std::endl ;
                valid = false ;
            }
            if( GM.one_interface( i ).geological_feature() == GME::STRATI
                && layers.size() > 2 ) {
                GEO::Logger::warn( "GeoModel" ) << " Invalid horizon: "
                    << GM.one_interface( i ).gme_id() << " is in the boundary of "
                    << layers.size() << " Layers: " ;
                for( index_t j = 0; j < layers.size(); ++j ) {
                    GEO::Logger::warn( "GeoModel" ) << layers[ j ] << " ; " ;
                }
                GEO::Logger::warn( "GeoModel" ) << std::endl ;
                valid = false ;
            }
        }
        return valid ;
    }


    /*!
     * @brief Implementation class for validity check on a GeoModel
     */
    class GeoModelValidityCheck {
    public:
        GeoModelValidityCheck( const GeoModel& geomodel ) :
            geomodel_(geomodel), valid_(true)
        {
            // Ensure that the model vertices are computed and up-to-date
            // Without them we cannot do anything        
            geomodel_.mesh.vertices.test_and_initialize() ;
            do_check_validity() ;
        }

        bool is_geomodel_valid()
        {
            return valid_ ;
        }

    private:        
        const GeoModel& geomodel() const
        {
            return geomodel_ ;
        }

        void test_global_element_access()
        {
            index_t nb_global = geomodel().nb_elements( GME::ALL_TYPES ) ;
            index_t sum_all_types = geomodel().nb_corners() + geomodel().nb_lines() 
                + geomodel().nb_surfaces() + geomodel().nb_regions() 
                + geomodel().nb_contacts() + geomodel().nb_interfaces()
                + geomodel().nb_layers() ;

            if( nb_global != sum_all_types ) {
                set_invalid_model() ;
            }
        }
        
        // Verify the validity of all GeoModelElements
        void test_model_elements_validity()
        {
            if( !are_geomodel_elements_valid( geomodel_ ) ) {
                set_invalid_model() ;
            }
        }

        // Verify the geological validity if the model has
        // interfaces and layers
        void test_geological_validity()
        {
            if( geomodel().nb_interfaces() > 0 && geomodel().nb_layers() > 0 ) {                
                if( !is_geomodel_geology_valid( geomodel() ) ) {
                    set_invalid_model() ;
                }
            }
        }
        // Check that the model has a finite extension 
        // The boundary of the universe region is a one connected component 
        //  manifold closed surface 
        void test_finite_extension()
        {
            if( !is_region_valid( geomodel().universe() ) ) {
                set_invalid_model() ;                
            }
        }
        // Check geometrical-connectivity consistency
        void test_geometry_connectivity_consistency()
        {
            // Powerful check on relationships between GeoModelElements
            // containing the same point of the model
            if( !check_model_points_validity( geomodel() ) ) {
                set_invalid_model() ;
            }
            // Check on the edges - Only for Surfaces - A boundary edge has to be in a Line
            /// @todo Add consistency test for facets on boundary of Regions 
            for( index_t i = 0; i < geomodel().nb_surfaces(); ++i ) {
                if( !surface_boundary_valid( geomodel().surface( i ) ) ) {
                    set_invalid_model() ;
                }
            }
        }

        /// @todo Check that all Line segments correspond to a Surface
        /// edge that is on the boundary
        // With the current tests, it is possible we miss this problem,
        // but I am not sure (JP - 08/2015)
        void do_check_validity()
        {           
            test_global_element_access() ; 
            test_model_elements_validity() ;            
            test_geological_validity() ;
            test_finite_extension() ;
            test_geometry_connectivity_consistency() ;          
            test_non_manifold_edges() ;
            test_facet_intersections() ;
        }
    
        // Creates a Mesh from the GeoModel and triangulates it
        void create_model_mesh()
        {
            GEO::Logger::instance()->set_quiet( true ) ;

            bool connect_facets = false ;
            build_mesh_from_geomodel( geomodel(),
                triangulated_global_model_mesh_, connect_facets ) ;
            GEO::mesh_repair( triangulated_global_model_mesh_, MESH_REPAIR_TRIANGULATE ) ;

            GEO::Logger::instance()->set_quiet( false ) ;
        }

        /*!
        * @brief Returns true if there are non-manifold edges that are
        *        not in any Line of the model
        * @note Connect the facets of the global mesh
        * @note This is a quite expensive test.
        */
        void test_non_manifold_edges()
        {   
            create_model_mesh() ;
            std::vector< index_t > non_manifold_edges ;
            connect_mesh_facets_except_on_mesh_edges( 
                triangulated_global_model_mesh_, non_manifold_edges ) ;

            if( !non_manifold_edges.empty() ) {
                GEO::Logger::warn( "GeoModel" ) << non_manifold_edges.size()/2
                    << "non-manifold edges " << std::endl ;
                debug_save_non_manifold_edges( geomodel(), non_manifold_edges ) ;                
                
                set_invalid_model() ;
            }
        }

        /*!
        * @brief Returns true if there intersections between facets
        * @details Operates on the global mesh
        * @note This is a very expensive test.
        */
        void test_facet_intersections()
        {
            index_t nb_intersections = detect_intersecting_facets(
                geomodel_, triangulated_global_model_mesh_ ) ;
            
            if( nb_intersections > 0 ) {
                GEO::Logger::warn( "GeoModel" ) << nb_intersections
                    << " facet intersections " << std::endl ;
                set_invalid_model();
            } 
        }

        void set_invalid_model()
        {
            valid_ = false ;
        }

    private:
        const GeoModel& geomodel_ ;
        bool valid_ ;

        // Intermediate stuff used to compute validity
        GEO::Mesh triangulated_global_model_mesh_ ;
    };

 
    bool is_geomodel_valid( const GeoModel& GM, bool check_surface_intersections )
    {
        /// @todo Set some options here
        GeoModelValidityCheck validity_checker( GM ) ;

        bool valid = validity_checker.is_geomodel_valid() ;

        // Feedback 
        if( valid ) {
            GEO::Logger::out( "GeoModel" ) << "Model " << GM.name() << " is valid "
                << std::endl << std::endl ;
        } else {
            GEO::Logger::warn( "GeoModel" ) << "Model " << GM.name()
                << " is invalid " << std::endl << std::endl ;
        }
        return valid ;
    }
     

} // namespace RINGMesh
