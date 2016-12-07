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

#include <ringmesh/geomodel/geo_model_validity.h>

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
#include <geogram/basic/command_line.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_intersection.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/triangle_intersection.h>

#include <geogram/points/colocate.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geogram_extension/geogram_extension.h>
#include <ringmesh/geogram_extension/geogram_mesh_repair.h>
#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/geomodel/geo_model_api.h>
#include <ringmesh/geomodel/geo_model_entity.h>
#include <ringmesh/geomodel/geo_model_mesh_entity.h>
#include <ringmesh/geomodel/geo_model_geological_entity.h>
#include <ringmesh/mesh/mesh.h>

/*!
 * @file ringmesh/geomodel/geo_model_validity.cpp
 * @brief Implementation of functions to check the validity of GeoModels
 * @author Jeanne Pellerin
 * @todo Refactor the functions - reorganize to have a cleaner code.
 */

namespace {

    using namespace RINGMesh ;
    using GEO::index_t ;
    using GEO::vec3 ;

    typedef GeoModelMeshEntity GMME ;
    typedef std::string EntityType ;

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
        const GEO::Mesh& M,
        index_t f1,
        index_t f2,
        GEO::vector< GEO::TriangleIsect >& sym )
    {
        ringmesh_assert( M.facets.nb_vertices( f1 ) == 3 ) ;
        ringmesh_assert( M.facets.nb_vertices( f2 ) == 3 ) ;
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

    bool is_edge_on_line( const Line& line, index_t v0, index_t v1 )
    {
        if( v0 > v1 ) {
            std::swap( v0, v1 ) ;
        }
        index_t delta_i = v1 - v0 ;

        if( delta_i == 1 ) {
            // There is an edge if their indices in the Line are i and i+1
            return true ;
        } else if( line.is_closed() && delta_i == line.nb_vertices() - 2 ) {
            // If the Line is closed we can also have 0; n-2 or n-1; 1
            return true ;
        } else {
            // The two points are on the same line but
            // do not define an edge
            return false ;
        }
    }

    /*!
     * @brief Returns the Line identification if the given points define
     *       an edge of one of the Line of the geomodel
     * @param[in] geomodel The GeoModel to consider
     * @param[in] v0 Index in the geomodel of the edge first point
     * @param[in] v1 Index in the geomodel of the edge second point
     */
    bool is_edge_on_line( const GeoModel& geomodel, index_t v0, index_t v1 )
    {
        std::vector< GMEVertex > v0_line_bme ;
        geomodel.mesh.vertices.gme_type_vertices( Line::type_name_static(), v0,
            v0_line_bme ) ;
        if( v0_line_bme.empty() ) {
            return false ;
        }
        std::vector< GMEVertex > v1_line_bme ;
        geomodel.mesh.vertices.gme_type_vertices( Line::type_name_static(), v1,
            v1_line_bme ) ;
        if( v1_line_bme.empty() ) {
            return false ;
        }

        bool found_line = false ;
        for( index_t i = 0; i < v0_line_bme.size(); ++i ) {
            index_t line0_id = v0_line_bme[i].gme_id.index ;
            for( index_t j = 0; j < v1_line_bme.size(); ++j ) {
                if( line0_id == v1_line_bme[j].gme_id.index ) {
                    if( !is_edge_on_line( geomodel.line( line0_id ),
                        v0_line_bme[i].v_id, v1_line_bme[j].v_id ) ) {
                        return false ;
                    }
                    found_line = true ;
                    break ;
                }
            }
        }
        return found_line ;
    }

    /*!
     * @brief Returns the Line identification if the given points define
     *       an edge of one of the Line of the geomodel
     */
    bool is_edge_on_line( const GeoModel& geomodel, const vec3& p0, const vec3& p1 )
    {
        // Get the ids in the geomodel of these 2 points
        index_t v0 = geomodel.mesh.vertices.index( p0 ) ;
        index_t v1 = geomodel.mesh.vertices.index( p1 ) ;
        ringmesh_assert( v0 != NO_ID && v1 != NO_ID ) ;

        return is_edge_on_line( geomodel, v0, v1 ) ;
    }

    /*!
     * @brief Returns true if the facets @param f1 and @param f2
     *        of the mesh @param M share an edge
     *        that is on one Line of the boundary geomodel @param BM
     * @pre The mesh M is triangulated
     *
     */
    bool facets_share_line_edge(
        const GEO::Mesh& M,
        const GeoModel& BM,
        index_t f1,
        index_t f2 )
    {
        ringmesh_assert( M.facets.nb_vertices( f1 ) == 3 ) ;
        ringmesh_assert( M.facets.nb_vertices( f2 ) == 3 ) ;

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
                        ringmesh_assert( v10 != NO_ID && v11 != NO_ID ) ;

                        index_t v20 = BM.mesh.vertices.index( p20 ) ;
                        index_t v21 = BM.mesh.vertices.index( p21 ) ;

                        if( v10 == v20 && v11 == v21
                            && is_edge_on_line( BM, p20, p21 ) ) {
                            return true ;
                        }
                        if( v10 == v21 && v11 == v20
                            && is_edge_on_line( BM, p20, p21 ) ) {
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
    bool facets_are_adjacent( const GEO::Mesh& M, index_t f1, index_t f2 )
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
            const GEO::Mesh& M,
            const GeoModel& BM,
            GEO::vector< index_t >& has_isect )
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
        const GEO::Mesh& M_ ;
        const GeoModel& BM_ ;
        GEO::vector< index_t >& has_intersection_ ;
        GEO::vector< GEO::TriangleIsect > sym_ ;
    } ;

    void save_mesh_locating_geomodel_inconsistencies(
        const GEO::Mesh& mesh,
        const std::ostringstream& file )
    {
        if( GEO::CmdLine::get_arg_bool( "in:validity_save" ) ) {
            GEO::mesh_save( mesh, file.str() ) ;
        }
    }

    /** \note Copied from geogram
     * \brief Detect intersecting facets in a TRIANGULATED mesh
     * \param[in] M the mesh
     * \return number of intersecting facets
     */
    index_t detect_intersecting_facets( const GeoModel& geomodel, GEO::Mesh& M )
    {
        geo_assert( M.vertices.dimension() >= 3 ) ;

        GEO::vector< index_t > has_intersection ;
        StoreIntersections action( M, geomodel, has_intersection ) ;
        GEO::MeshFacetsAABB AABB( M ) ;
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
            save_mesh_locating_geomodel_inconsistencies( mesh, file ) ;
            Logger::out( "I/O" ) << std::endl ;
        }
        return nb_intersections ;
    }

    /***************************************************************************/

    /*!
     * @brief Check if entity @param is of the @param geomodel is in the
     *        in_boundary vector of entity @param in.
     */
    bool is_in_in_boundary( const GeoModel& geomodel, gme_t is, gme_t in )
    {
        const GeoModelMeshEntity& E = geomodel.mesh_entity( in ) ;
        for( index_t i = 0; i < E.nb_in_boundary(); ++i ) {
            if( E.in_boundary_gme( i ) == is ) {
                return true ;
            }
        }
        return false ;
    }

    void save_invalid_points(
        const std::ostringstream& file,
        const GeoModel& M,
        const std::vector< bool >& valid )
    {
        GEO::Mesh point_mesh ;
        for( index_t i = 0; i < valid.size(); ++i ) {
            if( !valid[i] ) {
                const vec3& V = M.mesh.vertices.vertex( i ) ;
                point_mesh.vertices.create_vertex( V.data() ) ;
            }
        }
        save_mesh_locating_geomodel_inconsistencies( point_mesh, file ) ;
    }

    /*!
     * @brief Check the geometrical-topological consistency of the geomodel
     * @details Verification is based on the information stored by the unique
     *          vertices of the geomodel which validity must be checked beforehand
     * @todo Check that the geomodel vertices are consistent with the geomodel_vertex_ids
     *       stored at by the GMME
     * @todo Implementation for regions
     * @todo Split in smaller functions
     */
    bool check_model_points_validity( const GeoModel& M )
    {
        // For all the vertices of the geomodel 
        // We check that the entities in which they are are consistent 
        // to have a valid B-Rep geomodel
        std::vector< bool > valid( M.mesh.vertices.nb(), true ) ;
        for( index_t i = 0; i < M.mesh.vertices.nb(); ++i ) {
            bool valid_vertex = true ;

            // Get the mesh entities in which this vertex is            
            index_t corner = NO_ID ;
            std::vector< index_t > lines ;
            std::vector< index_t > surfaces ;
            std::vector< index_t > regions ;

            std::vector< GMEVertex > bmes ;
            M.mesh.vertices.gme_vertices( i, bmes ) ;

            for( index_t j = 0; j < bmes.size(); ++j ) {
                const std::string& T = bmes[j].gme_id.type ;
                index_t id = bmes[j].gme_id.index ;
                if( T == Region::type_name_static() ) {
                    regions.push_back( id ) ;
                } else if( T == Surface::type_name_static() ) {
                    surfaces.push_back( id ) ;
                } else if( T == Line::type_name_static() ) {
                    lines.push_back( id ) ;
                } else if( T == Corner::type_name_static() ) {
                    if( corner != NO_ID ) {
                        Logger::warn( "GeoModel" ) << " Vertex " << i
                            << " is in at least 2 Corners" << std::endl ;
                        valid_vertex = false ;
                    } else {
                        corner = id ;
                    }
                } else {
                    Logger::warn( "GeoModel" ) << " Vertex " << i
                        << " is in no Entity of the Model" << std::endl ;
                    valid_vertex = false ;
                    break ;
                }
            }

            if( valid_vertex ) {
                if( surfaces.empty() ) {
                    if( regions.size() != 1 ) {
                        Logger::warn( "GeoModel" ) << " Vertex " << i << " is in "
                            << regions.size() << " Regions: " ;
                        for( index_t j = 0; j < surfaces.size(); ++j ) {
                            Logger::warn( "GeoModel" ) << regions[j] << " ; " ;
                        }
                        Logger::warn( "GeoModel" ) << std::endl ;
                        valid_vertex = false ;
                    } /// @todo Implement the other conditions for Region point validity
                } else if( corner == NO_ID && lines.empty() ) {
                    // This is a point on one SURFACE and only one
                    if( surfaces.size() != 1 ) {
                        Logger::warn( "GeoModel" ) << " Vertex " << i << " is in "
                            << surfaces.size() << " Surfaces: " ;
                        for( index_t j = 0; j < surfaces.size(); ++j ) {
                            Logger::warn( "GeoModel" ) << surfaces[j] << " ; " ;
                        }
                        Logger::warn( "GeoModel" ) << std::endl ;
                        valid_vertex = false ;
                    }
                } else if( corner == NO_ID && !lines.empty() ) {
                    // This is a point on one LINE 
                    if( lines.size() != 1 ) {
                        Logger::warn( "GeoModel" ) << " Vertex " << i << " is in "
                            << lines.size() << " Lines " ;
                        for( index_t j = 0; j < lines.size(); ++j ) {
                            Logger::warn( "GeoModel" ) << lines[j] << " ; " ;
                        }
                        Logger::warn( "GeoModel" ) << std::endl ;
                        valid_vertex = false ;
                    } else {
                        // This point must also be in at least one SURFACE
                        if( surfaces.empty() ) {
                            Logger::warn( "GeoModel" ) << " Vertex " << i
                                << " is in a Line but in no Surface " << std::endl ;
                            valid_vertex = false ;
                        }
                        // Check that one point is no more than twice in a SURFACE
                        for( index_t k = 0; k < surfaces.size(); ++k ) {
                            index_t nb = static_cast< index_t >( std::count(
                                surfaces.begin(), surfaces.end(), surfaces[k] ) ) ;
                            if( nb > 2 ) {
                                Logger::warn( "GeoModel" ) << " Vertex " << i
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
                                    Logger::warn( "GeoModel" ) << " Vertex " << i
                                        << " appears " << nb << " times in Surface "
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
                                gme_t s_id( Surface::type_name_static(),
                                    surfaces[k] ) ;
                                gme_t l_id( Line::type_name_static(), lines[l] ) ;
                                if( !is_in_in_boundary( M, s_id, l_id ) ) {
                                    Logger::warn( "GeoModel" )
                                        << " Inconsistent Line-Surface connectivity "
                                        << " Vertex " << i << " shows that " << s_id
                                        << " must be in the boundary of " << l_id
                                        << std::endl ;
                                    valid_vertex = false ;
                                }
                            }
                        }
                    }
                } else if( corner != NO_ID ) {
                    // This is one point at a CORNER
                    // It must be in at least one LINE
                    if( lines.empty() ) {
                        Logger::warn( "GeoModel" ) << " Vertex " << i
                            << " is at a Corner but in no Line " << std::endl ;
                        valid_vertex = false ;
                    } else {
                        if( lines.size() < 2 ) {
                            Logger::warn( "GeoModel" ) << " Vertex " << i
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
                                    Logger::warn( "GeoModel" ) << " Vertex " << i
                                        << " is twice in Line " << lines[k]
                                        << std::endl ;
                                    valid_vertex = false ;
                                }
                            }
                            if( nb > 2 ) {
                                Logger::warn( "GeoModel" ) << " Vertex " << i
                                    << " appears " << nb << " times in Line "
                                    << lines[k] << std::endl ;
                                valid_vertex = false ;
                                break ;
                            }
                        }
                        // Check that all the lines are in in_boundary of this corner
                        for( index_t k = 0; k < lines.size(); ++k ) {
                            gme_t l_id( Line::type_name_static(), lines[k] ) ;
                            gme_t c_id( Corner::type_name_static(), corner ) ;
                            if( !is_in_in_boundary( M, l_id, c_id ) ) {
                                Logger::warn( "GeoModel" )
                                    << " Inconsistent Line-Corner connectivity "
                                    << " vertex " << i << " shows that " << l_id
                                    << " must be in the boundary of " << c_id
                                    << std::endl ;
                                valid_vertex = false ;
                            }
                        }
                    }
                    // It must also be in a least one surface ? perhaps 2
                    if( surfaces.empty() ) {
                        Logger::warn( "GeoModel" ) << " Vertex " << i
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
            file << validity_errors_directory << "/invalid_global_vertices.mesh" ;
            save_invalid_points( file, M, valid ) ;

            Logger::warn( "GeoModel" ) << nb_invalid << " invalid vertices "
                << std::endl << "Saved in file: " << file.str() << std::endl ;
            return false ;
        } else {
            return true ;
        }
    }

    void save_edges(
        const std::ostringstream& file,
        const GeoModel& M,
        const std::vector< index_t >& e )
    {
        GEO::Mesh edge_mesh ;
        index_t previous_vertex_id = NO_ID ;
        for( index_t i = 0; i < e.size(); ++i ) {
            index_t cur_vertex_id = edge_mesh.vertices.create_vertex(
                M.mesh.vertices.vertex( e[i] ).data() ) ;
            if( i % 2 == 0 ) {
                ringmesh_assert( previous_vertex_id == NO_ID ) ;
                previous_vertex_id = cur_vertex_id ;
                continue ;
            }
            ringmesh_assert( previous_vertex_id != NO_ID ) ;
            edge_mesh.edges.create_edge( previous_vertex_id, cur_vertex_id ) ;
            previous_vertex_id = NO_ID ;
        }
        save_mesh_locating_geomodel_inconsistencies( edge_mesh, file ) ;
    }

    void save_facets(
        const std::string& file,
        const Surface& surface,
        const std::vector< index_t >& facets )
    {
        GEO::Mesh mesh ;
        for( index_t f = 0; f < facets.size(); ++f ) {
            index_t cur_facet = facets[f] ;
            index_t nb_vertices_in_facet = surface.nb_mesh_element_vertices(
                cur_facet ) ;
            GEO::vector< index_t > vertices ;
            vertices.reserve( nb_vertices_in_facet ) ;
            for( index_t v = 0; v < nb_vertices_in_facet; v++ ) {
                index_t new_vertex = mesh.vertices.create_vertex(
                    surface.mesh_element_vertex( cur_facet, v ).data() ) ;
                vertices.push_back( new_vertex ) ;
            }
            mesh.facets.create_polygon( vertices ) ;
        }
        GEO::mesh_save( mesh, file ) ;
    }

    /*!
     * @brief Check boundary of a surface
     * @details All the edges on the boundary of a surface must be in a Line
     *          of the associated geomodel
     *          The Line boundaries must form a closed manifold line.
     */
    bool surface_boundary_valid( const Surface& S )
    {
        const GeoModelMeshVertices& geomodel_vertices = S.geomodel().mesh.vertices ;
        std::vector< index_t > invalid_corners ;
        for( index_t f = 0; f < S.nb_mesh_elements(); ++f ) {
            for( index_t v = 0; v < S.nb_mesh_element_vertices( f ); ++v ) {
                if( S.facet_adjacent_index( f, v ) == NO_ID
                    && !is_edge_on_line( S.geomodel(),
                        geomodel_vertices.geomodel_vertex_id( S.gme_id(), f, v ),
                        geomodel_vertices.geomodel_vertex_id( S.gme_id(), f,
                            S.next_facet_vertex_index( f, v ) ) ) ) {
                    invalid_corners.push_back(
                        geomodel_vertices.geomodel_vertex_id( S.gme_id(), f, v ) ) ;
                    invalid_corners.push_back(
                        geomodel_vertices.geomodel_vertex_id( S.gme_id(), f,
                            S.next_facet_vertex_index( f, v ) ) ) ;
                }
            }
        }
        if( !invalid_corners.empty() ) {
            std::ostringstream file ;
            file << validity_errors_directory << "/invalid_boundary_surface_"
                << S.index() << ".mesh" ;
            save_edges( file, S.geomodel(), invalid_corners ) ;

            Logger::warn( "GeoModel" ) << " Invalid surface boundary: "
                << invalid_corners.size() / 2 << " boundary edges of " << S.gme_id()
                << "  are in no line of the geomodel " << std::endl
                << " Saved in file: " << file.str() << std::endl ;
            return false ;
        } else {
            return true ;
        }
    }

    /*!
     * @brief Save in a .lin file the
     */
    void debug_save_non_manifold_edges(
        const GeoModel& geomodel,
        const std::vector< index_t >& edge_vertices )
    {
        std::ostringstream file_name(
            validity_errors_directory + "/non_manifold_edges.mesh" ) ;

        save_edges( file_name, geomodel, edge_vertices ) ;
    }

    bool is_surface_conformal_to_volume(
        const Surface& surface,
        const ColocaterANN& cell_facet_barycenter_ann )
    {
        std::vector< index_t > unconformal_facets ;
        for( index_t f = 0; f < surface.nb_mesh_elements(); f++ ) {
            vec3 center = surface.mesh_element_barycenter( f ) ;
            std::vector< index_t > result ;
            if( !cell_facet_barycenter_ann.get_neighbors( center, result,
                surface.geomodel().epsilon() ) ) {
                unconformal_facets.push_back( f ) ;
            }
        }
        if( !unconformal_facets.empty() ) {
            std::ostringstream file ;
            file << validity_errors_directory << "/unconformal_surface_"
                << surface.index() << ".mesh" ;
            save_facets( file.str(), surface, unconformal_facets ) ;

            Logger::warn( "GeoModel" ) << " Unconformal surface: "
                << unconformal_facets.size() << " facets of " << surface.gme_id()
                << " are unconformal with the geomodel cells " << std::endl
                << " Saved in file: " << file.str() << std::endl ;
            return false ;
        } else {
            return true ;
        }
    }

    /*!
     * @brief Implementation class for validity checks on a GeoModel
     */
    class GeoModelValidityCheck {
    public:
        GeoModelValidityCheck(
            const GeoModel& geomodel,
            bool check_surface_intersections )
            :
                geomodel_( geomodel ),
                valid_( true ),
                check_surface_intersections_( check_surface_intersections )
        {
            // Ensure that the geomodel vertices are computed and up-to-date
            // Without that we cannot do anything        
            geomodel_.mesh.vertices.test_and_initialize() ;
            geomodel_.mesh.cells.test_and_initialize() ;
        }

        bool is_geomodel_valid()
        {
            do_check_validity() ;
            return valid_ ;
        }

    private:
        void do_check_validity()
        {
            test_model_entities_validity() ;
            test_finite_extension() ;
            test_geometry_connectivity_consistency() ;
            test_non_manifold_edges() ;
            if( check_surface_intersections_ ) {
                test_facet_intersections() ;
            }
        }
        /*! 
         * @brief Verify the validity of all GeoModelEntities
         */
        void test_model_entities_validity()
        {
            if( !are_geomodel_meshed_entities_valid( geomodel_ ) ) {
                set_invalid_model() ;
            }
            if( !are_geomodel_geological_entities_valid( geomodel_ ) ) {
                set_invalid_model() ;
            }
        }
        /*!
         * @brief Check that the geomodel has a finite extension 
         * @details The boundary of the universe region is a one connected component 
         * manifold closed surface.
         */
        void test_finite_extension()
        {
            if( !geomodel_.universe().is_valid() ) {
                set_invalid_model() ;
            }
        }
        /*!
         * Check geometrical-connectivity consistency
         * @todo Add consistency test for facets on boundary of Regions 
         * @todo Check that all Line segments correspond to a Surface
         *  edge that is on the boundary.
         */
        void test_geometry_connectivity_consistency()
        {
            // Check relationships between GeoModelEntities
            // sharing the same point of the geomodel
            if( !check_model_points_validity( geomodel_ ) ) {
                set_invalid_model() ;
            }
            // Check on that Surface edges are in a Line
            for( index_t i = 0; i < geomodel_.nb_surfaces(); ++i ) {
                if( !surface_boundary_valid( geomodel_.surface( i ) ) ) {
                    set_invalid_model() ;
                }
            }
            if( geomodel_.mesh.cells.nb() > 0 ) {
                // Check the consistency between Surface facets and Region cell facets
                const ColocaterANN& ann =
                    geomodel_.mesh.cells.cell_facet_colocater() ;
                for( index_t i = 0; i < geomodel_.nb_surfaces(); ++i ) {
                    if( !is_surface_conformal_to_volume( geomodel_.surface( i ),
                        ann ) ) {
                        set_invalid_model() ;
                    }
                }
            }
        }
        /*! 
         * @brief Creates a Mesh from the GeoModel and triangulates it
         */
        void create_model_mesh()
        {
            bool logger_status = Logger::instance()->is_quiet() ;
            Logger::instance()->set_quiet( true ) ;

            bool connect_facets = false ;
            build_mesh_from_geomodel( geomodel_, triangulated_global_model_mesh_,
                connect_facets ) ;
            GEO::mesh_repair( triangulated_global_model_mesh_,
                GEO::MESH_REPAIR_TRIANGULATE ) ;

            Logger::instance()->set_quiet( logger_status ) ;
        }
        /*!
         * @brief Returns true if there are non-manifold edges that are
         *        not in any Line of the geomodel
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
                Logger::warn( "GeoModel" ) << non_manifold_edges.size() / 2
                    << "non-manifold edges " << std::endl ;
                debug_save_non_manifold_edges( geomodel_, non_manifold_edges ) ;

                set_invalid_model() ;
            }
        }
        /*!
         * @brief Returns true if there are intersections between facets
         * @details Operates on the global mesh
         * @note This is a very expensive test.
         */
        void test_facet_intersections()
        {
            index_t nb_intersections = detect_intersecting_facets( geomodel_,
                triangulated_global_model_mesh_ ) ;

            if( nb_intersections > 0 ) {
                Logger::warn( "GeoModel" ) << nb_intersections
                    << " facet intersections " << std::endl ;
                set_invalid_model() ;
            }
        }
        void set_invalid_model()
        {
            valid_ = false ;
        }

    private:
        const GeoModel& geomodel_ ;
        bool valid_ ;
        bool check_surface_intersections_ ;
        // Global mesh of the GeoModel used for some validity checks
        GEO::Mesh triangulated_global_model_mesh_ ;
    } ;

} // anonymous namespace

namespace RINGMesh {

    void set_validity_errors_directory( const std::string& directory )
    {
        // If trailing / or \ is not removed, the test fails on Windows
        std::string copy( directory ) ;
        if( *copy.rbegin() == '/' || *copy.rbegin() == '\\' ) {
            copy.erase( copy.end() - 1 ) ;
        }
        if( GEO::FileSystem::is_directory( copy ) ) {
            validity_errors_directory = copy + '/' ;
        }
    }

    bool are_geomodel_meshed_entities_valid( const GeoModel& geomodel )
    {
        const std::vector< EntityType >& meshed_types =
            EntityTypeManager::mesh_entity_types() ;
        index_t count_invalid = 0 ;
        for( index_t i = 0; i < meshed_types.size(); ++i ) {
            const EntityType& type = meshed_types[i] ;
            index_t nb_entities = geomodel.nb_mesh_entities( type ) ;
            for( index_t i = 0; i < nb_entities; ++i ) {
                const GeoModelEntity& E = geomodel.mesh_entity( type, i ) ;
                if( !E.is_valid() ) {
                    count_invalid++ ;
                }
            }
        }
        if( count_invalid != 0 ) {
            Logger::warn( "GeoModel" ) << count_invalid
                << " mesh entities of the geomodel are invalid " << std::endl ;
        }
        return count_invalid == 0 ;
    }

    bool are_geomodel_geological_entities_valid( const GeoModel& geomodel )
    {
        const std::vector< EntityType >& geological_types =
            geomodel.entity_type_manager().geological_entity_types() ;
        index_t count_invalid = 0 ;
        for( index_t i = 0; i < geological_types.size(); ++i ) {
            const EntityType& type = geological_types[i] ;
            index_t nb_entities = geomodel.nb_geological_entities( type ) ;
            for( index_t i = 0; i < nb_entities; ++i ) {
                const GeoModelEntity& E = geomodel.geological_entity( type, i ) ;
                if( !E.is_valid() ) {
                    count_invalid++ ;
                }
            }
        }
        if( count_invalid != 0 ) {
            Logger::warn( "GeoModel" ) << count_invalid
                << " geological entities of the geomodel are invalid " << std::endl ;
        }
        return count_invalid == 0 ;
    }

    bool is_geomodel_valid( const GeoModel& GM )
    {
        GeoModelValidityCheck validity_checker( GM,
            GEO::CmdLine::get_arg_bool( "in:intersection_check" ) ) ;

        bool valid = validity_checker.is_geomodel_valid() ;

        if( valid ) {
            Logger::out( "GeoModel" ) << "Model " << GM.name() << " is valid "
                << std::endl << std::endl ;
        } else {
            Logger::warn( "GeoModel" ) << "Model " << GM.name() << " is invalid "
                << std::endl << std::endl ;
        }
        return valid ;
    }

    bool check_volume_watertightness( const GeoModel& geomodel, const gme_t& gme_id )
    {
        std::vector< gme_t > volume_boundaries ;

        // Check if the given Entity is a MeshEntity or a GeologicalEntity
        // or the Universe and fill the volume_boundaries vector.
        if( gme_id.type == geomodel.universe().type_name() ) {
            index_t nb_boundaries = geomodel.universe().nb_boundaries() ;
            volume_boundaries.resize( nb_boundaries ) ;
            for( index_t b = 0; b < nb_boundaries; b++ ) {
                volume_boundaries[b] = geomodel.universe().boundary_gme( b ) ;
            }
        } else if( geomodel.entity_type_manager().is_mesh_entity_type(
            gme_id.type ) ) {
            index_t nb_boundaries = geomodel.mesh_entity( gme_id ).nb_boundaries() ;
            volume_boundaries.resize( nb_boundaries ) ;
            for( index_t b = 0; b < nb_boundaries; b++ ) {
                volume_boundaries[b] = geomodel.mesh_entity( gme_id ).boundary_gme(
                    b ) ;
            }
        } else {
            Logger::warn( "GeoModel" ) << "Checking for volume watertightness of "
                "a geological entity (" << gme_id << ") is not yet implemented."
                << std::endl ;
            ringmesh_assert_not_reached ;
        }

        if( volume_boundaries.empty() ) {
            Logger::warn( "GeoModel" ) << gme_id << " has no boundary Surface"
                << std::endl ;
            return false ;
        } else {
            GEO::Mesh mesh ;
            bool logger_status = Logger::instance()->is_quiet() ;

            Logger::instance()->set_quiet( true ) ;
            build_mesh_from_geomodel_mesh_entities( geomodel, volume_boundaries,
                mesh ) ;
            GEO::mesh_repair( mesh ) ;
            Logger::instance()->set_quiet( logger_status ) ;

            bool valid = true ;
            index_t nb_cc = GEO::mesh_nb_connected_components( mesh ) ;
            signed_index_t nb_b = GEO::mesh_nb_borders( mesh ) ;
            if( nb_cc != 1 ) {
                Logger::warn( "GeoModel" ) << " Surface boundary of " << gme_id
                    << " has " << nb_cc << " connected components " << std::endl ;
                valid = false ;
            }
            if( nb_b != 0 ) {
                Logger::warn( "GeoModel" ) << " Surface boundary of " << gme_id
                    << " has " << nb_b << " border connected components "
                    << std::endl ;
                valid = false ;
            }
            if( !valid ) {
                std::ostringstream file ;
                file << validity_errors_directory << "/boundary_surface_region_"
                    << gme_id.index << ".mesh" ;
                if( GEO::CmdLine::get_arg_bool( "in:validity_save" ) ) {
                    GEO::mesh_save( mesh, file.str() ) ;
                }
                return false ;
            } else {
                return true ;
            }
        }
    }

} // namespace RINGMesh
