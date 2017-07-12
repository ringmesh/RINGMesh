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

#include <ringmesh/geomodel/geomodel_validity.h>

#include <thread>

#include <geogram/mesh/triangle_intersection.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>

#include <ringmesh/mesh/geogram_mesh.h>
#include <ringmesh/mesh/geogram_mesh_builder.h>

/*!
 * @file ringmesh/geomodel/geomodel_validity.cpp
 * @brief Implementation of functions to check the validity of GeoModels
 * @author Jeanne Pellerin
 * @todo Refactor the functions - reorganize to have a cleaner code.
 */

namespace {
    using namespace RINGMesh;

    bool triangles_intersect(
        const GeoModel& geomodel,
        const GeoModelMeshPolygons& polygons,
        index_t triangle1,
        index_t triangle2 )
    {
        ringmesh_assert( polygons.nb_vertices( triangle1 ) == 3 );
        ringmesh_assert( polygons.nb_vertices( triangle2 ) == 3 );
        const GeoModelMeshVertices& vertices = geomodel.mesh.vertices;
        const vec3& p1 = vertices.vertex( polygons.vertex( triangle1, 0 ) );
        const vec3& p2 = vertices.vertex( polygons.vertex( triangle1, 1 ) );
        const vec3& p3 = vertices.vertex( polygons.vertex( triangle1, 2 ) );

        const vec3& q1 = vertices.vertex( polygons.vertex( triangle2, 0 ) );
        const vec3& q2 = vertices.vertex( polygons.vertex( triangle2, 1 ) );
        const vec3& q3 = vertices.vertex( polygons.vertex( triangle2, 2 ) );
        GEO::vector< GEO::TriangleIsect > sym;
        return triangles_intersections( p1, p2, p3, q1, q2, q3, sym );
    }

    bool triangle_quad_intersect(
        const GeoModel& geomodel,
        const GeoModelMeshPolygons& polygons,
        index_t triangle,
        index_t quad )
    {
        ringmesh_assert( polygons.nb_vertices( triangle ) == 3 );
        ringmesh_assert( polygons.nb_vertices( quad ) == 4 );
        const GeoModelMeshVertices& vertices = geomodel.mesh.vertices;
        const vec3& p1 = vertices.vertex( polygons.vertex( triangle, 0 ) );
        const vec3& p2 = vertices.vertex( polygons.vertex( triangle, 1 ) );
        const vec3& p3 = vertices.vertex( polygons.vertex( triangle, 2 ) );

        const vec3& q1 = vertices.vertex( polygons.vertex( quad, 0 ) );
        const vec3& q2 = vertices.vertex( polygons.vertex( quad, 1 ) );
        const vec3& q3 = vertices.vertex( polygons.vertex( quad, 2 ) );
        const vec3& q4 = vertices.vertex( polygons.vertex( quad, 3 ) );
        GEO::vector< GEO::TriangleIsect > sym;
        if( triangles_intersections( p1, p2, p3, q1, q2, q3, sym ) ) {
            return true;
        }
        if( triangles_intersections( p1, p2, p3, q1, q3, q4, sym ) ) {
            return true;
        }
        return false;
    }

    bool quad_quad_intersect(
        const GeoModel& geomodel,
        const GeoModelMeshPolygons& polygons,
        index_t quad1,
        index_t quad2 )
    {
        ringmesh_assert( polygons.nb_vertices( quad1 ) == 4 );
        ringmesh_assert( polygons.nb_vertices( quad2 ) == 4 );
        const GeoModelMeshVertices& vertices = geomodel.mesh.vertices;
        const vec3& p1 = vertices.vertex( polygons.vertex( quad1, 0 ) );
        const vec3& p2 = vertices.vertex( polygons.vertex( quad1, 1 ) );
        const vec3& p3 = vertices.vertex( polygons.vertex( quad1, 2 ) );
        const vec3& p4 = vertices.vertex( polygons.vertex( quad1, 3 ) );

        const vec3& q1 = vertices.vertex( polygons.vertex( quad2, 0 ) );
        const vec3& q2 = vertices.vertex( polygons.vertex( quad2, 1 ) );
        const vec3& q3 = vertices.vertex( polygons.vertex( quad2, 2 ) );
        const vec3& q4 = vertices.vertex( polygons.vertex( quad2, 3 ) );
        GEO::vector< GEO::TriangleIsect > sym;
        if( triangles_intersections( p1, p2, p3, q1, q2, q3, sym ) ) {
            return true;
        }
        if( triangles_intersections( p1, p2, p3, q1, q3, q4, sym ) ) {
            return true;
        }
        if( triangles_intersections( p1, p3, p4, q1, q2, q3, sym ) ) {
            return true;
        }
        if( triangles_intersections( p1, p3, p4, q1, q3, q4, sym ) ) {
            return true;
        }
        return false;
    }

    bool is_edge_on_line( const Line& line, index_t v0, index_t v1 )
    {
        if( v0 > v1 ) {
            std::swap( v0, v1 );
        }
        index_t delta_i = v1 - v0;

        if( delta_i == 1 ) {
            // There is an edge if their indices in the Line are i and i+1
            return true;
        } else if( line.is_closed() && delta_i == line.nb_vertices() - 2 ) {
            // If the Line is closed we can also have 0; n-2 or n-1; 1
            return true;
        } else {
            // The two points are on the same line but
            // do not define an edge
            return false;
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
        MeshEntityType line_type = Line::type_name_static();
        std::vector< GMEVertex > v0_line_bme =
            geomodel.mesh.vertices.gme_type_vertices( line_type, v0 );
        if( v0_line_bme.empty() ) {
            return false;
        }
        std::vector< GMEVertex > v1_line_bme =
            geomodel.mesh.vertices.gme_type_vertices( line_type, v1 );
        if( v1_line_bme.empty() ) {
            return false;
        }

        bool found_line = false;
        for( const GMEVertex& vertex0 : v0_line_bme ) {
            index_t line0_id = vertex0.gmme.index();
            for( const GMEVertex& vertex1 : v1_line_bme ) {
                if( line0_id == vertex1.gmme.index() ) {
                    if( !is_edge_on_line( geomodel.line( line0_id ), vertex0.v_index,
                        vertex1.v_index ) ) {
                        return false;
                    }
                    found_line = true;
                    break;
                }
            }
        }
        return found_line;
    }

    /*!
     * @brief Returns true if the polygons @param p1 and @param p2
     *        of the mesh @param polygons share an edge
     *        that is on one Line of the boundary geomodel @param geomodel
     *
     */
    bool polygons_share_line_edge(
        const GeoModel& geomodel,
        const GeoModelMeshPolygons& polygons,
        index_t p1,
        index_t p2 )
    {
        // Only test the edges on boundary
        for( index_t v1 = 0; v1 < polygons.nb_vertices( p1 ); v1++ ) {
            if( polygons.adjacent( p1, v1 ) != NO_ID ) {
                continue;
            }
            index_t v10 = polygons.vertex( p1, v1 );
            index_t v11 = polygons.vertex( p1,
                ( v1 + 1 ) % polygons.nb_vertices( p1 ) );
            for( index_t v2 = 0; v2 < polygons.nb_vertices( p2 ); v2++ ) {
                if( polygons.adjacent( p2, v2 ) != NO_ID ) {
                    continue;
                }
                index_t v20 = polygons.vertex( p2, v2 );
                index_t v21 = polygons.vertex( p2,
                    ( v2 + 1 ) % polygons.nb_vertices( p2 ) );

                if( ( v10 == v20 && v11 == v21 ) || ( v10 == v21 && v11 == v20 ) ) {
                    if( is_edge_on_line( geomodel, v20, v21 ) ) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    bool polygons_are_adjacent(
        const GeoModelMeshPolygons& polygons,
        index_t p1,
        index_t p2 )
    {
        if( p1 == p2 ) {
            return true;
        }
        for( index_t v = 0; v < polygons.nb_vertices( p1 ); v++ ) {
            if( polygons.adjacent( p1, v ) == p2 ) {
                return true;
            }
        }
        return false;
    }

    /*!
     * @brief Action class for storing intersections when traversing
     *  a AABBTree.
     */
    class StoreIntersections {
    public:
        /*!
         * @brief Constructs the StoreIntersections
         * @param[in] geomodel the geomodel
         * @param[out] has_isect the flag that indicates for each polygon
         *  whether it has intersections
         */
        StoreIntersections(
            const GeoModel& geomodel,
            std::vector< bool >& has_isect )
            :
                geomodel_( geomodel ),
                polygons_( geomodel.mesh.polygons ),
                has_intersection_( has_isect )
        {
            has_intersection_.assign( polygons_.nb(), 0 );
        }

        /*!
         * @brief Determines the intersections between two polygons
         * @details It is a callback for AABBTree traversal
         * @param[in] p1 index of the first polygon
         * @param[in] p2 index of the second polygon
         */
        void operator()( index_t p1, index_t p2 )
        {
            if( p1 == p2 || polygons_are_adjacent( polygons_, p1, p2 )
                || polygons_share_line_edge( geomodel_, polygons_, p1, p2 ) ) {
                return;
            }

            if( is_triangle( p1 ) ) {
                if( is_triangle( p2 ) ) {
                    if( triangles_intersect( geomodel_, polygons_, p1, p2 ) ) {
                        has_intersection_[p1] = 1;
                        has_intersection_[p2] = 1;
                    }
                } else if( is_quad( p2 ) ) {
                    if( triangle_quad_intersect( geomodel_, polygons_, p1, p2 ) ) {
                        has_intersection_[p1] = 1;
                        has_intersection_[p2] = 1;
                    }
                } else {
                    ringmesh_assert_not_reached;
                }
            } else if( is_quad( p1 ) ) {
                if( is_triangle( p2 ) ) {
                    if( triangle_quad_intersect( geomodel_, polygons_, p2, p1 ) ) {
                        has_intersection_[p1] = 1;
                        has_intersection_[p2] = 1;
                    }
                } else if( is_quad( p2 ) ) {
                    if( quad_quad_intersect( geomodel_, polygons_, p1, p2 ) ) {
                        has_intersection_[p1] = 1;
                        has_intersection_[p2] = 1;
                    }
                } else {
                    ringmesh_assert_not_reached;
                }
            } else {
                ringmesh_assert_not_reached;
            }
        }

        bool is_triangle( index_t p ) const
        {
            index_t index;
            return polygons_.type( p, index ) == PolygonType::TRIANGLE;
        }
        bool is_quad( index_t p ) const
        {
            index_t index;
            return polygons_.type( p, index ) == PolygonType::QUAD;
        }

    private:
        const GeoModel& geomodel_;
        const GeoModelMeshPolygons& polygons_;
        std::vector< bool >& has_intersection_;
    };

    void save_mesh_locating_geomodel_inconsistencies(
        const GEO::Mesh& mesh,
        const std::ostringstream& file )
    {
        if( GEO::CmdLine::get_arg_bool( "validity_save" ) ) {
            GEO::mesh_save( mesh, file.str() );
        }
    }

    /***************************************************************************/

    /*!
     * @brief Check if entity @param is of the @param geomodel is in the
     *        incident_entity vector of entity @param in.
     */
    bool is_in_incident_entity(
        const GeoModel& geomodel,
        const gmme_id& is,
        const gmme_id& in )
    {
        const GeoModelMeshEntity& E = geomodel.mesh_entity( in );
        for( index_t i = 0; i < E.nb_incident_entities(); ++i ) {
            if( E.incident_entity_gmme( i ) == is ) {
                return true;
            }
        }
        return false;
    }

    void save_invalid_points(
        const std::ostringstream& file,
        const GeoModel& geomodel,
        const std::vector< bool >& valid )
    {
        GEO::Mesh point_mesh;
        for( index_t i = 0; i < valid.size(); ++i ) {
            if( !valid[i] ) {
                const vec3& V = geomodel.mesh.vertices.vertex( i );
                point_mesh.vertices.create_vertex( V.data() );
            }
        }
        save_mesh_locating_geomodel_inconsistencies( point_mesh, file );
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
    bool check_model_points_validity( const GeoModel& geomodel )
    {
        // For all the vertices of the geomodel 
        // We check that the entities in which they are are consistent 
        // to have a valid B-Rep geomodel
        std::vector< bool > valid( geomodel.mesh.vertices.nb(), true );
        for( index_t i = 0; i < geomodel.mesh.vertices.nb(); ++i ) {
            bool valid_vertex = true;

            // Get the mesh entities in which this vertex is            
            index_t corner = NO_ID;
            std::vector< index_t > lines;
            std::vector< index_t > surfaces;
            std::vector< index_t > regions;

            const std::vector< GMEVertex >& bmes =
                geomodel.mesh.vertices.gme_vertices( i );

            for( const GMEVertex& vertex : bmes ) {
                const MeshEntityType& T = vertex.gmme.type();
                index_t id = vertex.gmme.index();
                if( T == Region::type_name_static() ) {
                    regions.push_back( id );
                } else if( T == Surface::type_name_static() ) {
                    surfaces.push_back( id );
                } else if( T == Line::type_name_static() ) {
                    lines.push_back( id );
                } else if( T == Corner::type_name_static() ) {
                    if( corner != NO_ID ) {
                        Logger::warn( "GeoModel", " Vertex ", i,
                            " is in at least 2 Corners" );
                        valid_vertex = false;
                    } else {
                        corner = id;
                    }
                } else {
                    Logger::warn( "GeoModel", " Vertex ", i,
                        " is in no Entity of the Model" );
                    valid_vertex = false;
                    break;
                }
            }

            if( valid_vertex ) {
                if( surfaces.empty() ) {
                    if( regions.size() != 1 ) {
                        std::ostringstream oss;
                        oss << " Vertex " << i << " is in " << regions.size()
                            << " Regions: ";
                        for( index_t region : regions ) {
                            oss << region << " ; ";
                        }
                        Logger::warn( "GeoModel", oss.str() );
                        valid_vertex = false;
                    } /// @todo Implement the other conditions for Region point validity
                } else if( corner == NO_ID && lines.empty() ) {
                    // This is a point on one SURFACE and only one
                    if( surfaces.size() != 1 ) {
                        std::ostringstream oss;
                        oss << " Vertex " << i << " is in " << surfaces.size()
                            << " Surfaces: ";
                        for( index_t surface : surfaces ) {
                            oss << surface << " ; ";
                        }
                        Logger::warn( "GeoModel", oss.str() );
                        valid_vertex = false;
                    }
                } else if( corner == NO_ID && !lines.empty() ) {
                    // This is a point on one LINE 
                    if( lines.size() != 1 ) {
                        std::ostringstream oss;
                        oss << " Vertex " << i << " is in " << lines.size()
                            << " Lines ";
                        for( index_t line : lines ) {
                            oss << line << " ; ";
                        }
                        Logger::warn( "GeoModel", oss.str() );
                        valid_vertex = false;
                    } else {
                        // This point must also be in at least one SURFACE
                        if( surfaces.empty() ) {
                            Logger::warn( "GeoModel", " Vertex ", i,
                                " is in a Line but in no Surface " );
                            valid_vertex = false;
                        }
                        // Check that one point is no more than twice in a SURFACE
                        for( index_t surface : surfaces ) {
                            index_t nb = static_cast< index_t >( std::count(
                                surfaces.begin(), surfaces.end(), surface ) );
                            if( nb > 2 ) {
                                Logger::warn( "GeoModel", " Vertex ", i, " is ", nb,
                                    " times in Surface ",
                                    geomodel.surface( surface ).gmme() );
                                valid_vertex = false;
                            } else if( nb == 2 ) {
                                // If a point is twice in a SURFACE, it must be
                                // on an internal boundary Line.
                                bool internal_boundary = false;
                                for( index_t line : lines ) {
                                    if( geomodel.line( line ).is_inside_border(
                                        geomodel.surface( surface ) ) ) {
                                        internal_boundary = true;
                                        break;
                                    }
                                }
                                if( !internal_boundary ) {
                                    Logger::warn( "GeoModel", " Vertex ", i,
                                        " appears ", nb, " times in Surface ",
                                        geomodel.surface( surface ).gmme() );
                                    valid_vertex = false;
                                }
                            }
                        }
                        // Check that all the surfaces are in incident_entity of all
                        // the lines 
                        for( index_t surface : surfaces ) {
                            for( index_t line : lines ) {
                                gmme_id s_id( Surface::type_name_static(), surface );
                                gmme_id l_id( Line::type_name_static(), line );
                                if( !is_in_incident_entity( geomodel, s_id, l_id ) ) {
                                    Logger::warn( "GeoModel",
                                        " Inconsistent Line-Surface connectivity ",
                                        " Vertex ", i, " shows that ", s_id,
                                        " must be in the boundary of ", l_id );
                                    valid_vertex = false;
                                }
                            }
                        }
                    }
                } else if( corner != NO_ID ) {
                    // This is one point at a CORNER
                    // It must be in at least one LINE
                    if( lines.empty() ) {
                        Logger::warn( "GeoModel", " Vertex ", i,
                            " is at a Corner but in no Line " );
                        valid_vertex = false;
                    } else {
                        if( lines.size() < 2 ) {
                            Logger::warn( "GeoModel", " Vertex ", i,
                                " is in at a Corner but in one Line only: ",
                                lines.front() );
                            valid_vertex = false;
                        }
                        // Check that a point is no more than twice in a LINE
                        for( index_t line : lines ) {
                            index_t nb = static_cast< index_t >( std::count(
                                lines.begin(), lines.end(), line ) );
                            if( nb == 2 ) {
                                // The line must be closed
                                if( !geomodel.line( line ).is_closed() ) {
                                    Logger::warn( "GeoModel", " Vertex ", i,
                                        " is twice in Line ", line );
                                    valid_vertex = false;
                                }
                            }
                            if( nb > 2 ) {
                                Logger::warn( "GeoModel", " Vertex ", i, " appears ",
                                    nb, " times in Line ", line );
                                valid_vertex = false;
                                break;
                            }
                        }
                        // Check that all the lines are in incident_entity of this corner
                        for( index_t line : lines ) {
                            gmme_id l_id( Line::type_name_static(), line );
                            gmme_id c_id( Corner::type_name_static(), corner );
                            if( !is_in_incident_entity( geomodel, l_id, c_id ) ) {
                                Logger::warn( "GeoModel",
                                    " Inconsistent Line-Corner connectivity ",
                                    " vertex ", i, " shows that ", l_id,
                                    " must be in the boundary of ", c_id );
                                valid_vertex = false;
                            }
                        }
                    }
                    // It must also be in a least one surface ? perhaps 2
                    if( surfaces.empty() ) {
                        Logger::warn( "GeoModel", " Vertex ", i,
                            " is at a Corner but in no Surface " );
                        valid_vertex = false;
                    }
                }
            }
            valid[i] = valid_vertex;
        }
        index_t nb_invalid = static_cast< index_t >( std::count( valid.begin(),
            valid.end(), false ) );

        if( nb_invalid > 0 ) {
            std::ostringstream file;
            file << validity_errors_directory << "/invalid_global_vertices.geogram";
            save_invalid_points( file, geomodel, valid );

            if( GEO::CmdLine::get_arg_bool( "validity_save" ) ) {
                Logger::warn( "GeoModel", nb_invalid, " invalid vertices" );
                Logger::warn( "GeoModel", "Saved in file: ", file.str() );
            }

            return false;
        } else {
            return true;
        }
    }

    void save_edges(
        const std::ostringstream& file,
        const GeoModel& geomodel,
        const std::vector< index_t >& e )
    {
        GEO::Mesh edge_mesh;
        index_t previous_vertex_id = NO_ID;
        for( index_t i = 0; i < e.size(); ++i ) {
            index_t cur_vertex_id = edge_mesh.vertices.create_vertex(
                geomodel.mesh.vertices.vertex( e[i] ).data() );
            if( i % 2 == 0 ) {
                ringmesh_assert( previous_vertex_id == NO_ID );
                previous_vertex_id = cur_vertex_id;
                continue;
            }
            ringmesh_assert( previous_vertex_id != NO_ID );
            edge_mesh.edges.create_edge( previous_vertex_id, cur_vertex_id );
            previous_vertex_id = NO_ID;
        }
        save_mesh_locating_geomodel_inconsistencies( edge_mesh, file );
    }

    void save_polygons(
        const std::string& file,
        const Surface& surface,
        const std::vector< index_t >& polygons )
    {
        GEO::Mesh mesh;
        for( index_t cur_polygon : polygons ) {
            index_t nb_vertices_in_polygon = surface.nb_mesh_element_vertices(
                cur_polygon );
            GEO::vector< index_t > vertices;
            vertices.reserve( nb_vertices_in_polygon );
            for( index_t v = 0; v < nb_vertices_in_polygon; v++ ) {
                index_t new_vertex = mesh.vertices.create_vertex(
                    surface.mesh_element_vertex( cur_polygon, v ).data() );
                vertices.push_back( new_vertex );
            }
            mesh.facets.create_polygon( vertices );
        }
        GEO::mesh_save( mesh, file );
    }

    /*!
     * @brief Check boundary of a surface
     * @details All the edges on the boundary of a surface must be in a Line
     *          of the associated geomodel
     *          The Line boundaries must form a closed manifold line.
     */
    bool surface_boundary_valid( const Surface& surface )
    {
        const GeoModelMeshVertices& geomodel_vertices =
            surface.geomodel().mesh.vertices;
        std::vector< index_t > invalid_corners;
        gmme_id S_id = surface.gmme();
        for( index_t p = 0; p < surface.nb_mesh_elements(); ++p ) {
            for( index_t v = 0; v < surface.nb_mesh_element_vertices( p ); ++v ) {
                if( surface.polygon_adjacent_index( p, v ) == NO_ID
                    && !is_edge_on_line( surface.geomodel(),
                        geomodel_vertices.geomodel_vertex_id( S_id, p, v ),
                        geomodel_vertices.geomodel_vertex_id( S_id, p,
                            surface.next_polygon_vertex_index( p, v ) ) ) ) {
                    invalid_corners.push_back(
                        geomodel_vertices.geomodel_vertex_id( S_id, p, v ) );
                    invalid_corners.push_back(
                        geomodel_vertices.geomodel_vertex_id( S_id, p,
                            surface.next_polygon_vertex_index( p, v ) ) );
                }
            }
        }
        if( !invalid_corners.empty() ) {
            std::ostringstream file;
            file << validity_errors_directory << "/invalid_boundary_surface_"
                << surface.index() << ".geogram";
            save_edges( file, surface.geomodel(), invalid_corners );

            if( GEO::CmdLine::get_arg_bool( "validity_save" ) ) {
                Logger::warn( "GeoModel", " Invalid surface boundary: ",
                    invalid_corners.size() / 2, " boundary edges of ", S_id,
                    "  are in no line of the geomodel " );
                Logger::warn( "GeoModel", " Saved in file: ", file.str() );
            }
            return false;
        } else {
            return true;
        }
    }

    /*!
     * @brief Save in a .lin file the
     */
    void debug_save_non_manifold_edges(
        const GeoModel& geomodel,
        const std::vector< index_t >& edge_indices,
        const std::vector< index_t >& non_manifold_edges )
    {
        GeogramLineMesh mesh;
        GeogramLineMeshBuilder builder;
        builder.set_mesh( mesh );
        index_t nb_edges = static_cast< index_t >( non_manifold_edges.size() );
        builder.create_vertices( 2 * nb_edges );
        builder.create_edges( nb_edges );
        const GeoModelMeshVertices& vertices = geomodel.mesh.vertices;
        for( index_t e = 0; e < non_manifold_edges.size(); e++ ) {
            index_t edge_id = non_manifold_edges[e];
            const vec3& v0 = vertices.vertex( edge_indices[edge_id] );
            const vec3& v1 = vertices.vertex( edge_indices[edge_id + 1] );
            builder.set_vertex( 2 * e, v0 );
            builder.set_vertex( 2 * e + 1, v1 );
            builder.set_edge_vertex( e, 0, 2 * e );
            builder.set_edge_vertex( e, 1, 2 * e + 1 );
        }
        mesh.save_mesh( validity_errors_directory + "/non_manifold_edges.geogram" );
    }

    bool is_surface_conformal_to_volume(
        const Surface& surface,
        const NNSearch& cell_facet_barycenter_nn_search )
    {
        std::vector< index_t > unconformal_polygons;
        for( index_t p = 0; p < surface.nb_mesh_elements(); p++ ) {
            vec3 center = surface.mesh_element_barycenter( p );
            std::vector< index_t > result =
                cell_facet_barycenter_nn_search.get_neighbors( center,
                    surface.geomodel().epsilon() );
            if( result.empty() ) {
                unconformal_polygons.push_back( p );
            }
        }
        if( !unconformal_polygons.empty() ) {
            std::ostringstream file;
            file << validity_errors_directory << "/unconformal_surface_"
                << surface.index() << ".geogram";
            save_polygons( file.str(), surface, unconformal_polygons );

            if( GEO::CmdLine::get_arg_bool( "validity_save" ) ) {
                Logger::warn( "GeoModel", " Unconformal surface: ",
                    unconformal_polygons.size(), " polygons of ", surface.gmme(),
                    " are unconformal with the geomodel cells " );
                Logger::warn( "GeoModel", " Saved in file: ", file.str() );
            }

            return false;
        } else {
            return true;
        }
    }

    void compute_border_edges(
        const GeoModel& geomodel,
        std::vector< index_t >& edge_indices )
    {
        const GeoModelMeshPolygons& polygons = geomodel.mesh.polygons;
        for( index_t s = 0; s < geomodel.nb_surfaces(); s++ ) {
            for( index_t p = 0; p < polygons.nb_polygons( s ); p++ ) {
                index_t polygon_id = polygons.polygon( s, p );
                for( index_t v = 0; v < polygons.nb_vertices( polygon_id ); v++ ) {
                    index_t adj = polygons.adjacent( polygon_id, v );
                    if( adj == NO_ID ) {
                        edge_indices.push_back( polygons.vertex( polygon_id, v ) );
                        index_t next_v = ( v + 1 )
                            % polygons.nb_vertices( polygon_id );
                        edge_indices.push_back(
                            polygons.vertex( polygon_id, next_v ) );
                    }
                }
            }
        }
    }

    void compute_border_edge_barycenters(
        const GeoModel& geomodel,
        const std::vector< index_t >& edge_indices,
        std::vector< vec3 >& edge_barycenters )
    {
        const GeoModelMeshVertices& vertices = geomodel.mesh.vertices;
        index_t nb_edges = static_cast< index_t >( edge_indices.size() / 2 );
        edge_barycenters.reserve( nb_edges );
        for( index_t e = 0; e < edge_indices.size(); e += 2 ) {
            const vec3& v0 = vertices.vertex( edge_indices[e] );
            const vec3& v1 = vertices.vertex( edge_indices[e + 1] );
            edge_barycenters.push_back( ( v0 + v1 ) * 0.5 );
        }
    }

    void compute_edge_on_lines(
        const GeoModel& geomodel,
        const std::vector< vec3 >& edge_barycenters,
        std::vector< bool >& edge_on_lines )
    {
        edge_on_lines.resize( edge_barycenters.size(), false );
        NNSearch nn( edge_barycenters );
        for( index_t l = 0; l < geomodel.nb_lines(); l++ ) {
            const Line& line = geomodel.line( l );
            for( index_t e = 0; e < line.nb_mesh_elements(); e++ ) {
                const vec3 query = line.mesh_element_barycenter( e );
                std::vector< index_t > results = nn.get_neighbors( query,
                    geomodel.epsilon() );
                for( index_t edge : results ) {
                    edge_on_lines[edge] = true;
                }
            }
        }
    }

    void compute_non_manifold_edges(
        const std::vector< bool >& edge_on_lines,
        std::vector< index_t >& non_manifold_edges )
    {
        for( index_t e = 0; e < edge_on_lines.size(); e++ ) {
            if( !edge_on_lines[e] ) {
                non_manifold_edges.push_back( e );
            }
        }
    }

    /*!
     * @brief Implementation class for validity checks on a GeoModel
     */
    class GeoModelValidityCheck {
    public:
        GeoModelValidityCheck(
            const GeoModel& geomodel,
            const ValidityCheckMode validity_check_mode )
            : geomodel_( geomodel ), valid_( true ), mode_( validity_check_mode )
        {
            if( mode_ == ValidityCheckMode::UNDEFINED ) {
                // If not defined, reset the check mode to the largest possible
                mode_ = ValidityCheckMode::ALL;
            }
            // Ensure that the geomodel vertices are computed and up-to-date
            // Without that we cannot do anything        
            geomodel_.mesh.vertices.test_and_initialize();
            geomodel_.mesh.polygons.test_and_initialize();
        }

        /*!
         * @brief Run the geomodel validity check in accordance to the defined mode.
         */
        bool is_geomodel_valid()
        {
            do_check_validity( mode_ );
            return valid_;
        }

    private:
        void do_check_validity( ValidityCheckMode mode )
        {
            std::vector< std::thread > threads;
            threads.reserve( 8 );
            if( mode == ValidityCheckMode::GEOMETRY
                || mode == ValidityCheckMode::ALL ) {
                threads.emplace_back(
                    &GeoModelValidityCheck::test_polygon_intersections, this );
            }
            if( mode != ValidityCheckMode::TOPOLOGY ) {
                // Add geometrical validity check
                threads.emplace_back(
                    &GeoModelValidityCheck::test_geomodel_mesh_entities_validity,
                    this );
                threads.emplace_back(
                    &GeoModelValidityCheck::test_region_surface_mesh_conformity,
                    this );
                threads.emplace_back(
                    &GeoModelValidityCheck::test_non_manifold_edges, this );

            }
            if( mode != ValidityCheckMode::GEOMETRY ) {
                // Add topological validity check
                threads.emplace_back(
                    &GeoModelValidityCheck::test_geomodel_connectivity_validity,
                    this );
                threads.emplace_back( &GeoModelValidityCheck::test_finite_extension,
                    this );
                threads.emplace_back(
                    &GeoModelValidityCheck::test_surface_line_mesh_conformity,
                    this );
                /// TODO: find a way to add this test for Model3d. See BC.
//                threads.emplace_back(
//                    &GeoModelValidityCheck::test_non_free_line_at_two_interfaces_intersection,
//                    this );
            }

            // Geological validity must always be checked
            threads.emplace_back(
                &GeoModelValidityCheck::test_geomodel_geological_validity, this );

            for( index_t i = 0; i < threads.size(); i++ ) {
                threads[i].join();
            }
        }

        /*! 
         * @brief Verify the validity of all GeoModelMeshEntities
         */
        void test_geomodel_mesh_entities_validity()
        {
            if( !are_geomodel_mesh_entities_mesh_valid( geomodel_ ) ) {
                set_invalid_model();
            }
        }

        /*!
         * @brief Verify the validity of all GeoModelEntities
         */
        void test_geomodel_connectivity_validity()
        {
            if( !are_geomodel_mesh_entities_connectivity_valid( geomodel_ ) ) {
                set_invalid_model();
            }
        }

        /*!
         * @brief Verify the validity of all GeoModelGeologicalEntities
         */
        void test_geomodel_geological_validity()
        {
            if( !are_geomodel_geological_entities_valid( geomodel_ ) ) {
                set_invalid_model();
            }
            if( !are_geomodel_mesh_entities_parent_valid( geomodel_ ) ) {
                set_invalid_model();
            }
        }
        /*!
         * @brief Check that the geomodel has a finite extension
         * @details The boundary of the universe region is a one connected component
         * manifold closed surface.
         */
        void test_finite_extension()
        {
            std::vector< index_t > voi_surfaces;
            std::tie( voi_surfaces, std::ignore ) = geomodel_.get_voi_surfaces();

            index_t nb_points_all_voi_surfaces = 0;
            for( index_t voi_surface_id : voi_surfaces ) {
                const Surface& cur_surface = geomodel_.surface( voi_surface_id );
                nb_points_all_voi_surfaces += cur_surface.nb_vertices();
            }

            std::vector< vec3 > all_points;
            all_points.reserve( nb_points_all_voi_surfaces );

            for( index_t voi_surface_id : voi_surfaces ) {
                const Surface& cur_surface = geomodel_.surface( voi_surface_id );
                for( index_t v_i = 0; v_i < cur_surface.nb_vertices(); ++v_i ) {
                    all_points.push_back( cur_surface.vertex( v_i ) );
                }
            }

            NNSearch make_unique_surf( all_points );
            std::vector< index_t > unique_id;
            std::vector< vec3 > facet_points;
            make_unique_surf.get_colocated_index_mapping( geomodel_.epsilon(),
                unique_id, facet_points );

            index_t offset_vertices = 0;
            std::vector< index_t > facet_indices;
            std::vector< index_t > facet_ptr;
            index_t count_facet_vertices = 0;
            facet_ptr.push_back( count_facet_vertices );
            for( index_t voi_surface_id : voi_surfaces ) {

                const Surface& cur_surf = geomodel_.surface( voi_surface_id );

                for( index_t facet_itr = 0; facet_itr < cur_surf.nb_mesh_elements();
                    ++facet_itr ) {
                    for( index_t point_i = 0;
                        point_i < cur_surf.nb_mesh_element_vertices( facet_itr );
                        ++point_i ) {

                        index_t index = cur_surf.mesh_element_vertex_index(
                            facet_itr, point_i );
                        facet_indices.push_back(
                            unique_id[index + offset_vertices] );

                    }
                    count_facet_vertices += cur_surf.nb_mesh_element_vertices(
                        facet_itr );
                    facet_ptr.push_back( count_facet_vertices );
                }

                offset_vertices += cur_surf.nb_vertices();
            }

            std::unique_ptr< SurfaceMesh > surface = SurfaceMesh::create_mesh();
            ringmesh_assert( surface != nullptr );
            std::unique_ptr< SurfaceMeshBuilder > builder =
                SurfaceMeshBuilder::create_builder( *surface );
            ringmesh_assert( builder != nullptr );
            index_t start = builder->create_vertices(
                static_cast< index_t >( facet_points.size() ) );
            ringmesh_unused( start );
            ringmesh_assert( start == 0 );
            for( index_t v_i = 0; v_i < facet_points.size(); ++v_i ) {
                builder->set_vertex( v_i, facet_points[v_i] );
            }
            builder->create_polygons( facet_indices, facet_ptr );

            for( index_t p = 0; p < surface->nb_polygons(); p++ ) {
                for( index_t v = 0; v < surface->nb_polygon_vertices( p );
                    v++ ) {
                    builder->set_polygon_adjacent( p, v, NO_ID );
                }
            }

            builder->connect_polygons();
            index_t nb_connected_components = NO_ID;
            std::tie( nb_connected_components, std::ignore ) =
                surface->get_connected_components();

            DEBUG(voi_surfaces.size());
            DEBUG(nb_connected_components);
            surface->save_mesh("toto.geogram");
            if( nb_connected_components != 1 ) {
                set_invalid_model();
            }
        }

        /*!
         * Check geometrical-connectivity consistency
         * @todo Check that all Line segments correspond to a Surface
         *  edge that is on the boundary.
         */
        void test_surface_line_mesh_conformity()
        {
            // Check relationships between GeoModelEntities
            // sharing the same point of the geomodel
            if( !check_model_points_validity( geomodel_ ) ) {
                set_invalid_model();
            }
            // Check on that Surface edges are in a Line
            for( index_t i = 0; i < geomodel_.nb_surfaces(); ++i ) {
                if( !surface_boundary_valid( geomodel_.surface( i ) ) ) {
                    set_invalid_model();
                }
            }
        }

        void test_region_surface_mesh_conformity()
        {
            if( geomodel_.mesh.cells.nb() > 0 ) {
                // Check the consistency between Surface polygons and Region cell facets
                const NNSearch& nn_search =
                    geomodel_.mesh.cells.cell_facet_nn_search();
                for( index_t i = 0; i < geomodel_.nb_surfaces(); ++i ) {
                    if( !is_surface_conformal_to_volume( geomodel_.surface( i ),
                        nn_search ) ) {
                        set_invalid_model();
                    }
                }
            }
        }

        void test_non_free_line_at_two_interfaces_intersection()
        {
            if( !geomodel_.entity_type_manager().geological_entity_manager.is_valid_type(
                Interface::type_name_static() ) ) {
                return;
            }
            for( index_t line_i = 0; line_i < geomodel_.nb_lines(); ++line_i ) {
                const Line& cur_line = geomodel_.line( line_i );
                if( cur_line.nb_incident_entities() == 1 ) {
                    continue;
                }

                const index_t first_interface_id =
                    cur_line.incident_entity( 0 ).parent_gmge(
                        Interface::type_name_static() ).index();
                ringmesh_assert( first_interface_id != NO_ID );
                bool at_least_two_different_interfaces = false;
                for( index_t in_boundary_i = 1;
                    in_boundary_i < cur_line.nb_incident_entities(); ++in_boundary_i ) {
                    const index_t cur_interface_id =
                        cur_line.incident_entity( in_boundary_i ).parent_gmge(
                            Interface::type_name_static() ).index();
                    ringmesh_assert( cur_interface_id != NO_ID );
                    if( cur_interface_id != first_interface_id ) {
                        at_least_two_different_interfaces = true;
                        break;
                    }
                }

                if( !at_least_two_different_interfaces ) {
                    Logger::warn( "GeoModel",
                        "All in boundaries (surfaces) of line ", line_i,
                        " are children of a same interface." );
                    set_invalid_model();
                }
            }
        }

        /*!
         * @brief Returns true if there are non-manifold edges that are
         *        not in any Line of the geomodel
         * @note Connect the polygons of the global mesh
         * @note This is a quite expensive test.
         */
        void test_non_manifold_edges()
        {
            std::vector< index_t > edge_indices;
            compute_border_edges( geomodel_, edge_indices );
            std::vector< vec3 > edge_barycenters;
            compute_border_edge_barycenters( geomodel_, edge_indices,
                edge_barycenters );
            std::vector< bool > edge_on_lines;
            compute_edge_on_lines( geomodel_, edge_barycenters, edge_on_lines );
            std::vector< index_t > non_manifold_edges;
            compute_non_manifold_edges( edge_on_lines, non_manifold_edges );

            if( !non_manifold_edges.empty() ) {
                Logger::warn( "GeoModel", non_manifold_edges.size(),
                    " non-manifold edges " );
                debug_save_non_manifold_edges( geomodel_, edge_indices,
                    non_manifold_edges );

                set_invalid_model();
            }
        }

        /*!
         * @brief Returns true if there are intersections between polygons
         * @details Operates on the global mesh
         * @note This is a very expensive test.
         */
        void test_polygon_intersections()
        {
            if( geomodel_.mesh.polygons.nb()
                == geomodel_.mesh.polygons.nb_triangle()
                    + geomodel_.mesh.polygons.nb_quad() ) {
                std::vector< bool > has_intersection;
                StoreIntersections action( geomodel_, has_intersection );
                const SurfaceAABBTree& AABB = geomodel_.mesh.polygons.aabb();
                AABB.compute_self_element_bbox_intersections( action );

                index_t nb_intersections = static_cast< index_t >( std::count(
                    has_intersection.begin(), has_intersection.end(), 1 ) );

                if( nb_intersections > 0 ) {
                    GEO::Mesh mesh;
                    for( index_t p = 0; p < has_intersection.size(); p++ ) {
                        if( !has_intersection[p] ) continue;
                        GEO::vector< index_t > vertices;
                        vertices.reserve( geomodel_.mesh.polygons.nb_vertices( p ) );
                        for( index_t v = 0;
                            v < geomodel_.mesh.polygons.nb_vertices( p ); v++ ) {
                            index_t id =
                                mesh.vertices.create_vertex(
                                    geomodel_.mesh.vertices.vertex(
                                        geomodel_.mesh.polygons.vertex( p, v ) ).data() );
                            vertices.push_back( id );
                        }
                        mesh.facets.create_polygon( vertices );
                    }
                    std::ostringstream file;
                    file << validity_errors_directory
                        << "/intersected_polygons.geogram";
                    save_mesh_locating_geomodel_inconsistencies( mesh, file );
                    Logger::out( "I/O" );
                    Logger::warn( "GeoModel", nb_intersections,
                        " polygon intersections " );
                    set_invalid_model();
                }
            } else {
                Logger::warn( "GeoModel",
                    "Polygonal intersection check not implemented yet" );
            }
        }

        void set_invalid_model()
        {
            valid_ = false;
        }

    private:
        const GeoModel& geomodel_;
        bool valid_;
        ValidityCheckMode mode_;
    };

}
// anonymous namespace

namespace RINGMesh {

    void set_validity_errors_directory( const std::string& directory )
    {
        // If trailing / or \ is not removed, the test fails on Windows
        std::string copy( directory );
        if( *copy.rbegin() == '/' || *copy.rbegin() == '\\' ) {
            copy.erase( copy.end() - 1 );
        }
        if( GEO::FileSystem::is_directory( copy ) ) {
            validity_errors_directory = copy + '/';
        }
    }

    bool are_geomodel_mesh_entities_mesh_valid( const GeoModel& geomodel )
    {
        const std::vector< MeshEntityType >& meshed_types =
            MeshEntityTypeManager::mesh_entity_types();
        index_t count_invalid = 0;
        for( const MeshEntityType& type : meshed_types ) {
            index_t nb_entities = geomodel.nb_mesh_entities( type );
            for( index_t i = 0; i < nb_entities; ++i ) {
                const GeoModelMeshEntity& E = geomodel.mesh_entity( type, i );
                if( !E.is_valid() ) {
                    count_invalid++;
                }
            }
        }
        if( count_invalid != 0 ) {
            Logger::warn( "GeoModel", count_invalid,
                " mesh entities of the geomodel have an invalid mesh." );
        }
        return count_invalid == 0;
    }

    bool are_geomodel_mesh_entities_connectivity_valid( const GeoModel& geomodel )
    {
        const std::vector< MeshEntityType >& meshed_types =
            MeshEntityTypeManager::mesh_entity_types();
        index_t count_invalid = 0;
        for( const MeshEntityType& type : meshed_types ) {
            index_t nb_entities = geomodel.nb_mesh_entities( type );
            for( index_t i = 0; i < nb_entities; ++i ) {
                const GeoModelMeshEntity& E = geomodel.mesh_entity( type, i );
                if( !E.is_connectivity_valid() ) {
                    count_invalid++;
                }
            }
        }
        if( count_invalid != 0 ) {
            Logger::warn( "GeoModel", count_invalid,
                " mesh entities of the geomodel have an invalid connectivity." );
        }
        return count_invalid == 0;
    }

    bool are_geomodel_geological_entities_valid( const GeoModel& geomodel )
    {
        const std::vector< GeologicalEntityType >& geological_types =
            geomodel.entity_type_manager().geological_entity_manager.geological_entity_types();
        index_t count_invalid = 0;
        for( const GeologicalEntityType& type : geological_types ) {
            index_t nb_entities = geomodel.nb_geological_entities( type );
            for( index_t i = 0; i < nb_entities; ++i ) {
                const GeoModelGeologicalEntity& E = geomodel.geological_entity( type,
                    i );
                if( !E.is_valid() ) {
                    count_invalid++;
                }
            }
        }
        if( count_invalid != 0 ) {
            Logger::warn( "GeoModel", count_invalid,
                " geological entities of the geomodel are invalid " );
        }
        return count_invalid == 0;
    }

    bool are_geomodel_mesh_entities_parent_valid( const GeoModel& geomodel )
    {
        const std::vector< MeshEntityType >& meshed_types =
            MeshEntityTypeManager::mesh_entity_types();
        index_t count_invalid = 0;
        for( const MeshEntityType& type : meshed_types ) {
            index_t nb_entities = geomodel.nb_mesh_entities( type );
            for( index_t i = 0; i < nb_entities; ++i ) {
                const GeoModelMeshEntity& E = geomodel.mesh_entity( type, i );
                if( !E.is_parent_connectivity_valid() ) {
                    count_invalid++;
                }
            }
        }
        if( count_invalid != 0 ) {
            Logger::warn( "GeoModel", count_invalid,
                " mesh entities of the geomodel have an invalid ",
                "parent connectivity (geological relationships)." );
        }
        return count_invalid == 0;
    }

    bool is_geomodel_valid(
        const GeoModel& geomodel,
        ValidityCheckMode validity_check_mode )
    {
        if( validity_check_mode == ValidityCheckMode::GEOMETRY
            && !GEO::CmdLine::get_arg_bool( "in:intersection_check" ) ) {
            validity_check_mode =
                ValidityCheckMode::GEOMETRY_EXCEPT_FACET_INTERSECTION;
        } else if( validity_check_mode == ValidityCheckMode::ALL
            && !GEO::CmdLine::get_arg_bool( "in:intersection_check" ) ) {
            validity_check_mode = ValidityCheckMode::ALL_EXCEPT_FACET_INTERSECTION;
        }

        GeoModelValidityCheck validity_checker( geomodel, validity_check_mode );

        bool valid = validity_checker.is_geomodel_valid();

        if( valid ) {
            Logger::out( "GeoModel", "Model ", geomodel.name(), " is valid " );
        } else {
            Logger::warn( "GeoModel", "Model ", geomodel.name(), " is invalid " );
            if( !GEO::CmdLine::get_arg_bool( "validity_save" ) ) {
                Logger::out( "Info", "To save geomodel invalidities in files ",
                    "(.geogram) set \"validity_save\" to true in the command line." );
            }
        }
        return valid;
    }

} // namespace RINGMesh
