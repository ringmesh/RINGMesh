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

    template< index_t DIMENSION >
    bool triangles_intersect(
        const GeoModel< DIMENSION >& geomodel,
        const GeoModelMeshPolygons< DIMENSION >& polygons,
        index_t triangle1,
        index_t triangle2 )
    {
        ringmesh_assert( polygons.nb_vertices( triangle1 ) == 3 );
        ringmesh_assert( polygons.nb_vertices( triangle2 ) == 3 );
        const GeoModelMeshVertices< DIMENSION >& vertices = geomodel.mesh.vertices;
        const vecn< DIMENSION >& p1 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( triangle1, 0 ) ) );
        const vecn< DIMENSION >& p2 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( triangle1, 1 ) ) );
        const vecn< DIMENSION >& p3 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( triangle1, 2 ) ) );

        const vecn< DIMENSION >& q1 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( triangle2, 0 ) ) );
        const vecn< DIMENSION >& q2 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( triangle2, 1 ) ) );
        const vecn< DIMENSION >& q3 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( triangle2, 2 ) ) );
        GEO::vector< GEO::TriangleIsect > sym;
        return triangles_intersections( p1, p2, p3, q1, q2, q3, sym );
    }

    template< index_t DIMENSION >
    bool triangle_quad_intersect(
        const GeoModel< DIMENSION >& geomodel,
        const GeoModelMeshPolygons< DIMENSION >& polygons,
        index_t triangle,
        index_t quad )
    {
        ringmesh_assert( polygons.nb_vertices( triangle ) == 3 );
        ringmesh_assert( polygons.nb_vertices( quad ) == 4 );
        const GeoModelMeshVertices< DIMENSION >& vertices = geomodel.mesh.vertices;
        const vecn< DIMENSION >& p1 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( triangle, 0 ) ) );
        const vecn< DIMENSION >& p2 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( triangle, 1 ) ) );
        const vecn< DIMENSION >& p3 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( triangle, 2 ) ) );

        const vecn< DIMENSION >& q1 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad, 0 ) ) );
        const vecn< DIMENSION >& q2 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad, 1 ) ) );
        const vecn< DIMENSION >& q3 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad, 2 ) ) );
        const vecn< DIMENSION >& q4 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad, 3 ) ) );
        GEO::vector< GEO::TriangleIsect > sym;
        if( triangles_intersections( p1, p2, p3, q1, q2, q3, sym ) ) {
            return true;
        }
        if( triangles_intersections( p1, p2, p3, q1, q3, q4, sym ) ) {
            return true;
        }
        return false;
    }

    template< index_t DIMENSION >
    bool quad_quad_intersect(
        const GeoModel< DIMENSION >& geomodel,
        const GeoModelMeshPolygons< DIMENSION >& polygons,
        index_t quad1,
        index_t quad2 )
    {
        ringmesh_assert( polygons.nb_vertices( quad1 ) == 4 );
        ringmesh_assert( polygons.nb_vertices( quad2 ) == 4 );
        const GeoModelMeshVertices< DIMENSION >& vertices = geomodel.mesh.vertices;
        const vecn< DIMENSION >& p1 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad1, 0 ) ) );
        const vecn< DIMENSION >& p2 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad1, 1 ) ) );
        const vecn< DIMENSION >& p3 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad1, 2 ) ) );
        const vecn< DIMENSION >& p4 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad1, 3 ) ) );

        const vecn< DIMENSION >& q1 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad2, 0 ) ) );
        const vecn< DIMENSION >& q2 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad2, 1 ) ) );
        const vecn< DIMENSION >& q3 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad2, 2 ) ) );
        const vecn< DIMENSION >& q4 = vertices.vertex(
            polygons.vertex( ElementLocalVertex( quad2, 3 ) ) );
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

    template< index_t DIMENSION >
    bool is_edge_on_line( const Line< DIMENSION >& line, index_t v0, index_t v1 )
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
    template< index_t DIMENSION >
    bool is_edge_on_line(
        const GeoModel< DIMENSION >& geomodel,
        index_t v0,
        index_t v1 )
    {
        MeshEntityType line_type = Line< DIMENSION >::type_name_static();
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
    template< index_t DIMENSION >
    bool polygons_share_line_edge(
        const GeoModel< DIMENSION >& geomodel,
        const GeoModelMeshPolygons< DIMENSION >& polygons,
        index_t p1,
        index_t p2 )
    {
        // Only test the edges on boundary
        for( index_t v1 : range( polygons.nb_vertices( p1 ) ) ) {
            if( polygons.adjacent( PolygonLocalEdge( p1, v1 ) ) != NO_ID ) {
                continue;
            }
            index_t v10 = polygons.vertex( ElementLocalVertex( p1, v1 ) );
            index_t v11 = polygons.vertex(
                ElementLocalVertex( p1, ( v1 + 1 ) % polygons.nb_vertices( p1 ) ) );
            for( index_t v2 : range( polygons.nb_vertices( p2 ) ) ) {
                if( polygons.adjacent( PolygonLocalEdge( p2, v2 ) ) != NO_ID ) {
                    continue;
                }
                index_t v20 = polygons.vertex( ElementLocalVertex( p2, v2 ) );
                index_t v21 = polygons.vertex(
                    ElementLocalVertex( p2,
                        ( v2 + 1 ) % polygons.nb_vertices( p2 ) ) );

                if( ( v10 == v20 && v11 == v21 ) || ( v10 == v21 && v11 == v20 ) ) {
                    if( is_edge_on_line( geomodel, v20, v21 ) ) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    template< index_t DIMENSION >
    bool polygons_are_adjacent(
        const GeoModelMeshPolygons< DIMENSION >& polygons,
        index_t p1,
        index_t p2 )
    {
        if( p1 == p2 ) {
            return true;
        }
        for( index_t v : range( polygons.nb_vertices( p1 ) ) ) {
            if( polygons.adjacent( PolygonLocalEdge( p1, v ) ) == p2 ) {
                return true;
            }
        }
        return false;
    }

    /*!
     * @brief Action class for storing intersections when traversing
     *  a AABBTree.
     */
    template< index_t DIMENSION >
    class StoreIntersections {
    public:
        /*!
         * @brief Constructs the StoreIntersections
         * @param[in] geomodel the geomodel
         * @param[out] has_isect the flag that indicates for each polygon
         *  whether it has intersections
         */
        StoreIntersections(
            const GeoModel< DIMENSION >& geomodel,
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
            PolygonType type;
            std::tie( type, std::ignore ) = polygons_.type( p );
            return type == PolygonType::TRIANGLE;
        }
        bool is_quad( index_t p ) const
        {
            PolygonType type;
            std::tie( type, std::ignore ) = polygons_.type( p );
            return type == PolygonType::QUAD;
        }

    private:
        const GeoModel< DIMENSION >& geomodel_;
        const GeoModelMeshPolygons< DIMENSION >& polygons_;
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
    template< index_t DIMENSION >
    bool is_boundary_entity(
        const GeoModel< DIMENSION >& geomodel,
        const gmme_id& entity,
        const gmme_id& boundary )
    {
        const GeoModelMeshEntity< DIMENSION >& E = geomodel.mesh_entity( entity );
        for( index_t i : range( E.nb_boundaries() ) ) {
            if( E.boundary_gmme( i ) == boundary ) {
                return true;
            }
        }
        return false;
    }

    template< index_t DIMENSION >
    void save_invalid_points(
        const std::ostringstream& file,
        const GeoModel< DIMENSION >& geomodel,
        const std::vector< bool >& valid )
    {
        GEO::Mesh point_mesh;
        for( index_t i : range( valid.size() ) ) {
            if( !valid[i] ) {
                const vecn< DIMENSION >& V = geomodel.mesh.vertices.vertex( i );
                point_mesh.vertices.create_vertex( V.data() );
            }
        }
        save_mesh_locating_geomodel_inconsistencies( point_mesh, file );
    }

    template< index_t DIMENSION >
    std::map< MeshEntityType, std::vector< index_t > > get_entities(
        const GeoModel< DIMENSION >& geomodel,
        index_t i )
    {
        std::map< MeshEntityType, std::vector< index_t > > entities;
        const std::vector< MeshEntityType >& types =
            geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();
        for( const MeshEntityType& type : types ) {
            entities[type];
        }

        const std::vector< GMEVertex >& bmes = geomodel.mesh.vertices.gme_vertices(
            i );
        for( const GMEVertex& vertex : bmes ) {
            const MeshEntityType& T = vertex.gmme.type();
            index_t id = vertex.gmme.index();
            entities[T].push_back( id );
        }
        return entities;
    }

    void print_error(
        const std::vector< index_t >& entities,
        const std::string& entity_name )
    {
        std::ostringstream oss;
        oss << " Vertex is in " << entities.size() << " " << entity_name << ": ";
        for( index_t entity : entities ) {
            oss << entity << " ; ";
        }
        Logger::warn( "GeoModel", oss.str() );
    }

    template< template< index_t > class ENTITY, index_t DIMENSION >
    bool is_vertex_valid(
        const GeoModel< DIMENSION >& geomodel,
        const std::map< MeshEntityType, std::vector< index_t > >& entities )
    {
        MeshEntityType type = ENTITY< DIMENSION >::type_name_static();
        MeshEntityType boundary_type =
            geomodel.entity_type_manager().mesh_entity_manager.boundary_entity_type(
                type );
        const std::vector< index_t >& type_entities = entities.find( type )->second;
        if( entities.find( boundary_type )->second.empty() ) {
            if( !type_entities.empty() ) {
                if( type_entities.size() != 1 ) {
                    print_error( type_entities, type.to_string() + "s" );
                    Logger::warn( "GeoModel", "It should be in only one ",
                        boundary_type );
                    return false;
                } else {
                    return true;
                }
            } else {
                return true;
            }

        } else {
            if( type_entities.empty() ) {
                Logger::warn( "GeoModel", " Vertex is in a ", boundary_type,
                    " but in no ", type );
                return false;
            } else {
                const std::vector< index_t >& boundary_entities = entities.find(
                    boundary_type )->second;
                // Check that one point is no more than twice in a SURFACE
                for( index_t entity : type_entities ) {
                    index_t nb = static_cast< index_t >( std::count(
                        type_entities.begin(), type_entities.end(), entity ) );
                    if( nb > 2 ) {
                        Logger::warn( "GeoModel", " Vertex is ", nb, " times in ",
                            geomodel.mesh_entity( type, entity ).gmme() );
                        return false;
                    } else if( nb == 2 ) {
                        // If a point is twice in a SURFACE, it must be
                        // on an internal boundary Line.
                        bool internal_boundary = false;
                        for( index_t line : boundary_entities ) {
                            if( geomodel.mesh_entity( boundary_type, line ).is_inside_border(
                                geomodel.mesh_entity( type, entity ) ) ) {
                                internal_boundary = true;
                                break;
                            }
                        }
                        if( !internal_boundary ) {
                            Logger::warn( "GeoModel", " Vertex appears ", nb,
                                " times in ",
                                geomodel.mesh_entity( type, entity ).gmme() );
                            return false;
                        }
                    }
                }
            }
            return true;
        }
    }

    template< index_t DIMENSION >
    bool is_region_vertex_valid(
        const GeoModel< DIMENSION >& geomodel,
        const std::map< MeshEntityType, std::vector< index_t > >& entities )
    {
        if( geomodel.nb_regions() > 0 && geomodel.region( 0 ).is_meshed() ) {
            return is_vertex_valid< Region >( geomodel, entities );
        }
        return true;
    }

    template< index_t DIMENSION >
    bool is_surface_vertex_valid(
        const GeoModel< DIMENSION >& geomodel,
        const std::map< MeshEntityType, std::vector< index_t > >& entities )
    {
        return is_vertex_valid< Surface >( geomodel, entities );
    }

    template< index_t DIMENSION >
    bool is_line_vertex_valid(
        const GeoModel< DIMENSION >& geomodel,
        const std::map< MeshEntityType, std::vector< index_t > >& entities )
    {
        const std::vector< index_t >& lines = entities.find(
            Line< DIMENSION >::type_name_static() )->second;
        if( entities.find( Corner< DIMENSION >::type_name_static() )->second.empty() ) {
            if( !lines.empty() ) {
                if( lines.size() != 1 ) {
                    print_error( lines, "Lines" );
                    Logger::warn( "GeoModel", "It should be in only one Line" );
                    return false;
                } else {
                    return true;
                }
            } else {
                return true;
            }
        } else {
            if( lines.size() < 2 ) {
                print_error( lines, "Lines" );
                Logger::warn( "GeoModel", "It should be in at least 2 Lines" );
                return false;
            } else {
                for( index_t line : lines ) {
                    index_t nb = static_cast< index_t >( std::count( lines.begin(),
                        lines.end(), line ) );
                    if( nb == 2 ) {
                        if( !geomodel.line( line ).is_closed() ) {
                            Logger::warn( "GeoModel", " Vertex"
                                " is twice in Line ", line );
                            return false;
                        }
                    } else if( nb > 2 ) {
                        Logger::warn( "GeoModel", " Vertex appears ", nb,
                            " times in Line ", line );
                        return false;
                    }
                }
                // Check that all the lines are in incident_entity of this corner
                gmme_id corner_id( Corner< DIMENSION >::type_name_static(),
                    entities.find( Corner< DIMENSION >::type_name_static() )->second.front() );
                for( index_t line : lines ) {
                    gmme_id line_id( Line< DIMENSION >::type_name_static(), line );
                    if( !is_boundary_entity( geomodel, line_id, corner_id ) ) {
                        Logger::warn( "GeoModel",
                            " Inconsistent Line-Corner connectivity ",
                            " vertex shows that ", line_id,
                            " must be in the boundary of ", corner_id );
                        return false;
                    }
                }
                return true;
            }
        }
    }

    template< index_t DIMENSION >
    bool is_corner_valid(
        const std::map< MeshEntityType, std::vector< index_t > >& entities )
    {
        const std::vector< index_t >& corners = entities.find(
            Corner< DIMENSION >::type_name_static() )->second;
        if( corners.size() > 1 ) {
            print_error( corners, "Corners" );
            Logger::warn( "GeoModel", "It should be in only one Corner" );
            return false;
        }
        return true;
    }

    template< index_t DIMENSION >
    bool is_geomodel_vertex_valid_base(
        const GeoModel< DIMENSION >& geomodel,
        std::map< MeshEntityType, std::vector< index_t > >& entities )
    {
        if( !is_corner_valid< DIMENSION >( entities ) ) {
            return false;
        }
        if( !is_line_vertex_valid< DIMENSION >( geomodel, entities ) ) {
            return false;
        }
        if( !is_surface_vertex_valid< DIMENSION >( geomodel, entities ) ) {
            return false;
        }

        return true;
    }

    template< index_t DIMENSION >
    bool is_geomodel_vertex_valid(
        const GeoModel< DIMENSION >& geomodel,
        index_t i );

    template< >
    bool is_geomodel_vertex_valid( const GeoModel3D& geomodel, index_t i )
    {
        // Get the mesh entities in which this vertex is
        std::map< MeshEntityType, std::vector< index_t > > entities = get_entities(
            geomodel, i );

        if( !is_geomodel_vertex_valid_base( geomodel, entities ) ) {
            return false;
        }
        if( !is_region_vertex_valid< 3 >( geomodel, entities ) ) {
            return false;
        }

        return true;
    }

    template< >
    bool is_geomodel_vertex_valid( const GeoModel2D& geomodel, index_t i )
    {
        // Get the mesh entities in which this vertex is
        std::map< MeshEntityType, std::vector< index_t > > entities = get_entities(
            geomodel, i );

        return is_geomodel_vertex_valid_base( geomodel, entities );
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
    template< index_t DIMENSION >
    bool check_model_points_validity( const GeoModel< DIMENSION >& geomodel )
    {
        // For all the vertices of the geomodel 
        // We check that the entities in which they are are consistent 
        // to have a valid B-Rep geomodel
        std::vector< bool > valid( geomodel.mesh.vertices.nb(), true );
        for( index_t i : range( geomodel.mesh.vertices.nb() ) ) {
            valid[i] = is_geomodel_vertex_valid( geomodel, i );
            if( !valid[i] ) {
                Logger::warn( "GeoModel", " Vertex ", i, " is not valid" );
            }
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

    template< index_t DIMENSION >
    void save_edges(
        const std::ostringstream& file,
        const GeoModel< DIMENSION >& geomodel,
        const std::vector< index_t >& e )
    {
        GEO::Mesh edge_mesh;
        index_t previous_vertex_id = NO_ID;
        for( index_t i : range( e.size() ) ) {
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

    template< index_t DIMENSION >
    void save_polygons(
        const std::string& file,
        const Surface< DIMENSION >& surface,
        const std::vector< index_t >& polygons )
    {
        GEO::Mesh mesh;
        for( index_t cur_polygon : polygons ) {
            index_t nb_vertices_in_polygon = surface.nb_mesh_element_vertices(
                cur_polygon );
            GEO::vector< index_t > vertices;
            vertices.reserve( nb_vertices_in_polygon );
            for( index_t v : range( nb_vertices_in_polygon ) ) {
                index_t new_vertex = mesh.vertices.create_vertex(
                    surface.mesh_element_vertex(
                        ElementLocalVertex( cur_polygon, v ) ).data() );
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
    template< index_t DIMENSION >
    bool surface_boundary_valid( const Surface< DIMENSION >& surface )
    {
        const GeoModelMeshVertices< DIMENSION >& geomodel_vertices =
            surface.geomodel().mesh.vertices;
        std::vector< index_t > invalid_corners;
        gmme_id S_id = surface.gmme();
        for( index_t p : range( surface.nb_mesh_elements() ) ) {
            for( index_t v : range( surface.nb_mesh_element_vertices( p ) ) ) {
                if( surface.polygon_adjacent_index( PolygonLocalEdge( p, v ) )
                    == NO_ID
                    && !is_edge_on_line( surface.geomodel(),
                        geomodel_vertices.geomodel_vertex_id( S_id,
                            ElementLocalVertex( p, v ) ),
                        geomodel_vertices.geomodel_vertex_id( S_id,
                            surface.low_level_mesh_storage().next_polygon_vertex(
                                ElementLocalVertex( p, v ) ) ) ) ) {
                    invalid_corners.push_back(
                        geomodel_vertices.geomodel_vertex_id( S_id,
                            ElementLocalVertex( p, v ) ) );
                    invalid_corners.push_back(
                        geomodel_vertices.geomodel_vertex_id( S_id,
                            surface.low_level_mesh_storage().next_polygon_vertex(
                                ElementLocalVertex( p, v ) ) ) );
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
    template< index_t DIMENSION >
    void debug_save_non_manifold_edges(
        const GeoModel< DIMENSION >& geomodel,
        const std::vector< index_t >& edge_indices,
        const std::vector< index_t >& non_manifold_edges )
    {
        GeogramLineMesh< DIMENSION > mesh;
        GeogramLineMeshBuilder< DIMENSION > builder( mesh );
        index_t nb_edges = static_cast< index_t >( non_manifold_edges.size() );
        builder.create_vertices( 2 * nb_edges );
        builder.create_edges( nb_edges );
        const GeoModelMeshVertices< DIMENSION >& vertices = geomodel.mesh.vertices;
        for( index_t e : range( non_manifold_edges.size() ) ) {
            index_t edge_id = non_manifold_edges[e];
            const vecn< DIMENSION >& v0 = vertices.vertex( edge_indices[edge_id] );
            const vecn< DIMENSION >& v1 = vertices.vertex(
                edge_indices[edge_id + 1] );
            builder.set_vertex( 2 * e, v0 );
            builder.set_vertex( 2 * e + 1, v1 );
            builder.set_edge_vertex( e, 0, 2 * e );
            builder.set_edge_vertex( e, 1, 2 * e + 1 );
        }
        mesh.save_mesh( validity_errors_directory + "/non_manifold_edges.geogram" );
    }

    template< index_t DIMENSION >
    bool is_surface_conformal_to_volume(
        const Surface< DIMENSION >& surface,
        const NNSearch< DIMENSION >& cell_facet_barycenter_nn_search )
    {
        std::vector< index_t > unconformal_polygons;
        for( index_t p : range( surface.nb_mesh_elements() ) ) {
            vecn< DIMENSION > center = surface.mesh_element_barycenter( p );
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

    template< index_t DIMENSION >
    std::vector< index_t > compute_border_edges(
        const GeoModel< DIMENSION >& geomodel )
    {
        std::vector< index_t > edge_indices;
        const GeoModelMeshPolygons< DIMENSION >& polygons = geomodel.mesh.polygons;
        for( const auto& surface : geomodel.surfaces() ) {
            for( index_t p : range( polygons.nb_polygons( surface.index() ) ) ) {
                index_t polygon_id = polygons.polygon( surface.index(), p );
                for( index_t v : range( polygons.nb_vertices( polygon_id ) ) ) {
                    index_t adj = polygons.adjacent(
                        PolygonLocalEdge( polygon_id, v ) );
                    if( adj == NO_ID ) {
                        edge_indices.push_back(
                            polygons.vertex( ElementLocalVertex( polygon_id, v ) ) );
                        index_t next_v = ( v + 1 )
                            % polygons.nb_vertices( polygon_id );
                        edge_indices.push_back(
                            polygons.vertex(
                                ElementLocalVertex( polygon_id, next_v ) ) );
                    }
                }
            }
        }
        return edge_indices;
    }

    template< index_t DIMENSION >
    std::vector< vecn< DIMENSION > > compute_border_edge_barycenters(
        const GeoModel< DIMENSION >& geomodel,
        const std::vector< index_t >& edge_indices )
    {
        const GeoModelMeshVertices< DIMENSION >& vertices = geomodel.mesh.vertices;
        index_t nb_edges = static_cast< index_t >( edge_indices.size() / 2 );
        std::vector< vecn< DIMENSION > > edge_barycenters;
        edge_barycenters.reserve( nb_edges );
        for( index_t e = 0; e < edge_indices.size(); e += 2 ) {
            const vecn< DIMENSION >& v0 = vertices.vertex( edge_indices[e] );
            const vecn< DIMENSION >& v1 = vertices.vertex( edge_indices[e + 1] );
            edge_barycenters.push_back( ( v0 + v1 ) * 0.5 );
        }
        return edge_barycenters;
    }

    template< index_t DIMENSION >
    std::vector< bool > compute_edge_on_lines(
        const GeoModel< DIMENSION >& geomodel,
        const std::vector< vecn< DIMENSION > >& edge_barycenters )
    {
        std::vector< bool > edge_on_lines( edge_barycenters.size(), false );
        NNSearch< DIMENSION > nn( edge_barycenters );
        for( const auto& line : geomodel.lines() ) {
            for( index_t e : range( line.nb_mesh_elements() ) ) {
                const vecn< DIMENSION > query = line.mesh_element_barycenter( e );
                std::vector< index_t > results = nn.get_neighbors( query,
                    geomodel.epsilon() );
                for( index_t edge : results ) {
                    edge_on_lines[edge] = true;
                }
            }
        }
        return edge_on_lines;
    }

    std::vector< index_t > compute_non_manifold_edges(
        const std::vector< bool >& edge_on_lines )
    {
        std::vector< index_t > non_manifold_edges;
        for( index_t e : range( edge_on_lines.size() ) ) {
            if( !edge_on_lines[e] ) {
                non_manifold_edges.push_back( e );
            }
        }
        return non_manifold_edges;
    }

    /*!
     * @brief Implementation class for validity checks on a GeoModel
     */
    template< index_t DIMENSION >
    class GeoModelValidityCheck {
    public:
        GeoModelValidityCheck(
            const GeoModel< DIMENSION >& geomodel,
            const ValidityCheckMode validity_check_mode )
            : geomodel_( geomodel ), valid_( true ), mode_( validity_check_mode )
        {
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
            do_check_validity();
            return valid_;
        }

    private:
        void add_base_checks( std::vector< std::thread >& threads )
        {
            if( enum_contains( mode_, ValidityCheckMode::FINITE_EXTENSION ) ) {
                threads.emplace_back( &GeoModelValidityCheck::test_finite_extension,
                    this );
            }
            if( enum_contains( mode_, ValidityCheckMode::GEOMODEL_CONNECTIVITY ) ) {
                threads.emplace_back(
                    &GeoModelValidityCheck::test_geomodel_connectivity_validity,
                    this );
            }
            if( enum_contains( mode_, ValidityCheckMode::GEOLOGICAL_ENTITIES ) ) {
                threads.emplace_back(
                    &GeoModelValidityCheck::test_geomodel_geological_validity,
                    this );

            }
            if( enum_contains( mode_,
                ValidityCheckMode::SURFACE_LINE_MESH_CONFORMITY ) ) {
                threads.emplace_back(
                    &GeoModelValidityCheck::test_surface_line_mesh_conformity,
                    this );
            }
            if( enum_contains( mode_, ValidityCheckMode::MESH_ENTITIES ) ) {
                threads.emplace_back(
                    &GeoModelValidityCheck::test_geomodel_mesh_entities_validity,
                    this );
                /// TODO: find a way to add this test for Model3d. See BC.
//                threads.emplace_back(
//                    &GeoModelValidityCheck::test_non_free_line_at_two_interfaces_intersection,
//                    this );
            }
            if( enum_contains( mode_, ValidityCheckMode::NON_MANIFOLD_EDGES ) ) {
                threads.emplace_back(
                    &GeoModelValidityCheck::test_non_manifold_edges, this );
            }
        }

        void add_checks( std::vector< std::thread >& threads )
        {
            add_base_checks( threads );
        }

        void do_check_validity()
        {
            std::vector< std::thread > threads;
            threads.reserve( 8 );
            add_checks( threads );

            for( std::thread& thread : threads ) {
                thread.join();
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
         * @todo Implement this check
         */
        void test_finite_extension()
        {
            if( !geomodel_.universe().is_valid() ) {
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
            for( const auto& surface : geomodel_.surfaces() ) {
                if( !surface_boundary_valid( surface ) ) {
                    set_invalid_model();
                }
            }
        }

        void test_region_surface_mesh_conformity()
        {
            if( geomodel_.mesh.cells.nb() > 0 ) {
                // Check the consistency between Surface polygons and Region cell facets
                const NNSearch< DIMENSION >& nn_search =
                    geomodel_.mesh.cells.cell_facet_nn_search();
                for( const auto& surface : geomodel_.surfaces() ) {
                    if( !is_surface_conformal_to_volume( surface, nn_search ) ) {
                        set_invalid_model();
                    }
                }
            }
        }

        void test_non_free_line_at_two_interfaces_intersection()
        {
            if( !geomodel_.entity_type_manager().geological_entity_manager.is_valid_type(
                Interface< DIMENSION >::type_name_static() ) ) {
                return;
            }
            for( const auto& line : geomodel_.lines() ) {
                if( line.nb_incident_entities() == 1 ) {
                    continue;
                }

                const index_t first_interface_id =
                    line.incident_entity( 0 ).parent_gmge(
                        Interface< DIMENSION >::type_name_static() ).index();
                ringmesh_assert( first_interface_id != NO_ID );
                bool at_least_two_different_interfaces = false;
                for( index_t in_boundary_i : range( 1, line.nb_incident_entities() ) ) {
                    const index_t cur_interface_id = line.incident_entity(
                        in_boundary_i ).parent_gmge(
                        Interface< DIMENSION >::type_name_static() ).index();
                    ringmesh_assert( cur_interface_id != NO_ID );
                    if( cur_interface_id != first_interface_id ) {
                        at_least_two_different_interfaces = true;
                        break;
                    }
                }

                if( !at_least_two_different_interfaces ) {
                    Logger::warn( "GeoModel",
                        "All in boundaries (surfaces) of line ", line.index(),
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
            std::vector< index_t > edge_indices = compute_border_edges( geomodel_ );
            std::vector< vecn< DIMENSION > > edge_barycenters =
                compute_border_edge_barycenters( geomodel_, edge_indices );
            std::vector< bool > edge_on_lines = compute_edge_on_lines( geomodel_,
                edge_barycenters );
            std::vector< index_t > non_manifold_edges = compute_non_manifold_edges(
                edge_on_lines );

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
                StoreIntersections< DIMENSION > action( geomodel_,
                    has_intersection );
                const SurfaceAABBTree< DIMENSION >& AABB =
                    geomodel_.mesh.polygons.aabb();
                AABB.compute_self_element_bbox_intersections( action );

                index_t nb_intersections = static_cast< index_t >( std::count(
                    has_intersection.begin(), has_intersection.end(), 1 ) );

                if( nb_intersections > 0 ) {
                    GEO::Mesh mesh;
                    for( index_t p : range( has_intersection.size() ) ) {
                        if( !has_intersection[p] ) continue;
                        GEO::vector< index_t > vertices;
                        vertices.reserve( geomodel_.mesh.polygons.nb_vertices( p ) );
                        for( index_t v : range(
                            geomodel_.mesh.polygons.nb_vertices( p ) ) ) {
                            index_t id = mesh.vertices.create_vertex(
                                geomodel_.mesh.vertices.vertex(
                                    geomodel_.mesh.polygons.vertex(
                                        ElementLocalVertex( p, v ) ) ).data() );
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
        const GeoModel< DIMENSION >& geomodel_;
        bool valid_;
        ValidityCheckMode mode_;
    };
    template< >
    void GeoModelValidityCheck< 3 >::add_checks(
        std::vector< std::thread >& threads )
    {
        if( enum_contains( mode_, ValidityCheckMode::POLYGON_INTERSECTIONS ) ) {
            threads.emplace_back( &GeoModelValidityCheck::test_polygon_intersections,
                this );
        }
        if( enum_contains( mode_,
            ValidityCheckMode::REGION_SURFACE_MESH_CONFORMITY ) ) {
            threads.emplace_back(
                &GeoModelValidityCheck::test_region_surface_mesh_conformity, this );
        }
        add_base_checks( threads );

    }

}

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

    template< index_t DIMENSION >
    bool are_geomodel_mesh_entities_mesh_valid(
        const GeoModel< DIMENSION >& geomodel )
    {
        const std::vector< MeshEntityType >& meshed_types =
            geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();
        index_t count_invalid = 0;
        for( const MeshEntityType& type : meshed_types ) {
            index_t nb_entities = geomodel.nb_mesh_entities( type );
            for( index_t i : range( nb_entities ) ) {
                const GeoModelMeshEntity< DIMENSION >& E = geomodel.mesh_entity(
                    type, i );
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

    template< index_t DIMENSION >
    bool are_geomodel_mesh_entities_connectivity_valid(
        const GeoModel< DIMENSION >& geomodel )
    {
        const std::vector< MeshEntityType >& meshed_types =
            geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();
        index_t count_invalid = 0;
        for( const MeshEntityType& type : meshed_types ) {
            index_t nb_entities = geomodel.nb_mesh_entities( type );
            for( index_t i : range( nb_entities ) ) {
                const GeoModelMeshEntity< DIMENSION >& E = geomodel.mesh_entity(
                    type, i );
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

    template< index_t DIMENSION >
    bool are_geomodel_geological_entities_valid(
        const GeoModel< DIMENSION >& geomodel )
    {
        const std::vector< GeologicalEntityType >& geological_types =
            geomodel.entity_type_manager().geological_entity_manager.geological_entity_types();
        index_t count_invalid = 0;
        for( const GeologicalEntityType& type : geological_types ) {
            index_t nb_entities = geomodel.nb_geological_entities( type );
            for( index_t i : range( nb_entities ) ) {
                const GeoModelGeologicalEntity< DIMENSION >& E =
                    geomodel.geological_entity( type, i );
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

    template< index_t DIMENSION >
    bool are_geomodel_mesh_entities_parent_valid(
        const GeoModel< DIMENSION >& geomodel )
    {
        const std::vector< MeshEntityType >& meshed_types =
            geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();
        index_t count_invalid = 0;
        for( const MeshEntityType& type : meshed_types ) {
            index_t nb_entities = geomodel.nb_mesh_entities( type );
            for( index_t i : range( nb_entities ) ) {
                const GeoModelMeshEntity< DIMENSION >& E = geomodel.mesh_entity(
                    type, i );
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

    template< index_t DIMENSION >
    bool is_geomodel_valid(
        const GeoModel< DIMENSION >& geomodel,
        ValidityCheckMode validity_check_mode )
    {
        if( !GEO::CmdLine::get_arg_bool( "in:intersection_check" ) ) {
            validity_check_mode = validity_check_mode
                ^ ValidityCheckMode::POLYGON_INTERSECTIONS;
        }

        GeoModelValidityCheck< DIMENSION > validity_checker( geomodel,
            validity_check_mode );

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

    template bool RINGMESH_API is_geomodel_valid< 2 >(
        const GeoModel2D&,
        ValidityCheckMode );
    template bool RINGMESH_API are_geomodel_mesh_entities_mesh_valid(
        const GeoModel2D& );
    template bool RINGMESH_API are_geomodel_mesh_entities_connectivity_valid(
        const GeoModel2D& );
    template bool RINGMESH_API are_geomodel_mesh_entities_parent_valid(
        const GeoModel2D& );
    template bool RINGMESH_API are_geomodel_geological_entities_valid(
        const GeoModel2D& );

    template bool RINGMESH_API is_geomodel_valid< 3 >(
        const GeoModel3D&,
        ValidityCheckMode );
    template bool RINGMESH_API are_geomodel_mesh_entities_mesh_valid(
        const GeoModel3D& );
    template bool RINGMESH_API are_geomodel_mesh_entities_connectivity_valid(
        const GeoModel3D& );
    template bool RINGMESH_API are_geomodel_mesh_entities_parent_valid(
        const GeoModel3D& );
    template bool RINGMESH_API are_geomodel_geological_entities_valid(
        const GeoModel3D& );

}
