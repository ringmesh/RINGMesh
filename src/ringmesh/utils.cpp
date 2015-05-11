/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#include <ringmesh/utils.h>
#include <ringmesh/boundary_model.h>
#include <ringmesh/boundary_model_element.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/numerics/predicates.h>

#include <geogram/mesh/mesh_intersection.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/triangle_intersection.h>


#include <iostream>
#include <sstream>
#include <stack>
#include <algorithm>
#include <cstring>

namespace RINGMesh {

    /*!
     * Computes the volume of a Mesh cell
     * @param[in] M the mesh
     * @param[in] c the cell index
     * @return the volume of the cell
     */
    double Geom::mesh_cell_volume( const GEO::Mesh& M, index_t c )
    {
        switch( M.cells.type( c ) ) {
            case GEO::MESH_TET:
                return GEO::Geom::tetra_volume(
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 1 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 3 ) ) ) ;
            case GEO::MESH_PYRAMID:
                return GEO::Geom::tetra_volume(
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 1 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 4 ) ) )
                    + GEO::Geom::tetra_volume(
                        GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                        GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                        GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 3 ) ),
                        GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 4 ) ) ) ;
            case GEO::MESH_PRISM:
            case GEO::MESH_HEX: {
                vec3 ori( 0, 0, 0 ) ;
                double volume = 0 ;
                for( index_t f = 0; f < M.cells.nb_facets( c ); f++ ) {
                    switch( M.cells.facet_nb_vertices( c, f ) ) {
                        case 3:
                            volume += GEO::Geom::tetra_signed_volume(
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 0 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 1 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 2 ) ), ori ) ;
                            break ;
                        case 4:
                            volume += GEO::Geom::tetra_signed_volume(
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 0 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 1 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 2 ) ), ori ) ;
                            volume += GEO::Geom::tetra_signed_volume(
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 0 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 2 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 3 ) ), ori ) ;
                            break ;
                        default:
                            ringmesh_assert_not_reached;
                            return 0 ;
                        }
                    }
                ringmesh_debug_assert( volume > 0 ) ;
                return volume ;
            }
            default:
                ringmesh_assert_not_reached;
                return 0 ;
            }
        }

        /*!
         * Computes the Mesh cell facet barycenter
         * @param[in] M the mesh
         * @param[in] cell the cell index
         * @param[in] f the facet index in the cell
         * @return the cell facet center
         */
    vec3 Geom::mesh_cell_facet_center( const GEO::Mesh& M, index_t cell, index_t f )
    {
        vec3 result( 0.0, 0.0, 0.0 ) ;
        double count = 0.0 ;
        for( index_t v = 0; v < M.cells.facet_nb_vertices( cell, f ); ++v ) {
            result += GEO::Geom::mesh_vertex( M,
                M.cells.facet_vertex( cell, f, v ) ) ;
            count += 1.0 ;
        }
        return ( 1.0 / count ) * result ;
    }

    /*!
     * Computes the Mesh cell facet normal
     * @param[in] M the mesh
     * @param[in] c the cell index
     * @param[in] f the facet index in the cell
     * @return the cell facet normal
     */
    vec3 Geom::mesh_cell_facet_normal( const GEO::Mesh& M, index_t c, index_t f )
    {
        const vec3& p1 = GEO::Geom::mesh_vertex( M,
            M.cells.facet_vertex( c, f, 0 ) ) ;
        const vec3& p2 = GEO::Geom::mesh_vertex( M,
            M.cells.facet_vertex( c, f, 1 ) ) ;
        const vec3& p3 = GEO::Geom::mesh_vertex( M,
            M.cells.facet_vertex( c, f, 2 ) ) ;
        return cross( p2 - p1, p3 - p1 ) ;
    }

    /*!
     * Computes the Mesh cell barycenter
     * @param[in] M the mesh
     * @param[in] c the cell index
     * @return the cell center
     */
    vec3 Geom::mesh_cell_center( const GEO::Mesh& M, index_t cell )
    {
        vec3 result( 0.0, 0.0, 0.0 ) ;
        double count = 0.0 ;
        for( index_t v = 0; v < M.cells.nb_vertices( cell ); ++v ) {
            result += GEO::Geom::mesh_vertex( M, M.cells.vertex( cell, v ) ) ;
            count += 1.0 ;
        }
        return ( 1.0 / count ) * result ;
    }

    MakeUnique::MakeUnique( const std::vector< vec3 >& points )
        : points_( points )
    {
        index_t nb_points = points_.size() ;
        indices_.resize( nb_points ) ;
        for( index_t i = 0; i < nb_points; i++ ) {
            indices_[i] = i ;
        }
    }

    /*!
     * Tests if a tetrahedron has an egde between two given points
     * @param[in] mesh the mesh
     * @param[in] t Tetrahedron index
     * @param[in] p0 First vertex index
     * @param[in] p1 Second vertex index
     * @param[out] edge Output edge index
     * @return The result of the test
     */
    bool Geom::has_edge(
        const GEO::Mesh& mesh,
        index_t t,
        index_t p0,
        index_t p1,
        index_t& edge )
    {
        for( uint8 e = 0; e < 6; e++ ) {
            index_t v0 = mesh.cells.edge_vertex( t, e, 0 ) ;
            index_t v1 = mesh.cells.edge_vertex( t, e, 1 ) ;
            if( ( p0 == v0 && p1 == v1 ) || ( p0 == v1 && p1 == v0 ) ) {
                edge = e ;
                return true ;
            }
        }
        return false ;
    }

    /*!
     * Gets all the next adjacent tetrahedra sharing an edge
     * @param[in] mesh the mesh
     * @param[in] t Starting tetrahedron index to test, should contain the edge
     * @param[in] prev Previous tetrahedron index
     * (if propagation arround the edge, prevent to go back were we came from)
     * @param[in] p0 First vertex index of the edge
     * @param[in] p1 Second vertex index of the edge
     * @return The edge index
     * \pre the mesh needs to be tetrahedralized
     */
    index_t Geom::next_arround_edge(
        const GEO::Mesh& mesh,
        index_t t,
        index_t prev,
        index_t p0,
        index_t p1 )
    {
        for( index_t adj = 0; adj < mesh.cells.nb_facets( t ); adj++ ) {
            index_t t_adj = mesh.cells.adjacent( t, adj ) ;
            if( t_adj == GEO::NO_CELL || t_adj == prev ) continue ;
            index_t edge ;
            if( has_edge( mesh, t_adj, p0, p1, edge ) ) {
                //todo handles any cell type
                return 6 * t_adj + edge ;
            }
        }
        return GEO::NO_CELL ;
    }

    /*!
     * Gets all the edge indices arround one edge
     * @param[in] mesh the mesh
     * @param[in] t First tetrahedron index to test, should include the edge
     * @param[in] p0 First vertex index of the edge
     * @param[in] p1 Second vertex index of the edge
     * @param[out] result Output list of edge indices
     */
    void Geom::edges_arround_edge(
        const GEO::Mesh& mesh,
        index_t t,
        index_t p0,
        index_t p1,
        std::vector< index_t >& result )
    {
        index_t prev = t ;
        int cur = t ;
        do {
            index_t info = next_arround_edge( mesh, cur, prev, p0, p1 ) ;
            if( info == GEO::NO_CELL ) return ;
            result.push_back( info ) ;
            prev = cur ;
            cur = info / 6 ;
        } while( cur != t ) ;
    }

    /*!
     * Get vertices when an edge is divide into \p nb_parts parts
     * @param[in] mesh the mesh
     * @param[in] edge the edge id in \p mesh
     * @param[in] nb_parts the number of edge division
     * @param[out] points the points which divide the edge
     */
    void Geom::divide_edge_in_parts(
        const GEO::Mesh& mesh,
        index_t edge,
        index_t nb_parts,
        std::vector< vec3 >& points )
    {
        points.resize( nb_parts - 1 ) ;
        double pond = 1. / nb_parts ;
        vec3 node0 = GEO::Geom::mesh_vertex(mesh,mesh.edges.vertex( edge, 0 ) ) ;
        vec3 node1 = GEO::Geom::mesh_vertex(mesh,mesh.edges.vertex( edge, 1 ) ) ;
        for( index_t i = 0; i < nb_parts - 1; i++ ) {
            for( index_t j = 0; j < 3; j++ ) {
                points[i][j] = ( i + 1 ) * pond * node1[j]
                    + ( 1. - ( i + 1 ) * pond ) * node0[j] ;
            }
        }

    }
    /*!
     * Gets the closest local vertex index in a mesh cell of a point
     * @param[in] mesh the mesh
     * @param[in] p the point to test
     * @param[in] t the cell index
     * @return the local vertex index
     */
    index_t Utils::get_nearest_vertex_index(
        const GEO::Mesh& mesh,
        const vec3& p,
        index_t t )
    {
        float64 dist = GEO::Numeric::max_float64() ;
        index_t result = -1 ;
        for( index_t v = 0; v < mesh.cells.nb_vertices( t ); v++ ) {
            float64 distance = length2(
                GEO::Geom::mesh_vertex( mesh, mesh.cells.vertex( t, v ) ) - p ) ;
            if( distance < dist ) {
                result = v ;
            }
        }
        return result ;
    }

    /*!
     * Tests if two adjacent facets have the same orientation
     * @param[in] mesh the mesh
     * @param[in[ f1 the first facet index
     * @param[in] c11 the corner index in the first facet
     * @param[in] f2 the second facet index
     * @return the result of the test
     */
    bool Utils::facets_have_same_orientation(
        const GEO::Mesh& mesh,
        index_t f1,
        index_t c11,
        index_t f2 )
    {
        index_t c12 = mesh.facets.next_corner_around_facet( f1, c11 ) ;
        index_t v11 = mesh.facet_corners.vertex( c11 ) ;
        index_t v12 = mesh.facet_corners.vertex( c12 ) ;
        for( index_t c21 = mesh.facets.corners_begin( f2 );
            c21 < mesh.facets.corners_end( f2 ); c21++ ) {
            index_t c22 = mesh.facets.next_corner_around_facet( f2, c21 ) ;
            index_t v21 = mesh.facet_corners.vertex( c21 ) ;
            index_t v22 = mesh.facet_corners.vertex( c22 ) ;
            if( v11 == v21 && v12 == v22 ) {
                return false ;
            }
            if( v11 == v22 && v12 == v21 ) {
                return true ;
            }
        }
        return true ;
    }

    /*!
     * Repair the consistency between a BoundaryModel region
     * and its volumetric Mesh. It repairs duplicated facets and facet orientation
     * @param[in] region the BoundaryModel region
     * @param[in] mesh the mesh to repair
     * @param[in] check_duplicated_facet the test of duplicated facets is optional
     */
    void Utils::check_and_repair_mesh_consistency(
        const BoundaryModelElement& region,
        GEO::Mesh& mesh,
        bool check_duplicated_facet )
    {
        if( mesh.facets.nb() == 0 ) return ;

        GEO::Attribute< index_t > attribute( mesh.facets.attributes(),
            surface_att_name ) ;

        /// 0 - Remove duplicated facets (optionnal)
        if( check_duplicated_facet ) {
            std::vector< vec3 > barycenters( mesh.facets.nb(), vec3( 0, 0, 0 ) ) ;
            for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                barycenters[f] = GEO::Geom::mesh_facet_center( mesh, f ) ;
            }

            MakeUnique unique( barycenters ) ;
            unique.unique() ;
            const std::vector< index_t > indices = unique.indices() ;
            GEO::vector< index_t > facet_to_remove( mesh.facets.nb(), 0 ) ;
            signed_index_t cur_id = 0 ;
            for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                if( cur_id == indices[f] ) {
                    cur_id++ ;
                } else {
                    facet_to_remove[f] = 1 ;
                }
            }
            mesh.facets.delete_elements( facet_to_remove ) ;

            if( GEO::Attribute< index_t >::is_defined( mesh.facets.attributes(),
                surface_att_name ) ) {
                GEO::Attribute< index_t > attribute( mesh.facets.attributes(),
                    surface_att_name ) ;
                index_t offset = 0 ;
                cur_id = 0 ;
                for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                    if( cur_id == indices[f] ) {
                        cur_id++ ;
                        attribute[f - offset] = attribute[f] ;
                    } else {
                        offset++ ;
                    }
                }
                attribute.redim( attribute.size() - offset ) ;
            }
            mesh.facets.connect() ;
        }

        /// 1 - Check facet adjacencies for non-manifold surfaces
        std::vector< index_t > temp ;
        temp.reserve( 6 ) ;
        std::vector< std::vector< index_t > > stars( mesh.vertices.nb(), temp ) ;
        for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
            for( index_t c = mesh.facets.corners_begin( f );
                c < mesh.facets.corners_end( f ); c++ ) {
                stars[mesh.facet_corners.vertex( c )].push_back( f ) ;
            }
        }
        for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
            index_t surface_id = attribute[f] ;
            for( index_t c = mesh.facets.corners_begin( f );
                c < mesh.facets.corners_end( f ); c++ ) {
                index_t f_adj = mesh.facet_corners.adjacent_facet( c ) ;
                if( f_adj != GEO::NO_FACET && attribute[f_adj] != surface_id ) {
                    f_adj = GEO::NO_FACET ;
                }
                if( f_adj == GEO::NO_FACET ) {
                    const std::vector< index_t >& star0 =
                        stars[mesh.facet_corners.vertex( c )] ;
                    const std::vector< index_t >& star1 =
                        stars[mesh.facet_corners.vertex(
                            mesh.facets.next_corner_around_facet( f, c ) )] ;
                    std::vector< index_t > intersect(
                        std::min( star0.size(), star1.size() ) ) ;
                    intersect.erase(
                        std::set_intersection( star0.begin(), star0.end(),
                            star1.begin(), star1.end(), intersect.begin() ),
                        intersect.end() ) ;
                    if( intersect.size() > 1 ) {
                        for( index_t i = 0; i < intersect.size(); i++ ) {
                            index_t cur_f = intersect[i] ;
                            if( cur_f != f && attribute[cur_f] == surface_id ) {
                                f_adj = cur_f ;
                            }
                        }
                    }
                }
                mesh.facet_corners.set_adjacent_facet( c, f_adj ) ;
            }
        }

        /// 2 - Reorient in the same direction using propagation
        std::vector< bool > facet_visited( mesh.facets.nb(), false ) ;
        for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
            if( facet_visited[f] ) continue ;
            index_t surface_id = attribute[f] ;
            std::stack< index_t > S ;
            S.push( f ) ;
            do {
                index_t cur_f = S.top() ;
                S.pop() ;
                if( facet_visited[cur_f] ) continue ;
                facet_visited[cur_f] = true ;
                for( index_t c = mesh.facets.corners_begin( cur_f );
                    c < mesh.facets.corners_end( cur_f ); c++ ) {
                    index_t f_adj = mesh.facet_corners.adjacent_facet( c ) ;
                    if( f_adj == GEO::NO_FACET || attribute[f_adj] != surface_id
                        || facet_visited[f_adj] ) continue ;
                    if( !facets_have_same_orientation( mesh, cur_f, c, f_adj ) ) {
                        mesh.facets.flip( f_adj ) ;
                    }
                    S.push( f_adj ) ;
                }
            } while( !S.empty() ) ;
        }

        /// 3 - Check for consistent orientation with BoundaryModel
        GEO::MeshFacetsAABB aabb( mesh ) ;
        std::vector< bool > flip_surface( region.model().nb_surfaces(), false ) ;
        bool flip_sthg = false ;
        for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
            const Surface& surface = dynamic_cast< const Surface& >( region.boundary(
                s ) ) ;
            vec3 barycenter = surface.facet_barycenter( 0 ) ;
            vec3 nearest_point ;
            float64 distance ;
            index_t f = aabb.nearest_facet( barycenter, nearest_point, distance ) ;
            ringmesh_debug_assert( surface.bme_id().index == attribute[f] ) ;

            vec3 ori_normal = surface.facet_normal( 0 ) ;
            vec3 test_normal = GEO::Geom::mesh_facet_normal( mesh, f ) ;
            if( dot( ori_normal, test_normal ) < 0 ) {
                flip_surface[surface.bme_id().index] = true ;
                flip_sthg = true ;
            }
        }
        if( flip_sthg ) {
            for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                index_t surface_id = attribute[f] ;
                if( flip_surface[surface_id] ) {
                    mesh.facets.flip( f ) ;
                }
            }
        }
    }

    void Utils::print_bounded_attributes( const GEO::Mesh& M )
    {
        {
            GEO::vector< std::string > names ;
            M.vertices.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size(), false ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    if( names[a] == "point" ) continue ;
                    is_bounded[a] = M.vertices.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on vertices:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.edges.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.edges.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on edges:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.facets.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.facets.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on facets:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.facet_corners.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] =
                        M.facet_corners.attributes().find_attribute_store( names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on facet_corners:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.cells.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.cells.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on cells:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.cell_corners.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.cell_corners.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on cell_corners:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.cell_facets.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.cell_facets.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on cell_facets:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
    }

    /*!
     * Compares the contains of two files
     * @param[in] f1 the first filename
     * @param[in] f2 the second filename
     * @return return True if the files are identical
     */
    bool Utils::compare_file( const std::string& f1, const std::string& f2 )
    {
        const unsigned int MAX_LINE_LEN = 65535 ;

        std::ifstream lFile( f1.c_str() ) ;
        std::ifstream rFile( f2.c_str() ) ;

        char* lBuffer = new char[MAX_LINE_LEN]() ;
        char* rBuffer = new char[MAX_LINE_LEN]() ;

        do {
            lFile.read( lBuffer, MAX_LINE_LEN ) ;
            rFile.read( rBuffer, MAX_LINE_LEN ) ;
            unsigned int numberOfRead = lFile.gcount() ;

            if( std::memcmp( lBuffer, rBuffer, numberOfRead ) != 0 ) {
                delete[] lBuffer ;
                delete[] rBuffer ;
                return false ;
            }
        } while( lFile.good() || rFile.good() ) ;
        delete[] lBuffer ;
        delete[] rBuffer ;
        return true ;
    }

    /*!
     * Computes the intersection(s) between a circle and a plane
     * @param[in] O_plane a point on the plane
     * @param[in] N_plane the normal of the plane
     * @param[in] O_circle the center of the circle
     * @param[in] N_circle the normal of the plane supporting the circle
     * @param[in] r the radius of the circle
     * @param[out] result the intersected points
     * @return returns true if there is at least one intersection
     */
    bool Math::circle_plane_intersection(
        const vec3& O_plane,
        const vec3& N_plane,
        const vec3& O_circle,
        const vec3& N_circle,
        float64 r,
        std::vector< vec3 >& result )
    {
        vec3 O_inter, D_inter ;
        if( !plan_plane_intersection( O_plane, N_plane, O_circle, N_circle, O_inter,
            D_inter ) ) {
            return false ;
        }

        // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
        // Locate one or two points that are on the circle and line.  If the
        // line is t*D+P, the circle center is C, and the circle radius is r,
        // then r^2 = |t*D+P-C|^2 = |D|^2*t^2 + 2*Dot(D,P-C)*t + |P-C|^2.  This
        // is a quadratic equation of the form:  a2*t^2 + 2*a1*t + a0 = 0.
        vec3 diff = O_inter - O_circle ;
        float64 a2 = D_inter.length2() ;
        float64 a1 = dot( diff, D_inter ) ;
        float64 a0 = diff.length2() - r * r ;

        float64 discr = a1 * a1 - a0 * a2 ;
        if( discr < 0.0 ) return false ;

        if( fabs( a2 ) < epsilon ) return false ;
        float64 inv = 1.0 / a2 ;
        if( discr < epsilon ) {
            result.push_back( vec3( O_inter - ( a1 * inv ) * D_inter ) ) ;
        } else {
            float64 root = sqrt( discr ) ;
            result.push_back( vec3( O_inter - ( ( a1 + root ) * inv ) * D_inter ) ) ;
            result.push_back( vec3( O_inter - ( ( a1 - root ) * inv ) * D_inter ) ) ;
        }
        return true ;
    }

    /*!
     * Computes the intersection between two planes
     * @param[in] O_P0 a point on the first plane
     * @param[in] N_P0 the normal of the frst plane
     * @param[in] O_P1 a point on the second plane
     * @param[in] N_P1 the normal of the second plane
     * @param[out] O_inter a point on the intersected line
     * @param[out] D_inter the direction of the intersected line
     * @return true is there is an intersection between the planes
     */
    bool Math::plan_plane_intersection(
        const vec3& O_P0,
        const vec3& N_P0,
        const vec3& O_P1,
        const vec3& N_P1,
        vec3& O_inter,
        vec3& D_inter )
    {
        // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
        // If N0 and N1 are parallel, either the planes are parallel and separated
        // or the same plane.  In both cases, 'false' is returned.  Otherwise,
        // the intersection line is
        //   L(t) = t*Cross(N0,N1)/|Cross(N0,N1)| + c0*N0 + c1*N1
        // for some coefficients c0 and c1 and for t any real number (the line
        // parameter).  Taking dot products with the normals,
        //   d0 = Dot(N0,L) = c0*Dot(N0,N0) + c1*Dot(N0,N1) = c0 + c1*d
        //   d1 = Dot(N1,L) = c0*Dot(N0,N1) + c1*Dot(N1,N1) = c0*d + c1
        // where d = Dot(N0,N1).  These are two equations in two unknowns.  The
        // solution is
        //   c0 = (d0 - d*d1)/det
        //   c1 = (d1 - d*d0)/det
        // where det = 1 - d^2.

        float64 d = dot( N_P0, N_P1 ) ;
        if( fabs( d - 1 ) < epsilon ) return false ;

        float64 invDet = 1.0 / ( 1.0 - d * d ) ;
        float64 const_P0 = dot( N_P0, O_P0 ) ;
        float64 const_P1 = dot( N_P1, O_P1 ) ;
        float64 c0 = ( const_P0 - d * const_P1 ) * invDet ;
        float64 c1 = ( const_P1 - d * const_P0 ) * invDet ;
        O_inter = c0 * N_P0 + c1 * N_P1 ;
        D_inter = cross( N_P0, N_P1 ) ;
        return true ;
    }

    /*!
     * Computes the intersection(s) between a circle and a triangle
     * @param[in] p0 the first vertex of the triangle
     * @param[in] p1 the second vertex of the triangle
     * @param[in] p2 the third vertex of the triangle
     * @param[in] O_circle the center of the circle
     * @param[in] N_circle the normal of the plane supporting the circle
     * @param[in] r the radius of the circle
     * @param[out] result the intersected points
     * @return returns true if there is at least one intersection
     */
    bool Math::circle_triangle_intersection(
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& O_circle,
        const vec3& N_circle,
        float64 r,
        std::vector< vec3 >& result )
    {
        vec3 N_triangle = normalize( cross( p1 - p0, p2 - p0 ) ) ;
        vec3 barycenter = ( p0 + p1 + p2 ) / 3 ;
        std::vector< vec3 > inter_circle_plane ;
        if( circle_plane_intersection( barycenter, N_triangle, O_circle, N_circle, r,
            inter_circle_plane ) ) {
            for( index_t i = 0; i < inter_circle_plane.size(); i++ ) {
                const vec3& p = inter_circle_plane[i] ;
                if( point_inside_triangle( p, p0, p1, p2 ) ) {
                    result.push_back( p ) ;
                }
            }
        }
        return !result.empty() ;
    }

    /*!
     * Computes the orthogonal projection of a point on a segment
     * @param[in] p the point to project
     * @param[in] p0 the first vertex of the segment
     * @param[in] p1 the second vertex of the segment
     * @param[out] new_p the projected point
     * @return returns true if the projection is possible
     */
    bool Math::point_segment_projection(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& new_p )
    {
        vec3 center = ( p0 + p1 ) * 0.5 ;
        vec3 diff = p - center ;
        vec3 edge = p1 - p0 ;
        float64 extent = 0.5 * edge.length() ;
        edge = normalize( edge ) ;
        float64 d = dot( edge, diff ) ;

        if( fabs( d ) <= extent ) {
            new_p = center + d * edge ;
            return true ;
        }
        return false ;
    }

    /*!
     * Computes the smallest distance between a point and a quad
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the quad
     * @param[in] p1 the second vertex of the quad
     * @param[in] p2 the third vertex of the quad
     * @param[in] p3 the fourth vertex of the quad
     * @param[out] nearest_p the closest point on the quad
     * @return the smallest distance
     */
    float64 Math::point_quad_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p )
    {
        const vec3 center( ( p0 + p1 + p2 + p3 ) / 4. ) ;
        vec3 edge0( p1 - p0 ) ;
        vec3 edge1( p3 - p0 ) ;
        vec3 axis[2] ;
        axis[0] = normalize( edge0 ) ;
        axis[1] = normalize( edge1 ) ;
        float64 extent[2] ;
        extent[0] = edge0.length() / 2. ;
        extent[1] = edge1.length() / 2. ;

        vec3 diff = center - p ;
        float64 b0 = dot( diff, axis[0] ) ;
        float64 b1 = dot( diff, axis[1] ) ;
        float64 s0 = -b0 ;
        float64 s1 = -b1 ;
        float64 sqrDistance = dot( diff, diff ) ;

        if( s0 < -extent[0] ) {
            s0 = -extent[0] ;
        } else if( s0 > extent[0] ) {
            s0 = extent[0] ;
        }
        sqrDistance += s0 * ( s0 + 2. * b0 ) ;

        if( s1 < -extent[1] ) {
            s1 = -extent[1] ;
        } else if( s1 > extent[1] ) {
            s1 = extent[1] ;
        }
        sqrDistance += s1 * ( s1 + 2. * b1 ) ;

        // Account for numerical round-off error.
        if( sqrDistance < 0 ) {
            sqrDistance = 0 ;
        }

        float64 distance = sqrt( sqrDistance ) ;
        nearest_p = center ;
        nearest_p += s0 * axis[0] ;
        nearest_p += s1 * axis[1] ;

        return distance ;
    }

    /*!
     * Computes the intersection of a segment and a triangle
     * @param[in] seg0 the first vertex of the segment
     * @param[in] seg1 the second vertex of the segment
     * @param[in] trgl0 the first vertex of the triangle
     * @param[in] trgl1 the second vertex of the triangle
     * @param[in] trgl2 the third vertex of the triangle
     * @param[out] result the intersected point
     * @return true is there is an intersection
     */
    bool Math::segment_triangle_intersection(
        const vec3& seg0,
        const vec3& seg1,
        const vec3& trgl0,
        const vec3& trgl1,
        const vec3& trgl2,
        vec3& result )
    {
        // http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
        // Compute the offset origin, edges, and normal.
        vec3 seg_center = ( seg0 + seg1 ) / 2 ;
        vec3 diff = seg_center - trgl0 ;
        vec3 edge1 = trgl1 - trgl0 ;
        vec3 edge2 = trgl2 - trgl0 ;
        vec3 normal = cross( edge1, edge2 ) ;

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = segment direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        vec3 D = normalize( seg1 - seg0 ) ;
        float64 DdN = dot( D, normal ) ;
        signed_index_t sign ;
        if( DdN > epsilon ) {
            sign = 1 ;
        } else if( DdN < -epsilon ) {
            sign = -1 ;
            DdN = -DdN ;
        } else {
            // Segment and triangle are parallel, call it a "no intersection"
            // even if the segment does intersect.
            return false ;
        }

        float64 DdQxE2 = sign * dot( D, cross( diff, edge2 ) ) ;
        if( DdQxE2 >= 0 ) {
            float64 DdE1xQ = sign * dot( D, cross( edge1, diff ) ) ;
            if( DdE1xQ >= 0 ) {
                if( DdQxE2 + DdE1xQ <= DdN ) {
                    // Line intersects triangle, check if segment does.
                    float64 QdN = -sign * dot( diff, normal ) ;
                    float64 extDdN = length( seg1 - seg0 ) * DdN / 2. ;
                    if( -extDdN <= QdN && QdN <= extDdN ) {
                        // Segment intersects triangle.
                        float64 inv = 1. / DdN ;
                        float64 seg_parameter = QdN * inv ;

                        result = seg_center + seg_parameter * D ;
                        return true ;
                    }
                    // else: |t| > extent, no intersection
                }
                // else: b1+b2 > 1, no intersection
            }
            // else: b2 < 0, no intersection
        }
        // else: b1 < 0, no intersection
        return false ;
    }

    /*!
     * Tests if a point is inside a triangle, more precisely if it is inside
     * a prism based on the triangle and its normal
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the triangle
     * @param[in] p1 the second vertex of the triangle
     * @param[in] p2 the third vertex of the triangle
     * @return returns true if the point is inside
     */
    bool Math::point_inside_triangle(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2 )
    {
        vec3 n = cross( p2 - p0, p1 - p0 ) ;
        vec3 q = p + n ;

        Sign s1 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p0.data(), p1.data() ) ) ;
        Sign s2 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p1.data(), p2.data() ) ) ;
        Sign s3 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p2.data(), p0.data() ) ) ;
        if( s1 == ZERO || s2 == ZERO || s3 == ZERO ) {
            return true ; // Arbitrary choice !!!!
        }

        return s1 == s2 && s2 == s3 ;
    }

    /*!
     * Tests if a point is inside a quad, more precisely if it is inside the box
     * based on the quad and its normal
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the quad
     * @param[in] p1 the second vertex of the quad
     * @param[in] p2 the third vertex of the quad
     * @param[in] p3 the fourth vertex of the quad
     * @return returns true if the point is inside
     */
    bool Math::point_inside_quad(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 )
    {
        vec3 n = cross( p2 - p0, p1 - p0 ) ;
        vec3 q = p + n ;

        Sign s1 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p0.data(), p1.data() ) ) ;
        Sign s2 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p1.data(), p2.data() ) ) ;
        Sign s3 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p2.data(), p3.data() ) ) ;
        Sign s4 = sign(
            GEO::PCK::orient_3d( p.data(), q.data(), p3.data(), p0.data() ) ) ;

        if( s1 == ZERO || s2 == ZERO || s3 == ZERO || s4 == ZERO ) {
            if( Utils::inexact_equal( p, p0 ) || Utils::inexact_equal( p, p1 )
                || Utils::inexact_equal( p, p2 ) || Utils::inexact_equal( p, p3 ) ) {
                return true ;
            }
            return false ; // Arbitrary choice !!!!
        }

        return s1 == s2 && s2 == s3 && s3 == s4 ;
    }

    /*!
     * Computes the distance between a point and a tetrahedron
     * @param[in] p the point
     * @param[in] p0 the first vertex of the tetrahedron
     * @param[in] p1 the second vertex of the tetrahedron
     * @param[in] p2 the third vertex of the tetrahedron
     * @param[in] p3 the fourth vertex of the tetrahedron
     * @param[out] nearest_p the nearest point on the tetrahedron
     * @return the distance between the point and the tetrahedron facets
     */
    float64 Math::point_tetra_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p )
    {
        vec3 vertices[4] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        float64 dist = big_float64 ;
        for( uint8 f = 0; f < tetra_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance = point_triangle_distance( p,
                vertices[tetra_descriptor.facet[f][0]],
                vertices[tetra_descriptor.facet[f][1]],
                vertices[tetra_descriptor.facet[f][2]], cur_p ) ;
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    /*!
     * Computes the distance between a point and a pyramid
     * @param[in] p the point
     * @param[in] p0 the first vertex of the pyramid
     * @param[in] p1 the second vertex of the pyramid
     * @param[in] p2 the third vertex of the pyramid
     * @param[in] p3 the fourth vertex of the pyramid
     * @param[in] p4 the fifth vertex of the pyramid
     * @param[out] nearest_p the nearest point on the pyramid
     * @return the distance between the point and the pyramid facets
     */
    float64 Math::point_pyramid_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        vec3& nearest_p )
    {
        vec3 vertices[5] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        float64 dist = big_float64 ;
        for( uint8 f = 0; f < pyramid_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance ;
            uint8 nb_vertices = pyramid_descriptor.nb_vertices_in_facet[f] ;
            if( nb_vertices == 3 ) {
                distance = point_triangle_distance( p,
                    vertices[pyramid_descriptor.facet[f][0]],
                    vertices[pyramid_descriptor.facet[f][1]],
                    vertices[pyramid_descriptor.facet[f][2]], cur_p ) ;
            } else if( nb_vertices == 4 ) {
                distance = point_quad_distance( p,
                    vertices[pyramid_descriptor.facet[f][0]],
                    vertices[pyramid_descriptor.facet[f][1]],
                    vertices[pyramid_descriptor.facet[f][2]],
                    vertices[pyramid_descriptor.facet[f][3]], cur_p ) ;
            } else {
                ringmesh_assert_not_reached;
            }
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    /*!
     * Computes the distance between a point and a prism
     * @param[in] p the point
     * @param[in] p0 the first vertex of the prism
     * @param[in] p1 the second vertex of the prism
     * @param[in] p2 the third vertex of the prism
     * @param[in] p3 the fourth vertex of the prism
     * @param[in] p4 the fifth vertex of the prism
     * @param[in] p5 the sixth vertex of the prism
     * @param[out] nearest_p the nearest point on the prism
     * @return the distance between the point and the prism facets
     */
    float64 Math::point_prism_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        vec3& nearest_p )
    {
        vec3 vertices[6] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        vertices[5] = p5 ;
        float64 dist = big_float64 ;
        for( uint8 f = 0; f < prism_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance ;
            uint8 nb_vertices = prism_descriptor.nb_vertices_in_facet[f] ;
            if( nb_vertices == 3 ) {
                distance = point_triangle_distance( p,
                    vertices[prism_descriptor.facet[f][0]],
                    vertices[prism_descriptor.facet[f][1]],
                    vertices[prism_descriptor.facet[f][2]], cur_p ) ;
            } else if( nb_vertices == 4 ) {
                distance = point_quad_distance( p,
                    vertices[prism_descriptor.facet[f][0]],
                    vertices[prism_descriptor.facet[f][1]],
                    vertices[prism_descriptor.facet[f][2]],
                    vertices[prism_descriptor.facet[f][3]], cur_p ) ;
            } else {
                ringmesh_assert_not_reached;
            }
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }
    /*!
     * Computes the distance between a point and a hexahedron
     * @param[in] p the point
     * @param[in] p0 the first vertex of the hexahedron
     * @param[in] p1 the second vertex of the hexahedron
     * @param[in] p2 the third vertex of the hexahedron
     * @param[in] p3 the fourth vertex of the hexahedron
     * @param[in] p4 the fifth vertex of the hexahedron
     * @param[in] p5 the sixth vertex of the hexahedron
     * @param[in] p6 the seventh vertex of the hexahedron
     * @param[in] p7 the heith vertex of the hexahedron
     * @param[out] nearest_p the nearest point on the hexahedron
     * @return the distance between the point and the hexahedron facets
     */
    float64 Math::point_hexa_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        const vec3& p6,
        const vec3& p7,
        vec3& nearest_p )
    {
        vec3 vertices[8] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        vertices[5] = p5 ;
        vertices[6] = p6 ;
        vertices[7] = p7 ;
        float64 dist = big_float64 ;
        for( uint8 f = 0; f < hexa_descriptor.nb_facets; f++ ) {
            vec3 cur_p ;
            float64 distance = point_quad_distance( p,
                vertices[hexa_descriptor.facet[f][0]],
                vertices[hexa_descriptor.facet[f][1]],
                vertices[hexa_descriptor.facet[f][2]],
                vertices[hexa_descriptor.facet[f][3]], cur_p ) ;
            if( distance < dist ) {
                dist = distance ;
                nearest_p = cur_p ;
            }
        }
        return dist ;
    }

    /*!
     * Tests if a point is inside a tetrahedron
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the tetrahedron
     * @param[in] p1 the second vertex of the tetrahedron
     * @param[in] p2 the third vertex of the tetrahedron
     * @param[in] p3 the fourth vertex of the tetrahedron
     * @return returns true if the point is inside the tetrahedron
     */
    bool Math::point_inside_tetra(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3 )
    {
        vec3 vertices[4] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        for( uint8 f = 0; f < tetra_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[tetra_descriptor.facet[f][1]]
                    - vertices[tetra_descriptor.facet[f][0]],
                vertices[tetra_descriptor.facet[f][2]]
                    - vertices[tetra_descriptor.facet[f][0]] ) ;
            vec3 n = p
                - ( ( vertices[tetra_descriptor.facet[f][0]]
                    + vertices[tetra_descriptor.facet[f][1]]
                    + vertices[tetra_descriptor.facet[f][2]] ) / 3. ) ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    /*!
     * Tests if a point is inside a pyramid
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the pyramid
     * @param[in] p1 the second vertex of the pyramid
     * @param[in] p2 the third vertex of the pyramid
     * @param[in] p3 the fourth vertex of the pyramid
     * @param[in] p4 the fifth vertex of the pyramid
     * @return returns true if the point is inside the pyramid
     */
    bool Math::point_inside_pyramid(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4 )
    {
        vec3 vertices[5] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        for( uint8 f = 0; f < pyramid_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[pyramid_descriptor.facet[f][1]]
                    - vertices[pyramid_descriptor.facet[f][0]],
                vertices[pyramid_descriptor.facet[f][2]]
                    - vertices[pyramid_descriptor.facet[f][0]] ) ;
            uint8 nb_vertices = pyramid_descriptor.nb_vertices_in_facet[f] ;
            vec3 barycenter ;
            if( nb_vertices == 3 )
                barycenter = ( ( vertices[pyramid_descriptor.facet[f][0]]
                    + vertices[pyramid_descriptor.facet[f][1]]
                    + vertices[pyramid_descriptor.facet[f][2]] ) / 3. ) ;
            else if( nb_vertices == 4 )
                barycenter = ( ( vertices[pyramid_descriptor.facet[f][0]]
                    + vertices[pyramid_descriptor.facet[f][1]]
                    + vertices[pyramid_descriptor.facet[f][2]]
                    + vertices[pyramid_descriptor.facet[f][3]] ) / 4. ) ;
            else
                ringmesh_assert_not_reached;
            vec3 n = p - barycenter ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    /*!
     * Tests if a point is inside a prism
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the prism
     * @param[in] p1 the second vertex of the prism
     * @param[in] p2 the third vertex of the prism
     * @param[in] p3 the fourth vertex of the prism
     * @param[in] p4 the fifth vertex of the prism
     * @param[in] p5 the sixth vertex of the prism
     * @return returns true if the point is inside the prism
     */
    bool Math::point_inside_prism(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5 )
    {
        vec3 vertices[6] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        vertices[5] = p5 ;
        for( uint8 f = 0; f < prism_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[prism_descriptor.facet[f][1]]
                    - vertices[prism_descriptor.facet[f][0]],
                vertices[prism_descriptor.facet[f][2]]
                    - vertices[prism_descriptor.facet[f][0]] ) ;
            uint8 nb_vertices = prism_descriptor.nb_vertices_in_facet[f] ;
            vec3 barycenter ;
            if( nb_vertices == 3 )
                barycenter = ( ( vertices[prism_descriptor.facet[f][0]]
                    + vertices[prism_descriptor.facet[f][1]]
                    + vertices[prism_descriptor.facet[f][2]] ) / 3. ) ;
            else if( nb_vertices == 4 )
                barycenter = ( ( vertices[prism_descriptor.facet[f][0]]
                    + vertices[prism_descriptor.facet[f][1]]
                    + vertices[prism_descriptor.facet[f][2]]
                    + vertices[prism_descriptor.facet[f][3]] ) / 4. ) ;
            else
                ringmesh_assert_not_reached;
            vec3 n = p - barycenter ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }
    /*!
     * Tests if a point is inside a hexahedron
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the hexahedron
     * @param[in] p1 the second vertex of the hexahedron
     * @param[in] p2 the third vertex of the hexahedron
     * @param[in] p3 the fourth vertex of the hexahedron
     * @param[in] p4 the fifth vertex of the hexahedron
     * @param[in] p5 the sixth vertex of the hexahedron
     * @param[in] p6 the seventh vertex of the hexahedron
     * @param[in] p7 the heigth vertex of the hexahedron
     * @return returns true if the point is inside the hexahedron
     */
    bool Math::point_inside_hexa(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        const vec3& p6,
        const vec3& p7 )
    {
        vec3 vertices[8] ;
        vertices[0] = p0 ;
        vertices[1] = p1 ;
        vertices[2] = p2 ;
        vertices[3] = p3 ;
        vertices[4] = p4 ;
        vertices[5] = p5 ;
        vertices[6] = p6 ;
        vertices[7] = p7 ;
        for( uint8 f = 0; f < hexa_descriptor.nb_facets; f++ ) {
            vec3 N = cross(
                vertices[hexa_descriptor.facet[f][1]]
                    - vertices[hexa_descriptor.facet[f][0]],
                vertices[hexa_descriptor.facet[f][2]]
                    - vertices[hexa_descriptor.facet[f][0]] ) ;
            vec3 barycenter = ( ( vertices[hexa_descriptor.facet[f][0]]
                + vertices[hexa_descriptor.facet[f][1]]
                + vertices[hexa_descriptor.facet[f][2]]
                + vertices[hexa_descriptor.facet[f][3]] ) / 4. ) ;
            vec3 n = p - barycenter ;
            if( dot( N, n ) > 0 ) return false ;
        }
        return true ;
    }

    /*!
     * Gets a vector of unique points (initial vector - colocated points)
     * @param[out] results the vector to fill
     */
    void MakeUnique::unique_points( std::vector< vec3 >& results ) const
    {
        results.reserve( indices_.size() ) ;
        signed_index_t offset = 0, cur_id = 0 ;
        for( index_t p = 0; p < indices_.size(); p++ ) {
            if( cur_id == indices_[p] ) {
                cur_id++ ;
                results.push_back( points_[indices_[p] + offset] ) ;
            } else {
                offset++ ;
            }
        }
    }

    /*!
     * Computes the unique database
     */
    void MakeUnique::unique()
    {
        ColocaterANN ann( points_ ) ;
        for( index_t i = 0; i < indices_.size(); i++ ) {
            if( indices_[i] != i ) continue ;
            std::vector< index_t > results ;
            ann.get_colocated( points_[i], results ) ;
            index_t id = *std::min_element( results.begin(), results.end() ) ;
            for( index_t j = 0; j < results.size(); j++ ) {
                if( id == results[j] ) continue ;
                indices_[results[j]] = id ;
            }
        }
        index_t offset = 0 ;
        for( index_t i = 0; i < indices_.size(); i++ ) {
            if( indices_[i] != i ) {
                indices_[i] = indices_[indices_[i]] ;
                offset++ ;
            } else {
                indices_[i] -= offset ;
            }
        }
    }
    /*!
     * Add edges to the initial vector
     * @param[in] points the edges to add
     */
    void MakeUnique::add_edges( const std::vector< Edge >& points )
    {
        signed_index_t offset = points_.size() ;
        points_.resize( offset + ( points.size() * 2 ) ) ;
        indices_.resize( offset + ( points.size() * 2 ) ) ;
        for( index_t p = 0; p < points.size(); p++ ) {
            points_[offset] = points[p].value( 0 ) ;
            indices_[offset] = offset ;
            offset++ ;
            points_[offset] = points[p].value( 1 ) ;
            indices_[offset] = offset ;
            offset++ ;
        }
    }
    /*!
     * Add points to the initial vector
     * @param[in] points the points to add
     */
    void MakeUnique::add_points( const std::vector< vec3 >& points )
    {
        signed_index_t offset = points_.size() ;
        points_.resize( offset + points.size() ) ;
        indices_.resize( offset + points.size() ) ;
        for( index_t p = 0; p < points.size(); p++, offset++ ) {
            points_[offset] = points[p] ;
            indices_[offset] = offset ;
        }
    }

    ColocaterANN::ColocaterANN()
        : ann_points_( nil )
    {
    }

    ColocaterANN::ColocaterANN( const Surface& mesh, const MeshLocation& location )
    {
        ann_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
        switch( location ) {
            case VERTICES: {
                index_t nb_vertices = mesh.nb_vertices() ;
                ann_points_ = new double[nb_vertices * 3] ;
                for( index_t i = 0; i < mesh.nb_vertices(); i++ ) {
                    index_t index_in_ann = 3 * i ;
                    ann_points_[index_in_ann] = mesh.vertex( i ).x ;
                    ann_points_[index_in_ann + 1] = mesh.vertex( i ).y ;
                    ann_points_[index_in_ann + 2] = mesh.vertex( i ).z ;
                }
                ann_tree_->set_points( nb_vertices, ann_points_ ) ;
                break ;
            }
            case FACETS: {
                index_t nb_vertices = mesh.nb_cells() ;
                ann_points_ = new double[nb_vertices * 3] ;
                for( index_t i = 0; i < mesh.nb_cells(); i++ ) {
                    vec3 center = mesh.facet_barycenter( i ) ;
                    index_t index_in_ann = 3 * i ;
                    ann_points_[index_in_ann] = center.x ;
                    ann_points_[index_in_ann + 1] = center.y ;
                    ann_points_[index_in_ann + 2] = center.z ;
                }
                ann_tree_->set_points( nb_vertices, ann_points_ ) ;
                break ;
            }
        }

    }

    ColocaterANN::ColocaterANN( const Line& mesh )
    {
        index_t nb_vertices = mesh.nb_vertices() ;
        ann_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
        ann_points_ = new double[nb_vertices * 3] ;
        for( index_t i = 0; i < mesh.nb_vertices(); i++ ) {
            index_t index_in_ann = 3 * i ;
            ann_points_[index_in_ann] = mesh.vertex( i ).x ;
            ann_points_[index_in_ann + 1] = mesh.vertex( i ).y ;
            ann_points_[index_in_ann + 2] = mesh.vertex( i ).z ;
        }
        ann_tree_->set_points( nb_vertices, ann_points_ ) ;
    }

    ColocaterANN::ColocaterANN( const GEO::Mesh& mesh, const MeshLocation& location )
    {
        ann_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
        switch( location ) {
            case VERTICES: {
                index_t nb_vertices = mesh.vertices.nb() ;
                ann_points_ = new double[nb_vertices * 3] ;
                for( index_t i = 0; i < mesh.vertices.nb(); i++ ) {
                    index_t index_in_ann = 3 * i ;
                    ann_points_[index_in_ann] = mesh.vertices.point_ptr( i )[0] ;
                    ann_points_[index_in_ann + 1] = mesh.vertices.point_ptr( i )[1] ;
                    ann_points_[index_in_ann + 2] = mesh.vertices.point_ptr( i )[2] ;
                }
                ann_tree_->set_points( nb_vertices, ann_points_ ) ;
                break ;
            }
            case FACETS: {
                index_t nb_vertices = mesh.facets.nb() ;
                ann_points_ = new double[nb_vertices * 3] ;
                for( index_t i = 0; i < mesh.facets.nb(); i++ ) {
                    vec3 center = GEO::Geom::mesh_facet_center( mesh, i ) ;
                    index_t index_in_ann = 3 * i ;
                    ann_points_[index_in_ann] = center.x ;
                    ann_points_[index_in_ann + 1] = center.y ;
                    ann_points_[index_in_ann + 2] = center.z ;
                }
                ann_tree_->set_points( nb_vertices, ann_points_ ) ;
                break ;
            }
            case CELLS: {
                index_t nb_vertices = mesh.cells.nb() ;
                ann_points_ = new double[nb_vertices * 3] ;
                for( index_t i = 0; i < mesh.cells.nb(); i++ ) {
                    vec3 center = Geom::mesh_cell_center( mesh, i ) ;
                    index_t index_in_ann = 3 * i ;
                    ann_points_[index_in_ann] = center.x ;
                    ann_points_[index_in_ann + 1] = center.y ;
                    ann_points_[index_in_ann + 2] = center.z ;
                }
                ann_tree_->set_points( nb_vertices, ann_points_ ) ;
                break ;
            }
        }
    }

    ColocaterANN::ColocaterANN( const std::vector< vec3 >& vertices )
    {
        index_t nb_vertices = vertices.size() ;
        ann_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
        ann_points_ = new double[nb_vertices * 3] ;
        for( index_t i = 0; i < nb_vertices; i++ ) {
            index_t index_in_ann = 3 * i ;
            ann_points_[index_in_ann] = vertices[i].x ;
            ann_points_[index_in_ann + 1] = vertices[i].y ;
            ann_points_[index_in_ann + 2] = vertices[i].z ;
        }
        ann_tree_->set_points( nb_vertices, ann_points_ ) ;
    }

    ColocaterANN::ColocaterANN( float64* vertices, index_t nb_vertices )
    {
        ann_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
        ann_points_ = new double[nb_vertices * 3] ;
        for( index_t i = 0; i < nb_vertices * 3; i++ ) {
            ann_points_[i] = vertices[i] ;

        }
        ann_tree_->set_points( nb_vertices, ann_points_ ) ;
    }

    ColocaterANN::ColocaterANN( const std::vector< Edge >& edges )
    {
        index_t nb_vertices = edges.size() ;
        ann_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
        ann_points_ = new double[nb_vertices * 3] ;
        for( index_t i = 0; i < nb_vertices; i++ ) {
            vec3 barycenter( ( edges[i].value( 0 ) + edges[i].value( 1 ) ) / 2.0 ) ;
            index_t index_in_ann = 3 * i ;
            ann_points_[index_in_ann] = barycenter.x ;
            ann_points_[index_in_ann + 1] = barycenter.y ;
            ann_points_[index_in_ann + 2] = barycenter.z ;
        }
        ann_tree_->set_points( nb_vertices, ann_points_ ) ;
    }

    void ColocaterANN::set_points( const std::vector< vec3 >& vertices )
    {
        index_t nb_vertices = vertices.size() ;
        ann_tree_ = GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
        ann_points_ = new double[nb_vertices * 3] ;
        for( index_t i = 0; i < nb_vertices; i++ ) {
            index_t index_in_ann = 3 * i ;
            ann_points_[index_in_ann] = vertices[i].x ;
            ann_points_[index_in_ann + 1] = vertices[i].y ;
            ann_points_[index_in_ann + 2] = vertices[i].z ;
        }
        ann_tree_->set_points( nb_vertices, ann_points_ ) ;
    }

    /*!
     * Compute the colocated point(s) of a given point
     * @param[in] v the point to test
     * @param[out] result the colocated point indices
     * @return return true if there is at least one intersections
     */
    bool ColocaterANN::get_colocated(
        const vec3& v,
        std::vector< index_t >& result ) const
    {
        index_t nb_neighbors = std::min( index_t( 5 ), ann_tree_->nb_points() ) ;
        result.clear() ;
        std::vector< index_t > neighbors ;
        index_t cur_neighbor = 0 ;
        index_t prev_neighbor = 0 ;
        do {
            prev_neighbor = cur_neighbor ;
            cur_neighbor += nb_neighbors ;
            neighbors.resize( cur_neighbor ) ;
            double* dist = (double*) alloca( sizeof(double) * cur_neighbor ) ;
            nb_neighbors = get_neighbors( v, cur_neighbor, neighbors, dist ) ;
            for( index_t i = prev_neighbor; i < cur_neighbor; ++i ) {
                if( dist[i] > epsilon_sq ) {
                    break ;
                }
                result.push_back( neighbors[i] ) ;
            }
        } while( result.size() == cur_neighbor ) ;

        return !result.empty() ;
    }

    /*!
     * Gets the neighboring points of a given one sorted by increasing distance
     * @param[in] v the point to test
     * @param[in] nb_neighbors the number of neighbors to return
     * @param[out] result the neighboring points
     * @param[out] dist the distance between each neigbhor and the point \p v
     * @return the number of neighbors returned (can be less than \p nb_neighbors
     * if there is not enough points)
     */
    index_t ColocaterANN::get_neighbors(
        const vec3& v,
        index_t nb_neighbors,
        std::vector< index_t >& result,
        double* dist ) const
    {
        if( ann_tree_->nb_points() == 0 ) return 0 ;
        if( !dist ) {
            dist = (double*) alloca( sizeof(double) * nb_neighbors ) ;
        }
        nb_neighbors = std::min( nb_neighbors, ann_tree_->nb_points() ) ;
        result.resize( nb_neighbors ) ;
        ann_tree_->get_nearest_neighbors( nb_neighbors, v.data(), &result[0],
            dist ) ;
        return nb_neighbors ;
    }

    SortTriangleAroundEdge::TriangleToSort::TriangleToSort(
        index_t index,
        index_t surface_index,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2 )
        :
            index_( index ),
            surface_index_( surface_index ),
            N_(),
            B_A_(),
            angle_( -99999 ),
            side_( false )
    {
        ringmesh_assert( p0 != p1 ) ;
        ringmesh_assert( p0 != p2 ) ;
        ringmesh_assert( p1 != p2 ) ;

        vec3 e1 = normalize( p1 - p0 ) ;
        vec3 e2 = normalize( p2 - p0 ) ;

        N_ = normalize( cross( e1, e2 ) ) ;
        ringmesh_assert( dot( N_, e1 ) < epsilon ) ;

        vec3 B = 0.5 * p1 + 0.5 * p0 ;
        vec3 p2B = p2 - B ;
        B_A_ = normalize( p2B - dot( p2B, e1 ) * e1 ) ;

        ringmesh_assert( dot( B_A_, e1 ) < epsilon ) ;
        ringmesh_assert( B_A_.length() > epsilon ) ;
    }

    vec3 SortTriangleAroundEdge::rotate(
        const vec3& axis,
        double angle,
        const vec3& V )
    {
        vec3 q = axis ;
        if( q.length() > 0 ) {
            double s = 1.0 / q.length() ;
            q[0] *= s ;
            q[1] *= s ;
            q[2] *= s ;
        }
        q *= sinf( 0.5 * angle ) ;

        float quat[4] = { q[0], q[1], q[2], cosf( 0.5 * angle ) } ;

        double m[4][4] ;

        m[0][0] = 1 - 2.0 * ( quat[1] * quat[1] + quat[2] * quat[2] ) ;
        m[0][1] = 2.0 * ( quat[0] * quat[1] + quat[2] * quat[3] ) ;
        m[0][2] = 2.0 * ( quat[2] * quat[0] - quat[1] * quat[3] ) ;
        m[0][3] = 0.0 ;

        m[1][0] = 2.0 * ( quat[0] * quat[1] - quat[2] * quat[3] ) ;
        m[1][1] = 1 - 2.0 * ( quat[2] * quat[2] + quat[0] * quat[0] ) ;
        m[1][2] = 2.0 * ( quat[1] * quat[2] + quat[0] * quat[3] ) ;
        m[1][3] = 0.0 ;

        m[2][0] = 2.0 * ( quat[2] * quat[0] + quat[1] * quat[3] ) ;
        m[2][1] = 2.0 * ( quat[1] * quat[2] - quat[0] * quat[3] ) ;
        m[2][2] = 1 - 2.0 * ( quat[1] * quat[1] + quat[0] * quat[0] ) ;
        m[2][3] = 0.0 ;

        m[3][0] = 0.0 ;
        m[3][1] = 0.0 ;
        m[3][2] = 0.0 ;
        m[3][3] = 1.0 ;

        double x = V[0] * m[0][0] + V[1] * m[1][0] + V[2] * m[2][0] + m[3][0] ;
        double y = V[0] * m[0][1] + V[1] * m[1][1] + V[2] * m[2][1] + m[3][1] ;
        double z = V[0] * m[0][2] + V[1] * m[1][2] + V[2] * m[2][2] + m[3][2] ;
        double w = V[0] * m[0][3] + V[1] * m[1][3] + V[2] * m[2][3] + m[3][3] ;
        return vec3( x / w, y / w, z / w ) ;
    }

    void SortTriangleAroundEdge::sort()
    {
        ringmesh_assert( triangles_.size() > 0 ) ;

        std::pair< index_t, bool > default_pair( index_t( -1 ), false ) ;
        sorted_triangles_.resize( 2 * triangles_.size(), default_pair ) ;

        // If there is only one Triangle to sort - nothing to do
        if( triangles_.size() == 1 ) {
            sorted_triangles_[0] = std::pair< index_t, bool >(
                triangles_[0].surface_index_, true ) ;
            sorted_triangles_[1] = std::pair< index_t, bool >(
                triangles_[0].surface_index_, false ) ;
            return ;
        }

        // Initialization
        // We start on the plus (true) side of the first Triangle            
        sorted_triangles_[0] = std::pair< index_t, bool >(
            triangles_[0].surface_index_, true ) ;

        // Reference vectors with wich angles will be computed
        vec3 N_ref = triangles_[0].N_ ;
        vec3 B_A_ref = triangles_[0].B_A_ ;
        vec3 Ax_ref = normalize( cross( B_A_ref, N_ref ) ) ;

        // The minus (false) side of the start triangle will the last one encountered
        triangles_[0].angle_ = 2 * M_PI ;
        triangles_[0].side_ = false ;

        for( index_t i = 1; i < triangles_.size(); ++i ) {
            TriangleToSort& cur = triangles_[i] ;
            // Compute the angle RADIANS between the reference and the current
            // triangle 
            double cos = dot( B_A_ref, cur.B_A_ ) ;
            // Remove invalid values
            if( cos < -1 )
                cos = -1 ;
            else if( cos > 1 ) cos = 1 ;
            cur.angle_ = std::acos( cos ) ;
            // Put the angle between PI and 2PI if necessary
            if( dot( cross( B_A_ref, cur.B_A_ ), Ax_ref ) < 0. ) {
                cur.angle_ = 2 * M_PI - cur.angle_ ;
            }

            // Get the side of the surface first encountered
            // when rotating in the N_ref direction
            vec3 N_rotate = rotate( Ax_ref, -cur.angle_, cur.N_ ) ;
            cur.side_ = dot( N_rotate, N_ref ) > 0 ? false : true ;
        }

        // Sort the Surfaces according to the angle
        std::sort( triangles_.begin(), triangles_.end() ) ;

        // Fill the sorted surfaces adding the side
        index_t it = 1 ;
        for( index_t i = 0; i < triangles_.size(); ++i ) {
            TriangleToSort& cur = triangles_[i] ;
            if( triangles_[i].index_ == 0 ) { // The last to add
                ringmesh_assert( i == triangles_.size() - 1 ) ;
                sorted_triangles_[it].first = cur.surface_index_ ;
                sorted_triangles_[it].second = cur.side_ ;
            } else {
                sorted_triangles_[it].first = cur.surface_index_ ;
                sorted_triangles_[it].second = cur.side_ ;
                sorted_triangles_[it + 1].first = cur.surface_index_ ;
                sorted_triangles_[it + 1].second = !cur.side_ ;
                it += 2 ;
            }
        }
        // All the surfaces must have been sorted
        ringmesh_assert(
            std::count( sorted_triangles_.begin(), sorted_triangles_.end(),
                default_pair ) == 0 ) ;
    }

    const std::pair< index_t, bool >& SortTriangleAroundEdge::next(
        const std::pair< index_t, bool >& in ) const
    {
        for( index_t i = 0; i < sorted_triangles_.size(); ++i ) {
            if( sorted_triangles_[i] == in ) {
                if( i == sorted_triangles_.size() - 1 )
                    return sorted_triangles_[sorted_triangles_.size() - 2] ;
                if( i == 0 ) return sorted_triangles_[1] ;

                if( sorted_triangles_[i + 1].first == sorted_triangles_[i].first ) {
                    // The next has the same surface id, check its sign
                    if( sorted_triangles_[i + 1].second
                        != sorted_triangles_[i].second ) {
                        return sorted_triangles_[i - 1] ;
                    } else {
                        // Sign is the same
                        return sorted_triangles_[i + 1] ;
                    }
                } else {
                    ringmesh_assert(
                        sorted_triangles_[i - 1].first
                            == sorted_triangles_[i].first ) ;
                    if( sorted_triangles_[i - 1].second
                        != sorted_triangles_[i].second ) {
                        return sorted_triangles_[i + 1] ;
                    } else {
                        return sorted_triangles_[i - 1] ;
                    }
                }
            }
        }
        ringmesh_assert_not_reached;
    }


    /**********************************************************/

    // code below is copied and modified from geogram\mesh\mesh_intersection.cpp
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
    */
namespace {

    using namespace GEO;

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

    /**  // Modified  JP  
    * \brief Tests whether two facets are adjacent
    * \details Two facets are adjacents if they share an edge

    * \param[in] M the mesh
    * \param[in] f1 index of the first facet
    * \param[in] f2 index of the second facet
    * \return true if facets \p f1 and \p f2 share an edge, false
    *  otherwise
    */
    bool facets_are_adjacent( const Surface& S, index_t f1, index_t f2 )
    {
        if( f1 == f2 ) {
            return true;
        }
        for( index_t v = 0; v < S.nb_vertices_in_facet( f1 ); ++v ) {
            if( S.adjacent( f1, v ) == f2 ) {
                return true;
            }
            else if( S.adjacent( f1, v ) == NO_ID )
            {
                index_t p0 = S.model_vertex_id( f1, v ) ;
                index_t p1 = S.model_vertex_id( f1, S.next_in_facet( f1, v ) ); 
                // Check if the edge on border is not the same
                // than an edge on the border in f2
                for( index_t v2 = 0; v2 < S.nb_vertices_in_facet( f2 ); ++v2 ) {
                    if( S.adjacent( f2, v2 ) == NO_ID ) {
                        index_t p02 = S.model_vertex_id( f2, v ) ;
                        index_t p12 = S.model_vertex_id( f2, S.next_in_facet( f2, v ) );
                        if( p0 == p02 && p1 == p12 ) {
                            return true ;
                        } else if( p0 = p12 && p1 == p02 ) {
                            return true ;
                        }
                    }
                }
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
            const Surface& S, vector<index_t>& has_isect
            ) :
            S_( S ),
            has_intersection_( has_isect )
        {
            has_intersection_.assign( S.mesh().facets.nb(), 0 );
        }

        /**
        * \brief Determines the intersections between two facets
        * \details It is a callback for AABBTree traversal
        * \param[in] f1 index of the first facet
        * \param[in] f2 index of the second facet
        */
        void operator() ( index_t f1, index_t f2 )
        {
            if(
                !facets_are_adjacent( S_, f1, f2 ) &&
                f1 != f2 &&
                triangles_intersect( S_.mesh(), f1, f2, sym_ )
                ) {
                has_intersection_[ f1 ] = 1;
                has_intersection_[ f2 ] = 1;
            }


        }

    private:
        const Surface& S_;
        vector<index_t>& has_intersection_;
        vector<TriangleIsect> sym_;
    };

    
}

    /**
    * \brief Detect intersecting facets in a mesh TRIANGULATED !!
    * \param[in] M the mesh
    * \return
    */
    index_t detect_intersecting_facets( const Surface& S )
    {
        GEO::Mesh& M = S.mesh() ;
        geo_assert( M.vertices.dimension() >= 3 );

        vector<index_t> has_intersection;
        StoreIntersections action( S, has_intersection );
        MeshFacetsAABB AABB( M );
        AABB.compute_facet_bbox_intersections( action );

        return std::count( has_intersection.begin(), has_intersection.end(), 0 ) ;
    }
}
