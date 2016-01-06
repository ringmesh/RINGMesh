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
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_builder.h>
#include <ringmesh/geometry.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/tetra_gen.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/mesh/mesh_geometry.h>

#include <iostream>
#include <iomanip>

namespace {
    using namespace RINGMesh ;

    /*!
    * @brief Total number of facets in the geomodel Surfaces
    */
    index_t nb_facets( const GeoModel& geomodel )
    {
        index_t result = 0 ;
        for( index_t i = 0; i < geomodel.nb_surfaces(); ++i ) {
            result += geomodel.surface( i ).nb_cells() ;
        }
        return result ;
    }

    index_t nb_cells( const GeoModel& geomodel )
    {
        index_t nb_cells = 0 ;
        for( index_t i = 0; i < geomodel.nb_regions(); ++i ) {
            nb_cells += geomodel.region( i ).nb_cells() ;
        }
        return nb_cells ;
    }


    index_t nb_edges( const GeoModel& geomodel )
    {
        index_t nb_edges = 0 ;
        for( index_t i = 0; i < geomodel.nb_lines(); ++i ) {
            nb_edges += geomodel.line( i ).nb_cells() ;
        }
        return nb_edges ;
    }

    index_t nb_facet_elements(
        const GeoModel& geomodel,
        index_t& nb_triangles,
        index_t& nb_quads,
        index_t& nb_polygons )
    {
        for( index_t s = 0; s < geomodel.nb_surfaces(); s++ ) {
            const Surface& surface = geomodel.surface( s ) ;
            for( index_t f = 0; f < surface.nb_cells(); f++ ) {
                index_t nb_vertices = surface.nb_vertices_in_facet( f ) ;
                switch( nb_vertices ) {
                    case 3:
                        nb_triangles++ ;
                        break ;
                    case 4:
                        nb_quads++ ;
                        break ;
                    default:
                        nb_polygons++ ;
                        break ;
                }
            }
        }
        return nb_triangles + nb_quads + nb_polygons ;
    }

    index_t nb_cell_elements(
        const GeoModel& geomodel,
        index_t& nb_tet,
        index_t& nb_pyramids,
        index_t& nb_prisms,
        index_t& nb_hex,
        index_t& nb_poly )
    {
        for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
            const Region& region = geomodel.region( r ) ;
            for( index_t c = 0; c < region.nb_cells(); c++ ) {
                index_t nb_vertices = region.nb_vertices_in_cell( c ) ;
                switch( nb_vertices ) {
                    case 4:
                        nb_tet++ ;
                        break ;
                    case 5:
                        nb_pyramids++ ;
                        break ;
                    case 6:
                        nb_prisms++ ;
                        break ;
                    case 8:
                        nb_hex++ ;
                        break ;
                    default:
                        nb_poly++ ;
                        break ;
                }
            }
        }
        return nb_tet + nb_pyramids + nb_prisms + nb_hex + nb_poly ;
    }

    double cell_volume(
        const GeoModel& geomodel,
        double& tet_volume,
        double& pyramid_volume,
        double& prism_volume,
        double& hex_volume,
        double& poly_volume )
    {
        for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
            const Region& region = geomodel.region( r ) ;
            for( index_t c = 0; c < region.nb_cells(); c++ ) {
                index_t nb_vertices = region.nb_vertices_in_cell( c ) ;
                switch( nb_vertices ) {
                    case 4:
                        tet_volume += RINGMesh::mesh_cell_volume( region.mesh(), c ) ;
                        break ;
                    case 5:
                        pyramid_volume += RINGMesh::mesh_cell_volume( region.mesh(), c ) ;
                        break ;
                    case 6:
                        prism_volume += RINGMesh::mesh_cell_volume( region.mesh(), c ) ;
                        break ;
                    case 8:
                        hex_volume += RINGMesh::mesh_cell_volume( region.mesh(), c ) ;
                        break ;
                    default:
                        poly_volume += RINGMesh::mesh_cell_volume( region.mesh(), c ) ;
                        break ;
                }
            }
        }
        return tet_volume + pyramid_volume + prism_volume + hex_volume + poly_volume ;
    }

    void print_nb_cell_stat(
        const index_t& nb_cell,
        const index_t& nb_cell_total,
        const std::string& cell_type )
    {
        GEO::Logger::out( "GeoModel" ) << "* " << nb_cell << " " << cell_type << " ("
            << nb_cell * 100 / nb_cell_total << "%)\n" ;
    }

    void print_cell_volume_stat(
        const double& cell_volume,
        const double& cell_volume_total,
        const std::string& cell_type )
    {
        GEO::Logger::out( "GeoModel" )
            << "* " << cell_type << " volume " << cell_volume << " ("
            << static_cast< index_t >( cell_volume * 100 / cell_volume_total + 0.5 ) << "%)\n" ;
    }

}

namespace RINGMesh {
    void print_geomodel( const GeoModel& geomodel )
    {
        GEO::Logger::out( "GeoModel" ) << "Model " << geomodel.name() << " has\n"
            << std::setw( 10 ) << std::left
            << geomodel.mesh.vertices.nb() << " vertices\n"
            << std::setw( 10 ) << std::left
            << nb_facets( geomodel ) << " facets\n"
            << std::endl ;

        for( index_t t = GME::CORNER; t < GME::NO_TYPE; ++t ) {
            GME::TYPE T = static_cast<GME::TYPE>( t ) ;
            GEO::Logger::out( "GeoModel" ) << std::setw( 10 ) << std::left
                << geomodel.nb_elements( T ) << " " << GME::type_name( T )
                << std::endl ;
        }
    }

    void print_geomodel_nb_elements( const GeoModel& geomodel )
    {
        GEO::Logger::out( "GeoModel" ) << "Model " << geomodel.name() << " is made of\n"
            << std::setw( 10 ) << std::left
            << geomodel.mesh.vertices.nb() << " vertices\n"
            << std::setw( 10 ) << std::left
            << nb_edges( geomodel ) << " edges\n" ;

        index_t nb_triangles = 0 ;
        index_t nb_quads = 0 ;
        index_t nb_polygons = 0 ;
        index_t nb_facets =  nb_facet_elements( geomodel, nb_triangles, nb_quads, nb_polygons ) ;
        GEO::Logger::out( "GeoModel" ) << std::setw( 10 ) << std::left
            << nb_facets << " facets\n" ;
        if( nb_triangles > 0 ) {
            print_nb_cell_stat( nb_triangles, nb_facets, "triangles" ) ;
        }
        if( nb_quads > 0 ) {
            print_nb_cell_stat( nb_quads, nb_facets, "quads" ) ;
        }
        if( nb_polygons > 0 ) {
            print_nb_cell_stat( nb_polygons, nb_facets, "polygons" ) ;
        }

        index_t nb_tet = 0 ;
        index_t nb_pyramids = 0 ;
        index_t nb_prisms = 0 ;
        index_t nb_hex = 0 ;
        index_t nb_poly = 0 ;
        index_t nb_cells = nb_cell_elements( geomodel, nb_tet, nb_pyramids,
            nb_prisms, nb_hex, nb_poly ) ;
        GEO::Logger::out( "GeoModel" ) << std::setw( 10 ) << std::left
            << nb_cells << " cells\n" ;
        if( nb_tet > 0 ) {
            print_nb_cell_stat( nb_tet, nb_cells, "tet" ) ;
        }
        if( nb_pyramids > 0 ) {
            print_nb_cell_stat( nb_pyramids, nb_cells, "pyramids" ) ;
        }
        if( nb_prisms > 0 ) {
            print_nb_cell_stat( nb_prisms, nb_cells, "prisms" ) ;
        }
        if( nb_hex > 0 ) {
            print_nb_cell_stat( nb_hex, nb_cells, "hex" ) ;
        }
        if( nb_poly > 0 ) {
            print_nb_cell_stat( nb_poly, nb_cells, "polyhedra" ) ;
        }
        GEO::Logger::out( "GeoModel" ) << std::endl ;
    }

    void print_geomodel_volume( const GeoModel& geomodel )
    {
        double tet_volume = 0;
        double pyramid_volume = 0;
        double prism_volume = 0;
        double hex_volume = 0;
        double poly_volume = 0 ;
        double volume = cell_volume( geomodel, tet_volume, pyramid_volume,
            prism_volume, hex_volume, poly_volume ) ;
        GEO::Logger::out( "GeoModel" ) << "Model " << geomodel.name()
            << " has a volume of " << volume << "\n" ;
        if( tet_volume > 0 ) {
            print_cell_volume_stat( tet_volume, volume, "tet" ) ;
        }
        if( pyramid_volume > 0 ) {
            print_cell_volume_stat( pyramid_volume, volume, "pyramid" ) ;
        }
        if( prism_volume > 0 ) {
            print_cell_volume_stat( prism_volume, volume, "prism" ) ;
        }
        if( hex_volume > 0 ) {
            print_cell_volume_stat( hex_volume, volume, "hex" ) ;
        }
        if( poly_volume > 0 ) {
            print_cell_volume_stat( poly_volume, volume, "polyhedron" ) ;
        }
        GEO::Logger::out( "GeoModel" ) << std::endl ;
    }

    bool are_geomodel_surface_meshes_simplicial( const GeoModel& geomodel )
    {
        for( index_t i = 0; i != geomodel.nb_surfaces(); ++i ) {
            if( !geomodel.surface( i ).is_simplicial() ) {
                return false ;
            }
        }
        return true ;
    }

    bool are_geomodel_region_meshes_simplicial( const GeoModel& geomodel )
    {
        for( index_t i = 0; i != geomodel.nb_regions(); ++i ) {
            if( !geomodel.region( i ).is_simplicial() ) {
                return false ;
            }
        }
        return true ;
    }


	///@todo A class encapsulating the copy from a GeoModel to a Mesh ?
    /// See what has been done in GeoModelMeshBuilder	

    void add_geomodel_vertices_to_mesh( const GeoModel& geomodel, GEO::Mesh& M )
    {
        index_t nbv = geomodel.mesh.vertices.nb() ;
        M.vertices.create_vertices( nbv ) ;

        // We need to copy the point one after another since we do not have access
        // to the storage of the geomodel.vertices
        // I do not want to provide this access [JP]
        for( index_t v = 0; v < nbv; ++v ) {
            M.vertices.point( v ) = geomodel.mesh.vertices.vertex( v ) ;
        }

        GEO::Attribute<index_t> corner_attribute( M.vertices.attributes(), "region" ) ;
        corner_attribute.fill( NO_ID ) ;
        for( index_t i = 0; i < geomodel.nb_corners(); ++i ) {
            index_t vertex_index = geomodel.corner( i ).model_vertex_id() ;
            corner_attribute[ vertex_index ] = i ;
        }
        corner_attribute.unbind() ;
    }

    void add_line_edges_to_mesh( const Line& line, GEO::Mesh& M )
    {
        index_t from = M.edges.create_edges( line.nb_cells() ) ;
        for( index_t i = 0; i < line.nb_cells(); ++i ) {
            index_t v0 = line.model_vertex_id( i, 0 ) ;
            index_t v1 = line.model_vertex_id( i, 1 ) ;
            M.edges.set_vertex( from + i, 0, v0 ) ;
            M.edges.set_vertex( from + i, 1, v1  );
        }
    }

    void create_and_fill_line_index_attribute( const GeoModel& geomodel,
                                               const std::string& attribute_name,
                                               GEO::Mesh& M )
    {
        GEO::Attribute<index_t> line_attribute( M.edges.attributes(), "region" ) ;
        line_attribute.fill( NO_ID ) ;
        index_t edge_counter = 0 ;
        for( index_t i = 0; i < geomodel.nb_lines(); ++i ) {
            index_t nb_line_edges = geomodel.line( i ).nb_cells() ;
            index_t line_edges_start = edge_counter ;
            index_t line_edges_end = edge_counter + nb_line_edges ;

            for( index_t e = line_edges_start; e != line_edges_end; ++e ) {
                line_attribute[ e ] = i ;
            }
            edge_counter += nb_line_edges ;
        }
        line_attribute.unbind() ;
    }

    void add_geomodel_line_edges_to_mesh( const GeoModel& geomodel, GEO::Mesh& M )
    {
        for( index_t i = 0; i < geomodel.nb_lines(); ++i ) {
            const Line& line( geomodel.line( i ) ) ;
            add_line_edges_to_mesh( line, M ) ;
        }
        create_and_fill_line_index_attribute( geomodel, "region", M ) ;
    }

    void add_surface_facets_to_mesh( const Surface& surface, GEO::Mesh& M )
    {
        for( index_t j = 0; j < surface.nb_cells(); ++j ) {
            index_t nbv = surface.nb_vertices_in_facet( j ) ;
            GEO::vector< index_t > ids( nbv ) ;
            for( index_t v = 0; v < nbv; ++v ) {
                ids[ v ] = surface.model_vertex_id( j, v ) ;
            }
            M.facets.create_polygon( ids ) ;
        }
    }

    void add_surfaces_triangles_to_mesh( const GeoModel& geomodel, GEO::Mesh& M )
    {
        GEO::vector< index_t > triangles( 3 * nb_facets( geomodel ) ) ;
        index_t triangle_index = 0 ;
        for( index_t i = 0; i < geomodel.nb_surfaces(); ++i ) {
            const Surface& S = geomodel.surface( i ) ;         
            index_t nb_surface_triangles = S.nb_cells() ;
            for( index_t j = 0 ; j != nb_surface_triangles; ++j ) {
                triangles[ 3*triangle_index ] = S.model_vertex_id( j, 0 ) ;
                triangles[ 3*triangle_index + 1 ] = S.model_vertex_id( j, 1 ) ;
                triangles[ 3*triangle_index + 2 ] = S.model_vertex_id( j, 2 );
                ++triangle_index ;
            }
        }
        M.facets.assign_triangle_mesh( triangles, true ) ;
    }

    void create_and_fill_surface_index_attribute( const GeoModel& geomodel,
                                                  const std::string& attribute_name,
                                                  GEO::Mesh& M )
    {
        GEO::Attribute<index_t> surface_attribute( M.facets.attributes(), attribute_name ) ;
        surface_attribute.fill( NO_ID ) ;
        index_t facet_counter = 0 ;
        for( index_t i = 0; i < geomodel.nb_surfaces(); ++i ) {
            index_t nb_surface_facets = geomodel.surface( i ).nb_cells() ;
            index_t surface_facet_start = facet_counter ;
            index_t surface_facet_end = facet_counter + nb_surface_facets ;

            for( index_t f = surface_facet_start; f != surface_facet_end; ++f ) {
                surface_attribute[ f ] = i ;
            }
            facet_counter += nb_surface_facets ;
        }
        surface_attribute.unbind() ;
    }

    void add_geomodel_surface_facets_to_mesh( const GeoModel& geomodel, GEO::Mesh& M )
    {
        if( are_geomodel_surface_meshes_simplicial( geomodel ) ) {
            add_surfaces_triangles_to_mesh( geomodel, M ) ;
        } else {
            for( index_t i = 0; i < geomodel.nb_surfaces(); ++i ) {
                const Surface& S = geomodel.surface( i ) ; 
                add_surface_facets_to_mesh( S, M ) ;
            }
        }
        create_and_fill_surface_index_attribute( geomodel, "region", M ) ;
    }

    void add_regions_tets_to_mesh( const GeoModel& geomodel, GEO::Mesh& M )
    {
        GEO::vector< index_t > tets( 4*nb_cells( geomodel ) ) ;
        index_t tet_index = 0 ;
        for( index_t i = 0; i < geomodel.nb_regions(); ++i ) {
            const Region& region = geomodel.region( i ) ;
            index_t nb_region_tets = region.nb_cells() ;
            for( index_t j = 0; j < nb_region_tets; ++j ) {
                tets[ 4*tet_index ] = region.model_vertex_id( j, 0 ) ;
                tets[ 4*tet_index + 1 ] = region.model_vertex_id( j, 1 ) ;
                tets[ 4*tet_index + 2 ] = region.model_vertex_id( j, 2 ) ;
                tets[ 4*tet_index + 3 ] = region.model_vertex_id( j, 3 ) ;
                ++tet_index ;
            }
        }
        M.cells.assign_tet_mesh( tets, true ) ;
    }

    void create_and_fill_region_index_attribute( const GeoModel& geomodel,
                                                 const std::string& attribute_name,
                                                 GEO::Mesh& M )
    {
        GEO::Attribute<index_t> region_attribute( M.cells.attributes(), attribute_name ) ;
        region_attribute.fill( NO_ID ) ;
        index_t cell_counter = 0 ;
        for( index_t i = 0; i < geomodel.nb_regions(); ++i ) {
            index_t nb_region_cells = geomodel.region( i ).nb_cells() ;
            index_t region_cells_start = cell_counter ;
            index_t region_cells_end = cell_counter + nb_region_cells ;

            for( index_t j = region_cells_start; j != region_cells_end; ++j ) {
                region_attribute[ j ] = i ;
            }
            cell_counter += nb_region_cells ;
        }
        region_attribute.unbind() ;
    }

    /*!
    * @pre Regions meshes are all tetrahedral !
    * @todo to implement for other types of cells
    */
    void add_geomodel_region_tets_to_mesh( const GeoModel& geomodel, GEO::Mesh& M )
    {
        if( are_geomodel_region_meshes_simplicial( geomodel ) ) {
            add_regions_tets_to_mesh( geomodel, M ) ;
        } else {
            ringmesh_assert( false ) ;
        }
        create_and_fill_region_index_attribute( geomodel, "region", M ) ;
    }


    void build_mesh_from_geomodel( const GeoModel& geomodel, GEO::Mesh& M )
    {
        // Keep the attributes when clearing the mesh, otherwise we crash
        M.clear( true ) ;

        add_geomodel_vertices_to_mesh( geomodel, M ) ;
        add_geomodel_line_edges_to_mesh( geomodel, M ) ;
        add_geomodel_surface_facets_to_mesh( geomodel, M ) ;
        add_geomodel_region_tets_to_mesh( geomodel, M ) ;
    }


    double model_element_size( const GeoModelElement& E )
    {
        double result = 0. ;
        if( E.nb_children() ) {
            // Sum up the size of children elements
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                result += model_element_size( E.child( i ) ) ;
            }
            return result ;
        } else {
            switch( E.type() ) {
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

        switch( E.type() ) {
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
                return GEO::Geom::distance( L.vertex( c, 0 ), L.vertex( c, 1 ) ) ;
            }
        }
        ringmesh_assert_not_reached ;
        return result ;
    }

    vec3 model_element_center( const GeoModelElement& E )
    {
        vec3 result( 0., 0., 0. ) ;
        index_t nb_vertices = 0 ;

        if( GeoModelElement::has_mesh( E.type() ) ) {
            // @todo Improve efficiency, overload the functions to avoid
            // casting each time
            const GeoModelMeshElement& M =
                dynamic_cast< const GeoModelMeshElement& >( E ) ;
            for( index_t v = 0; v < M.nb_vertices(); v++ ) {
                result += M.vertex( v ) ;
            }
            return result / static_cast< double >( M.nb_vertices() ) ;
        }
        else if( E.nb_children() > 0 ) {
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                const GeoModelMeshElement& F =
                    dynamic_cast< const GeoModelMeshElement& >( E.child( i ) ) ;
                nb_vertices += F.nb_vertices() ;
                result += model_element_center( F ) * F.nb_vertices() ;
            }
            return result / static_cast< double >( nb_vertices ) ;
        } else {
            return result ;
        }
    }

    vec3 model_element_cell_center( const GeoModelMeshElement& E, index_t c )
    {
        vec3 result( 0., 0., 0. ) ;       
        const GEO::Mesh& mesh = E.mesh() ;
        switch( E.type() ) {
            case GeoModelElement::REGION: {                
                return RINGMesh::mesh_cell_center( mesh, c ) ;
            }
            case GeoModelElement::SURFACE: {             
                return GEO::Geom::mesh_facet_center( mesh, c ) ;
            }
            case GeoModelElement::LINE: {
                index_t v0 = mesh.edges.vertex( c, 0 ) ;
                index_t v1 = mesh.edges.vertex( c, 1 ) ;
                return 0.5 * ( mesh.vertices.point( v0 ) + mesh.vertices.point( v1 ) ) ;
            }
            case GeoModelElement::CORNER: {
                return mesh.vertices.point(0) ;
            }
        }
        ringmesh_assert_not_reached ;
        return result ;
    }

    void translate( GeoModel& M, const vec3& translation_vector )
    {
        for( index_t v = 0; v < M.mesh.vertices.nb(); ++v ) {
            // Coordinates are not directly modified to 
            // update the matching vertices in geomodel entities
            vec3 p = M.mesh.vertices.vertex( v ) ;
            for( index_t i = 0; i < 3; i++ ) {
                p[i] += translation_vector[i] ;
            }
            M.mesh.vertices.update_point( v, p ) ;
        }
    }

    void rotate( GeoModel& M, const vec3& origin, const vec3& axis, double theta, bool degrees )
    {
        if( length( axis ) < epsilon ) {
            GEO::Logger::err( "GeoModel" )
                << "Rotation around an epsilon length axis is impossible"
                << std::endl ;
            return ;
        }
        GEO::Matrix< double, 4 > rot_mat ;
        rotation_matrix_about_arbitrary_axis( origin, axis, theta, degrees, rot_mat ) ;

        for( index_t v = 0; v < M.mesh.vertices.nb(); ++v ) {
            const vec3& p = M.mesh.vertices.vertex( v ) ;
            double old[4] = { p[0], p[1], p[2], 1. } ;
            double new_p[4] = { 0, 0, 0, 1. } ;
            GEO::mult( rot_mat, old, new_p ) ;
            /*! @todo You need an epsilon tolerance here [JP] */
            ringmesh_debug_assert( new_p[3] == 1. ) ;

            M.mesh.vertices.update_point( v, vec3( new_p[0], new_p[1], new_p[2] ) ) ;
        }
    }



    /*!
     * @brief Generate a point that lies strictly a Region defined by its Surface boundaries.
     * @details Returnsthe midpoint of barycenter of the first facet of the first surface on 
     * the region boundary and the closest point in the other surfaces.
     * @warning Incomplete implementation.
     */
    vec3 generate_point_in_region( const Region& region )
    {
        /// @todo To implement for bubbles
        ringmesh_assert( region.nb_boundaries() > 1 ) ;
        
        const GeoModel& geomodel = region.model() ;

        const Surface& first_boundary_surface = geomodel.surface( region.boundary_gme( 0 ).index ) ; 
        double facet_area = first_boundary_surface.facet_area( 0 ) ; 
        vec3 barycenter = first_boundary_surface.facet_barycenter( 0 ) ;                
        /// @todo Check that this is the right condition to have a correct enough barycenter
        ringmesh_assert( facet_area > epsilon) ;

        double minimum_distance = DBL_MAX ;
        vec3 nearest_point ;        
        for( index_t i = 1; i != region.nb_boundaries(); ++i ) {
            const Surface& S = geomodel.surface( region.boundary_gme(i).index ) ;                        
            SurfaceTools tool_on_surface( S ) ;
            double distance = DBL_MAX ;
            vec3 point ; 
            tool_on_surface.aabb().nearest_facet( barycenter, point, distance ) ;

            if( distance < minimum_distance) {
                minimum_distance = distance ;
                nearest_point = point ;
            }            
        } 
        /// @todo Change implementation to use second triangle if that one failed, and futher surfaces
        ringmesh_assert( minimum_distance > epsilon ) ;
        return 0.5*( barycenter + nearest_point ) ;
    }

    void get_one_point_per_geomodel_regions( const GeoModel& geomodel, 
                                             std::vector< vec3 >& one_point_one_region )
    {
        one_point_one_region.resize( geomodel.nb_regions() ) ;
        for( index_t i = 0; i != geomodel.nb_regions(); ++i ) {
            vec3 point = generate_point_in_region( geomodel.region(i) ) ;
            one_point_one_region[i] = point ; 
        }    
    }

#ifdef RINGMESH_WITH_TETGEN

    void tetgen_tetrahedralize_geomodel_regions( GeoModel& geomodel ) 
    {
        GEO::Mesh mesh ;
        build_mesh_from_geomodel( geomodel, mesh ) ;
        
        std::vector< vec3 > points_in_regions ;
        get_one_point_per_geomodel_regions( geomodel, points_in_regions ) ;
       
        TetgenMesher mesher ;
        mesher.tetrahedralize( mesh, points_in_regions, "QpO0YA", mesh ) ; 

        GeoModelBuilderMesh builder ( geomodel, mesh, "", "region" ) ;
        builder.build_regions() ;

        // Force recomputation of global mesh vertices - otherwise we crash sooner or later
        // because of model_vertex_id crazy sharing [JP]
        geomodel.mesh.vertices.clear() ;
        geomodel.mesh.vertices.test_and_initialize() ;
    }
#endif

    void tetrahedralize( GeoModel& M, const std::string& method, index_t region_id, bool add_steiner_points )
    {
        /* @todo Review: Maybe rethink these functions
         *       to have a function that can mesh a region of a geomodel
         *       taking only one vector of points [JP]
         */
        std::vector< std::vector< vec3 > > internal_vertices( M.nb_regions() ) ;
        tetrahedralize( M, method, region_id, add_steiner_points,
            internal_vertices ) ;
    }

    void tetrahedralize( GeoModel& M, const std::string& method, index_t region_id, bool add_steiner_points,
                         const std::vector< std::vector< vec3 > >& internal_vertices )
    {
        if( region_id == NO_ID ) {
            GEO::Logger::out( "Info" ) << "Using " << method << std::endl ;
            GEO::ProgressTask progress( "Compute", M.nb_regions() ) ;
            for( index_t i = 0; i < M.nb_regions(); i++ ) {
                tetrahedralize( M, method, i, add_steiner_points, internal_vertices ) ;                
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

        // The GeoModelMesh should be updated, just erase everything
        // and it will be re-computed during its next access.
        M.mesh.vertices.clear() ;
    }

}
