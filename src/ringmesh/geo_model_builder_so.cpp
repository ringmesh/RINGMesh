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

#include <ringmesh/geo_model_builder_so.h>

#include <iostream>
#include <iomanip>

#include <geogram/basic/logger.h>

#include <ringmesh/geometry.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/utils.h>

namespace RINGMesh {

    bool GeoModelBuilderTSolid::load_file()
    {
        time_t t1, t2, t3, t4, t5, t6, t7, t8 ;
        time ( &t1 ) ;
        GME::gme_t cur_region ;

        // First : count the number of vertex and tetras
        // in each region
        std::vector< index_t > nb_elements_per_region ;
        read_number_of_mesh_elements( nb_elements_per_region ) ;
        print_number_of_mesh_elements( nb_elements_per_region ) ;

        // Region vertices
        std::vector< vec3 > region_vertices ;
        // Region tetraedron corners
        std::vector< index_t > tetra_corners ;
        // Vector which maps the indices of vertices from Gocad .so file
        // to the local (in region) indices of vertices
        std::vector< index_t > gocad_vertices2region_vertices ;
        // Vector which maps the indices of vertices from Gocad .so file
        // to the index of the region they belong to
        std::vector< index_t > gocad_vertices2region_id ;

        GME::gme_t current_interface ;
        GME::gme_t current_surface ;
        std::vector< index_t > cur_surf_facets ;
        std::vector< index_t > cur_surf_facet_ptr (1, 0) ;

        index_t nb_vertex_properties = 0 ;
        index_t nb_cell_properties = 0 ;
        std::vector< std::string > vertex_property_names ;
        std::vector< std::string > cell_property_names ;

        time ( &t2 ) ;
        std::cerr << "t2 - t1 : " << t2 - t1 << std::endl ;

        // Then : Reading .so file
        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if( in_.nb_fields() > 0 ) {
                if( in_.field_matches( 0, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                    read_and_set_gocad_coordinates_system() ;
                } else if( in_.field_matches( 0, "PROPERTIES" ) ) {
                    nb_vertex_properties = in_.nb_fields() - 1 ;
                } else if( in_.field_matches( 0, "PROPERTY_CLASS_HEADER" ) ) {
                    add_new_property( vertex_property_names,
                        model_.mesh.vertex_attribute_manager() ) ;
                } else if( in_.field_matches( 0, "TETRA_PROPERTIES" ) ) {
                    nb_cell_properties = in_.nb_fields() - 1 ;
                } else if( in_.field_matches( 0, "TETRA_PROPERTY_CLASS_HEADER" ) ) {
                    add_new_property( cell_property_names,
                        model_.mesh.cell_attribute_manager() ) ;
                } else if( in_.field_matches( 0, "TVOLUME" ) ) {
                    if( region_vertices.size() > 0 ) {
                        set_region_geometry( cur_region.index,
                            region_vertices,
                            tetra_corners ) ;
                        region_vertices.clear() ;
                        tetra_corners.clear() ;
                    }
                    cur_region = create_region() ;
                } else if( in_.field_matches( 0, "VRTX" ) ||
                    in_.field_matches( 0, "PVRTX" ) ) {
                    gocad_vertices2region_vertices.push_back(
                        region_vertices.size() ) ;
                    gocad_vertices2region_id.push_back( cur_region.index ) ;
                    vec3 vertex ;
                    read_vertex_coordinates( vertex ) ;
                    region_vertices.push_back( vertex ) ;
                } else if( in_.field_matches( 0, "ATOM" ) ||
                    in_.field_matches( 0, "PATOM" ) ) {
                    const index_t referring_vertex = in_.field_as_double( 2 ) - 1 ;
                    index_t referred_vertex_local_id =
                            gocad_vertices2region_vertices[referring_vertex] ;
                    index_t referred_vertex_region_id =
                            gocad_vertices2region_id[referring_vertex] ;

                    if( referred_vertex_region_id < cur_region.index ) {
                        // If the atom referred to a vertex of another region,
                        // acting like for a vertex
                        gocad_vertices2region_vertices.push_back(
                            region_vertices.size() );
                        gocad_vertices2region_id.push_back( cur_region.index ) ;
                        region_vertices.push_back(
                                model_.region( referred_vertex_region_id ).vertex(
                                        referred_vertex_local_id ) ) ;
                    } else {
                        // If the atom referred to an atom of the same region
                        gocad_vertices2region_vertices.push_back(
                            referred_vertex_local_id ) ;
                        gocad_vertices2region_id.push_back(
                            referred_vertex_region_id ) ;
                    }
                } else if( in_.field_matches( 0, "TETRA" ) ) {
                    // Reading and create a tetra
                    std::vector< index_t > corners(4) ;
                    read_tetraedra( gocad_vertices2region_vertices, corners ) ;
                    tetra_corners.insert( tetra_corners.end(),
                        corners.begin(),
                        corners.end() ) ;
                } else if( in_.field_matches( 0, "name:" ) ) {
                    // GeoModel name is set to the TSolid name.
                    set_model_name( in_.field( 1 ) ) ;
                } else if( in_.field_matches( 0, "MODEL" ) ) {
                    if( region_vertices.size() > 0 ) {
                        set_region_geometry ( cur_region.index,
                            region_vertices,
                            tetra_corners ) ;
                        region_vertices.clear() ;
                        tetra_corners.clear() ;
                    }
                } else if( in_.field_matches( 0, "SURFACE" ) ) {
                    current_interface = create_element( GME::INTERFACE ) ;
                    set_element_name( current_interface, in_.field( 1 ) ) ;
                } else if( in_.field_matches( 0, "TFACE" ) ) {
                    // Compute the surface
                    if( cur_surf_facets.size() > 0 ) {
                        build_surface( current_surface.index,
                                cur_surf_facets,
                                cur_surf_facet_ptr,
                                gocad_vertices2region_id,
                                gocad_vertices2region_vertices ) ;
                    }
                    // Create a new surface
                    current_surface = create_element( GME::SURFACE ) ;
                    set_element_parent( current_surface, current_interface ) ;
                    add_element_child( current_interface, current_surface ) ;
                } else if( in_.field_matches( 0, "TRGL" ) ) {
                    cur_surf_facets.push_back( in_.field_as_uint( 1 ) - 1 ) ;
                    cur_surf_facets.push_back( in_.field_as_uint( 2 ) - 1 ) ;
                    cur_surf_facets.push_back( in_.field_as_uint( 3 ) - 1 ) ;
                    cur_surf_facet_ptr.push_back( cur_surf_facets.size() ) ;
                } else if( in_.field_matches( 0, "MODEL_REGION" ) ) {
                    // Compute the last surface
                    if( cur_surf_facets.size() > 0 ) {
                        build_surface( current_surface.index,
                                cur_surf_facets,
                                cur_surf_facet_ptr,
                                gocad_vertices2region_id,
                                gocad_vertices2region_vertices ) ;
                    }
                }
            }
        }

        time ( &t3 ) ;
        std::cerr << "t3 - t2 : " << t3 - t2 << std::endl ;

        // Compute internal borders (by removing adjacencies on
        // triangle edges common to at least two surfaces)
        compute_internal_borders() ;

        time ( &t4 ) ;
        std::cerr << "t4 - t3 : " << t4 - t3 << std::endl ;

        // Build GeoModel Lines and Corners from the surfaces
        model_.mesh.vertices.test_and_initialize() ;
        build_lines_and_corners_from_surfaces() ;

        time ( &t5 ) ;
        std::cerr << "t5 - t4 : " << t5 - t4 << std::endl ;

        // Regions boundaries
        compute_boundaries_of_geomodel_regions() ;

        time ( &t6 ) ;
        std::cerr << "t6 - t5 : " << t6 - t5 << std::endl ;

        // Universe boundaries
        compute_universe_boundaries() ;

        time ( &t7 ) ;
        std::cerr << "t7 - t6 : " << t7 - t6 << std::endl ;

        // Contacts building
        build_contacts() ;

        time ( &t8 ) ;
        std::cerr << "t8 - t7 : " << t8 - t7 << std::endl ;

        std::cerr << "Total : " << t8 - t1 << std::endl ;

        return true ;

    }

    void GeoModelBuilderTSolid::read_and_set_gocad_coordinates_system()
    {
        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if( in_.field_matches( 0, "END_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                return ;
            } else if( in_.field_matches( 0, "NAME" ) ) {
                // Useless for the moment
            } else if( in_.field_matches( 0, "PROJECTION" ) ) {
                // Useless for the moment
            } else if( in_.field_matches( 0, "DATUM" ) ) {
                // Useless for the moment
            } else if( in_.field_matches( 0, "AXIS_NAME" ) ) {
                set_gocad_coordinates_system_axis_name() ;
            } else if( in_.field_matches( 0, "AXIS_UNIT" ) ) {
                set_gocad_coordinates_system_axis_unit() ;
            } else if( in_.field_matches( 0, "ZPOSITIVE" ) ) {
                set_gocad_coordinates_system_z_sign() ;
            }
        }
    }

    void GeoModelBuilderTSolid::set_gocad_coordinates_system_axis_name()
    {
        gocad_coordinates_system_axis_name_.push_back( in_.field(1) ) ;
        gocad_coordinates_system_axis_name_.push_back( in_.field(2) ) ;
        gocad_coordinates_system_axis_name_.push_back( in_.field(3) ) ;
    }

    void GeoModelBuilderTSolid::set_gocad_coordinates_system_axis_unit()
    {
        gocad_coordinates_system_axis_unit_.push_back( in_.field(1) ) ;
        gocad_coordinates_system_axis_unit_.push_back( in_.field(2) ) ;
        gocad_coordinates_system_axis_unit_.push_back( in_.field(3) ) ;
    }

    void GeoModelBuilderTSolid::set_gocad_coordinates_system_z_sign()
    {
        if( in_.field_matches( 1, "Elevation" ) ) {
            z_sign_ = 1 ;
        } else if( in_.field_matches( 1, "Depth" ) ) {
            z_sign_ = -1 ;
        } else {
            ringmesh_assert_not_reached ;
        }
    }

    void GeoModelBuilderTSolid::read_number_of_mesh_elements(
            std::vector< index_t >& nb_elements_par_region )
    {
        nb_elements_par_region.clear() ;

        // Define a new LineInput counting number of elements
        GEO::LineInput line_input( filename_ ) ;

        // Initialize counters
        index_t cur_region = NO_ID ;
        index_t nb_vertices_in_region = 0 ;
        index_t nb_tetras_in_region = 0 ;
        index_t nb_surfaces_in_bmodel = 0 ;
        index_t nb_triangles_in_bmodel = 0 ;

        // Reading file
        while( !line_input.eof() && line_input.get_line() ) {
            line_input.get_fields() ;
            if( line_input.nb_fields() > 0 ) {
                if( line_input.field_matches( 0, "TVOLUME" ) ||
                        line_input.field_matches( 0, "MODEL" ) ) {
                    if( cur_region != NO_ID ) {
                        nb_elements_par_region.push_back( nb_vertices_in_region ) ;
                        nb_elements_par_region.push_back( nb_tetras_in_region ) ;
                        nb_vertices_in_region = 0 ;
                        nb_tetras_in_region = 0 ;
                    }
                    ++cur_region ;
                } else if( line_input.field_matches( 0, "VRTX" ) ||
                           line_input.field_matches( 0, "PVRTX" ) ||
                           line_input.field_matches( 0, "ATOM" ) ||
                           line_input.field_matches( 0, "PATOM" ) ) {
                    ++nb_vertices_in_region ;
                } else if( line_input.field_matches( 0, "TETRA" ) ) {
                    ++nb_tetras_in_region ;
                }
            }
        }
    }

    void GeoModelBuilderTSolid::print_number_of_mesh_elements(
            const std::vector< index_t >& nb_elements_per_region ) const
    {
        const index_t nb_regions = 0.5 * nb_elements_per_region.size() ;
        GEO::Logger::out( "Mesh" )
            << "Mesh has " << nb_regions << " regions "
            << std::endl ;
        for( index_t i = 0 ; i < nb_regions ; ++i ) {
            GEO::Logger::out( "Mesh" )
                << "Region " << i << " has"
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements_per_region.at( 2*i ) << " vertices "
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements_per_region.at( 2*i + 1 ) << " tetras "
                << std::endl ;
        }
    }

    void GeoModelBuilderTSolid::add_new_property(
            std::vector< std::string >& property_names,
            GEO::AttributesManager& attribute_manager )
    {
        property_names.push_back( in_.field(1) ) ;
        /// @todo All the property types are double.
        /// Change in order to have the good type for each property.
        GEO::Attribute< double > property( attribute_manager, in_.field(1) ) ;
    }

    GME::gme_t GeoModelBuilderTSolid::create_region()
    {
        GME::gme_t cur_region = create_element( GME::REGION ) ;
        set_element_name( cur_region, in_.field( 1 ) ) ;
        return cur_region ;
    }

    void GeoModelBuilderTSolid::read_vertex_coordinates( vec3& vertex )
    {
        vertex.x = in_.field_as_double( 2 ) ;
        vertex.y = in_.field_as_double( 3 ) ;
        vertex.z = z_sign_ * in_.field_as_double( 4 ) ;
    }

    void GeoModelBuilderTSolid::read_tetraedra(
            const std::vector< index_t >& gocad_vertices2region_vertices,
            std::vector< index_t >& corners_id )
    {
        ringmesh_debug_assert( corners_id.size() == 4 ) ;
        corners_id[0] =
            gocad_vertices2region_vertices[ in_.field_as_uint( 1 ) - 1 ] ;
        corners_id[1] =
            gocad_vertices2region_vertices[ in_.field_as_uint( 2 ) - 1 ] ;
        corners_id[2] =
            gocad_vertices2region_vertices[ in_.field_as_uint( 3 ) - 1 ] ;
        corners_id[3] =
            gocad_vertices2region_vertices[ in_.field_as_uint( 4 ) - 1 ] ;
    }

    void GeoModelBuilderTSolid::compute_boundaries_of_geomodel_regions()
    {
        std::vector< ColocaterANN* > reg_anns( model_.nb_regions(), nil ) ;
        for( index_t r = 0 ; r < model_.nb_regions() ; ++r ) {
            const Region& reg = model_.region( r ) ;
            const index_t nb_cells = reg.nb_cells() ;
            std::vector< vec3 > cell_facet_centers ;
            cell_facet_centers.reserve( 4 * nb_cells ) ;
            for( index_t c = 0 ; c < nb_cells ; ++c ) {
                cell_facet_centers.push_back(
                    mesh_cell_facet_center( reg.mesh(), c, 0 ) ) ;
                cell_facet_centers.push_back(
                    mesh_cell_facet_center( reg.mesh(), c, 1 ) ) ;
                cell_facet_centers.push_back(
                    mesh_cell_facet_center( reg.mesh(), c, 2 ) ) ;
                cell_facet_centers.push_back(
                    mesh_cell_facet_center( reg.mesh(), c, 3 ) ) ;
            }
            reg_anns[r] = new ColocaterANN( cell_facet_centers, true ) ;
        }
        for( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            const Surface& surf = model_.surface( s ) ;
            vec3 first_facet_center = surf.facet_barycenter( 0 ) ;
            vec3 first_facet_normal = surf.facet_normal( 0 ) ;
            for( index_t r = 0 ; r < model_.nb_regions() ; ++r ) {
                std::vector< index_t > result ;
                reg_anns[r]->get_colocated( first_facet_center, result ) ;
                switch( result.size() ) {
                    case 0 :
                        break ;
                    case 1 :
                        {
                            index_t local_facet_id = result[0] % 4 ;
                            index_t cell_id =
                                0.25 * ( result[0] - local_facet_id ) ;
                            vec3 cell_facet_normal =
                                mesh_cell_facet_normal( model_.region(r).mesh(),
                                    cell_id,
                                    local_facet_id ) ;
                            bool side =
                                dot( first_facet_normal, cell_facet_normal ) > 0 ;
                            add_element_boundary(
                                GME::gme_t( GME::REGION, r ),
                                GME::gme_t( GME::SURFACE, s ),
                                side ) ;
                            add_element_in_boundary(
                                GME::gme_t( GME::SURFACE, s ),
                                GME::gme_t( GME::REGION, r ) ) ;
                            break ;
                        }
                    case 2 :
                        {
                            add_element_boundary(
                                GME::gme_t( GME::REGION, r ),
                                GME::gme_t( GME::SURFACE, s ),
                                true ) ;
                            add_element_in_boundary(
                                GME::gme_t( GME::SURFACE, s ),
                                GME::gme_t( GME::REGION, r ) ) ;
                            add_element_boundary(
                                GME::gme_t( GME::REGION, r ),
                                GME::gme_t( GME::SURFACE, s ),
                                false ) ;
                            add_element_in_boundary(
                                GME::gme_t( GME::SURFACE, s ),
                                GME::gme_t( GME::REGION, r ) ) ;
                            break ;
                        }
                    default :
                        ringmesh_assert_not_reached ;
                }
            }
        }
        for( index_t r = 0 ; r < model_.nb_regions() ; ++r ) {
            delete reg_anns[r] ;
        }
    }

    void GeoModelBuilderTSolid::compute_universe_boundaries()
    {
        std::vector< bool > surf_side_minus( model_.nb_surfaces(), false ) ;
        std::vector< bool > surf_side_plus( model_.nb_surfaces(), false ) ;
        for( index_t r = 0 ; r < model_.nb_regions() ; ++r ) {
            for( index_t s = 0 ; s < model_.region(r).nb_boundaries() ; ++s ) {
                if( model_.region(r).side(s) ) {
                    surf_side_plus[ model_.region(r).boundary(s).index() ] = true ;
                } else if( !model_.region(r).side(s) ) {
                    surf_side_minus[ model_.region(r).boundary(s).index() ] = true ;
                } else {
                    ringmesh_assert_not_reached
                }
            }
        }
        for( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            if( surf_side_minus[s] && !surf_side_plus[s] ) {
                add_element_boundary(
                    GME::gme_t( GME::REGION, NO_ID ),
                    GME::gme_t( GME::SURFACE, s ),
                    false ) ;
            } else if( !surf_side_minus[s] && surf_side_plus[s] ) {
                add_element_boundary(
                    GME::gme_t( GME::REGION, NO_ID ),
                    GME::gme_t( GME::SURFACE, s ),
                    true ) ;
            }
        }
    }

    void GeoModelBuilderTSolid::build_surface(
            index_t surface_id,
            std::vector< index_t >& facet_corners,
            std::vector< index_t >& facet_ptr,
            const std::vector< index_t >& gocad_vertices2region_id,
            const std::vector< index_t >& gocad_vertices2region_vertices )
    {
        std::vector< vec3 > cur_surf_points ;
        std::vector< index_t > cur_surf_facets ;
        std::vector< index_t > gocad_vertices2cur_surf_points(
            gocad_vertices2region_vertices.size(), NO_ID ) ;
        for( index_t co = 0 ; co < facet_corners.size() ; ++co ) {
            const index_t corner_gocad_id = facet_corners[ co ] ;
            if( gocad_vertices2cur_surf_points[ corner_gocad_id ] == NO_ID ) {
                // First time this facet corner is met in facet_corners
                const index_t corner_local_id =
                    gocad_vertices2region_vertices[ corner_gocad_id ] ;
                const index_t corner_region =
                    gocad_vertices2region_id[ corner_gocad_id ] ;
                const vec3 point =
                    model_.region( corner_region ).vertex( corner_local_id ) ;
                cur_surf_facets.push_back( cur_surf_points.size() ) ;
                gocad_vertices2cur_surf_points[ corner_gocad_id ] =
                    cur_surf_points.size() ;
                cur_surf_points.push_back( point ) ;
            } else {
                // If this facet corner have already been met in facet_corners
                cur_surf_facets.push_back(
                    gocad_vertices2cur_surf_points[ corner_gocad_id ] ) ;
            }
        }
        set_surface_geometry( surface_id,
            cur_surf_points,
            cur_surf_facets,
            facet_ptr ) ;
        cur_surf_points.clear() ;
        facet_corners.clear() ;
        facet_ptr.clear() ;
        facet_ptr.push_back( 0 ) ;
    }

    void GeoModelBuilderTSolid::compute_internal_borders()
    {
        std::vector< ColocaterANN* > anns( model_.nb_surfaces(), nil ) ;
        std::vector< Box3d > boxes( model_.nb_surfaces() ) ;
        for( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            const Surface& S = model_.surface(s) ;
            for( index_t p = 0; p < S.nb_vertices(); p++ ) {
                boxes[s].add_point( S.vertex( p ) ) ;
            }
            std::vector< vec3 > facet_edge_barycenters ;
            for( index_t f = 0 ; f < S.nb_cells() ; ++f ) {
                for( index_t e = 0 ; e < 3 ; ++e ) {
                    if(S.is_on_border(f,e)) {
                        const vec3 barycenter =
                            0.5 * ( S.vertex( f, e ) + S.vertex( f, ( e+1 )%3 ) ) ;
                        facet_edge_barycenters.push_back( barycenter ) ;
                    }
                }
            }
            anns[s] = new ColocaterANN( facet_edge_barycenters, true ) ;
        }
        for( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            const Surface& S = model_.surface(s) ;
            for( index_t f = 0 ; f < S.nb_cells() ; ++f ) {
                for( index_t e = 0 ; e < 3 ; ++e ) {
                   if( !S.is_on_border(f,e) ) {
                       const vec3 barycenter =
                           0.5 * ( S.vertex(f, e) + S.vertex(f, (e+1)%3 ) ) ;
                       std::vector< index_t > result ;
                       index_t tested_surf = 0 ;
                       while( result.empty() && tested_surf < anns.size() ) {
                           if( boxes[tested_surf].contains( barycenter ) ) {
                               anns[tested_surf]->
                               get_colocated( barycenter, result ) ;
                           }
                           ++tested_surf ;
                       }
                       if( !result.empty() ) {
                           S.mesh().facets.set_adjacent( f,e, GEO::NO_FACET ) ;
                       }
                   }
                }
            }
        }
        for( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            delete anns[s] ;
        }
    }
}
