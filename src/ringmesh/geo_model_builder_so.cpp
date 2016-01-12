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

#include <algorithm>
#include <iostream>
#include <iomanip>

#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>

#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/geometry.h>
#include <ringmesh/io.h>
#include <ringmesh/utils.h>

namespace RINGMesh {


    // Sorry but you have to rewrite your code without accessing the Mesh of any of the 
    // GeoModelElements and do not either access the GeoModelElements themselves
    bool GeoModelBuilderTSolid::load_file()
    {
        ///@todo Split this function into smaller functions when it will works

        // mmmhhh You should split it now. It will be much easier to debug [Jeanne]
        index_t last_tetra = 0 ;

        bool has_model_in_file = false ;

        GME::gme_t cur_region ;

        // First : count the number of vertex and tetras
        // in each region
        ///@todo:reserve space before and THEN assign
        std::vector< index_t > nb_mesh_elements_per_region =
                read_number_of_mesh_elements() ;
        print_number_of_mesh_elements( nb_mesh_elements_per_region ) ;

        ///@todo Link these vectors with the results above
        // Region vertices
        std::vector < vec3 > region_vertices ;
        // Region tetraedron corners
        std::vector< index_t > tetra_corners ;
        // Vector which transforms the indices of vertices from Gocad .so file
        // to the local (region) indices of vertices
        std::vector< index_t > gocad_vertices2region_vertices ;
        // Vector which transforms the indices of vertices from Gocad .so file
        // to the index of the region they belong to
        std::vector< index_t > gocad_vertices2region_id ;


        GME::gme_t current_interface ;
        GME::gme_t current_surface ;
        std::vector< index_t > cur_surf_facets ;
        std::vector< index_t > cur_surf_facet_ptr (1, 0) ;

        index_t nb_vertex_properties = 0 ;
        index_t nb_cell_properties = 0 ;
        std::vector < std::string > vertex_property_names ;
        std::vector < std::string > cell_property_names ;

        time_t start_reading_vol, end_reading_vol, end_reading_model, end_building_model ;
        time( &start_reading_vol ) ;
        std::cout << "Reading volume..." << std::endl ;

        // Then : Reading .so file
        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if( in_.nb_fields() > 0 ) {
                if( in_.field_matches( 0, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                    read_and_set_gocad_coordinates_system() ;
                } else if ( in_.field_matches( 0, "PROPERTIES" ) ) {
                    nb_vertex_properties = in_.nb_fields() - 1 ;
                } else if ( in_.field_matches( 0, "PROPERTY_CLASS_HEADER" ) ) {
                    add_new_property( vertex_property_names,
                            model_.mesh.vertex_attribute_manager() ) ;
                } else if ( in_.field_matches( 0, "TETRA_PROPERTIES" ) ) {
                    nb_cell_properties = in_.nb_fields() - 1 ;
                } else if ( in_.field_matches( 0, "TETRA_PROPERTY_CLASS_HEADER" ) ) {
                    add_new_property( cell_property_names,
                            model_.mesh.cell_attribute_manager() ) ;
                } else if( in_.field_matches( 0, "TVOLUME" ) ) {
                    if ( region_vertices.size() > 0 ) {
                        set_region_geometry ( cur_region.index, region_vertices, tetra_corners ) ;
                        region_vertices.clear() ;
                        tetra_corners.clear() ;
                    }
                    cur_region = create_region() ;
                } else if( in_.field_matches( 0, "VRTX" ) || in_.field_matches( 0, "PVRTX" ) ) {
                    gocad_vertices2region_vertices.push_back( region_vertices.size() );
                    gocad_vertices2region_id.push_back( cur_region.index ) ;
                    vec3 vertex ;
                    read_vertex_coordinates( vertex ) ;
                    region_vertices.push_back( vertex ) ;
                } else if( in_.field_matches( 0, "ATOM" ) || in_.field_matches( 0, "PATOM" ) ) {
                    index_t referring_vertex = in_.field_as_double( 2 ) - 1 ;

                    index_t referred_vertex_local_id =
                            gocad_vertices2region_vertices[referring_vertex] ;
                    index_t referred_vertex_region_id =
                            gocad_vertices2region_id[referring_vertex] ;

                    if ( referred_vertex_region_id < cur_region.index ) {
                        // If the atom referred to a vertex of another region,
                        // acting like for a vertex
                        gocad_vertices2region_vertices.push_back( region_vertices.size() );
                        gocad_vertices2region_id.push_back( cur_region.index ) ;
                        region_vertices.push_back(
                                model_.region( referred_vertex_region_id ).vertex(
                                        referred_vertex_local_id ) ) ;
                    } else {
                        // If the atom referred to an atom of the same region
                        gocad_vertices2region_vertices.push_back( referred_vertex_local_id ) ;
                        gocad_vertices2region_id.push_back( referred_vertex_region_id ) ;
                    }
                } else if( in_.field_matches( 0, "TETRA" ) ) {
                    // Reading and create a tetra
                    std::vector< index_t > corners(4) ;
                    read_tetraedra( gocad_vertices2region_vertices, corners ) ;
                    tetra_corners.insert( tetra_corners.end(), corners.begin(), corners.end() ) ;
                } else if( in_.field_matches( 0, "name:" ) ) {
                    // GeoModel name is set to the TSolid name.
                    set_model_name( in_.field( 1 ) ) ;
                } else if( in_.field_matches( 0, "MODEL" ) ) {
                    if ( region_vertices.size() > 0 ) {
                        set_region_geometry ( cur_region.index, region_vertices, tetra_corners ) ;
                        region_vertices.clear() ;
                        tetra_corners.clear() ;
                    }
                    time( &end_reading_vol ) ;
                    std::cout << "Timing : " << difftime( end_reading_vol, start_reading_vol ) << " seconds." << std::endl ;
                    std::cout << "Reading BRep model info..." << std::endl ;
                    has_model_in_file = true ;
                } else if( in_.field_matches( 0, "SURFACE" ) ) {
                    current_interface = create_element( GME::INTERFACE );
                    set_element_name( current_interface, in_.field( 1 ) ) ;
                } else if( in_.field_matches( 0, "TFACE" ) ) {
                    // Compute the surface
                    if ( cur_surf_facets.size() > 0 ) {
                        set_surface_geometry( current_surface.index, cur_surf_facets, cur_surf_facet_ptr ) ;
                        cur_surf_facets.clear() ;
                        cur_surf_facet_ptr.clear() ;
                        cur_surf_facet_ptr.push_back( 0 ) ;
                    }
                    // Create a new surface
                    current_surface = create_element( GME::SURFACE ) ;
                    set_element_parent( current_surface, current_interface ) ;
                    add_element_child( current_interface, current_surface ) ;
                } else if( in_.field_matches( 0, "TRGL" ) ) {
                    ///@todo think to reserve space before
                    cur_surf_facets.push_back( gocad_vertices2region_vertices[ in_.field_as_uint( 1 ) - 1 ] ) ;
                    cur_surf_facets.push_back( gocad_vertices2region_vertices[ in_.field_as_uint( 2 ) - 1 ] ) ;
                    cur_surf_facets.push_back( gocad_vertices2region_vertices[ in_.field_as_uint( 3 ) - 1 ] ) ;
                    cur_surf_facet_ptr.push_back( cur_surf_facets.size() ) ;
                } else if( in_.field_matches( 0, "MODEL_REGION" ) ) {
                    // Compute the last surface
                    if ( cur_surf_facets.size() > 0 ) {
                        // Good :) [Jeanne]
                        set_surface_geometry( current_surface.index, cur_surf_facets, cur_surf_facet_ptr ) ;
                        cur_surf_facets.clear() ;
                        cur_surf_facet_ptr.clear() ;
                        cur_surf_facet_ptr.push_back( 0 ) ;
                    }
                }
            }
        }
        time( &end_reading_model ) ;
        std::cout << "Timing : " << difftime( end_reading_model, end_reading_vol ) << " seconds." << std::endl ;
        std::cout << "Building Model..." << std::endl ;
        time_t step1, step2, step3 ;

        // Build GeoModel Lines and Corners from the surfaces
        // TODO : What are you doing here? Check all the nice functions in GeoModelBuilder [Jeanne]
        model_.mesh.vertices.test_and_initialize() ;
        time( &step1 ) ;
        std::cout << "Timing step 1 : " << difftime( step1, end_reading_model ) << " seconds." << std::endl ;

        // Compute internal borders (remove adjacencies)
        std::vector< ColocaterANN* > anns( model_.nb_surfaces(), nil ) ;
        std::vector< Box3d > boxes( model_.nb_surfaces() ) ;
        for ( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            const Surface& S = model_.surface(s) ;
            for( index_t p = 0; p < S.nb_vertices(); p++ ) {
                boxes[s].add_point( S.vertex( p ) ) ;
            }
            std::vector < vec3 > facet_edge_barycenters ;
            for ( index_t f = 0 ; f < S.nb_cells() ; ++f ) {
                for ( index_t e = 0 ; e < 3 ; ++e ) {
                    if (S.is_on_border(f,e)) {
                        facet_edge_barycenters.push_back( ( S.vertex(f, e) + S.vertex(f, (e+1)%3 ) ) * 0.5 );
                    }
                }
            }
            anns[s] = new ColocaterANN( facet_edge_barycenters, true ) ;
        }

//        index_t nb_calls = 0 ;
        // TODO : What is this ????? [Jeanne]
        // It is not at all the job of this function to take care of tasks like this one
        // All the functions to do that are implemented in Mesh or somewhere else
        // 2nd remark  DO NOT EVER use geometry to compute combinatorial things [Jeanne]
        for ( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            const Surface& S = model_.surface(s) ;
//            std::cout << "Surface " << s << std::endl ;
            for ( index_t f = 0 ; f < S.nb_cells() ; ++f ) {
                for ( index_t e = 0 ; e < 3 ; ++e ) {
                   if ( !S.is_on_border(f,e) ) {
                       vec3 barycenter = ( S.vertex(f, e) + S.vertex(f, (e+1)%3 ) ) * 0.5 ;
                       std::vector< index_t > result ;
                       index_t tested_surf = 0 ;
                       while ( result.empty() && tested_surf < anns.size() ) {
                           if ( boxes[tested_surf].contains( barycenter ) ) {
                               anns[tested_surf]->get_colocated(barycenter, result) ;
//                               ++nb_calls ;
                           }
                           ++tested_surf ;
                       }
                       if ( !result.empty() ) {
                           S.mesh().facets.set_adjacent( f,e, GEO::NO_FACET ) ;
                       }
                   }
                }
            }
        }
//        std::cout << "Number of calls of get_colocated function is " << nb_calls << std::endl ;

        for ( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            delete anns[s];
        }

        build_lines_and_corners_from_surfaces() ;

        time( &step2 ) ;
        std::cout << "Timing step 2 : " << difftime( step2, step1 ) << " seconds." << std::endl ;

        // Regions boundaries
        compute_boundaries_of_geomodel_regions() ;

        // Universe boundaries
        compute_universe_boundaries() ;

        // Contacts building
        build_contacts() ;

        return true ;

    }

    void GeoModelBuilderTSolid::read_and_set_gocad_coordinates_system()
    {
        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if ( in_.field_matches( 0, "END_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                return ;
            } else if ( in_.field_matches( 0, "NAME" ) ) {
                // Useless for the moment
            } else if ( in_.field_matches( 0, "PROJECTION" ) ) {
                // Useless for the moment
            } else if ( in_.field_matches( 0, "DATUM" ) ) {
                // Useless for the moment
            } else if ( in_.field_matches( 0, "AXIS_NAME" ) ) {
                gocad_coordinates_system_axis_name_.push_back( in_.field(1) ) ;
                gocad_coordinates_system_axis_name_.push_back( in_.field(2) ) ;
                gocad_coordinates_system_axis_name_.push_back( in_.field(3) ) ;
            } else if ( in_.field_matches( 0, "AXIS_UNIT" ) ) {
                gocad_coordinates_system_axis_unit_.push_back( in_.field(1) ) ;
                gocad_coordinates_system_axis_unit_.push_back( in_.field(2) ) ;
                gocad_coordinates_system_axis_unit_.push_back( in_.field(3) ) ;
            } else if ( in_.field_matches( 0, "ZPOSITIVE" ) ) {
                if( in_.field_matches( 1, "Elevation" ) ) {
                    z_sign_ = 1 ;
                } else if( in_.field_matches( 1, "Depth" ) ) {
                    z_sign_ = -1 ;
                } else {
                    ringmesh_assert_not_reached ;
                }
            }
        }
    }

    // TODO: Returning a vector is never a great idea, except if you're sure that it size is
    // limited [Jeanne]
    // In c++11 the vector is moved, but we are copying it.
    std::vector< index_t > GeoModelBuilderTSolid::read_number_of_mesh_elements()
    {
        GEO::LineInput line_input ( filename_ ) ;

        std::vector< index_t > nb_vertices_and_tets_per_region ;

        index_t cur_region = 0 ;
        index_t nb_vertices_in_region = 0 ;
        index_t nb_tetras_in_region = 0 ;
        index_t nb_surfaces_in_bmodel = 0 ;
        index_t nb_triangles_in_bmodel = 0 ;

        while( !line_input.eof() && line_input.get_line() ) {
            line_input.get_fields() ;
            if( line_input.nb_fields() > 0 ) {
                if( line_input.field_matches( 0, "TVOLUME" ) ||
                        line_input.field_matches( 0, "MODEL" ) ) {
                    if ( cur_region ){ // What is this test on a index_T ? [Jeanne]
                        nb_vertices_and_tets_per_region.push_back( nb_vertices_in_region ) ;
                        nb_vertices_and_tets_per_region.push_back( nb_tetras_in_region ) ;
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
        return nb_vertices_and_tets_per_region ;
    }

    void GeoModelBuilderTSolid::print_number_of_mesh_elements(
            const std::vector< index_t >& nb_elements_per_region ) const
    {
        index_t nb_regions = nb_elements_per_region.size()/2 ;
        GEO::Logger::out( "Mesh" )
            << "Mesh has " << nb_regions << " regions "
            << std::endl ;
        for (int i = 0 ; i < nb_regions ; ++i) {
            GEO::Logger::out( "Mesh" )
                << "Region " << i << " has"
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements_per_region.at( 2 * i ) << " vertices "
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements_per_region.at( 2 * i + 1 ) << " tetras "
                << std::endl ;
        }
    }

    void GeoModelBuilderTSolid::add_new_property(
            std::vector < std::string >& property_names,
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
            std::vector< index_t >& gocad_vertices2region_vertices,
            std::vector< index_t >& corners_id )
    {
        ringmesh_debug_assert( corners_id.size() == 4 ) ;
        corners_id[0] = gocad_vertices2region_vertices[ in_.field_as_uint( 1 ) - 1 ] ;
        corners_id[1] = gocad_vertices2region_vertices[ in_.field_as_uint( 2 ) - 1 ] ;
        corners_id[2] = gocad_vertices2region_vertices[ in_.field_as_uint( 3 ) - 1 ] ;
        corners_id[3] = gocad_vertices2region_vertices[ in_.field_as_uint( 4 ) - 1 ] ;
    }

    void GeoModelBuilderTSolid::compute_boundaries_of_geomodel_regions()
    {
        ///@todo Find how to accelerate this part which take more than 80% of the time (think to ColocaterANN)
        for ( index_t r = 0 ; r < model_.nb_regions() ; ++r ) {
            index_t id_reg = model_.region(r).index() ;
            for( index_t c = 0; c < model_.mesh.cells.nb_tet( id_reg ); ++c ) {
                for ( index_t f = 0; f < 4 ; ++f ) {
                    index_t facet = NO_ID ;
                    bool side = false ;
                    if (model_.mesh.cells.is_cell_facet_on_surface( model_.mesh.cells.tet( id_reg, c ), f, facet, side )) {
                        index_t surface = model_.mesh.facets.surface( facet ) ;
                        bool surface_in_boundary = false ;
                        bool surface_in_boundary_side = false ;
                        index_t b = NO_ID ;
                        while ( !(surface_in_boundary && side == surface_in_boundary_side )
                                && ++b < model_.region(r).nb_boundaries() ) {
                            if ( model_.region(r).boundary(b).gme_id() == model_.surface( surface ).gme_id() ) {
                                surface_in_boundary = true ;
                                surface_in_boundary_side = model_.region(r).side(b) ;
                            }
                        }
                        if ( !surface_in_boundary ) {
                            add_element_boundary(
                                GME::gme_t( GME::REGION, id_reg ),
                                GME::gme_t( GME::SURFACE, surface ),
                                side ) ;
                            add_element_in_boundary(
                                GME::gme_t( GME::SURFACE, surface ),
                                GME::gme_t( GME::REGION, id_reg ) ) ;
                        } else if (surface_in_boundary && side != surface_in_boundary_side ) {
                            // Case in which both sides of the surface are in boundaries of the region.
                            add_element_boundary(
                                GME::gme_t( GME::REGION, id_reg ),
                                GME::gme_t( GME::SURFACE, surface ),
                                side ) ;
                        }
                    }
                }
            }
        }
    }

    void GeoModelBuilderTSolid::compute_universe_boundaries()
    {
        std::vector< bool > surf_side_minus ( model_.nb_surfaces(), false ) ;
        std::vector< bool > surf_side_plus ( model_.nb_surfaces(), false ) ;
        for ( index_t r = 0 ; r < model_.nb_regions() ; ++r ) {
            for ( index_t s = 0 ; s < model_.region(r).nb_boundaries() ; ++s ) {
                if ( model_.region(r).side(s) ) {
                    surf_side_plus[model_.region(r).boundary(s).index()] = true ;
                } else if ( !model_.region(r).side(s) ) {
                    surf_side_minus[model_.region(r).boundary(s).index()] = true ;
                } else {
                    ringmesh_assert_not_reached
                }
            }
        }
        for ( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            if ( surf_side_minus[s] && !surf_side_plus[s] ) {
                add_element_boundary(
                    GME::gme_t( GME::REGION, NO_ID ),
                    GME::gme_t( GME::SURFACE, s ),
                    false ) ;
            } else if ( !surf_side_minus[s] && surf_side_plus[s] ) {
                add_element_boundary(
                    GME::gme_t( GME::REGION, NO_ID ),
                    GME::gme_t( GME::SURFACE, s ),
                    true ) ;
            }
        }
    }
}
