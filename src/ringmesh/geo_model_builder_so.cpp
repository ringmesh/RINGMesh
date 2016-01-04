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
#include <ringmesh/io.h>

namespace RINGMesh {

    bool GeoModelBuilderTSolid::load_file()
    {
        ///@todo Split this function into smaller functions when it will works
        index_t last_tetra = 0 ;

        bool has_model_in_file = false ;

        GME::gme_t cur_region = model_.universe().gme_id() ;
        GEO::Mesh* mesh = &(model_.universe().mesh()) ;

        // First : count the number of vertex and tetras
        // in each region
        ///@todo:reserve space before and THEN assign
        std::vector< index_t > nb_mesh_elements_per_region =
                read_number_of_mesh_elements() ;
        print_number_of_mesh_elements( nb_mesh_elements_per_region ) ;

        ///@todo Link these vector with the results above
        std::vector< index_t > vertices_id_in_region ;
        index_t nb_non_dupl_vrtx = 0 ;
        std::vector< index_t > map_gocad2gmm_vertices ;

        GME::gme_t current_interface ;
        GME::gme_t current_surface ;
        std::vector< index_t > cur_surf_facets ;
        std::vector< index_t > cur_surf_facet_ptr (1, 0) ;

        // Then : Reading .so file
        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if( in_.nb_fields() > 0 ) {
                if( in_.field_matches( 0, "TVOLUME" ) ) {
                    // Create a region in the GeoModel.
                    // Its mesh will be update while reading file.
                    cur_region = create_element( GME::REGION ) ;
                    set_element_name( cur_region, in_.field( 1 ) ) ;
                    mesh = &( model_.region(cur_region.index).mesh() ) ;
                } else if( in_.field_matches( 0, "VRTX" ) || in_.field_matches( 0, "PVRTX" ) ) {
                    // Add a vertex to the mesh of the current region.
                    double coord[3] = { in_.field_as_double( 2 ),
                                        in_.field_as_double( 3 ),
                                        z_sign_ * in_.field_as_double( 4 ) } ;
                    vertices_id_in_region.push_back( mesh->vertices.create_vertex( coord ) ) ;
                    map_gocad2gmm_vertices.push_back( nb_non_dupl_vrtx );
                    ++nb_non_dupl_vrtx ;
                } else if( in_.field_matches( 0, "ATOM" ) || in_.field_matches( 0, "PATOM" ) ) {
                    index_t referring_vertex = in_.field_as_double( 2 ) - 1 ;
                    index_t referring_vertex_region_id = NO_ID ;
                    index_t cum_nb_vertex = 0 ;
                    while ( cum_nb_vertex < referring_vertex ) {
                        ++referring_vertex_region_id ;
                        cum_nb_vertex += model_.region( referring_vertex_region_id ).mesh().vertices.nb() ;
                    }
                    if ( referring_vertex_region_id == NO_ID ) {
                        referring_vertex_region_id = 0 ;
                    }
                    ringmesh_debug_assert( referring_vertex_region_id !=  NO_ID ) ;
                    double* coord = model_.region( referring_vertex_region_id ).mesh().vertices.point_ptr( vertices_id_in_region[referring_vertex] ) ;
                    vertices_id_in_region.push_back( mesh->vertices.create_vertex( coord ) ) ;
                    map_gocad2gmm_vertices.push_back( map_gocad2gmm_vertices[referring_vertex] ) ;
                } else if( in_.field_matches( 0, "TETRA" ) ) {
                    // Reading and create a tetra
                    std::vector< index_t > vertices (4) ;
                    vertices[0] = vertices_id_in_region[ in_.field_as_uint( 1 ) - 1 ] ;
                    vertices[1] = vertices_id_in_region[ in_.field_as_uint( 2 ) - 1 ] ;
                    vertices[2] = vertices_id_in_region[ in_.field_as_uint( 3 ) - 1 ] ;
                    vertices[3] = vertices_id_in_region[ in_.field_as_uint( 4 ) - 1 ] ;
                    last_tetra = mesh->cells.create_tet( vertices[0], vertices[1], vertices[2], vertices[3] ) ;
                    ///@todo Improve following to avoid creation of the object "attribute_region"
                    GEO::Attribute <index_t> attribute_region ( mesh->cells.attributes(), "region" ) ;
                    attribute_region[ last_tetra ] = cur_region.index ;
                } else if( in_.field_matches( 0, "#" ) && in_.field_matches( 1, "CTETRA" ) ) {/*
                    // Read information about tetra faces
                    for ( index_t f = 0 ; f < 4 ; ++f ) {
                        if ( !in_.field_matches( f + 3, "none" ) ) {
                            // If a tetra face match with a triangle of an surface
                            std::string read_interface = in_.field( f + 3 ) ;
                            std::string interface_name = read_interface.substr(1) ;
                            bool sign_plus = (read_interface[0] == '+') ? true : false ;

                            // 1 - Find or create the interface
                            index_t id_model_interface = NO_ID ;
                            for ( index_t interf = 0 ; interf < model_.nb_interfaces() ; ++interf ) {
                                if ( model_.one_interface( interf ).name() == interface_name ) {
                                    id_model_interface = interf ;
                                    break ;
                                }
                            }
                            if ( id_model_interface == NO_ID ) {
                                GME::gme_t created_interface = create_element( GME::INTERFACE );
                                set_element_name( created_interface, interface_name ) ;
                                id_model_interface = created_interface.index ;
                            }

                            // 2 - Check if a child of the interface (a surface) is in boundary of the region.
                            GME::gme_t model_surface ( GME::SURFACE, NO_ID ) ;
                            if ( model_.one_interface( id_model_interface ).nb_children() == 0 ) {
                                // The interface has no children
                                model_surface = create_element( GME::SURFACE ) ;
                                set_element_parent( model_surface, model_.one_interface( id_model_interface ).gme_id() ) ;
                                add_element_child( model_.one_interface( id_model_interface ).gme_id(), model_surface ) ;
                                add_element_boundary( cur_region, model_surface, sign_plus ) ;
                            } else {
                                // If the interface has at least one child,
                                // check if a child is a boundary of the current region.
                                for (index_t b = 0 ; b < model_.region( cur_region.index ).nb_boundaries() ; ++b ) {
                                    if ( model_.region( cur_region.index ).boundary(b).parent()
                                            == model_.one_interface( id_model_interface ) ) {
                                        model_surface = model_.region( cur_region.index ).boundary(b).gme_id() ;
                                        break ;
                                    }
                                }
                                if ( model_surface.index == NO_ID ) {
                                    // Interface child(ren) is/are not a boundary of the region
                                    model_surface = create_element( GME::SURFACE ) ;
                                    set_element_parent( model_surface, model_.one_interface( id_model_interface ).gme_id() ) ;
                                    add_element_child( model_.one_interface( id_model_interface ).gme_id(), model_surface ) ;
                                    add_element_boundary( cur_region, model_surface, sign_plus ) ;
                                }
                            }

                            // 3 - Creation of the facet
                            index_t facet = mesh->facets.create_triangle( mesh->cells.tet_facet_vertex(last_tetra, f, 0),
                                    mesh->cells.tet_facet_vertex(last_tetra, f, 1),
                                    mesh->cells.tet_facet_vertex(last_tetra, f, 2) ) ;

                            // 4 - Set attribute
                            GEO::Attribute <index_t> attribute_interface ( mesh->facets.attributes(), "interface" ) ;
                            attribute_interface[ facet ] = id_model_interface ;

                            // 5 - Add the triangle in the surface
                            ///@todo use set_surface_geometry... When all the facets are known
                        }
                    }
                */} else if( in_.field_matches( 0, "name:" ) ) {
                    // GeoModel name is set to the TSolid name.
                    set_model_name( in_.field( 1 ) ) ;
                } else if( in_.field_matches( 0, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                    read_GCS() ;
                } else if( in_.field_matches( 0, "MODEL" ) ) {
                    has_model_in_file = true ;
                    model_.mesh.vertices.test_and_initialize() ;
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
                    cur_surf_facets.push_back( map_gocad2gmm_vertices[ in_.field_as_uint( 1 ) - 1 ] ) ;
                    cur_surf_facets.push_back( map_gocad2gmm_vertices[ in_.field_as_uint( 2 ) - 1 ] ) ;
                    cur_surf_facets.push_back( map_gocad2gmm_vertices[ in_.field_as_uint( 3 ) - 1 ] ) ;
                    cur_surf_facet_ptr.push_back( cur_surf_facets.size() ) ;
                } else if( in_.field_matches( 0, "MODEL_REGION" ) ) {
                    // Compute the last surface
                      if ( cur_surf_facets.size() > 0 ) {
                          set_surface_geometry( current_surface.index, cur_surf_facets, cur_surf_facet_ptr ) ;
                          cur_surf_facets.clear() ;
                          cur_surf_facet_ptr.clear() ;
                          cur_surf_facet_ptr.push_back( 0 ) ;
                      }
                }
            }
        }

        // Build GeoModel Lines and Corners from the surfaces
        model_.mesh.vertices.test_and_initialize() ;

        // Determine internal borders
        std::vector< std::vector < index_t > > vertex_surf_states( model_.mesh.vertices.nb() ) ;
        for( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            const Surface& S = model_.surface(s) ;
            for ( index_t f = 0 ; f < S.nb_cells() ; ++f ) {
                // All the cells of the surface are triangles
                index_t v0 = S.model_vertex_id_at_facet_corner( 3 * f ) ;
                index_t v1 = S.model_vertex_id_at_facet_corner( 3 * f + 1 ) ;
                index_t v2 = S.model_vertex_id_at_facet_corner( 3 * f + 2 ) ;
                if ( vertex_surf_states[v0].size() == 0 ) {
                    std::vector < GMEVertex > gmes = model_.mesh.vertices.gme_vertices(v0) ;
                    for ( index_t gme = 0 ; gme < gmes.size() ; ++gme ) {
                        if ( gmes[gme].gme_id.type == GME::SURFACE ) {
                            vertex_surf_states[v0].push_back( gmes[gme].gme_id.index ) ;
                        }
                    }
                }
                if ( vertex_surf_states[v1].size() == 0 ) {
                    std::vector < GMEVertex > gmes = model_.mesh.vertices.gme_vertices(v1) ;
                    for ( index_t gme = 0 ; gme < gmes.size() ; ++gme ) {
                        if ( gmes[gme].gme_id.type == GME::SURFACE ) {
                            vertex_surf_states[v1].push_back( gmes[gme].gme_id.index ) ;
                        }
                    }
                }
                if ( vertex_surf_states[v2].size() == 0 ) {
                    std::vector < GMEVertex > gmes = model_.mesh.vertices.gme_vertices(v2) ;
                    for ( index_t gme = 0 ; gme < gmes.size() ; ++gme ) {
                        if ( gmes[gme].gme_id.type == GME::SURFACE ) {
                            vertex_surf_states[v2].push_back( gmes[gme].gme_id.index ) ;
                        }
                    }
                }

                if ( !S.is_on_border(f,0) &&
                        vertex_surf_states[v0].size() > 1 &&
                        vertex_surf_states[v1].size() > 1 ) {
                    bool internal_border = false ;
                    if ( vertex_surf_states[v0].size() == vertex_surf_states[v1].size()  &&
                            vertex_surf_states[v0] == vertex_surf_states[v1] ) {
                        internal_border = true ;
                    } else if ( vertex_surf_states[v0].size() < vertex_surf_states[v1].size() ) {
                        bool match = true ;
                        for ( index_t el = 0 ; el < vertex_surf_states[v0].size() ; ++el ) {
                            if ( std::find( vertex_surf_states[v0].begin(),
                                    vertex_surf_states[v0].end(),
                                    vertex_surf_states[v1][el] )
                                    == vertex_surf_states[v0].end() ) {
                                match = false ;
                                break ;
                            }
                        }
                        if ( match ) {
                            internal_border = true ;
                        }
                    } else {
                        bool match = true ;
                        for ( index_t el = 0 ; el < vertex_surf_states[v1].size() ; ++el ) {
                            if ( std::find( vertex_surf_states[v1].begin(),
                                    vertex_surf_states[v1].end(),
                                    vertex_surf_states[v0][el] )
                                    == vertex_surf_states[v1].end() ) {
                                match = false ;
                            }
                        }
                        if ( match ) {
                            internal_border = true ;
                        }
                    }
                    if ( internal_border ) {
                        S.mesh().facets.set_adjacent( f,0, GEO::NO_FACET ) ;
                    }
                }
                if ( !S.is_on_border(f,1) &&
                        vertex_surf_states[v1].size() > 1 &&
                        vertex_surf_states[v2].size() > 1 ) {
                    bool internal_border = false ;
                    if ( vertex_surf_states[v1].size() == vertex_surf_states[v2].size()  &&
                            vertex_surf_states[v1] == vertex_surf_states[v2] ) {
                        internal_border = true ;
                    } else if ( vertex_surf_states[v1].size() < vertex_surf_states[v2].size() ) {
                        bool match = true ;
                        for ( index_t el = 0 ; el < vertex_surf_states[v1].size() ; ++el ) {
                            if ( std::find( vertex_surf_states[v1].begin(),
                                    vertex_surf_states[v1].end(),
                                    vertex_surf_states[v2][el] )
                                    == vertex_surf_states[v1].end() ) {
                                match = false ;
                                break ;
                            }
                        }
                        if ( match ) {
                            internal_border = true ;
                        }
                    } else {
                        bool match = true ;
                        for ( index_t el = 0 ; el < vertex_surf_states[v2].size() ; ++el ) {
                            if ( std::find( vertex_surf_states[v2].begin(),
                                    vertex_surf_states[v2].end(),
                                    vertex_surf_states[v1][el] )
                                    == vertex_surf_states[v2].end() ) {
                                match = false ;
                            }
                        }
                        if ( match ) {
                            internal_border = true ;
                        }
                    }
                    if ( internal_border ) {
                        S.mesh().facets.set_adjacent( f,1, GEO::NO_FACET ) ;
                    }
                }
                if ( !S.is_on_border(f,2) &&
                        vertex_surf_states[v2].size() > 1 &&
                        vertex_surf_states[v0].size() > 1 ) {
                    bool internal_border = false ;
                    if ( vertex_surf_states[v2].size() == vertex_surf_states[v0].size()  &&
                            vertex_surf_states[v2] == vertex_surf_states[v0] ) {
                        internal_border = true ;
                    } else if ( vertex_surf_states[v2].size() < vertex_surf_states[v0].size() ) {
                        bool match = true ;
                        for ( index_t el = 0 ; el < vertex_surf_states[v2].size() ; ++el ) {
                            if ( std::find( vertex_surf_states[v2].begin(),
                                    vertex_surf_states[v2].end(),
                                    vertex_surf_states[v0][el] )
                                    == vertex_surf_states[v2].end() ) {
                                match = false ;
                                break ;
                            }
                        }
                        if ( match ) {
                            internal_border = true ;
                        }
                    } else {
                        bool match = true ;
                        for ( index_t el = 0 ; el < vertex_surf_states[v0].size() ; ++el ) {
                            if ( std::find( vertex_surf_states[v0].begin(),
                                    vertex_surf_states[v0].end(),
                                    vertex_surf_states[v2][el] )
                                    == vertex_surf_states[v0].end() ) {
                                match = false ;
                            }
                        }
                        if ( match ) {
                            internal_border = true ;
                        }
                    }
                    if ( internal_border ) {
                        S.mesh().facets.set_adjacent( f,2, GEO::NO_FACET ) ;
                    }
                }
            }
        }
        build_lines_and_corners_from_surfaces() ;

        // Regions boundaries
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
                        while ( !(surface_in_boundary && side == surface_in_boundary_side)
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
                            // Case if both sides of the surface are in boundaries of the region.
                            add_element_boundary(
                                GME::gme_t( GME::REGION, id_reg ),
                                GME::gme_t( GME::SURFACE, surface ),
                                side ) ;
                        }
                    }
                }
            }
        }

        // Universe boundaries
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



//        for ( index_t v = 0 ; v < model_.mesh.vertices.nb() ; ++v ) {
//            std::cout << "b vertex " << v << std::endl ;
//            const std::vector< GMEVertex >& gmes = model_.mesh.vertices.gme_vertices(v) ;
//            for ( index_t g = 0 ; g < gmes.size() ; ++g ) {
//                std::cout << "gme " << g << " type : " << gmes[g].gme_id.type << std::endl ;
//            }
//        }
//        for ( index_t v = 0 ; v < model_.mesh.vertices.nb() ; ++v ) {
//            for ( index_t g = 0 ; g < model_.mesh.vertices.gme_vertices(v).size() ; ++g ) {
//                if (model_.mesh.vertices.gme_vertices(v)[g].gme_id.type == GME::REGION ) {
//                    model_.mesh.vertices.gme_vertices(v).erase(model_.mesh.vertices.gme_vertices(v).begin()+g ) ;
//                    --g ;
//                }
//            }
//        }
//        for ( index_t v = 0 ; v < model_.mesh.vertices.nb() ; ++v ) {
//            std::cout << "a vertex " << v << std::endl ;
//            const std::vector< GMEVertex >& gmes = model_.mesh.vertices.gme_vertices(v) ;
//            for ( index_t g = 0 ; g < gmes.size() ; ++g ) {
//                std::cout << "gme " << g << " type : " << gmes[g].gme_id.type << std::endl ;
//            }
//        }



        return true ;

    }
    void GeoModelBuilderTSolid::read_GCS()
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
                GCS_axis_name_.push_back( in_.field(1) ) ;
                GCS_axis_name_.push_back( in_.field(2) ) ;
                GCS_axis_name_.push_back( in_.field(3) ) ;
            } else if ( in_.field_matches( 0, "AXIS_UNIT" ) ) {
                GCS_axis_unit_.push_back( in_.field(1) ) ;
                GCS_axis_unit_.push_back( in_.field(2) ) ;
                GCS_axis_unit_.push_back( in_.field(3) ) ;
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
    std::vector< index_t > GeoModelBuilderTSolid::read_number_of_mesh_elements()
    {
        GEO::LineInput lineInput_count ( filename_ ) ;

        std::vector< index_t > nb_VRTX_and_TETRA_per_region ;

        index_t cur_region = 0 ;
        index_t nb_vertices_in_region = 0 ;
        index_t nb_tetras_in_region = 0 ;
        index_t nb_interfaces_in_bmodel = 0 ;
        index_t nb_surfaces_in_bmodel = 0 ;
        index_t nb_triangles_in_bmodel = 0 ;

        while( !lineInput_count.eof() && lineInput_count.get_line() ) {
            lineInput_count.get_fields() ;
            if( lineInput_count.nb_fields() > 0 ) {
                if( lineInput_count.field_matches( 0, "TVOLUME" ) ||
                    lineInput_count.field_matches( 0, "MODEL" ) ) {
                    if ( cur_region ){
                        nb_VRTX_and_TETRA_per_region.push_back( nb_vertices_in_region ) ;
                        nb_VRTX_and_TETRA_per_region.push_back( nb_tetras_in_region ) ;
                        nb_vertices_in_region = 0 ;
                        nb_tetras_in_region = 0 ;
                    }
                    ++cur_region ;
                } else if( lineInput_count.field_matches( 0, "VRTX" ) ||
                           lineInput_count.field_matches( 0, "PVRTX" ) ||
                           lineInput_count.field_matches( 0, "ATOM" ) ||
                           lineInput_count.field_matches( 0, "PATOM" ) ) {
                    ++nb_vertices_in_region ;
                } else if( lineInput_count.field_matches( 0, "TETRA" ) ) {
                    ++nb_tetras_in_region ;
                }
            }
        }
//        std::cout << in_.eof() << std::endl ;
//        std::cout << lineInput_count.eof() << std::endl ;
//        std::cout << in_.line_number() << std::endl ;
//        std::cout << lineInput_count.line_number() << std::endl ;
        return nb_VRTX_and_TETRA_per_region ;
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
//        if ( nb_elements.size() == 2 * nb_regions + 3 ){
//            GEO::Logger::out( "Mesh" )
//                << "Model has" << std::endl
//                << std::setw( 10 ) << std::left
//                << nb_elements.at( 2 * nb_regions ) << " interfaces"
//                << std::endl
//                << std::setw( 10 ) << std::left
//                << nb_elements.at( 2 * nb_regions + 1 ) << " surfaces "
//                << std::endl
//                << std::setw( 10 ) << std::left
//                << nb_elements.at( 2 * nb_regions + 2 ) << " triangles "
//                << std::endl ;
//        } else {
//            GEO::Logger::out( "Mesh" )
//                << "File does not contain the model" << std::endl ;
//        }
    }
    ///@todo comment
    void GeoModelBuilderTSolid::mesh_regions(
            const std::vector< std::string >& region_names,
            const std::vector< vec3 >& vertices,
            const std::vector< index_t >& ptr_regions_first_vertex,
            const std::vector< index_t >& tetras_vertices,
            const std::vector< index_t >& ptr_regions_first_tetra )
    {

        for ( index_t r = 0 ; r < region_names.size() ; ++r ) {
            GME::gme_t region = create_element( GME::REGION ) ;
            set_element_name( region, region_names[r] ) ;

            GEO::Mesh& cur_region_mesh = model_.region( region.index ).mesh() ;


            // Set vertices into the region
            std::vector< vec3 > cur_region_vertices ( ptr_regions_first_vertex[ r + 1 ]
                                        - ptr_regions_first_vertex[r] ) ;
            for ( index_t v = ptr_regions_first_vertex[r] ;
                  v < ptr_regions_first_vertex[ r + 1 ] ;
                  ++v ) {
                cur_region_vertices[ v - ptr_regions_first_vertex[r] ] = vertices[v] ;
            }
            set_element_vertices( region, cur_region_vertices, false ) ;

            // Set tetras into the region
            for ( index_t i = ptr_regions_first_tetra[r] ;
                  i < ptr_regions_first_tetra[ r + 1 ] ;
                  i = i + 4 ) {
                ///@todo assert because difference on index_t
                cur_region_mesh.cells.create_tet(tetras_vertices[i] - 1 - ptr_regions_first_vertex[r],
                        tetras_vertices[ i + 1 ] - 1 - ptr_regions_first_vertex[r],
                        tetras_vertices[ i + 2 ] - 1 - ptr_regions_first_vertex[r],
                        tetras_vertices[ i + 3 ] - 1 - ptr_regions_first_vertex[r] ) ;

            }
            GEO::Attribute<index_t> attribute_region ( cur_region_mesh.cells.attributes(), "region" ) ;
            GeoModelBuilderMesh buildermesh( model_, model_.universe().mesh(), "", "region" ) ;
            buildermesh.create_and_build_regions() ;

        }
    }
//    void GeoModelBuilderTSolid::build_boundary_model()
//    {
//
//    }
}
