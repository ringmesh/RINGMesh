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
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/geo_model_api.h>
#include <ringmesh/geogram_extension.h>

namespace RINGMesh {

//    bool GeoModelBuilderTSolid::load_file()
//    {
//        std::vector< std::string > region_names ;
//        std::vector< vec3 > vertices ;
//        std::vector< index_t > ptr_regions_first_vertex ;
//        std::vector< index_t > tetras_vertices ;
//        std::vector< index_t > ptr_regions_first_tetra ;
//
//        bool has_model_in_file = false ;
//
//        // First file reading : count the number of vertex and tetras
//        // in each region
//        std::vector< index_t > nb_mesh_elements_per_region =
//                read_number_of_mesh_elements() ;
//        print_number_of_mesh_elements( nb_mesh_elements_per_region ) ;
//
//        ///@todo use two sub functions : one for reading and one for assigning
//        // Reading .so file
//        while( !in_.eof() && in_.get_line() ) {
//            in_.get_fields() ;
//            if( in_.nb_fields() > 0 ) {
//                if( in_.field_matches( 0, "TVOLUME" ) ) {
//                    region_names.push_back( in_.field( 1 ) ) ;
//                    ptr_regions_first_vertex.push_back( vertices.size() ) ;
//                    ptr_regions_first_tetra.push_back( tetras_vertices.size() ) ;
//                } else if( in_.field_matches( 0, "VRTX" ) || in_.field_matches( 0, "PVRTX" ) ) {
//                    vec3 vertex ( in_.field_as_double( 2 ),
//                                  in_.field_as_double( 3 ),
//                                  z_sign_ * in_.field_as_double( 4 ) ) ;
//                    vertices.push_back( vertex ) ;
//                } else if( in_.field_matches( 0, "ATOM" ) || in_.field_matches( 0, "PATOM" ) ) {
//                    vec3 vertex ( vertices[ in_.field_as_uint( 2 ) - 1 ] ) ;
//                    vertices.push_back( vertex ) ;
//                } else if( in_.field_matches( 0, "TETRA" ) ) {
//                    tetras_vertices.push_back( in_.field_as_uint( 1 ) ) ;
//                    tetras_vertices.push_back( in_.field_as_uint( 2 ) ) ;
//                    tetras_vertices.push_back( in_.field_as_uint( 3 ) ) ;
//                    tetras_vertices.push_back( in_.field_as_uint( 4 ) ) ;
//                } else if( in_.field_matches( 0, "name:" ) ) {
//                    // GeoModel name is set to the TSolid name.
//                    set_model_name( in_.field( 1 ) ) ;
//                    std::cout << "passed here and name is : " << model_.name() << std::endl ;
//                } else if( in_.field_matches( 0, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) ) {
//                    read_GCS() ;
//                } else if( in_.field_matches( 0, "MODEL" ) ) {
//                    has_model_in_file = true ;
//
//                    set_element_vertices( model_.universe().gme_id(), vertices, true ) ;
//                    for(index_t i = 0 ; i < tetras_vertices.size() ; i = i + 4 ){
//                        model_.universe().mesh().cells.create_tet(tetras_vertices[i],
//                                                tetras_vertices[ i + 1 ],
//                                                tetras_vertices[ i + 2 ],
//                                                tetras_vertices[ i + 3 ]) ;
//                    }
//                    model_.universe().mesh().show_stats() ;
//                    std::cout << "nb vert univ " << model_.universe().nb_vertices() << std::endl;
//                    std::cout << "nb cell univ " << model_.universe().nb_cells() << std::endl;
//                    std::cout << "create attribute " << std::endl;
//                    GEO::Attribute<index_t> attribute_region ( model_.universe().mesh().cells.attributes(), "region1" ) ;
//                    std::cout << "fill attribute " << std::endl;
//                    attribute_region.fill(1) ;
//                    std::cout << "create GMBMesh " << std::endl;
//                    model_.universe().mesh().show_stats() ;
//                    std::cout << "nb facet" << model_.universe().mesh().facets.nb() << std::endl ;
//                    GeoModelBuilderMesh::prepare_surface_mesh_from_connected_components( model_.universe().mesh(), "facet_attribute" ) ;
//                    GeoModelBuilderMesh builder( model_, model_.universe().mesh(), "", "region1" ) ;
////                    std::cout << "remove colocate vertices " << std::endl;
////                    repair_colocate_vertices( model_.universe().mesh(), epsilon ) ;
//                    std::cout << "call create and build regions " << std::endl;
//                    builder.create_and_build_regions() ;
//
//
//                    // TO DECOMMENT
////                    // Mesh the regions with read vertices and tetras
////                    ptr_regions_first_vertex.push_back( vertices.size() ) ;
////                    ptr_regions_first_tetra.push_back( tetras_vertices.size() ) ;
////                    mesh_regions( region_names,
////                                  vertices,
////                                  ptr_regions_first_vertex,
////                                  tetras_vertices,
////                                  ptr_regions_first_tetra ) ;
////                    GeoModelBuilderMesh::prepare_surface_mesh_from_connected_components(model_.region(0).mesh(), "") ;
//                    GeoModelBuilderMesh buildermesh( model_, model_.universe().mesh(), "", "" ) ;
////                    buildermesh.build_regions() ;
//                    model_.universe().mesh().show_stats();
//                    std::cout << model_.region(0).mesh().cells.attributes().nb() << std::endl ;
//                } else if( in_.field_matches( 0, "SURFACE" ) ) {
//                    // Create an Interface
//                } else if( in_.field_matches( 0, "TFACE" ) ) {
//                    // Create a Surface
//                } else if( in_.field_matches( 0, "KEYVERTICES" ) ) {
//                    // Check normals of triangles with propagation
//                } else if( in_.field_matches( 0, "MODEL_REGION" ) ) {
//                    // Don't know what to do
//                }
//            }
//        }
//        if ( !has_model_in_file ) {
//            ptr_regions_first_vertex.push_back( vertices.size() ) ;
//            ptr_regions_first_tetra.push_back( tetras_vertices.size() ) ;
//            mesh_regions( region_names,
//                          vertices,
//                          ptr_regions_first_vertex,
//                          tetras_vertices,
//                          ptr_regions_first_tetra ) ;
//            build_boundary_model() ;
//        }
//        return true ;
//    }
    bool GeoModelBuilderTSolid::load_file()
        {
            ///@todo Split this function into smaller functions when it will works
            std::vector< std::string > region_names ;
            index_t nb_read_vertices = 0 ;
            index_t last_tetra = 0 ;
            std::vector< std::string > surf_names ;

            bool has_model_in_file = false ;

            GEO::Mesh mesh ;
            GEO::Attribute<index_t> attribute_region ( mesh.cells.attributes(), "region" ) ;
            GEO::Attribute<index_t> attribute_surf ( mesh.cell_facets.attributes(), "surf" ) ;

            // First file reading : count the number of vertex and tetras
            // in each region
            std::vector< index_t > nb_mesh_elements_per_region =
                    read_number_of_mesh_elements() ;
            print_number_of_mesh_elements( nb_mesh_elements_per_region ) ;
            ///@todo use the number of vertex and THEN assign

            // Reading .so file
            while( !in_.eof() && in_.get_line() ) {
                in_.get_fields() ;
                if( in_.nb_fields() > 0 ) {
                    if( in_.field_matches( 0, "TVOLUME" ) ) {
                        region_names.push_back( in_.field( 1 ) ) ;
                    } else if( in_.field_matches( 0, "VRTX" ) || in_.field_matches( 0, "PVRTX" ) ) {
                        double coord[3] = { in_.field_as_double( 2 ),
                                            in_.field_as_double( 3 ),
                                            in_.field_as_double( 4 ) } ;
                        mesh.vertices.create_vertex( coord ) ;
                    } else if( in_.field_matches( 0, "ATOM" ) || in_.field_matches( 0, "PATOM" ) ) {
                        double* coord = mesh.vertices.point_ptr( in_.field_as_uint( 2 ) - 1 ) ;
                        mesh.vertices.create_vertex( coord ) ;
                    } else if( in_.field_matches( 0, "TETRA" ) ) {
                        // Reading a tetra and build the facets
                        std::vector< index_t > vertices ;
                        vertices.push_back( in_.field_as_uint( 1 ) - 1 ) ;
                        vertices.push_back( in_.field_as_uint( 2 ) - 1 ) ;
                        vertices.push_back( in_.field_as_uint( 3 ) - 1 ) ;
                        vertices.push_back( in_.field_as_uint( 4 ) - 1 ) ;
                        last_tetra = mesh.cells.create_tet(vertices[0],vertices[1],vertices[2],vertices[3]) ;
                        attribute_region[last_tetra] = region_names.size() - 1 ;
                    } else if( in_.field_matches( 0, "#" ) && in_.field_matches( 1, "CTETRA" ) ) {
                        // Read information about tetra facets
                        for ( index_t f = 0 ; f < 4 ; ++f ) {
                            if ( !in_.field_matches( f + 3, "none" ) ) {
                                std::string read_surf = in_.field( f + 3 ) ;
                                read_surf = read_surf.substr(1) ;
                                index_t surf_id = 0 ;
                                while ( surf_id < surf_names.size()
                                        && surf_names[surf_id] != read_surf ) {
                                    ++surf_id ;
                                }
                                if ( surf_id == surf_names.size() ) {
                                    // Case: the surface is not in the vector
                                    surf_names.push_back( read_surf ) ;
                                }
                                index_t cell_facet = mesh.cells.facet( last_tetra, f ) ;
                                attribute_surf[ cell_facet ] = surf_id ;
                            }
                        }
                    } else if( in_.field_matches( 0, "name:" ) ) {
                        // GeoModel name is set to the TSolid name.
                        set_model_name( in_.field( 1 ) ) ;
                    } else if( in_.field_matches( 0, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                        read_GCS() ;
                    }
                }
            }
            ///@todo try to remove colocated vertices
            mesh.show_stats() ;
            repair_colocate_vertices( mesh, epsilon ) ;
            mesh.show_stats() ;
            mesh.cells.connect() ;
            mesh.show_stats() ;
            mesh.cells.compute_borders() ;
            mesh.show_stats() ;
            GEO::Attribute< index_t > attribute_facet_surf (mesh.facets.attributes(), "facet_surf") ;
            for ( index_t f = 0 ; f < mesh.facets.nb() ; ++f ) {
                index_t v0 = mesh.facets.vertex(f,0) ;
                index_t v1 = mesh.facets.vertex(f,1) ;
                index_t v2 = mesh.facets.vertex(f,2) ;
//                std::cout << f << " : " << v0 << " " << v1 << " " << v2 << " " << std::endl ;
//                bool find_cell_facet = false ;
                index_t local_facet_index = GEO::NO_FACET ;
                index_t tetra_index = 0 ;
                while ( local_facet_index == GEO::NO_FACET ) {
                    local_facet_index = mesh.cells.find_tet_facet( tetra_index++, v0, v1, v2 ) ;
                }
//                std::cout << "find tet,lf = " << tetra_index - 1 << "," << local_facet_index << std::endl ;
                attribute_facet_surf[f] = attribute_surf[mesh.cells.facet( tetra_index - 1, local_facet_index ) ] ;

            }
//            for ( index_t f = 0 ; f < mesh.facets.nb() ; ++f ) {
//                std::cout << f << " : " << attribute_facet_surf[f] << std::endl ;
//            }
            mesh.show_stats() ;
            std::cout << "====================" << std::endl ;
            std::cout << "Step 1" << std::endl ;
//            GeoModelBuilderMesh::prepare_surface_mesh_from_connected_components(mesh, "facet_surf") ;
            mesh.show_stats() ;
            print_model( model_ ) ;
            std::cout << "====================" << std::endl ;
            std::cout << "Step 2" << std::endl ;
            GeoModelBuilderMesh buildermesh( model_, mesh, "facet_surf", "region" ) ;
            std::cout << "GMBM nb_surf_attri_values = " << buildermesh.nb_surface_attribute_values() << std::endl ;
            std::cout << "GMBM nb_reg_attri_values = " << buildermesh.nb_region_attribute_values() << std::endl ;
//            GeoModelBuilderMesh buildermesh2( model_, mesh, "created_surf", "region" ) ;
//            std::cout << "GMBM2 nb_surf_attri_values = " << buildermesh2.nb_surface_attribute_values() << std::endl ;
            mesh.show_stats() ;
            print_model( model_ ) ;
            std::cout << "====================" << std::endl ;
            std::cout << "Step 3" << std::endl ;
            buildermesh.create_and_build_surfaces() ;
            mesh.show_stats() ;
            print_model( model_ ) ;
            std::cout << "====================" << std::endl ;
            std::cout << "Step 4" << std::endl ;
            mesh.show_stats() ;
            print_model( model_ ) ;
            buildermesh.build_model_from_surfaces() ;
            mesh.show_stats() ;
            print_model( model_ ) ;
            std::cout << "====================" << std::endl ;
            std::cout << "Step 5" << std::endl ;
            buildermesh.create_and_build_regions() ;
            print_model(model_) ;
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
                    ringmesh_assert_not_reached;
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
    void GeoModelBuilderTSolid::build_boundary_model()
    {

    }
}
