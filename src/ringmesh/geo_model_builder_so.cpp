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

namespace RINGMesh {

    bool GeoModelBuilderTSolid::load_file()
    {
        std::vector< std::string > region_names ;
        std::vector< vec3 > vertices ;
        std::vector< index_t > ptr_regions_first_vertex ;
        std::vector< index_t > tetras_vertices ;
        std::vector< index_t > ptr_regions_first_tetra ;

        bool has_model_in_file = false ;

        ///@todo use two sub functions : one for reading and one for assigning
        // Reading .so file
        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if( in_.nb_fields() > 0 ) {
                if( in_.field_matches( 0, "TVOLUME" ) ) {
                    region_names.push_back( in_.field( 1 ) ) ;
                    ptr_regions_first_vertex.push_back( vertices.size() ) ;
                    ptr_regions_first_tetra.push_back( tetras_vertices.size() ) ;
                } else if( in_.field_matches( 0, "VRTX" ) || in_.field_matches( 0, "PVRTX" ) ) {
                    vec3 vertex ( in_.field_as_double( 2 ),
                                  in_.field_as_double( 3 ),
                                  z_sign_ * in_.field_as_double( 4 ) ) ;
                    vertices.push_back( vertex ) ;
                } else if( in_.field_matches( 0, "ATOM" ) || in_.field_matches( 0, "PATOM" ) ) {
                    vec3 vertex ( vertices[ in_.field_as_uint( 2 ) - 1 ] ) ;
                    vertices.push_back( vertex ) ;
                } else if( in_.field_matches( 0, "TETRA" ) ) {
                    tetras_vertices.push_back( in_.field_as_uint( 1 ) ) ;
                    tetras_vertices.push_back( in_.field_as_uint( 2 ) ) ;
                    tetras_vertices.push_back( in_.field_as_uint( 3 ) ) ;
                    tetras_vertices.push_back( in_.field_as_uint( 4 ) ) ;
                } else if( in_.field_matches( 0, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                    read_GCS() ;
                } else if( in_.field_matches( 0, "MODEL" ) ) {
                    has_model_in_file = true ;
                    // Mesh the regions with read vertices and tetras
                    ptr_regions_first_vertex.push_back( vertices.size() ) ;
                    ptr_regions_first_tetra.push_back( tetras_vertices.size() ) ;
                    mesh_regions( region_names,
                                  vertices,
                                  ptr_regions_first_vertex,
                                  tetras_vertices,
                                  ptr_regions_first_tetra ) ;
                } else if( in_.field_matches( 0, "SURFACE" ) ) {
                    // Create an Interface
                } else if( in_.field_matches( 0, "TFACE" ) ) {
                    // Create a Surface
                } else if( in_.field_matches( 0, "KEYVERTICES" ) ) {
                    // Check normals of triangles with propagation
                } else if( in_.field_matches( 0, "MODEL_REGION" ) ) {
                    // Don't know what to do
                }
            }
        }
        if ( !has_model_in_file ) {
            ptr_regions_first_vertex.push_back( vertices.size() ) ;
            ptr_regions_first_tetra.push_back( tetras_vertices.size() ) ;
            mesh_regions( region_names,
                          vertices,
                          ptr_regions_first_vertex,
                          tetras_vertices,
                          ptr_regions_first_tetra ) ;
            build_boundary_model() ;
        }
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
    ///@todo comment
    std::vector< std::string > GeoModelBuilderTSolid::read_number_of_mesh_elements(
            const index_t nb_regions,
            std::vector< index_t >& nb_elements_per_region )
    {
        std::vector< std::string > regions_name (nb_regions) ;

        // Check if the size of the given vector is enough
        // for all elements of all the regions.
        if ( nb_elements_per_region.size() < 2 * nb_regions ){
            GEO::Logger::err( "Mesh" )
                << "Not enough place for saving the number "
                        "of elements in each region" << std::endl ;
            return regions_name ;
        }

        index_t cur_region = -1 ;
        index_t nb_vertices_in_region = 0 ;
        index_t nb_tetras_in_region = 0 ;
        index_t nb_interfaces_in_bmodel = 0 ;
        index_t nb_surfaces_in_bmodel = 0 ;
        index_t nb_triangles_in_bmodel = 0 ;

        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if( in_.nb_fields() > 0 ) {
                if( in_.field_matches( 0, "TVOLUME" ) ) {
                    regions_name.at( ++cur_region ) = in_.field( 1 ) ;
                    if ( cur_region ){
                        nb_elements_per_region[ 2 * (cur_region - 1) ] = nb_vertices_in_region ;
                        nb_elements_per_region.at( 2 * (cur_region - 1) + 1 ) = nb_tetras_in_region ;
                        nb_vertices_in_region = 0 ;
                        nb_tetras_in_region = 0 ;
                    }
                } else if( in_.field_matches( 0, "VRTX" ) || in_.field_matches( 0, "PVRTX" )
                    || in_.field_matches( 0, "ATOM" ) || in_.field_matches( 0, "PATOM" ) ) {
                    ++nb_vertices_in_region ;
                } else if( in_.field_matches( 0, "TETRA" ) ) {
                    ++nb_tetras_in_region ;
                } else if( in_.field_matches( 0, "MODEL" ) ) {
                    nb_elements_per_region.at( 2 * cur_region ) = nb_vertices_in_region ;
                    nb_elements_per_region.at( 2 * cur_region + 1 ) = nb_tetras_in_region ;
                    nb_vertices_in_region = 0 ;
                    nb_tetras_in_region = 0 ;
                    nb_elements_per_region.resize( nb_elements_per_region.size() + 3 );
                } else if( in_.field_matches( 0, "SURFACE" ) ) {
                    ++nb_interfaces_in_bmodel ;
                } else if( in_.field_matches( 0, "TFACE" ) ) {
                    ++nb_surfaces_in_bmodel ;
                } else if( in_.field_matches( 0, "TRGL" ) ) {
                    ++nb_triangles_in_bmodel ;
                } else if( in_.field_matches( 0, "END" ) ) {
                    nb_elements_per_region.at( 2 * nb_regions ) = nb_interfaces_in_bmodel ;
                    nb_elements_per_region.at( 2 * nb_regions + 1 ) = nb_surfaces_in_bmodel ;
                    nb_elements_per_region.at( 2 * nb_regions + 2 ) = nb_triangles_in_bmodel ;
                }
            }
        }
        return regions_name ;
    }
    void GeoModelBuilderTSolid::print_number_of_mesh_elements(
            const index_t nb_regions,
            const std::vector< std::string >& regions_name,
            const std::vector< index_t >& nb_elements) const
    {

        GEO::Logger::out( "Mesh" )
            << "Mesh has " << nb_regions << " regions "
            << std::endl ;
        for (int i = 0 ; i < nb_regions ; ++i) {
            GEO::Logger::out( "Mesh" )
                << "Region " << regions_name.at(i) << " has"
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements.at( 2 * i ) << " vertices "
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements.at( 2 * i + 1 ) << " tetras "
                << std::endl ;
        }
        if ( nb_elements.size() == 2 * nb_regions + 3 ){
            GEO::Logger::out( "Mesh" )
                << "Model has" << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements.at( 2 * nb_regions ) << " interfaces"
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements.at( 2 * nb_regions + 1 ) << " surfaces "
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements.at( 2 * nb_regions + 2 ) << " triangles "
                << std::endl ;
        } else {
            GEO::Logger::out( "Mesh" )
                << "File does not contain the model" << std::endl ;
        }
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
        }
    }
    void GeoModelBuilderTSolid::build_boundary_model()
    {

    }
}
