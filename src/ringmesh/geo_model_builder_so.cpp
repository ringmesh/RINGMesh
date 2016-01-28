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
 *
 *
 *
 *
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

namespace {
    using namespace RINGMesh ;

    void count_nb_vertices_and_tetras_per_region(
        const std::string& filename,
        std::vector< index_t >& nb_elements_par_region,
        std::vector< index_t >& gocad_vertices2region_vertices,
        std::vector< index_t >& gocad_vertices2region_id )
    {
        nb_elements_par_region.clear() ;

        // Define a new LineInput counting number of elements
        GEO::LineInput line_input( filename ) ;

        // Initialize counters
        index_t cur_region = NO_ID ;
        index_t nb_vertices_in_region = 0 ;
        index_t nb_tetras_in_region = 0 ;
        index_t nb_surfaces_in_bmodel = 0 ;
        index_t nb_triangles_in_bmodel = 0 ;

        // Total (for the whole model) counter
        index_t nb_vertices_in_model = 0 ;

        // Reading file
        while( !line_input.eof() && line_input.get_line() ) {
            line_input.get_fields() ;
            if( line_input.nb_fields() > 0 ) {
                if( line_input.field_matches( 0, "TVOLUME" ) ||
                        line_input.field_matches( 0, "MODEL" ) ) {
                    if( cur_region != NO_ID ) {
                        nb_elements_par_region.push_back(
                            nb_vertices_in_region ) ;
                        nb_elements_par_region.push_back( nb_tetras_in_region ) ;
                        nb_vertices_in_model += nb_vertices_in_region ;
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
        gocad_vertices2region_id.reserve( nb_vertices_in_model ) ;
        gocad_vertices2region_vertices.reserve( nb_vertices_in_model ) ;
    }


    void print_nb_vertices_and_tetras_per_region(
        const std::vector< index_t >& nb_elements_per_region )
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

    void add_new_property(
        const std::string& prop_name,
        std::vector< std::string >& property_names,
        GEO::AttributesManager& attribute_manager )
    {
        property_names.push_back( prop_name ) ;
        /// @todo All the property types are double.
        /// Change in order to have the good type for each property.
        GEO::Attribute< double > property( attribute_manager, prop_name ) ;
    }

    index_t initialize_region(
        const std::string& region_name,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        GME::gme_t cur_region = geomodel_builder.create_element( GME::REGION ) ;
        geomodel_builder.set_element_name( cur_region, region_name ) ;
        return cur_region.index ;
    }

    void build_region(
        const index_t region_id,
        const index_t nb_vertices_in_next_region,
        const index_t nb_tetras_in_next_region,
        GeoModelBuilderTSolid& geomodel_builder,
        std::vector < vec3 >& region_vertices,
        std::vector < index_t >& tetra_corners )
    {
        geomodel_builder.set_region_geometry(
            region_id,
            region_vertices,
            tetra_corners ) ;

    }

    void reinitialize_region_vertices_and_teras(
        const index_t nb_vertices_in_next_region,
        const index_t nb_tetras_in_next_region,
        std::vector < vec3 >& region_vertices,
        std::vector < index_t >& tetra_corners )
    {
        region_vertices.clear() ;
        region_vertices.reserve( nb_vertices_in_next_region ) ;
        tetra_corners.clear() ;
        tetra_corners.reserve( nb_tetras_in_next_region ) ;
    }

    void read_vertex_coordinates(
        const GEO::LineInput& in,
        const int z_sign,
        vec3& vertex )
    {
        vertex.x = in.field_as_double( 2 ) ;
        vertex.y = in.field_as_double( 3 ) ;
        vertex.z = z_sign * in.field_as_double( 4 ) ;
    }

    void read_triangle(
        const GEO::LineInput& in,
        std::vector< index_t >& cur_surf_facets )
    {
        cur_surf_facets.push_back( in.field_as_uint( 1 ) - 1 ) ;
        cur_surf_facets.push_back( in.field_as_uint( 2 ) - 1 ) ;
        cur_surf_facets.push_back( in.field_as_uint( 3 ) - 1 ) ;
    }

    void read_tetraedra(
        const GEO::LineInput& in,
        const std::vector < index_t >& gocad_vertices2region_vertices,
        std::vector< index_t >& corners_id )
    {
        ringmesh_debug_assert( corners_id.size() == 4 ) ;
        corners_id[0] =
            gocad_vertices2region_vertices[ in.field_as_uint( 1 ) - 1 ] ;
        corners_id[1] =
            gocad_vertices2region_vertices[ in.field_as_uint( 2 ) - 1 ] ;
        corners_id[2] =
            gocad_vertices2region_vertices[ in.field_as_uint( 3 ) - 1 ] ;
        corners_id[3] =
            gocad_vertices2region_vertices[ in.field_as_uint( 4 ) - 1 ] ;
    }

    void compute_region_cell_facet_centers(
        const GeoModel& geomodel,
        const index_t region_id,
        std::vector< vec3 >& cell_facet_centers )
    {
        const Region& region = geomodel.region( region_id ) ;
        const index_t nb_cells = region.nb_cells() ;
        cell_facet_centers.reserve( 4 * nb_cells ) ;
        for( index_t c = 0 ; c < nb_cells ; ++c ) {
            for( index_t f = 0 ; f <= 3 ; ++f ) {
                cell_facet_centers.push_back(
                    mesh_cell_facet_center( region.mesh(), c, f ) ) ;
            }
        }
    }

    void compute_cell_facet_centers_region_anns(
        const GeoModel& geomodel,
        std::vector< ColocaterANN* >& region_anns )
    {
        for( index_t r = 0 ; r < geomodel.nb_regions() ; ++r ) {
            std::vector< vec3 > cell_facet_centers ;
            compute_region_cell_facet_centers( geomodel, r, cell_facet_centers ) ;
            region_anns[r] = new ColocaterANN( cell_facet_centers, true ) ;
        }
    }

    index_t are_surface_sides_region_boundaries(
        const Surface& surface,
        const ColocaterANN& region_ann,
        std::vector< index_t >& colocated_cell_facet_centers )
    {
        vec3 first_facet_center = surface.facet_barycenter( 0 ) ;
        region_ann.get_colocated( first_facet_center,
            colocated_cell_facet_centers ) ;
        return colocated_cell_facet_centers.size() ;
    }

    bool determine_surface_side_to_add(
        const GeoModel& geomodel,
        const index_t region_id,
        const index_t surface_id,
        const index_t cell_facet_center_id )
    {
        index_t local_facet_id = cell_facet_center_id % 4 ;
        index_t cell_id =
            0.25 * ( cell_facet_center_id - local_facet_id ) ;
        vec3 cell_facet_normal =
            mesh_cell_facet_normal(
                geomodel.region( region_id ).mesh(),
                cell_id,
                local_facet_id ) ;
        vec3 first_facet_normal =
            geomodel.surface( surface_id ).facet_normal( 0 ) ;
        return dot( first_facet_normal, cell_facet_normal ) > 0 ;
    }

    void fill_region_and_surface_boundaries_links(
        const index_t region_id,
        const index_t surface_id,
        const bool surf_side,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        geomodel_builder.add_element_boundary(
            GME::gme_t( GME::REGION, region_id ),
            GME::gme_t( GME::SURFACE, surface_id ),
            surf_side ) ;
        geomodel_builder.add_element_in_boundary(
            GME::gme_t( GME::SURFACE, surface_id ),
            GME::gme_t( GME::REGION, region_id ) ) ;
    }

    void add_both_surface_sides_to_region_boundaries(
        const index_t region_id,
        const index_t surface_id,
        GeoModelBuilderTSolid& geomodel_builder )

    {
        fill_region_and_surface_boundaries_links(
            region_id,
            surface_id,
            true,
            geomodel_builder ) ;
        fill_region_and_surface_boundaries_links(
            region_id,
            surface_id,
            false,
            geomodel_builder ) ;
    }

    void add_one_surface_side_to_region_boundaries(
        const index_t region_id,
        const index_t surface_id,
        const index_t cell_facet_center_id,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        bool side = determine_surface_side_to_add(
            geomodel_builder.model(),
            region_id,
            surface_id,
            cell_facet_center_id ) ;
        fill_region_and_surface_boundaries_links(
            region_id,
            surface_id,
            side,
            geomodel_builder ) ;
    }

    void add_surface_sides_to_region_boundaries(
        const index_t surface_id,
        const index_t region_id,
        const std::vector< index_t >& colocated_cell_facet_centers,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        switch( colocated_cell_facet_centers.size() ) {
            case 1 :
                add_one_surface_side_to_region_boundaries(
                    region_id,
                    surface_id,
                    colocated_cell_facet_centers[0],
                    geomodel_builder ) ;
                break ;
            case 2 :
                add_both_surface_sides_to_region_boundaries(
                    region_id,
                    surface_id,
                    geomodel_builder ) ;
                break ;
            default :
                ringmesh_assert_not_reached ;
        }
    }

    void add_surface_to_region_boundaries(
        const index_t surface_id,
        const std::vector< ColocaterANN* >& region_anns,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        const GeoModel& model = geomodel_builder.model() ;
        index_t cur_region = 0 ;
        index_t nb_added_surf_sides = 0 ;
        // Maximum 2 regions could be bounded by a single surface
        while ( cur_region < model.nb_regions() && nb_added_surf_sides < 2 ) {
            std::vector< index_t > colocated_cell_facet_centers ;
            index_t nb_surf_sides_are_boundary =
                are_surface_sides_region_boundaries(
                    model.surface( surface_id ),
                    *region_anns[cur_region],
                    colocated_cell_facet_centers ) ;
            if ( nb_surf_sides_are_boundary > 0 ) {
                add_surface_sides_to_region_boundaries(
                    surface_id,
                    cur_region,
                    colocated_cell_facet_centers,
                    geomodel_builder ) ;
                nb_added_surf_sides += nb_surf_sides_are_boundary ;
            }
            ++cur_region ;
        }
    }

    void compute_boundaries_of_geomodel_regions(
        GeoModelBuilderTSolid& geomodel_builder )
    {
        const GeoModel& model = geomodel_builder.model() ;
        std::vector< ColocaterANN* > reg_anns( model.nb_regions(), nil ) ;
        compute_cell_facet_centers_region_anns( model, reg_anns ) ;
        for( index_t s = 0 ; s < model.nb_surfaces() ; ++s ) {
            add_surface_to_region_boundaries( s, reg_anns, geomodel_builder ) ;
        }
        for( index_t r = 0 ; r < model.nb_regions() ; ++r ) {
            delete reg_anns[r] ;
        }
    }

    void add_surfaces_to_universe_boundaries(
        const std::vector< bool >& surf_minus_side,
        const std::vector< bool >& surf_plus_side,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        const GeoModel& model = geomodel_builder.model() ;
        for( index_t s = 0 ; s < model.nb_surfaces() ; ++s ) {
            if( surf_minus_side[s] && !surf_plus_side[s] ) {
                geomodel_builder.add_element_boundary(
                    GME::gme_t( GME::REGION, NO_ID ),
                    GME::gme_t( GME::SURFACE, s ),
                    false ) ;
            } else if( !surf_minus_side[s] && surf_plus_side[s] ) {
                geomodel_builder.add_element_boundary(
                    GME::gme_t( GME::REGION, NO_ID ),
                    GME::gme_t( GME::SURFACE, s ),
                    true ) ;
            }
        }
    }

    void determine_if_surface_sides_bound_regions(
        const GeoModel& geomodel,
        std::vector< bool >& surf_minus_side,
        std::vector< bool >& surf_plus_side )
    {
        for( index_t r = 0 ; r < geomodel.nb_regions() ; ++r ) {
            for( index_t s = 0 ; s < geomodel.region(r).nb_boundaries() ; ++s ) {
                if( geomodel.region(r).side(s) ) {
                    surf_plus_side[ geomodel.region(r).boundary(s).index() ] = true ;
                } else if( !geomodel.region(r).side(s) ) {
                    surf_minus_side[ geomodel.region(r).boundary(s).index() ] = true ;
                } else {
                    ringmesh_assert_not_reached
                }
            }
        }
    }

    void compute_universe_boundaries( GeoModelBuilderTSolid& geomodel_builder )
    {
        const GeoModel& model = geomodel_builder.model() ;
        // The universe boundaries are the surfaces with only one side in all
        // the boundaries of the other regions
        std::vector< bool > surf_minus_side( model.nb_surfaces(), false ) ;
        std::vector< bool > surf_plus_side( model.nb_surfaces(), false ) ;
        determine_if_surface_sides_bound_regions(
            model,
            surf_minus_side,
            surf_plus_side ) ;
        add_surfaces_to_universe_boundaries(
            surf_minus_side,
            surf_plus_side,
            geomodel_builder ) ;
    }

    void set_gocad_coordinate_system_axis_name(
        GEO::LineInput& in,
        std::vector< std::string >& gocad_coordinates_system_axis_name )
    {
        gocad_coordinates_system_axis_name.push_back( in.field(1) ) ;
        gocad_coordinates_system_axis_name.push_back( in.field(2) ) ;
        gocad_coordinates_system_axis_name.push_back( in.field(3) ) ;
    }

    void set_gocad_coordinate_system_axis_unit(
        GEO::LineInput& in,
        std::vector< std::string >& gocad_coordinates_system_axis_unit )
    {
        gocad_coordinates_system_axis_unit.push_back( in.field(1) ) ;
        gocad_coordinates_system_axis_unit.push_back( in.field(2) ) ;
        gocad_coordinates_system_axis_unit.push_back( in.field(3) ) ;
    }

    int set_gocad_coordinate_system_z_sign(
        GEO::LineInput& in )
    {
        if( in.field_matches( 1, "Elevation" ) ) {
            return 1 ;
        } else if( in.field_matches( 1, "Depth" ) ) {
            return -1 ;
        } else {
            ringmesh_assert_not_reached ;
        }
    }

} // anonymous namespace

namespace RINGMesh {

    bool GeoModelBuilderTSolid::load_file()
    {
        index_t cur_region ;

        // First : count the number of vertex and tetras
        // in each region
        std::vector< index_t > nb_elements_per_region ;
        count_nb_vertices_and_tetras_per_region(
            filename_,
            nb_elements_per_region,
            gocad_vertices2region_vertices_,
            gocad_vertices2region_id_) ;

        // Region vertices
        std::vector< vec3 > region_vertices ;
        region_vertices.reserve( nb_elements_per_region[0] ) ;
        // Region tetrahedron corners
        std::vector< index_t > tetra_corners ;
        tetra_corners.reserve( 4 * nb_elements_per_region[0] ) ;

        GME::gme_t current_interface ;
        GME::gme_t current_surface ;
        std::vector< index_t > cur_surf_facets ;
        std::vector< index_t > cur_surf_facet_ptr (1, 0) ;

        index_t nb_vertex_properties = 0 ;
        index_t nb_cell_properties = 0 ;
        std::vector< std::string > vertex_property_names ;
        std::vector< std::string > cell_property_names ;

        // Then : Reading .so file
        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if( in_.nb_fields() > 0 ) {
                if( in_.field_matches( 0, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                    read_and_set_gocad_coordinate_system() ;
                } else if( in_.field_matches( 0, "PROPERTIES" ) ) {
                    nb_vertex_properties = in_.nb_fields() - 1 ;
                } else if( in_.field_matches( 0, "PROPERTY_CLASS_HEADER" ) ) {
                    add_new_property( in_.field( 1 ),
                        vertex_property_names,
                        model_.mesh.vertex_attribute_manager() ) ;
                } else if( in_.field_matches( 0, "TETRA_PROPERTIES" ) ) {
                    nb_cell_properties = in_.nb_fields() - 1 ;
                } else if( in_.field_matches( 0, "TETRA_PROPERTY_CLASS_HEADER" ) ) {
                    add_new_property( in_.field( 1 ),
                        cell_property_names,
                        model_.mesh.cell_attribute_manager() ) ;
                } else if( in_.field_matches( 0, "TVOLUME" ) ) {
                    if( region_vertices.size() > 0 ) {
                        set_region_geometry(
                            cur_region,
                            region_vertices,
                            tetra_corners ) ;
                        reinitialize_region_vertices_and_teras(
                            nb_elements_per_region[ 2*cur_region ],
                            nb_elements_per_region[ 2*cur_region + 1 ],
                            region_vertices,
                            tetra_corners ) ;
//                        build_region(
//                            cur_region,
//                            nb_elements_per_region[ 2*cur_region ],
//                            nb_elements_per_region[ 2*cur_region + 1 ],
//                            region_vertices,
//                            tetra_corners ) ;
                    }
                    cur_region = initialize_region( in_.field( 1 ), *this ) ;
                } else if( in_.field_matches( 0, "VRTX" ) ||
                    in_.field_matches( 0, "PVRTX" ) ) {
                    read_and_add_vertex_to_region_vertices(
                        cur_region,
                        region_vertices ) ;
                } else if( in_.field_matches( 0, "ATOM" ) ||
                    in_.field_matches( 0, "PATOM" ) ) {
                    read_and_add_atom_to_region_vertices(
                        cur_region,
                        region_vertices ) ;
                } else if( in_.field_matches( 0, "TETRA" ) ) {
                    // Reading and create a tetra
                    std::vector< index_t > corners(4) ;
                    read_tetraedra(
                        in_, gocad_vertices2region_vertices_, corners ) ;
                    tetra_corners.insert( tetra_corners.end(),
                        corners.begin(),
                        corners.end() ) ;
                } else if( in_.field_matches( 0, "name:" ) ) {
                    // GeoModel name is set to the TSolid name.
                    set_model_name( in_.field( 1 ) ) ;
                } else if( in_.field_matches( 0, "MODEL" ) ) {
                    if( region_vertices.size() > 0 ) {
                        set_region_geometry(
                            cur_region,
                            region_vertices,
                            tetra_corners ) ;
                        reinitialize_region_vertices_and_teras(
                            0,
                            0,
                            region_vertices,
                            tetra_corners ) ;
//                        build_region(
//                            cur_region ,
//                            0,
//                            0,
//                            region_vertices,
//                            tetra_corners ) ;
                    }
                } else if( in_.field_matches( 0, "SURFACE" ) ) {
                    current_interface = create_element( GME::INTERFACE ) ;
                    set_element_name( current_interface, in_.field( 1 ) ) ;
                } else if( in_.field_matches( 0, "TFACE" ) ) {
                    // Compute the surface
                    if( cur_surf_facets.size() > 0 ) {
                        build_surface( current_surface.index,
                                cur_surf_facets,
                                cur_surf_facet_ptr ) ;
                    }
                    // Create a new surface
                    current_surface = create_element( GME::SURFACE ) ;
                    set_element_parent( current_surface, current_interface ) ;
                    add_element_child( current_interface, current_surface ) ;
                } else if( in_.field_matches( 0, "TRGL" ) ) {
                    read_triangle( in_, cur_surf_facets ) ;
                    cur_surf_facet_ptr.push_back( cur_surf_facets.size() ) ;
                } else if( in_.field_matches( 0, "MODEL_REGION" ) ) {
                    // Compute the last surface
                    if( cur_surf_facets.size() > 0 ) {
                        build_surface( current_surface.index,
                                cur_surf_facets,
                                cur_surf_facet_ptr ) ;
                    }
                }
            }
        }

        // Compute internal borders (by removing adjacencies on
        // triangle edges common to at least two surfaces)
        compute_surfaces_internal_borders() ;

        // Build GeoModel Lines and Corners from the surfaces
        model_.mesh.vertices.test_and_initialize() ;
        build_lines_and_corners_from_surfaces() ;

        // Regions boundaries
        compute_boundaries_of_geomodel_regions( *this ) ;

        // Universe boundaries
        compute_universe_boundaries( *this ) ;

        // Contacts building
        build_contacts() ;

        return true ;

    }

    void GeoModelBuilderTSolid::read_and_set_gocad_coordinate_system()
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
                set_gocad_coordinate_system_axis_name(
                    in_, gocad_coordinates_system_axis_name_ ) ;
            } else if( in_.field_matches( 0, "AXIS_UNIT" ) ) {
                set_gocad_coordinate_system_axis_unit(
                    in_, gocad_coordinates_system_axis_unit_ ) ;
            } else if( in_.field_matches( 0, "ZPOSITIVE" ) ) {
                z_sign_ = set_gocad_coordinate_system_z_sign( in_ ) ;
            }
        }
    }

    void GeoModelBuilderTSolid::read_and_add_vertex_to_region_vertices(
        const index_t region_id,
        std::vector < vec3 >& region_vertices )
    {
        gocad_vertices2region_vertices_.push_back(
            region_vertices.size() ) ;
        gocad_vertices2region_id_.push_back( region_id ) ;
        vec3 vertex ;
        read_vertex_coordinates( in_, z_sign_, vertex ) ;
        region_vertices.push_back( vertex ) ;
    }

    void GeoModelBuilderTSolid::read_and_add_atom_to_region_vertices(
        const index_t region_id,
        std::vector < vec3 >& region_vertices )
    {
        const index_t referring_vertex = in_.field_as_double( 2 ) - 1 ;
        const index_t referred_vertex_local_id =
                gocad_vertices2region_vertices_[referring_vertex] ;
        const index_t referred_vertex_region_id =
                gocad_vertices2region_id_[referring_vertex] ;
        if( referred_vertex_region_id < region_id ) {
            // If the atom referred to a vertex of another region,
            // acting like for a vertex
            gocad_vertices2region_vertices_.push_back(
                region_vertices.size() );
            gocad_vertices2region_id_.push_back( region_id ) ;
            region_vertices.push_back(
                    model_.region( referred_vertex_region_id ).vertex(
                            referred_vertex_local_id ) ) ;
        } else {
            // If the atom referred to an atom of the same region
            gocad_vertices2region_vertices_.push_back(
                referred_vertex_local_id ) ;
            gocad_vertices2region_id_.push_back(
                referred_vertex_region_id ) ;
        }
    }

    void GeoModelBuilderTSolid::build_surface(
        const index_t surface_id,
        std::vector< index_t >& facet_corners,
        std::vector< index_t >& facet_ptr )
    {
        std::vector< vec3 > cur_surf_points ;
        std::vector< index_t > cur_surf_facets ;
        get_surface_points_and_facets_from_gocad_indices(
            facet_corners,
            cur_surf_points,
            cur_surf_facets ) ;
        set_surface_geometry( surface_id,
            cur_surf_points,
            cur_surf_facets,
            facet_ptr ) ;
        cur_surf_points.clear() ;
        facet_corners.clear() ;
        facet_ptr.clear() ;
        facet_ptr.push_back( 0 ) ;
    }

    void GeoModelBuilderTSolid::get_surface_points_and_facets_from_gocad_indices(
        const std::vector< index_t >& facet_corners,
        std::vector< vec3 >& cur_surf_points,
        std::vector< index_t >& cur_surf_facets ) const
    {
        std::vector< index_t > gocad_vertices2cur_surf_points(
            gocad_vertices2region_vertices_.size(), NO_ID ) ;
        for( index_t co = 0 ; co < facet_corners.size() ; ++co ) {
            const index_t corner_gocad_id = facet_corners[ co ] ;
            get_surface_point_and_facet_from_gocad_index(
                corner_gocad_id,
                gocad_vertices2cur_surf_points,
                cur_surf_points,
                cur_surf_facets ) ;

        }
    }
    void GeoModelBuilderTSolid::get_surface_point_and_facet_from_gocad_index(
        const index_t vertex_gocad_id,
        std::vector< index_t >& gocad_vertices2cur_surf_points,
        std::vector< vec3 >& cur_surf_points,
        std::vector< index_t >& cur_surf_facets ) const
    {
        if( gocad_vertices2cur_surf_points[ vertex_gocad_id ] == NO_ID ) {
            // First time this facet corner is met in facet_corners
            vec3 point ;
            get_point_from_gocad_id( vertex_gocad_id, point ) ;
            cur_surf_facets.push_back( cur_surf_points.size() ) ;
            gocad_vertices2cur_surf_points[ vertex_gocad_id ] =
                cur_surf_points.size() ;
            cur_surf_points.push_back( point ) ;
        } else {
            // If this facet corner has already been met in facet_corners
            cur_surf_facets.push_back(
                gocad_vertices2cur_surf_points[ vertex_gocad_id ] ) ;
        }
    }

    void GeoModelBuilderTSolid::get_point_from_gocad_id(
        const index_t point_gocad_id,
        vec3& point ) const
    {
        const index_t point_local_id =
            gocad_vertices2region_vertices_[ point_gocad_id ] ;
        const index_t corner_region =
            gocad_vertices2region_id_[ point_gocad_id ] ;
        point =
            model_.region( corner_region ).vertex( point_local_id ) ;
    }


    void GeoModelBuilderTSolid::compute_surfaces_internal_borders()
    {
        std::vector< ColocaterANN* > anns( model_.nb_surfaces(), nil ) ;
        std::vector< Box3d > boxes( model_.nb_surfaces() ) ;
        compute_facet_edge_centers_anns_and_surface_boxes( anns, boxes ) ;
        for( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            compute_surface_internal_borders( s, anns, boxes ) ;
        }
        for( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            delete anns[s] ;
        }
    }

    void GeoModelBuilderTSolid::compute_surface_internal_borders(
        const index_t surface_id,
        const std::vector< ColocaterANN* >& surface_anns,
        const std::vector< Box3d >& surface_boxes )
    {
        const Surface& S = model_.surface( surface_id ) ;
        for( index_t f = 0 ; f < S.nb_cells() ; ++f ) {
            for( index_t e = 0 ; e < 3 ; ++e ) {
               if( !S.is_on_border(f,e) ) {
                   bool internal_border = is_edge_in_several_surfaces(
                       surface_id, f, e, surface_anns, surface_boxes ) ;
                   if( internal_border ) {
                       S.mesh().facets.set_adjacent( f, e, GEO::NO_FACET ) ;
                   }
               }
            }
        }
    }

    bool GeoModelBuilderTSolid::is_edge_in_several_surfaces(
        const index_t surface_id,
        const index_t facet,
        const index_t edge,
        const std::vector< ColocaterANN* >& surface_anns,
        const std::vector< Box3d >& surface_boxes ) const
    {
        const Surface& S = model_.surface( surface_id ) ;
        const vec3 barycenter = GEO::Geom::barycenter(
            S.vertex( facet, edge ),
            S.vertex( facet, ( edge+1 ) % 3 ) ) ;
       std::vector< index_t > result ;
       index_t tested_surf = 0 ;
       while( result.empty() &&
           tested_surf < surface_anns.size() ) {
           if( surface_boxes[tested_surf].contains( barycenter ) ) {
               surface_anns[tested_surf]->
                   get_colocated( barycenter, result ) ;
           }
           ++tested_surf ;
       }
       return !result.empty() ;
    }

    void GeoModelBuilderTSolid::compute_facet_edge_centers_anns_and_surface_boxes(
        std::vector< ColocaterANN* >& surface_anns,
        std::vector< Box3d >& surface_boxes ) const
    {
        for( index_t s = 0 ; s < model_.nb_surfaces() ; ++s ) {
            const Surface& S = model_.surface(s) ;
            for( index_t p = 0; p < S.nb_vertices(); p++ ) {
                surface_boxes[s].add_point( S.vertex( p ) ) ;
            }
            std::vector< vec3 > border_edge_barycenters ;
            get_surface_border_edge_barycenters( s, border_edge_barycenters ) ;
            surface_anns[s] = new ColocaterANN( border_edge_barycenters, true ) ;
        }
    }

    void GeoModelBuilderTSolid::get_surface_border_edge_barycenters(
        const index_t surface_id,
        std::vector< vec3 >& border_edge_barycenters ) const
    {
        const Surface& S = model_.surface( surface_id ) ;
        for( index_t f = 0 ; f < S.nb_cells() ; ++f ) {
            for( index_t e = 0 ; e < 3 ; ++e ) {
                if(S.is_on_border(f,e)) {
                    const vec3 barycenter = GEO::Geom::barycenter(
                        S.vertex( f, e ),
                        S.vertex( f, ( e+1 ) % 3 ) ) ;
                    border_edge_barycenters.push_back( barycenter ) ;
                }
            }
        }
    }
} // RINGMesh namespace
