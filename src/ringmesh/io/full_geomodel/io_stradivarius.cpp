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
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */
namespace {

    class StradivariusBuilder final: public GeoModelBuilderFile< 2 > {
    public:
        StradivariusBuilder( GeoModel< 2 >& geomodel, std::string filename )
            : GeoModelBuilderFile< 2 >( geomodel, std::move( filename ) )
        {

        }
        static const index_t SHIFT = 1;
        static const index_t SEABED = 0;

    private:
        void load_file() override
        {
            GEO::LineInput file( filename_ );
            while( !file.eof() && file.get_line() ) {
                file.get_fields();
                if( file.field_matches( 0, "//" ) ) {
                    continue;
                } else if( file.field_matches( 0, "liste" ) ) {
                    if( file.nb_fields() != 3 ) {
                        continue;
                    }
                    if( file.field_matches( 2, "milieux" ) ) {
                        load_surfaces( file );
                    } else if( file.field_matches( 2, "horizons" ) ) {
                        load_lines( file );
                    } else if( file.field_matches( 2, "points" ) ) {
                        load_points( file );
                    }
                }
            }
            build_lines();
            build_corners_from_lines();
        }

        void build_lines()
        {
            for( index_t l : range( line_used_.size() ) ) {
                if( !line_used_[l] ) {
                    continue;
                }
                const auto& indices = lines_vertices_[l];
                std::vector< vec2 > vertices;
                vertices.reserve( indices.size() );
                for( index_t index : indices ) {
                    vertices.push_back( points_[index] );
                }
                geometry.set_line( line_mapping_[l], vertices );
            }
        }

        void load_points( GEO::LineInput& file )
        {
            file.get_line();
            file.get_fields();
            index_t nb_points = file.field_as_uint( 0 );
            points_.resize( nb_points );
            for( index_t p : range( nb_points ) ) {
                file.get_line();
                file.get_fields();
                points_[p] = {file.field_as_double( 1 ), file.field_as_double( 2 )};
            }
        }

        void load_lines( GEO::LineInput& file )
        {
            index_t nb_lines;
            std::tie( nb_lines, line_mapping_ ) = compute_line_mapping();
            lines_vertices_.reserve( nb_lines );
            topology.create_mesh_entities( Line < 2 > ::type_name_static(),
                nb_lines );
            file.get_line();
            file.get_fields();
            index_t nb_lines_in_file = file.field_as_uint( 0 );
            for( index_t l : range( nb_lines_in_file ) ) {
                import_line( file, line_mapping_[l], !line_used_[l] );
            }
        }

        std::tuple< index_t, std::vector< index_t > > compute_line_mapping()
        {
            std::vector< index_t > mapping( line_used_.size(), NO_ID );
            index_t offset = 0;
            for( index_t l : range( line_used_.size() ) ) {
                if( line_used_[l] ) {
                    mapping[l] = offset++;
                }
            }
            return std::make_tuple( offset, mapping );
        }

        void import_line( GEO::LineInput& file, index_t line_id, bool skiped_line )
        {
            gmme_id line( Line< 2 >::type_name_static(), line_id );
            import_line_topology( file, line, skiped_line );
            lines_vertices_.push_back( import_line_geometry( file, line, skiped_line ) );
        }

        void import_line_topology( GEO::LineInput& file, const gmme_id& line, bool skiped_line )
        {
            file.get_line();
            file.get_fields();
            if( skiped_line ) {
                return;
            }
            index_t left_surface = index_t( file.field_as_int( 0 ) );
            if( left_surface != NO_ID && left_surface != SEABED ) {
                topology.add_mesh_entity_boundary_relation(
                    {   Surface < 2 > ::type_name_static(), left_surface - SHIFT}, line,
                    true );
            }
            index_t right_surface = index_t( file.field_as_int( 1 ) );
            if( right_surface != NO_ID && right_surface != SEABED ) {
                topology.add_mesh_entity_boundary_relation(
                    {   Surface < 2 > ::type_name_static(), right_surface - SHIFT}, line,
                    true );
            }
            info.set_mesh_entity_name( line, file.field( 4 ) );
        }

        std::vector< index_t > import_line_geometry(
            GEO::LineInput& file,
            const gmme_id& line,
            bool skiped_line )
        {
            index_t nb_vertices = file.field_as_uint( 2 );
            std::vector< index_t > vertices( nb_vertices );
            for( index_t v : range( nb_vertices ) ) {
                file.get_line();
                file.get_fields();
                vertices[v] = file.field_as_uint( 0 );
            }

            if( skiped_line ) {
                return std::vector< index_t >();
            } else {
                return vertices;
            }
        }

        void load_surfaces( GEO::LineInput& file )
        {
            index_t nb_surfaces = create_surfaces( file );
            for( index_t s : range( nb_surfaces ) ) {
                import_surface( file, s );
            }
        }

        void import_surface( GEO::LineInput& file, index_t surface_id )
        {
            file.get_line();
            file.get_fields();
            index_t nb_boundaries = file.field_as_uint( 1 );
            for( index_t b : range( nb_boundaries ) ) {
                ringmesh_unused( b );
                file.get_line();
                file.get_fields();
                index_t line_id = file.field_as_uint( 0 );
                mark_line_as_used( line_id );
            }
        }

        void mark_line_as_used( index_t line_id )
        {
            if( line_id >= line_used_.size() ) {
                line_used_.resize( line_id + 1, false );
            }
            line_used_[line_id] = true;
        }

        index_t create_surfaces( GEO::LineInput& file )
        {
            file.get_line();
            file.get_fields();
            index_t nb_surfaces = file.field_as_uint( 0 );
            topology.create_mesh_entities( Surface < 2 > ::type_name_static(),
                nb_surfaces );
            return nb_surfaces;
        }

    private:
        std::vector< bool > line_used_;
        std::vector< index_t > line_mapping_;
        std::vector< std::vector< index_t > > lines_vertices_;
        std::vector< vec2 > points_;
    };

    /*!
     * @brief Export for the GMSH format 2.2 which is described here:
     * http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format
     * NB : Mesh entities are also exported
     */
    class StradivariusIOHandler final: public GeoModelIOHandler< 2 > {
    public:
        void load( const std::string& filename, GeoModel< 2 >& geomodel ) final
        {
            StradivariusBuilder builder( geomodel, filename );
            builder.build_geomodel();
        }
        void save( const GeoModel< 2 >& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );
            check_stradivarius_validity( geomodel );
            save_header( out );
            save_milieux( out, geomodel );
            save_horizons( out, geomodel );
            save_faults( out, geomodel );
            save_vertices( out, geomodel );
        }

    private:
        void check_stradivarius_validity( const GeoModel< 2 >& geomodel )
        {
            for( index_t v : range( geomodel.mesh.vertices.nb() ) ) {
                model_box_.add_point( geomodel.mesh.vertices.vertex( v ) );
            }
        }
        void save_header( std::ofstream& out )
        {
            out << "// File wrote with RINGMesh" << std::endl;
        }

        void save_milieux( std::ofstream& out, const GeoModel< 2 >& geomodel )
        {
            out << "liste des milieux" << std::endl;
            out << geomodel.nb_surfaces() << std::endl;
            for( const auto& surface : geomodel.surfaces() ) {
                save_milieu( out, surface );
            }
        }

        void save_milieu( std::ofstream& out, const Surface< 2 >& surface )
        {
            save_milieu_index( out, surface );
            save_milieu_boundaries( out, surface );
        }

        void save_milieu_index( std::ofstream& out, const Surface< 2 >& surface )
        {
            out << surface.index() + shift_for_seabed << " "
                << surface.nb_boundaries() << std::endl;
        }

        void save_milieu_boundaries(
            std::ofstream& out,
            const Surface< 2 >& surface )
        {
            for( index_t boundary : range( surface.nb_boundaries() ) ) {
                out << surface.boundary( boundary ).index() << std::endl;
            }
        }

        /*!
         * @brief Here horizons is in term of stradivarius nomenclature
         */
        void save_horizons( std::ofstream& out, const GeoModel< 2 >& geomodel )
        {
            out << "liste des horizons" << std::endl;
            out << geomodel.nb_lines() + nb_horizons_for_seabed << " "
                << "// milG milD nbrePoint   numeroHorizon   nom d'horizon"
                << std::endl;
            ringmesh_assert( additional_points_for_seabed.size() == 2 );
            for( const auto& line : geomodel.lines() ) {
                save_horizon( out, line );
            }
            save_seabed( out, geomodel );

        }
        void save_horizon( std::ofstream& out, const Line< 2 >& line )
        {
            if( line.nb_incident_entities() == 2 ) {
                save_classical_horizon( out, line );

            } else if( line.nb_incident_entities() == 1 ) {
                save_horizon_as_boundary( out, line );
            } else {
                ringmesh_assert_not_reached;
            }
        }

        void save_classical_horizon( std::ofstream& out, const Line< 2 >& line )
        {
            const Surface< 2 > &incident_surface_1 = line.incident_entity( 0 );
            index_t left = NO_ID;
            index_t right = NO_ID;

            index_t boundary_id = 0;
            while( incident_surface_1.boundary_gmme( boundary_id ) != line.gmme() ) {
                boundary_id++;
            }
            if( incident_surface_1.side( boundary_id ) ) {
                left = 0;
                right = 1;
            } else {
                left = 1;
                right = 0;
            }
            save_horizon_topology( out,
                line.incident_entity( left ).index() + shift_for_seabed,
                line.incident_entity( right ).index() + shift_for_seabed,
                line.nb_vertices(), line.index(),
                "Line_" + GEO::String::to_string( line.index() ) );
            save_horizons_points( out, line );

        }

        void save_horizon_as_boundary( std::ofstream& out, const Line< 2 >& line )
        {
            vec2 v0 = line.vertex( 0 );
            vec2 v1 = line.vertex( 1 );
            if( v0.x == v1.x ) {
                if( v0.y < v1.y ) {
                    if( v0.x == model_box_.min().x ) {
                        //left boundary
                        left_is_outside( out, line );
                    } else if( v0.x == model_box_.max().x ) {
                        // right boundary
                        right_is_outside( out, line );
                    }
                } else {
                    if( v0.x == model_box_.min().x ) {
                        //left boundary
                        right_is_outside( out, line );
                    } else if( v0.x == model_box_.max().x ) {
                        // right boundary
                        left_is_outside( out, line );
                    }
                }
            } else if( v0.y == v1.y && v0.y == model_box_.min().y ) {
                // bottom boundary
                if( v0.x > v1.x ) {
                    left_is_outside( out, line );
                } else {
                    right_is_outside( out, line );
                }
            } else {
                DEBUG( line.index() );
                //top surface
                if( v0.x > v1.x ) {

                    right_is_outside( out, line );
                } else {
                    left_is_outside( out, line );
                }
            }
            save_horizons_points( out, line );
        }

        void left_is_outside( std::ofstream& out, const Line< 2 >& line )
        {
            save_horizon_topology( out, outside_region,
                line.incident_entity( 0 ).index() + shift_for_seabed,
                line.nb_vertices(), line.index(),
                "Line_" + GEO::String::to_string( line.index() ) );
        }

        void right_is_outside( std::ofstream& out, const Line< 2 >& line )
        {
            save_horizon_topology( out,
                line.incident_entity( 0 ).index() + shift_for_seabed, outside_region,
                line.nb_vertices(), line.index(),
                "Line_" + GEO::String::to_string( line.index() ) );
        }

        void save_horizons_points( std::ofstream& out, const Line< 2 >& line )
        {
            for( index_t vertex : range( line.nb_vertices() ) ) {
                out
                    << line.geomodel().mesh.vertices.geomodel_vertex_id( line.gmme(),
                        vertex ) << std::endl;
            }
        }

        void save_seabed( std::ofstream& out, const GeoModel< 2 >& geomodel )
        {
            // goal : find the maximum z of the model, add +200
            // to have the seabed
            for( index_t v : range( geomodel.mesh.vertices.nb() ) ) {
                model_box_.add_point( geomodel.mesh.vertices.vertex( v ) );
            }
            right_up_point_ = vec2( model_box_.max()[0],
                model_box_.max()[1] + default_sea_depth );
            left_up_point_ = vec2( model_box_.min()[0],
                model_box_.max()[1] + default_sea_depth );

            index_t nb_lines = geomodel.nb_lines();
            index_t nb_points = geomodel.mesh.vertices.nb();

            index_t right_up_point_of_the_model_wo_seabed =
                vertex_id_of_the_max_height_corner_at_constant_x( right_up_point_.x,
                    model_box_.min().y, geomodel );

            index_t left_up_point_of_the_model_wo_seabed =
                vertex_id_of_the_max_height_corner_at_constant_x( left_up_point_.x,
                    model_box_.min().y, geomodel );

            // Top of the seabed
            save_horizon_topology( out, outside_region, seabed_region,
                nb_points_for_seabed_boundaries, nb_lines, "SB_TOP" );
            out << nb_points << std::endl;
            out << nb_points + 1 << std::endl;

            // Left side of the seabed
            save_horizon_topology( out, outside_region, seabed_region,
                nb_points_for_seabed_boundaries, nb_lines + 1, "SB_LEFT" );
            out << left_up_point_of_the_model_wo_seabed << std::endl;
            out << nb_points << std::endl;

            // Right side of the seabed
            save_horizon_topology( out, outside_region, seabed_region,
                nb_points_for_seabed_boundaries, nb_lines + 2, "SB_RIGHT" );
            out << nb_points << std::endl;
            out << right_up_point_of_the_model_wo_seabed + 1 << std::endl;
        }

        index_t vertex_id_of_the_max_height_corner_at_constant_x(
            double x,
            double min_y,
            const GeoModel< 2 >& geomodel )
        {
            index_t vertex_index = NO_ID;
            double cur_y = min_y;
            for( const auto& corner : geomodel.corners() ) {
                const vec2& corner_vertex = corner.vertex( 0 );
                if( corner_vertex.x == x ) {
                    if( corner_vertex.y > cur_y ) {
                        cur_y = corner_vertex.y;
                        vertex_index = geomodel.mesh.vertices.geomodel_vertex_id(
                            corner.gmme(), 0 );
                    }
                }
            }
            ringmesh_assert( vertex_index != NO_ID );
            return vertex_index;
        }

        void save_horizon_topology(
            std::ofstream& out,
            int left_id,
            int right_id,
            index_t nb_vertices,
            index_t line_id,
            const std::string& name )
        {
            out << left_id << " " << right_id << " " << nb_vertices << " " << line_id
                << " " << name << std::endl;
        }

        void save_faults( std::ofstream& out, const GeoModel< 2 >& geomodel )
        {
            out << "liste des failles" << std::endl;
            out << 0 << std::endl;
        }

        void save_vertices( std::ofstream& out, const GeoModel< 2 >& geomodel )
        {
            out << "liste des points" << std::endl;
            out << geomodel.mesh.vertices.nb() + 2 << std::endl;
            for( index_t v : range( geomodel.mesh.vertices.nb() ) ) {
                out << v << " " << geomodel.mesh.vertices.vertex( v ) << " ";
                save_uncertainties( out );
                out << std::endl;
            }
            save_seabed_point( out, geomodel.mesh.vertices.nb(), left_up_point_ );
            save_seabed_point( out, geomodel.mesh.vertices.nb() + 1,
                right_up_point_ );
        }

        void save_seabed_point(
            std::ofstream& out,
            index_t vertex_id,
            const vec2& point )
        {
            out << vertex_id << " " << point.x << " " << -point.y;
            save_uncertainties( out );
            out << std::endl;
        }

        void save_uncertainties( std::ofstream& out )
        {
            out << " " << uncertainty_x << " " << uncertainty_y << " "
                << l_uncertainty_x << " " << l_uncertainty_y;
        }

    private:
        static const index_t nb_horizons_for_seabed = 3;
        static const double default_sea_depth;
        /// Sea milieu index is 0
        static const index_t seabed_region = 0;
        static const index_t nb_points_for_seabed_boundaries = 2;

        static const index_t shift_for_seabed = 1;

        static const int outside_region = -1;
        static const double uncertainty_x;
        static const double uncertainty_y;
        static const double l_uncertainty_x;
        static const double l_uncertainty_y;

        vec2 left_up_point_;
        vec2 right_up_point_;

        Box< 2 > model_box_;

    };

    const double StradivariusIOHandler::default_sea_depth = 200;

    const double StradivariusIOHandler::uncertainty_x = -1.0000;
    const double StradivariusIOHandler::uncertainty_y = -1.0000;
    const double StradivariusIOHandler::l_uncertainty_x = .0000;
    const double StradivariusIOHandler::l_uncertainty_y = .0000;

}
