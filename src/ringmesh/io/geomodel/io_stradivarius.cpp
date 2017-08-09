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

    /*!
     * @brief Import for the Stradivarius .model format
     */
    class StradivariusBuilder final: public GeoModelBuilderFile< 2 > {
    public:
        StradivariusBuilder( GeoModel2D& geomodel, std::string filename )
            : GeoModelBuilderFile< 2 >( geomodel, std::move( filename ) )
        {
        }

    private:
        void load_file() override
        {
            load_points();
            load_interfaces();
            load_media();

            horizon_m0_.insert(
                gmme_id( Surface2D::type_name_static(), index_t( 0 ) ) );
            removal.remove_mesh_entities( horizon_m0_ );
            build_corners_from_lines();
        }

        void load_points()
        {
            GEO::LineInput file( filename_ );
            move_to( file, "liste des points\n" );
            file.get_line();
            file.get_fields();
            index_t nb_points = file.field_as_uint( 0 );
            points_.resize( nb_points );
            for( index_t p : range( nb_points ) ) {
                file.get_line();
                file.get_fields();
                points_[p] = {file.field_as_double( 1 ), - file.field_as_double( 2 )};
            }
        }

        void load_interfaces() {
            GEO::LineInput file( filename_ );
            move_to( file, "liste des horizons\n" );
            file.get_line();
            file.get_fields();
            index_t nb_horizons = file.field_as_uint( 0 );
            topology.create_mesh_entities( Line2D::type_name_static(), nb_horizons );

            for (index_t horizon_id : range( nb_horizons )) {
                file.get_line();
                file.get_fields();

                gmme_id horizon( Line2D::type_name_static(), horizon_id );
                int medium_1 {file.field_as_int( 0 )};
                int medium_2 {file.field_as_int( 1 )};
                index_t nb_points = file.field_as_uint( 2 );
                info.set_mesh_entity_name(horizon, file.field( 4 ));

                std::vector<vec2> vertices(nb_points);
                for (index_t point_i : range( nb_points )) {
                    file.get_line();
                    file.get_fields();
                    index_t point_id = file.field_as_uint( 0 );
                    vertices[point_i] = points_[point_id];
                }
                geometry.set_line(horizon_id, vertices);

                if (((medium_1 == 0 && (medium_2 == -1))) || ((medium_1 == -1 && (medium_2 == 0)))) {
                    horizon_m0_.insert(horizon);
                }
            }
        }

        void load_media() {
            GEO::LineInput file( filename_ );
            move_to(file, "liste des milieux\n");
            file.get_line();
            file.get_fields();
            index_t nb_media = file.field_as_uint( 0 );
            topology.create_mesh_entities( Surface2D::type_name_static(), nb_media + 1 );

            for( const gmme_id& horizon : horizon_m0_ ) {
                topology.add_mesh_entity_boundary_relation(
                    {   Surface2D::type_name_static(), index_t( 0 )}, horizon,
                    true );
            }

            for (index_t milieu_i : range( nb_media )) {
                file.get_line();
                file.get_fields();
                index_t nb_interfaces = file.field_as_uint( 1 );
                for (index_t interface_i : range( nb_interfaces )) {
                    ringmesh_unused( interface_i );
                    file.get_line();
                    file.get_fields();
                    index_t interface_id = file.field_as_uint( 0 );
                    gmme_id horizon( Line2D::type_name_static(), interface_id );
                    topology.add_mesh_entity_boundary_relation(
                        {   Surface2D ::type_name_static(), milieu_i + 1}, horizon,
                        true );
                }
            }
        }

        bool move_to(GEO::LineInput& file, std::string selector) {
            while (!file.eof() && file.get_line()) {
                if (selector.compare(file.current_line()) == 0) {
                    return true;
                }
            }
            return false;
        }

    private:
        std::vector< vec2 > points_;
        std::set< gmme_id > horizon_m0_;
    };

    /*!
     * @brief Export for the Stradivarius .model format
     */
    class StradivariusIOHandler final: public GeoModelIOHandler< 2 > {
    public:
        void load( const std::string& filename, GeoModel2D& geomodel ) final
        {
            StradivariusBuilder builder( geomodel, filename );
            builder.build_geomodel();
        }
        void save( const GeoModel2D& geomodel, const std::string& filename ) final
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
        void check_stradivarius_validity( const GeoModel2D& geomodel )
        {
            for( index_t v : range( geomodel.mesh.vertices.nb() ) ) {
                model_box_.add_point( geomodel.mesh.vertices.vertex( v ) );
            }
        }
        void save_header( std::ofstream& out )
        {
            out << "// File wrote with RINGMesh" << std::endl;
        }

        void save_milieux( std::ofstream& out, const GeoModel2D& geomodel )
        {
            out << "liste des milieux" << std::endl;
            out << geomodel.nb_surfaces() << std::endl;
            for( const auto& surface : geomodel.surfaces() ) {
                save_milieu( out, surface );
            }
        }

        void save_milieu( std::ofstream& out, const Surface2D& surface )
        {
            save_milieu_index( out, surface );
            save_milieu_boundaries( out, surface );
        }

        void save_milieu_index( std::ofstream& out, const Surface2D& surface )
        {
            out << surface.index() + shift_for_seabed << " "
                << surface.nb_boundaries() << std::endl;
        }

        void save_milieu_boundaries( std::ofstream& out, const Surface2D& surface )
        {
            for( index_t boundary : range( surface.nb_boundaries() ) ) {
                out << surface.boundary( boundary ).index() << std::endl;
            }
        }

        /*!
         * @brief Here horizons is in term of stradivarius nomenclature
         */
        void save_horizons( std::ofstream& out, const GeoModel2D& geomodel )
        {
            out << "liste des horizons" << std::endl;
            out << geomodel.nb_lines() + nb_horizons_for_seabed << " "
                << "// milG milD nbrePoint   numeroHorizon   nom d'horizon"
                << std::endl;
            for( const auto& line : geomodel.lines() ) {
                save_horizon( out, line );
            }
            save_seabed( out, geomodel );

        }

        void save_horizon( std::ofstream& out, const Line2D& line )
        {
            if( line.nb_incident_entities() == 2 ) {
                save_classical_horizon( out, line );

            } else if( line.nb_incident_entities() == 1 ) {
                save_horizon_as_boundary( out, line );
            } else {
                ringmesh_assert_not_reached;
            }
        }

        void save_classical_horizon( std::ofstream& out, const Line2D& line )
        {
            const Surface2D &incident_surface_1 = line.incident_entity( 0 );
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

        void save_horizon_as_boundary( std::ofstream& out, const Line2D& line )
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

        void left_is_outside( std::ofstream& out, const Line2D& line )
        {
            save_horizon_topology( out, outside_region,
                line.incident_entity( 0 ).index() + shift_for_seabed,
                line.nb_vertices(), line.index(),
                "Line_" + GEO::String::to_string( line.index() ) );
        }

        void right_is_outside( std::ofstream& out, const Line2D& line )
        {
            save_horizon_topology( out,
                line.incident_entity( 0 ).index() + shift_for_seabed, outside_region,
                line.nb_vertices(), line.index(),
                "Line_" + GEO::String::to_string( line.index() ) );
        }

        void save_horizons_points( std::ofstream& out, const Line2D& line )
        {
            for( index_t vertex : range( line.nb_vertices() ) ) {
                out
                    << line.geomodel().mesh.vertices.geomodel_vertex_id( line.gmme(),
                        vertex ) << std::endl;
            }
        }

        void save_seabed( std::ofstream& out, const GeoModel2D& geomodel )
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
            const GeoModel2D& geomodel )
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

        void save_faults( std::ofstream& out, const GeoModel2D& geomodel )
        {
            out << "liste des failles" << std::endl;
            out << 0 << std::endl;
        }

        void save_vertices( std::ofstream& out, const GeoModel2D& geomodel )
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

        Box2D model_box_;

    };

    const double StradivariusIOHandler::default_sea_depth = 200;

    const double StradivariusIOHandler::uncertainty_x { -1. };
    const double StradivariusIOHandler::uncertainty_y { -1. };
    const double StradivariusIOHandler::l_uncertainty_x { 0. };
    const double StradivariusIOHandler::l_uncertainty_y { 0. };
}

