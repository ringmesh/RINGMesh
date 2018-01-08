/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*!
 * @author Nicolas Mastio
 */

namespace
{
    /*!
     * @brief Import for the Stradivarius .model format
     */
    class StradivariusBuilder final : public GeoModelBuilderFile< 2 >
    {
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

            horizon_m0_.emplace( Surface2D::type_name_static(), index_t( 0 ) );
            remove.remove_mesh_entities( horizon_m0_ );
            build_corners_from_lines();
        }

        void load_points()
        {
            GEO::LineInput file{ filename() };
            move_to( file, "liste des points\n" );
            file.get_line();
            file.get_fields();
            auto nb_points = file.field_as_uint( 0 );
            points_.resize( nb_points );
            for( auto p : range( nb_points ) )
            {
                file.get_line();
                file.get_fields();
                // convert coordinates from depth to elevation
                points_[p] = { file.field_as_double( 1 ),
                    -file.field_as_double( 2 ) };
            }
        }

        void load_interfaces()
        {
            GEO::LineInput file{ filename() };
            move_to( file, "liste des horizons\n" );
            file.get_line();
            file.get_fields();
            auto nb_horizons = file.field_as_uint( 0 );
            topology.create_mesh_entities(
                Line2D::type_name_static(), nb_horizons );

            for( auto horizon_id : range( nb_horizons ) )
            {
                file.get_line();
                file.get_fields();

                gmme_id horizon{ Line2D::type_name_static(), horizon_id };
                int medium_1{ file.field_as_int( 0 ) };
                int medium_2{ file.field_as_int( 1 ) };
                auto nb_points = file.field_as_uint( 2 );
                info.set_mesh_entity_name( horizon, file.field( 4 ) );

                std::vector< vec2 > vertices( nb_points );
                for( auto point_i : range( nb_points ) )
                {
                    file.get_line();
                    file.get_fields();
                    auto point_id = file.field_as_uint( 0 );
                    vertices[point_i] = points_[point_id];
                }
                geometry.set_line( horizon_id, vertices );

                if( ( ( medium_1 == 0 && ( medium_2 == -1 ) ) )
                    || ( ( medium_1 == -1 && ( medium_2 == 0 ) ) ) )
                {
                    horizon_m0_.insert( horizon );
                }
            }
        }

        void load_media()
        {
            GEO::LineInput file{ filename() };
            move_to( file, "liste des milieux\n" );
            file.get_line();
            file.get_fields();
            auto nb_media = file.field_as_uint( 0 );
            topology.create_mesh_entities(
                Surface2D::type_name_static(), nb_media + 1 );

            for( const auto& horizon : horizon_m0_ )
            {
                topology.add_surface_line_boundary_relation(
                    0, horizon.index(), true );
            }

            for( auto milieu_i : range( nb_media ) )
            {
                file.get_line();
                file.get_fields();
                auto nb_interfaces = file.field_as_uint( 1 );
                for( auto interface_i : range( nb_interfaces ) )
                {
                    ringmesh_unused( interface_i );
                    file.get_line();
                    file.get_fields();
                    auto interface_id = file.field_as_uint( 0 );
                    gmme_id horizon{ Line2D::type_name_static(), interface_id };
                    topology.add_surface_line_boundary_relation(
                        milieu_i + 1, horizon.index(), true );
                }
            }
        }

        bool move_to( GEO::LineInput& file, std::string selector )
        {
            while( !file.eof() && file.get_line() )
            {
                if( selector.compare( file.current_line() ) == 0 )
                {
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
    class StradivariusIOHandler final : public GeoModelInputHandler2D
    {
    public:
        void load( const std::string& filename, GeoModel2D& geomodel ) final
        {
            StradivariusBuilder builder{ geomodel, filename };
            builder.build_geomodel();
        }
    };
}
