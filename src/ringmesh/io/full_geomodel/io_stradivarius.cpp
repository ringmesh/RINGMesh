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
        using GeoModelBuilderFile< 2 >::GeoModelBuilderFile;
        static const index_t SHIFT = 1;
        static const index_t SEABED = 0;

    private:
        void load_file() override
        {
            GEO::LineInput file( filename_ );
            while( !file.eof() && file.get_line() ) {
                file.get_fields();
                DEBUG( file.current_line() );
                if( file.field_matches( 0, "//" ) ) {
                    continue;
                } else if( file.field_matches( 0, "liste" ) ) {
                    if( file.nb_fields() != 3 ) {
                        continue;
                    }
                    if( file.field_matches( 2, "milieux" ) ) {
                        DEBUG( "MILIEU" );
                        load_surfaces( file );
                    } else if( file.field_matches( 2, "horizons" ) ) {
                        DEBUG( "HORIZON" );
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
            if( right_surface != NO_ID && right_surface != SEABED  ) {
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
                DEBUG( s );
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

            DEBUG( geomodel.mesh.vertices.nb() );

            for( index_t v = 0; v < geomodel.mesh.vertices.nb(); v++ ) {
                DEBUG( geomodel.mesh.vertices.vertex( v ) );
            }
        }
        void save( const GeoModel< 2 >& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

        }

    private:
        /*!
         * @brief Find the gmsh type on an element using
         * the number of vertices and the mesh entity index
         * in which the element belong
         */
        index_t find_gmsh_element_type(
            index_t nb_vertices,
            index_t mesh_entity_type_index )
        {
            return element_type[nb_vertices + mesh_entity_type_index];
        }

        /*!
         * @brief Find the gmsh local vertex index of an element using
         * the number of vertices and the mesh entity index
         * in which the element belong
         */
        index_t find_gmsh_element_local_vertex_id(
            index_t nb_vertices,
            index_t mesh_entity_type_index,
            index_t local_vertex_index )
        {
            return vertices_in_elements[nb_vertices + mesh_entity_type_index][local_vertex_index];
        }

        /*!
         * @brief Count all the elements in GMSH
         * an Element can be :
         * - A Corner
         * - An Edge
         * - A Triangle
         * - A Quad
         * - A Tetrahedron
         * - A Pyramid
         * - A Prism
         * - An Hexaedron
         */
        index_t count_elements( const GeoModel< 2 >& geomodel )
        {
            index_t nb_elements = 0;
            const std::vector< MeshEntityType >& gmme_types =
                geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();
            for( MeshEntityType cur_mesh_entity_type : gmme_types ) {
                for( index_t index_of_gmme_of_the_current_type : range(
                    geomodel.nb_mesh_entities( cur_mesh_entity_type ) ) ) {
                    gmme_id cur_gmme_id = gmme_id( cur_mesh_entity_type,
                        index_of_gmme_of_the_current_type );
                    const GeoModelMeshEntity< 2 >& cur_gmme = geomodel.mesh_entity(
                        cur_gmme_id );
                    nb_elements += cur_gmme.nb_mesh_elements();
                }
            }
            return nb_elements;
        }
    };
}
