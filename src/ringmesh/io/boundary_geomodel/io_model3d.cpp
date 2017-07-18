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
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
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
     * @brief Total number of polygons in the Surfaces of a geomodel
     */
    index_t nb_polygons( const GeoModel< 3 >& geomodel )
    {
        index_t result = 0;
        for( index_t i : range( geomodel.nb_surfaces() ) ) {
            result += geomodel.surface( i ).nb_mesh_elements();
        }
        return result;
    }

    /*!
     * @brief Write a region information in a stream
     * @details Used by function to save the Model in a .ml file
     *
     * @param[in] count Region index in the file
     * @param[in] region The region to save
     * @param[in,out] out The file output stream
     */
    void save_region( index_t count, const Region< 3 >& region, std::ostream& out )
    {
        out << "REGION " << count << "  " << region.name() << " " << std::endl;
        index_t it = 0;

        for( index_t i : range( region.nb_boundaries() ) ) {
            out << "  ";
            if( region.side( i ) ) {
                out << "+";
            } else {
                out << "-";
            }
            out << region.boundary( i ).index() + 1;
            it++;
            if( it == 5 ) {
                out << std::endl;
                it = 0;
            }
        }
        out << "  0" << std::endl;
    }

    void save_universe(
        index_t count,
        const Universe< 3 >& universe,
        std::ostream& out )
    {
        out << "REGION " << count << "  " << universe.type_name() << " "
            << std::endl;
        index_t it = 0;

        for( index_t i : range( universe.nb_boundaries() ) ) {
            out << "  ";
            if( universe.side( i ) ) {
                out << "+";
            } else {
                out << "-";
            }
            out << universe.boundary_gmme( i ).index() + 1;
            it++;
            if( it == 5 ) {
                out << std::endl;
                it = 0;
            }
        }
        out << "  0" << std::endl;
    }

    /*!
     * @brief Write information for on layer in a stream
     * @details Used by function to save the Model in a .ml file
     *
     * @param[in] count Index of the layer in the file
     * @param[in] offset Offset of region indices in the file
     * @param[in] layer The layer to write
     * @param[in,out] out The output file stream
     */
    void save_layer(
        index_t count,
        index_t offset,
        const GeoModelGeologicalEntity< 3 >& layer,
        std::ostream& out )
    {
        out << "LAYER " << layer.name() << " " << std::endl;
        index_t it = 0;

        for( index_t i : range( layer.nb_children() ) ) {
            out << "  " << layer.child_gmme( i ).index() + offset + 1;
            it++;
            if( it == 5 ) {
                out << std::endl;
                it = 0;
            }
        }
        out << "  0" << std::endl;
    }

    /*!
     * @brief Write basic header for Gocad coordinate system.
     * @param[in,out] out Output .ml file stream
     */
    void save_coordinate_system( std::ostream& out )
    {
        out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl << "NAME Default"
            << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
            << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl << "ZPOSITIVE Elevation"
            << std::endl << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl;
    }

    /*!
     * @brief Check if the geomodel can be saved in a skua-gocad .ml file
     * @details It assumes that the geomodel is valid and verifies that:
     *   - all Interfaces have a name and geological feature
     *   - all Surfaces are in an Interface
     *   - all Surfaces are triangulated
     *   - all Regions have a name
     */
    bool check_gocad_validity( const GeoModel< 3 >& geomodel )
    {
        index_t nb_interfaces = geomodel.nb_geological_entities(
            Interface < 3 > ::type_name_static() );
        if( nb_interfaces == 0 ) {
            Logger::err( "", " The GeoModel ", geomodel.name(),
                " has no Interface" );
            return false;
        }
        for( index_t i : range( nb_interfaces ) ) {
            const GeoModelGeologicalEntity< 3 >& E = geomodel.geological_entity(
                Interface < 3 > ::type_name_static(), i );
            if( !E.has_geological_feature() ) {
                Logger::err( "", E.gmge(), " has no geological feature" );
                return false;
            }
        }
        for( index_t s : range( geomodel.nb_surfaces() ) ) {
            const Surface< 3 >& S = geomodel.surface( s );
            if( !S.has_parent() ) {
                Logger::err( "", S.gmme(),
                    " does not belong to any Interface of the geomodel" );
                return false;
            }
            if( !S.is_simplicial() ) {
                Logger::err( "", S.gmme(), " is not triangulated " );
                return false;
            }
        }
        return true;
    }

    /*! Brute force inefficient but I am debugging !!!! */
    bool has_surface_edge(
        const Surface< 3 >& surface,
        index_t v0_in,
        index_t v1_in )
    {
        for( index_t i : range( surface.nb_mesh_elements() ) ) {
            for( index_t j : range( surface.nb_mesh_element_vertices( i ) ) ) {
                index_t v0 = surface.mesh_element_vertex_index(
                    ElementLocalVertex( i, j ) );
                index_t v1 = surface.mesh_element_vertex_index(
                    surface.low_level_mesh_storage().next_polygon_vertex(
                        ElementLocalVertex( i, j ) ) );
                if( ( v0 == v0_in && v1 == v1_in )
                    || ( v0 == v1_in && v1 == v0_in ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    /*!
     * @brief Save the geomodel in a .ml file if it can
     * @param[in] geomodel the geomodel to save
     * @param[in,out] out Output file stream
     */
    void save_gocad_model3d( const GeoModel< 3 >& geomodel, std::ostream& out )
    {
        if( !check_gocad_validity( geomodel ) ) {
            throw RINGMeshException( "I/O",
                " The GeoModel " + geomodel.name()
                    + +" cannot be saved in .ml format" );
        }
        out.precision( 16 );

        // Gocad Model3d headers
        out << "GOCAD Model3d 1" << std::endl << "HEADER {" << std::endl << "name: "
            << geomodel.name() << std::endl << "}" << std::endl;

        save_coordinate_system( out );

        // Gocad::TSurf = RINGMesh::Interface
        index_t nb_interfaces = geomodel.nb_geological_entities(
            Interface < 3 > ::type_name_static() );
        for( index_t i : range( nb_interfaces ) ) {
            out << "TSURF "
                << geomodel.geological_entity( Interface < 3 > ::type_name_static(),
                    i ).name() << std::endl;
        }

        index_t count = 1;

        // Gocad::TFace = RINGMesh::Surface
        for( index_t s : range( geomodel.nb_surfaces() ) ) {
            const Surface< 3 >& cur_surface = geomodel.surface( s );
            const gmge_id& parent_interface = cur_surface.parent_gmge(
                Interface < 3 > ::type_name_static() );
            if( !parent_interface.is_defined() ) {
                throw RINGMeshException( "I/O",
                    "Failed to save GeoModel" " in .ml Gocad format "
                        "because Surface " + std::to_string( s )
                        + " has no Interface parent)" );
            }
            const GeoModelGeologicalEntity< 3 >::GEOL_FEATURE& cur_geol_feature =
                geomodel.geological_entity( parent_interface ).geological_feature();

            out << "TFACE " << count << "  ";
            out << GeoModelGeologicalEntity < 3 > ::geol_name( cur_geol_feature );
            out << " "
                << cur_surface.parent( Interface < 3 > ::type_name_static() ).name()
                << std::endl;

            // Print the key polygon which is the first three
            // vertices of the first polygon
            out << "  "
                << cur_surface.mesh_element_vertex( ElementLocalVertex( 0, 0 ) )
                << std::endl;
            out << "  "
                << cur_surface.mesh_element_vertex( ElementLocalVertex( 0, 1 ) )
                << std::endl;
            out << "  "
                << cur_surface.mesh_element_vertex( ElementLocalVertex( 0, 2 ) )
                << std::endl;

            ++count;
        }
        // Universe
        index_t offset_layer = count;
        save_universe( count, geomodel.universe(), out );
        ++count;
        // Regions
        for( index_t i : range( geomodel.nb_regions() ) ) {
            save_region( count, geomodel.region( i ), out );
            ++count;
        }
        // Layers
        if( geomodel.entity_type_manager().geological_entity_manager.is_valid_type(
            Layer < 3 > ::type_name_static() ) ) {
            index_t nb_layers = geomodel.nb_geological_entities(
                Layer < 3 > ::type_name_static() );
            for( index_t i : range( nb_layers ) ) {
                save_layer( count, offset_layer,
                    geomodel.geological_entity( Layer < 3 > ::type_name_static(),
                        i ), out );
                ++count;
            }
        }
        out << "END" << std::endl;

        const GeoModelMeshVertices< 3 >& geomodel_vertices = geomodel.mesh.vertices;
        // Save the geometry of the Surfaces, Interface per Interface
        for( index_t i : range( nb_interfaces ) ) {
            const GeoModelGeologicalEntity< 3 >& tsurf = geomodel.geological_entity(
                Interface < 3 > ::type_name_static(), i );
            // TSurf beginning header
            out << "GOCAD TSurf 1" << std::endl << "HEADER {" << std::endl << "name:"
                << tsurf.name() << std::endl << "name_in_model_list:" << tsurf.name()
                << std::endl << "}" << std::endl;
            save_coordinate_system( out );

            out << "GEOLOGICAL_FEATURE " << tsurf.name() << std::endl
                << "GEOLOGICAL_TYPE ";
            out << GeoModelGeologicalEntity < 3
                > ::geol_name( tsurf.geological_feature() );
            out << std::endl;
            out << "PROPERTY_CLASS_HEADER Z {" << std::endl << "is_z:on" << std::endl
                << "}" << std::endl;

            index_t vertex_count = 1;
            // TFace vertex index = Surface vertex index + offset
            index_t offset = vertex_count;

            // To collect Corners(BStones) indexes
            // and boundary (Line) first and second vertex indexes
            std::set< index_t > corners;
            std::set< std::pair< index_t, index_t > > lineindices;
            for( index_t j : range( tsurf.nb_children() ) ) {
                offset = vertex_count;
                const Surface< 3 >& surface =
                    dynamic_cast< const Surface< 3 >& >( tsurf.child( j ) );

                out << "TFACE" << std::endl;
                for( index_t k : range( surface.nb_vertices() ) ) {
                    out << "VRTX " << vertex_count << " " << surface.vertex( k )
                        << std::endl;
                    vertex_count++;
                }
                for( index_t k : range( surface.nb_mesh_elements() ) ) {
                    out << "TRGL "
                        << surface.mesh_element_vertex_index(
                            ElementLocalVertex( k, 0 ) ) + offset << " "
                        << surface.mesh_element_vertex_index(
                            ElementLocalVertex( k, 1 ) ) + offset << " "
                        << surface.mesh_element_vertex_index(
                            ElementLocalVertex( k, 2 ) ) + offset << std::endl;
                }
                for( index_t k : range( surface.nb_boundaries() ) ) {
                    const Line< 3 >& line = surface.boundary( k );
                    index_t v0_model_id = geomodel_vertices.geomodel_vertex_id(
                        line.gmme(), 0 );
                    index_t v1_model_id = geomodel_vertices.geomodel_vertex_id(
                        line.gmme(), 1 );

                    std::vector< index_t > v0_surface_ids =
                        geomodel_vertices.mesh_entity_vertex_id( surface.gmme(),
                            v0_model_id );
                    std::vector< index_t > v1_surface_ids =
                        geomodel_vertices.mesh_entity_vertex_id( surface.gmme(),
                            v1_model_id );

                    if( !surface.has_inside_border() ) {
                        index_t v0 = v0_surface_ids[0];
                        index_t v1 = v1_surface_ids[0];
                        v0 += offset;
                        v1 += offset;

                        lineindices.insert(
                            std::pair< index_t, index_t >( v0, v1 ) );
                        corners.insert( v0 );
                    } else {
                        // We need to get the right pair of v0 - v1  (not crossing the inside boundary)
                        // corner and a border
                        int count = 0;
                        bool to_break = false;
                        for( index_t v0 : v0_surface_ids ) {
                            for( index_t v1 : v1_surface_ids ) {
                                if( has_surface_edge( surface, v0, v1 ) ) {
                                    lineindices.insert(
                                        std::pair< index_t, index_t >( v0 + offset,
                                            v1 + offset ) );
                                    count++;
                                }
                                if( !line.is_inside_border( surface )
                                    && count == 1 ) {
                                    to_break = true;
                                    break;
                                } else if( count == 2 ) {
                                    to_break = true;
                                    break;
                                }
                            }
                            if( to_break ) {
                                corners.insert( v0 + offset );
                                break;
                            }
                        }
                    }
                    // Set a BSTONE at the line other extremity
                    const gmme_id& c1_id = line.boundary_gmme( 1 );
                    std::vector< index_t > gme_vertices =
                        geomodel_vertices.mesh_entity_vertex_id( surface.gmme(),
                            geomodel_vertices.geomodel_vertex_id( c1_id ) );
                    corners.insert( gme_vertices.front() + offset );
                }
            }
            // Add the remaining bstones that are not already in bstones
            for( auto it( corners.begin() ); it != corners.end(); ++it ) {
                out << "BSTONE " << *it << std::endl;
            }
            for( auto it( lineindices.begin() ); it != lineindices.end(); ++it ) {
                out << "BORDER " << vertex_count << " " << it->first << " "
                    << it->second << std::endl;
                vertex_count++;
            }
            out << "END" << std::endl;
        }
    }

    class MLIOHandler final: public GeoModelIOHandler< 3 > {
    public:
        /*! Load a .ml (Gocad file)
         * @pre Filename is valid
         */
        void load( const std::string& filename, GeoModel< 3 >& geomodel ) final
        {
            std::ifstream input( filename.c_str() );
            if( !input ) {
                throw RINGMeshException( "I/O", "Failed to open file " + filename );
            }
            GeoModelBuilderML builder( geomodel, filename );
            builder.build_geomodel();
        }

        void save( const GeoModel< 3 >& geomodel, const std::string& filename ) final
        {

            std::ofstream out( filename.c_str() );
            save_gocad_model3d( geomodel, out );
        }
    };
}
