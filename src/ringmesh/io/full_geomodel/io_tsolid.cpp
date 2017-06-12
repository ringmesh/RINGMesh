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

    class TSolidIOHandler final : public GeoModelIOHandler {
    public:
        TSolidIOHandler()
            : GeoModelIOHandler(), gocad_no_data_value_( -99999 )
        {
        }
        virtual bool load( const std::string& filename, GeoModel& geomodel ) final
        {
            std::ifstream input( filename.c_str() );
            if( input ) {
                GeoModelBuilderTSolid builder( geomodel, filename );

                time_t start_load, end_load;
                time( &start_load );

                builder.build_geomodel();
                print_geomodel( geomodel );
                // Check boundary geomodel validity
                bool is_valid = is_geomodel_valid( geomodel );

                time( &end_load );

                Logger::out( "I/O", " Loaded geomodel ", geomodel.name(), " from " );
                Logger::out( "I/O", filename, " timing: ",
                    difftime( end_load, start_load ), "sec" );
                return is_valid;
            } else {
                throw RINGMeshException( "I/O",
                    "Failed loading geomodel from file " + filename );
                return false;
            }
        }
        virtual void save( const GeoModel& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

            fill_top_header( geomodel, out );

            const GeoModelMesh& mesh = geomodel.mesh;
            //mesh.set_duplicate_mode( GeoModelMeshCells::ALL ) ;

            // att_v_double_names and att_c_double_names contain all the
            // attribute names found the regions. Only vertex and cell attributes.
            // If an attribute is defined in a region and not in another, in
            // the region where the attribute is not defined no data values
            // are stored.
            std::vector< std::string > att_v_double_names;
            std::vector< index_t > vertex_attr_dims;
            fill_vertex_attribute_header( geomodel, out, att_v_double_names,
                vertex_attr_dims );
            std::vector< std::string > att_c_double_names;
            std::vector< index_t > cell_attr_dims;
            fill_cell_attribute_header( geomodel, out, att_c_double_names,
                cell_attr_dims );

            std::vector< bool > vertex_exported( mesh.vertices.nb(), false );
            std::vector< bool > atom_exported( mesh.cells.nb_duplicated_vertices(),
                false );
            std::vector< index_t > vertex_exported_id( mesh.vertices.nb(), NO_ID );
            std::vector< index_t > atom_exported_id(
                mesh.cells.nb_duplicated_vertices(), NO_ID );
            index_t nb_vertices_exported = 1;
            for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                const RINGMesh::Region& region = geomodel.region( r );
                export_one_region( region, out, nb_vertices_exported,
                    vertex_exported, atom_exported, vertex_exported_id,
                    atom_exported_id, att_v_double_names, vertex_attr_dims,
                    att_c_double_names, cell_attr_dims );
            }

            export_model( geomodel, out, vertex_exported_id );
            export_model_region( geomodel, out );

            out << "END" << std::endl;
        }
    private:
        void fill_top_header( const GeoModel& geomodel, std::ofstream& out ) const
        {
            // Print Model3d headers
            out << "GOCAD TSolid 1" << std::endl << "HEADER {" << std::endl
                << "name:" << geomodel.name() << std::endl << "}" << std::endl;

            out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl << "NAME Default"
                << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
                << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl
                << "ZPOSITIVE Elevation" << std::endl
                << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl;
        }
        void fill_vertex_attribute_header(
            const GeoModel& geomodel,
            std::ofstream& out,
            std::vector< std::string >& att_v_double_names,
            std::vector< index_t >& vertex_attr_dims ) const
        {
            for( index_t reg_i = 0; reg_i < geomodel.nb_regions(); ++reg_i ) {
                const Region& cur_reg = geomodel.region( reg_i );
                GEO::AttributesManager& reg_vertex_attr_mgr =
                    cur_reg.vertex_attribute_manager();
                GEO::vector< std::string > att_v_names;
                reg_vertex_attr_mgr.list_attribute_names( att_v_names );
                ringmesh_assert( att_v_names.size() == reg_vertex_attr_mgr.nb() );
                for( const std::string& cur_att_v_name : att_v_names ) {

                    if( cur_att_v_name == "point" ) {
                        continue;
                    }

                    if( std::find( att_v_double_names.begin(),
                        att_v_double_names.end(), cur_att_v_name )
                        != att_v_double_names.end() ) {
                        continue;
                    }

                    const GEO::AttributeStore* attr_store =
                        reg_vertex_attr_mgr.find_attribute_store( cur_att_v_name );
                    ringmesh_assert( attr_store != nullptr );

                    if( !GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to(
                        attr_store ) ) {
                        continue;
                    }

                    att_v_double_names.push_back( cur_att_v_name );
                    index_t cur_dim = attr_store->dimension();
                    vertex_attr_dims.push_back( cur_dim );
                }
            }

            const index_t nb_att_v_double_names = att_v_double_names.size();
            if( !att_v_double_names.empty() ) {
                ringmesh_assert( nb_att_v_double_names > 0 );
                out << "PROPERTIES";
                for( const std::string& cur_att_v_double_name : att_v_double_names ) {
                    out << " " << cur_att_v_double_name;
                }
                out << std::endl;
                out << "PROP_LEGAL_RANGES";
                for( index_t i = 0; i < nb_att_v_double_names; ++i ) {
                    out << " **none**  **none**";
                }
                out << std::endl;
                out << "NO_DATA_VALUES";
                for( index_t i = 0; i < nb_att_v_double_names; ++i ) {
                    out << " " << GEO::String::to_string( gocad_no_data_value_ );
                }
                out << std::endl;
                out << "READ_ONLY";
                for( index_t i = 0; i < nb_att_v_double_names; ++i ) {
                    out << " 1";
                }
                out << std::endl;
                out << "PROPERTY_CLASSES";
                for( const std::string& cur_att_v_double_name : att_v_double_names ) {
                    out << " " << cur_att_v_double_name;
                }
                out << std::endl;
                out << "PROPERTY_KINDS";
                for( index_t i = 0; i < nb_att_v_double_names; ++i ) {
                    out << " \"Real Number\"";
                }
                out << std::endl;
                out << "PROPERTY_SUBCLASSES";
                for( index_t i = 0; i < nb_att_v_double_names; ++i ) {
                    out << " QUANTITY Float";
                }
                out << std::endl;
                out << "ESIZES";
                for( const index_t& cur_v_att_dim : vertex_attr_dims ) {
                    out << " " << GEO::String::to_string( cur_v_att_dim );
                }
                out << std::endl;
                out << "UNITS";
                for( index_t i = 0; i < nb_att_v_double_names; ++i ) {
                    out << " unitless";
                }
                out << std::endl;
                for( const std::string& cur_att_v_double_name : att_v_double_names ) {
                    out << "PROPERTY_CLASS_HEADER " << cur_att_v_double_name << " {"
                        << std::endl;
                    out << "kind: Real Number" << std::endl;
                    out << "unit: unitless" << std::endl;
                    out << "}" << std::endl;
                }
            }
        }

        void fill_cell_attribute_header(
            const GeoModel& geomodel,
            std::ofstream& out,
            std::vector< std::string >& att_c_double_names,
            std::vector< index_t >& cell_attr_dims ) const
        {
            for( index_t reg_i = 0; reg_i < geomodel.nb_regions(); ++reg_i ) {
                const Region& cur_reg = geomodel.region( reg_i );
                GEO::AttributesManager& reg_cell_attr_mgr =
                    cur_reg.cell_attribute_manager();
                GEO::vector< std::string > att_c_names;
                reg_cell_attr_mgr.list_attribute_names( att_c_names );
                ringmesh_assert( att_c_names.size() == reg_cell_attr_mgr.nb() );
                for( const std::string& cur_att_c_name : att_c_names ) {

                    if( std::find( att_c_double_names.begin(),
                        att_c_double_names.end(), cur_att_c_name )
                        != att_c_double_names.end() ) {
                        continue;
                    }

                    const GEO::AttributeStore* attr_store =
                        reg_cell_attr_mgr.find_attribute_store( cur_att_c_name );
                    ringmesh_assert( attr_store != nullptr );

                    if( !GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to(
                        attr_store ) ) {
                        continue;
                    }

                    att_c_double_names.push_back( cur_att_c_name );
                    index_t cur_dim = attr_store->dimension();
                    cell_attr_dims.push_back( cur_dim );
                }
            }

            const index_t nb_att_c_double_names = att_c_double_names.size();
            if( !att_c_double_names.empty() ) {
                out << "TETRA_PROPERTIES";
                for( const std::string& cur_att_c_double_name : att_c_double_names ) {
                    out << " " << cur_att_c_double_name;
                }
                out << std::endl;
                out << "TETRA_PROP_LEGAL_RANGES";
                for( index_t i = 0; i < nb_att_c_double_names; ++i ) {
                    out << " **none**  **none**";
                }
                out << std::endl;
                out << "TETRA_NO_DATA_VALUES";
                for( index_t i = 0; i < nb_att_c_double_names; ++i ) {
                    out << " " << GEO::String::to_string( gocad_no_data_value_ );
                }
                out << std::endl;
                out << "READ_ONLY";
                for( index_t i = 0; i < nb_att_c_double_names; ++i ) {
                    out << " 1";
                }
                out << std::endl;
                out << "TETRA_PROPERTY_CLASSES";
                for( const std::string& cur_att_c_double_name : att_c_double_names ) {
                    out << " " << cur_att_c_double_name;
                }
                out << std::endl;
                out << "TETRA_PROPERTY_KINDS";
                for( index_t i = 0; i < nb_att_c_double_names; ++i ) {
                    out << " \"Real Number\"";
                }
                out << std::endl;
                out << "TETRA_PROPERTY_SUBCLASSES";
                for( index_t i = 0; i < nb_att_c_double_names; ++i ) {
                    out << " QUANTITY Float";
                }
                out << std::endl;
                out << "TETRA_ESIZES";
                for( const index_t& cur_cell_attr_dim : cell_attr_dims ) {
                    out << " " << GEO::String::to_string( cur_cell_attr_dim );
                }
                out << std::endl;
                out << "TETRA_UNITS";
                for( index_t i = 0; i < nb_att_c_double_names; ++i ) {
                    out << " unitless";
                }
                out << std::endl;
                for( const std::string& cur_att_c_double_name : att_c_double_names ) {
                    out << "TETRA_PROPERTY_CLASS_HEADER " << cur_att_c_double_name
                        << " {" << std::endl;
                    out << "kind: Real Number" << std::endl;
                    out << "unit: unitless" << std::endl;
                    out << "}" << std::endl;
                }
            }
        }

        void export_one_region(
            const RINGMesh::Region& region,
            std::ofstream& out,
            index_t& nb_vertices_exported,
            std::vector< bool >& vertex_exported,
            std::vector< bool >& atom_exported,
            std::vector< index_t >& vertex_exported_id,
            std::vector< index_t >& atom_exported_id,
            const std::vector< std::string >& att_v_double_names,
            const std::vector< index_t >& vertex_attr_dims,
            const std::vector< std::string >& att_c_double_names,
            const std::vector< index_t >& cell_attr_dims ) const
        {
            out << "TVOLUME " << region.name() << std::endl;
            export_region_vertices( region, out, nb_vertices_exported,
                vertex_exported, vertex_exported_id, att_v_double_names,
                vertex_attr_dims );
            export_tetrahedra( region, out, vertex_exported_id, atom_exported_id,
                att_c_double_names, cell_attr_dims );
        }

        void export_region_vertices(
            const RINGMesh::Region& region,
            std::ofstream& out,
            index_t& nb_vertices_exported,
            std::vector< bool >& vertex_exported,
            std::vector< index_t >& vertex_exported_id,
            const std::vector< std::string >& att_v_double_names,
            const std::vector< index_t >& vertex_attr_dims ) const
        {
            GEO::AttributesManager& reg_vertex_attr_mgr =
                region.vertex_attribute_manager();
            const GeoModelMesh& mesh = region.geomodel().mesh;
            // Export not duplicated vertices
            for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
                index_t cell = mesh.cells.cell( region.gmme().index(), c );
                vec3 cell_center = mesh.cells.barycenter( cell );
                for( index_t v = 0; v < mesh.cells.nb_vertices( cell ); v++ ) {
                    index_t atom_id;
                    if( !mesh.cells.is_corner_duplicated( cell, v, atom_id ) ) {
                        index_t vertex_id = mesh.cells.vertex( cell, v );
                        if( vertex_exported[vertex_id] ) continue;
                        vertex_exported[vertex_id] = true;
                        vertex_exported_id[vertex_id] = nb_vertices_exported;
                        // PVRTX keyword must be used instead of VRTX keyword because
                        // properties are not read by Gocad if it is VRTX keyword.
                        out << "PVRTX " << nb_vertices_exported++ << " "
                            << mesh.vertices.vertex( vertex_id );

                        /// Export of vertex attributes
                        index_t vertex_id_in_reg = NO_ID;
                        bool vertex_id_in_reg_found = false;
                        /// As we export the non duplicated vertices,
                        /// gme_vertices should be at size 1 (to check),
                        /// so the loop should not be needed but I [BC]
                        /// keep it for now since duplicated nodes is not
                        /// operational yet...
                        const std::vector< GMEVertex >& gme_vertices =
                            mesh.vertices.gme_vertices( vertex_id );
                        for( const GMEVertex& cur_gme_vertex : gme_vertices ) {
                            if( cur_gme_vertex.gmme != region.gmme() ) {
                                continue;
                            }
                            std::vector< index_t > cells_around_vertex =
                                region.cells_around_vertex( cur_gme_vertex.v_index,
                                    NO_ID );
                            /// WARNING: the cell id in the region corresponding
                            /// to the cell id in the GMM "cell" is not the
                            /// variable "c" (in the for loop over the region cells).
                            for( const index_t& cur_cell_around_vertex : cells_around_vertex ) {
                                vec3 center = region.mesh_element_barycenter(
                                    cur_cell_around_vertex );
                                if( ( center - cell_center ).length()
                                    < region.geomodel().epsilon() ) {
                                    vertex_id_in_reg = cur_gme_vertex.v_index;
                                    vertex_id_in_reg_found = true;
                                    break;
                                }
                            }
                            if( vertex_id_in_reg_found ) {
                                break;
                            }
                        }
                        for( index_t attr_dbl_itr = 0;
                            attr_dbl_itr < att_v_double_names.size();
                            ++attr_dbl_itr ) {
                            if( reg_vertex_attr_mgr.is_defined(
                                att_v_double_names[attr_dbl_itr] ) ) {
                                const GEO::AttributeStore* attr_store =
                                    reg_vertex_attr_mgr.find_attribute_store(
                                        att_v_double_names[attr_dbl_itr] );
                                ringmesh_assert( attr_store != nullptr );
                                ringmesh_assert(
                                    attr_store->dimension()
                                        == vertex_attr_dims[attr_dbl_itr] );
                                ringmesh_assert(
                                    GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to(
                                        attr_store ) );
                                GEO::ReadOnlyScalarAttributeAdapter cur_attr(
                                    reg_vertex_attr_mgr,
                                    att_v_double_names[attr_dbl_itr] );
                                for( index_t dim_itr = 0;
                                    dim_itr < vertex_attr_dims[attr_dbl_itr];
                                    ++dim_itr ) {
                                    out << " "
                                        << cur_attr[vertex_id_in_reg
                                            * vertex_attr_dims[attr_dbl_itr]
                                            + dim_itr];
                                }
                            } else {
                                for( index_t dim_itr = 0;
                                    dim_itr < vertex_attr_dims[attr_dbl_itr];
                                    ++dim_itr ) {
                                    out << " "
                                        << GEO::String::to_string(
                                            gocad_no_data_value_ );
                                }
                            }
                        }
                        out << std::endl;
                    }
                }
            }
        }
        void export_tetrahedra(
            const RINGMesh::Region& region,
            std::ofstream& out,
            std::vector< index_t >& vertex_exported_id,
            std::vector< index_t >& atom_exported_id,
            const std::vector< std::string >& att_c_double_names,
            const std::vector< index_t >& cell_attr_dims ) const
        {
            // Export duplicated vertices
            /*for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
             index_t cell = mesh.cells.cell( r, c ) ;
             for( index_t v = 0; v < mesh.cells.nb_vertices( cell ); v++ ) {
             index_t atom_id ;
             if( mesh.cells.is_corner_duplicated( cell, v, atom_id ) ) {
             if( atom_exported[atom_id] ) continue ;
             atom_exported[atom_id] = true ;
             atom_exported_id[atom_id] = nb_vertices_exported ;
             index_t vertex_id = mesh.cells.vertex( cell, v ) ;
             out << "ATOM " << nb_vertices_exported++ << " "
             << vertex_exported_id[vertex_id] << std::endl ;
             }
             }
             }*/

            // Mark if a boundary is ending in the region
            std::map< index_t, index_t > sides;
            for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                if( sides.count( region.boundary_gmme( s ).index() ) > 0 )
                    // a surface is encountered twice, it is ending in the region
                    sides[region.boundary_gmme( s ).index()] = 2;
                else
                    sides[region.boundary_gmme( s ).index()] = region.side( s );
            }

            const GeoModelMesh& mesh = region.geomodel().mesh;
            GEO::AttributesManager& reg_cell_attr_mgr =
                region.cell_attribute_manager();
            for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
                out << "TETRA";
                index_t cell = mesh.cells.cell( region.gmme().index(), c );
                for( index_t v = 0;
                    v < region.geomodel().mesh.cells.nb_vertices( cell ); v++ ) {
                    index_t atom_id;
                    if( !mesh.cells.is_corner_duplicated( cell, v, atom_id ) ) {
                        index_t vertex_id = mesh.cells.vertex( cell, v );
                        out << " " << vertex_exported_id[vertex_id];
                    } else {
                        out << " " << atom_exported_id[atom_id];
                    }
                }

                /// Export cell attributes
                vec3 center = mesh.cells.barycenter( cell );
                const std::vector< index_t > c_in_reg =
                    region.cell_nn_search().get_neighbors( center,
                        region.geomodel().epsilon() );
                ringmesh_assert( c_in_reg.size() == 1 );
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    if( reg_cell_attr_mgr.is_defined(
                        att_c_double_names[attr_dbl_itr] ) ) {
                        const GEO::AttributeStore* attr_store =
                            reg_cell_attr_mgr.find_attribute_store(
                                att_c_double_names[attr_dbl_itr] );
                        ringmesh_assert( attr_store != nullptr );
                        ringmesh_assert(
                            attr_store->dimension()
                                == cell_attr_dims[attr_dbl_itr] );
                        ringmesh_assert(
                            GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to(
                                attr_store ) );
                        GEO::ReadOnlyScalarAttributeAdapter cur_attr(
                            reg_cell_attr_mgr, att_c_double_names[attr_dbl_itr] );
                        for( index_t dim_itr = 0;
                            dim_itr < cell_attr_dims[attr_dbl_itr]; ++dim_itr ) {
                            out << " "
                                << cur_attr[c_in_reg[0]
                                    * cell_attr_dims[attr_dbl_itr] + dim_itr];
                        }
                    } else {
                        for( index_t dim_itr = 0;
                            dim_itr < cell_attr_dims[attr_dbl_itr]; ++dim_itr ) {
                            out << " "
                                << GEO::String::to_string( gocad_no_data_value_ );
                        }
                    }
                }
                out << std::endl;
                out << "# CTETRA " << region.name();
                for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                    out << " ";
                    index_t polygon = NO_ID;
                    bool side;
                    if( mesh.cells.is_cell_facet_on_surface( c, f, polygon,
                        side ) ) {
                        index_t surface_id = mesh.polygons.surface( polygon );
                        side ? out << "+" : out << "-";
                        out
                            << region.geomodel().surface( surface_id ).parent( 0 ).name();
                    } else {
                        out << "none";
                    }
                }
                out << std::endl;
            }
        }
        void export_model(
            const GeoModel& geomodel,
            std::ofstream& out,
            std::vector< index_t >& vertex_exported_id ) const
        {
            out << "MODEL" << std::endl;
            int tface_count = 1;

            const GeoModelMeshPolygons& polygons = geomodel.mesh.polygons;
            for( index_t i = 0;
                i < geomodel.nb_geological_entities( Interface::type_name_static() );
                i++ ) {
                const RINGMesh::GeoModelGeologicalEntity& interf =
                    geomodel.geological_entity( Interface::type_name_static(), i );
                out << "SURFACE " << interf.name() << std::endl;
                for( index_t s = 0; s < interf.nb_children(); s++ ) {
                    out << "TFACE " << tface_count++ << std::endl;
                    index_t surface_id = interf.child_gmme( s ).index();
                    out << "KEYVERTICES";
                    index_t key_polygon_id = polygons.polygon( surface_id, 0 );
                    for( index_t v = 0; v < polygons.nb_vertices( key_polygon_id );
                        v++ ) {
                        out << " "
                            << vertex_exported_id[polygons.vertex( key_polygon_id,
                                v )];
                    }
                    out << std::endl;
                    for( index_t p = 0; p < polygons.nb_polygons( surface_id );
                        p++ ) {
                        index_t polygon_id = polygons.polygon( surface_id, p );
                        out << "TRGL";
                        for( index_t v = 0; v < polygons.nb_vertices( polygon_id );
                            v++ ) {
                            out << " "
                                << vertex_exported_id[polygons.vertex( polygon_id,
                                    v )];
                        }
                        out << std::endl;
                    }
                }
            }
        }
        void export_model_region(
            const GeoModel& geomodel,
            std::ofstream& out ) const
        {
            for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                const RINGMesh::Region& region = geomodel.region( r );
                out << "MODEL_REGION " << region.name() << " ";
                region.side( 0 ) ? out << "+" : out << "-";
                out << region.boundary_gmme( 0 ).index() + 1 << std::endl;
            }
        }
    private:
        const double gocad_no_data_value_;
    };

}
