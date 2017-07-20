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
            :
                GeoModelIOHandler(),
                out_(),
                numeric_like_vertex_attribute_names_(),
                vertex_attribute_dimensions_(),
                numeric_like_cell_attribute_names_(),
                cell_attribute_dimensions_(),
                vertex_exported_(),
                atom_exported_(),
                vertex_exported_id_(),
                atom_exported_id_()
        {
        }
        virtual void load( const std::string& filename, GeoModel& geomodel ) final
        {
            std::ifstream input( filename.c_str() );
            if( !input ) {
                throw RINGMeshException( "I/O",
                    "Failed loading geomodel from file " + filename );
            }
            GeoModelBuilderTSolid builder( geomodel, filename );
            builder.build_geomodel();
        }
        virtual void save( const GeoModel& geomodel, const std::string& filename ) final
        {
            ringmesh_assert( !out_.is_open() );
            out_.open( filename.c_str() );
            ringmesh_assert( out_.is_open() );
            out_.precision( 16 );

            fill_top_header( geomodel );

            const GeoModelMesh& mesh = geomodel.mesh;
            //mesh.set_duplicate_mode( GeoModelMeshCells::ALL ) ;

            fill_vertex_attribute_header( geomodel );
            fill_cell_attribute_header( geomodel );

            vertex_exported_.resize( mesh.vertices.nb(), false );
            atom_exported_.resize( mesh.cells.nb_duplicated_vertices(), false );
            vertex_exported_id_.resize( mesh.vertices.nb(), NO_ID );
            atom_exported_id_.resize( mesh.cells.nb_duplicated_vertices(), NO_ID );
            for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                const RINGMesh::Region& region = geomodel.region( r );
                export_one_region( region );
            }

            export_model( geomodel );
            export_model_region( geomodel );

            out_ << "END" << std::endl;
        }
    private:
        const GeoModelMesh& geomodel_mesh( const Region& region ) const
        {
            return region.geomodel().mesh;
        }
        void fill_top_header( const GeoModel& geomodel )
        {
            ringmesh_assert( out_.is_open() );
            // Print Model3d headers
            out_ << "GOCAD TSolid 1" << std::endl << "HEADER {" << std::endl
                << "name:" << geomodel.name() << std::endl << "}" << std::endl;

            out_ << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl << "NAME Default"
                << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
                << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl
                << "ZPOSITIVE Elevation" << std::endl
                << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl;
        }
        void fill_vertex_attribute_header( const GeoModel& geomodel )
        {
            ringmesh_assert( out_.is_open() );
            std::vector< bool > is_integer_like_attribute;
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

                    if( contains( numeric_like_vertex_attribute_names_,
                        cur_att_v_name ) ) {
                        continue;
                    }

                    const GEO::AttributeStore* attr_store =
                        reg_vertex_attr_mgr.find_attribute_store( cur_att_v_name );
                    ringmesh_assert( attr_store != nullptr );

                    if( !GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to(
                        attr_store ) ) {
                        continue;
                    }

                    numeric_like_vertex_attribute_names_.push_back( cur_att_v_name );
                    index_t cur_dim = attr_store->dimension();
                    vertex_attribute_dimensions_.push_back( cur_dim );

                    const GEO::ReadOnlyScalarAttributeAdapter adapter(
                        reg_vertex_attr_mgr, cur_att_v_name );
                    ringmesh_assert(
                        adapter.element_type()
                            != GEO::ReadOnlyScalarAttributeAdapter::ET_NONE );
                    is_integer_like_attribute.push_back(
                        adapter.element_type()
                            < GEO::ReadOnlyScalarAttributeAdapter::ET_FLOAT32 );
                }
            }

            const index_t nb_numeric_like_vertex_attribute_names =
                numeric_like_vertex_attribute_names_.size();
            if( !numeric_like_vertex_attribute_names_.empty() ) {
                ringmesh_assert( nb_numeric_like_vertex_attribute_names > 0 );
                out_ << "PROPERTIES";
                for( const std::string& cur_num_like_v_att_name : numeric_like_vertex_attribute_names_ ) {
                    out_ << " " << cur_num_like_v_att_name;
                }
                out_ << std::endl;
                out_ << "PROP_LEGAL_RANGES";
                for( index_t i = 0; i < nb_numeric_like_vertex_attribute_names;
                    ++i ) {
                    out_ << " **none**  **none**";
                }
                out_ << std::endl;
                out_ << "NO_DATA_VALUES";
                write_no_data_value( nb_numeric_like_vertex_attribute_names );
                out_ << std::endl;
                out_ << "READ_ONLY";
                for( index_t i = 0; i < nb_numeric_like_vertex_attribute_names;
                    ++i ) {
                    out_ << " 1";
                }
                out_ << std::endl;
                out_ << "PROPERTY_CLASSES";
                for( const std::string& cur_num_like_v_att_name : numeric_like_vertex_attribute_names_ ) {
                    out_ << " " << cur_num_like_v_att_name;
                }
                out_ << std::endl;
                out_ << "PROPERTY_KINDS";
                for( index_t i = 0; i < nb_numeric_like_vertex_attribute_names;
                    ++i ) {
                    if( is_integer_like_attribute[i] ) {
                        out_ << " \"Number\"";
                    } else {
                        out_ << " \"Real Number\"";
                    }
                }
                out_ << std::endl;
                out_ << "PROPERTY_SUBCLASSES";
                for( index_t i = 0; i < nb_numeric_like_vertex_attribute_names;
                    ++i ) {
                    out_ << " QUANTITY Float";
                }
                out_ << std::endl;
                out_ << "ESIZES";
                for( const index_t& cur_v_att_dim : vertex_attribute_dimensions_ ) {
                    out_ << " " << GEO::String::to_string( cur_v_att_dim );
                }
                out_ << std::endl;
                out_ << "UNITS";
                for( index_t i = 0; i < nb_numeric_like_vertex_attribute_names;
                    ++i ) {
                    out_ << " unitless";
                }
                out_ << std::endl;
                for( index_t i = 0; i < nb_numeric_like_vertex_attribute_names;
                    ++i ) {
                    out_ << "PROPERTY_CLASS_HEADER "
                        << numeric_like_vertex_attribute_names_[i] << " {"
                        << std::endl;
                    if( is_integer_like_attribute[i] ) {
                        out_ << "kind: Number" << std::endl;
                    } else {
                        out_ << "kind: Real Number" << std::endl;
                    }
                    out_ << "unit: unitless" << std::endl;
                    out_ << "}" << std::endl;
                }
            }
        }

        void fill_cell_attribute_header( const GeoModel& geomodel )
        {
            ringmesh_assert( out_.is_open() );
            std::vector< bool > is_integer_like_attribute;
            for( index_t reg_i = 0; reg_i < geomodel.nb_regions(); ++reg_i ) {
                const Region& cur_reg = geomodel.region( reg_i );
                GEO::AttributesManager& reg_cell_attr_mgr =
                    cur_reg.cell_attribute_manager();
                GEO::vector< std::string > att_c_names;
                reg_cell_attr_mgr.list_attribute_names( att_c_names );
                ringmesh_assert( att_c_names.size() == reg_cell_attr_mgr.nb() );
                for( const std::string& cur_att_c_name : att_c_names ) {

                    if( contains( numeric_like_cell_attribute_names_,
                        cur_att_c_name ) ) {
                        continue;
                    }

                    const GEO::AttributeStore* attr_store =
                        reg_cell_attr_mgr.find_attribute_store( cur_att_c_name );
                    ringmesh_assert( attr_store != nullptr );

                    if( !GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to(
                        attr_store ) ) {
                        continue;
                    }

                    numeric_like_cell_attribute_names_.push_back( cur_att_c_name );
                    index_t cur_dim = attr_store->dimension();
                    cell_attribute_dimensions_.push_back( cur_dim );

                    const GEO::ReadOnlyScalarAttributeAdapter adapter(
                        reg_cell_attr_mgr, cur_att_c_name );
                    ringmesh_assert(
                        adapter.element_type()
                            != GEO::ReadOnlyScalarAttributeAdapter::ET_NONE );
                    is_integer_like_attribute.push_back(
                        adapter.element_type()
                            < GEO::ReadOnlyScalarAttributeAdapter::ET_FLOAT32 );
                }
            }

            const index_t nb_numeric_like_cell_attribute_names =
                numeric_like_cell_attribute_names_.size();
            if( !numeric_like_cell_attribute_names_.empty() ) {
                out_ << "TETRA_PROPERTIES";
                for( const std::string& cur_num_like_c_att_name : numeric_like_cell_attribute_names_ ) {
                    out_ << " " << cur_num_like_c_att_name;
                }
                out_ << std::endl;
                out_ << "TETRA_PROP_LEGAL_RANGES";
                for( index_t i = 0; i < nb_numeric_like_cell_attribute_names; ++i ) {
                    out_ << " **none**  **none**";
                }
                out_ << std::endl;
                out_ << "TETRA_NO_DATA_VALUES";
                write_no_data_value( nb_numeric_like_cell_attribute_names );
                out_ << std::endl;
                out_ << "READ_ONLY";
                for( index_t i = 0; i < nb_numeric_like_cell_attribute_names; ++i ) {
                    out_ << " 1";
                }
                out_ << std::endl;
                out_ << "TETRA_PROPERTY_CLASSES";
                for( const std::string& cur_num_like_c_att_name : numeric_like_cell_attribute_names_ ) {
                    out_ << " " << cur_num_like_c_att_name;
                }
                out_ << std::endl;
                out_ << "TETRA_PROPERTY_KINDS";
                for( index_t i = 0; i < nb_numeric_like_cell_attribute_names; ++i ) {
                    if( is_integer_like_attribute[i] ) {
                        out_ << " \"Number\"";
                    } else {
                        out_ << " \"Real Number\"";
                    }
                }
                out_ << std::endl;
                out_ << "TETRA_PROPERTY_SUBCLASSES";
                for( index_t i = 0; i < nb_numeric_like_cell_attribute_names; ++i ) {
                    out_ << " QUANTITY Float";
                }
                out_ << std::endl;
                out_ << "TETRA_ESIZES";
                for( const index_t& cur_cell_attr_dim : cell_attribute_dimensions_ ) {
                    out_ << " " << GEO::String::to_string( cur_cell_attr_dim );
                }
                out_ << std::endl;
                out_ << "TETRA_UNITS";
                for( index_t i = 0; i < nb_numeric_like_cell_attribute_names; ++i ) {
                    out_ << " unitless";
                }
                out_ << std::endl;
                for( index_t i = 0; i < nb_numeric_like_cell_attribute_names; ++i ) {
                    out_ << "TETRA_PROPERTY_CLASS_HEADER "
                        << numeric_like_cell_attribute_names_[i] << " {"
                        << std::endl;
                    if( is_integer_like_attribute[i] ) {
                        out_ << "kind: Number" << std::endl;
                    } else {
                        out_ << "kind: Real Number" << std::endl;
                    }
                    out_ << "unit: unitless" << std::endl;
                    out_ << "}" << std::endl;
                }
            }
        }

        void export_one_region( const RINGMesh::Region& region )
        {
            ringmesh_assert( out_.is_open() );
            out_ << "TVOLUME " << region.name() << std::endl;
            export_region_vertices( region );
            export_tetrahedra( region );
        }

        void export_region_vertices( const RINGMesh::Region& region )
        {
            // Export not duplicated vertices
            for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
                export_region_cell_vertices( region, c );
            }
        }

        void export_region_cell_vertices( const RINGMesh::Region& region, index_t c )
        {
            const GeoModelMesh& mesh = geomodel_mesh( region );
            index_t cell = mesh.cells.cell( region.gmme().index(), c );
            vec3 cell_center = mesh.cells.barycenter( cell );
            for( index_t v = 0; v < mesh.cells.nb_vertices( cell ); v++ ) {
                export_region_cell_vertex( region, cell, v, cell_center );
            }
        }

        void export_region_cell_vertex(
            const RINGMesh::Region& region,
            index_t cell,
            index_t v,
            const vec3& cell_center )
        {
            ringmesh_assert( out_.is_open() );
            const GeoModelMesh& mesh = geomodel_mesh( region );
            index_t atom_id;
            if( !mesh.cells.is_corner_duplicated( cell, v, atom_id ) ) {
                index_t vertex_id = mesh.cells.vertex( cell, v );
                if( vertex_exported_[vertex_id] ) return;
                vertex_exported_[vertex_id] = true;
                vertex_exported_id_[vertex_id] = nb_vertices_exported_;
                // PVRTX keyword must be used instead of VRTX keyword because
                // properties are not read by Gocad if it is VRTX keyword.
                out_ << "PVRTX " << nb_vertices_exported_++ << " "
                    << mesh.vertices.vertex( vertex_id );

                /// Export of vertex attributes
                export_region_cell_vertex_attributes( region, vertex_id,
                    cell_center );
                out_ << std::endl;
            }
        }

        void export_region_cell_vertex_attributes(
            const RINGMesh::Region& region,
            index_t vertex_id,
            const vec3& cell_center )
        {
            ringmesh_assert( out_.is_open() );
            GEO::AttributesManager& reg_vertex_attr_mgr =
                region.vertex_attribute_manager();
            index_t vertex_id_in_reg = find_gmm_cell_in_gm_region( region, vertex_id,
                cell_center );
            ringmesh_assert( vertex_id_in_reg != NO_ID );

            for( index_t attr_dbl_itr = 0;
                attr_dbl_itr < numeric_like_vertex_attribute_names_.size();
                ++attr_dbl_itr ) {
                if( reg_vertex_attr_mgr.is_defined(
                    numeric_like_vertex_attribute_names_[attr_dbl_itr] ) ) {
                    const GEO::AttributeStore* attr_store =
                        reg_vertex_attr_mgr.find_attribute_store(
                            numeric_like_vertex_attribute_names_[attr_dbl_itr] );
                    ringmesh_assert( attr_store != nullptr );
                    ringmesh_assert(
                        attr_store->dimension()
                            == vertex_attribute_dimensions_[attr_dbl_itr] );
                    ringmesh_assert(
                        GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to(
                            attr_store ) );
                    GEO::ReadOnlyScalarAttributeAdapter cur_attr(
                        reg_vertex_attr_mgr,
                        numeric_like_vertex_attribute_names_[attr_dbl_itr] );
                    for( index_t dim_itr = 0;
                        dim_itr < vertex_attribute_dimensions_[attr_dbl_itr];
                        ++dim_itr ) {
                        out_ << " "
                            << cur_attr[vertex_id_in_reg
                                * vertex_attribute_dimensions_[attr_dbl_itr]
                                + dim_itr];
                    }
                } else {
                    write_no_data_value(
                        vertex_attribute_dimensions_[attr_dbl_itr] );
                }
            }
        }

        index_t find_gmm_cell_in_gm_region(
            const Region& region,
            index_t vertex_id,
            const vec3& cell_center ) const
        {
            const GeoModelMesh& mesh = geomodel_mesh( region );
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
                    region.cells_around_vertex( cur_gme_vertex.v_index, NO_ID );
                /// WARNING: the cell id in the region corresponding
                /// to the cell id in the GMM "cell" is not the
                /// variable "c" (in the for loop over the region cells).
                for( const index_t& cur_cell_around_vertex : cells_around_vertex ) {
                    vec3 center = region.mesh_element_barycenter(
                        cur_cell_around_vertex );
                    if( ( center - cell_center ).length()
                        < region.geomodel().epsilon() ) {
                        return cur_gme_vertex.v_index;
                    }
                }
            }
            ringmesh_assert_not_reached;
            return NO_ID;
        }

        void export_tetrahedra( const RINGMesh::Region& region )
        {
            // Export duplicated vertices
            export_duplicated_vertices();

            // Mark if a boundary is ending in the region
            mark_boundary_ending_in_region( region );

            const GeoModelMesh& mesh = geomodel_mesh( region );
            for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
                out_ << "TETRA";
                index_t cell = mesh.cells.cell( region.gmme().index(), c );
                export_tetra( region, cell );
                out_ << std::endl;
                out_ << "# CTETRA " << region.name();
                export_ctetra( region, c );
                out_ << std::endl;
            }
        }

        void export_duplicated_vertices() const
        {
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
        }

        /// TODO it seems that this code is not used for export since no writtings
        /// in out ofstream.
        void mark_boundary_ending_in_region( const Region& region )
        {
            std::map< index_t, index_t > sides;
            for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                if( sides.count( region.boundary_gmme( s ).index() ) > 0 ) {
                    // a surface is encountered twice, it is ending in the region
                    sides[region.boundary_gmme( s ).index()] = 2;
                } else {
                    sides[region.boundary_gmme( s ).index()] = region.side( s );
                }
            }
        }

        void export_tetra( const Region& region, index_t cell )
        {
            export_tetra_coordinates( region, cell );
            /// Export cell attributes
            export_tetra_attributes( region, cell );
        }

        void export_tetra_coordinates( const Region& region, index_t cell )
        {
            ringmesh_assert( out_.is_open() );
            const GeoModelMesh& mesh = geomodel_mesh( region );
            for( index_t v = 0; v < region.geomodel().mesh.cells.nb_vertices( cell );
                v++ ) {
                index_t atom_id;
                if( !mesh.cells.is_corner_duplicated( cell, v, atom_id ) ) {
                    index_t vertex_id = mesh.cells.vertex( cell, v );
                    out_ << " " << vertex_exported_id_[vertex_id];
                } else {
                    out_ << " " << atom_exported_id_[atom_id];
                }
            }
        }

        void export_tetra_attributes( const Region& region, index_t cell )
        {
            ringmesh_assert( out_.is_open() );
            GEO::AttributesManager& reg_cell_attr_mgr =
                region.cell_attribute_manager();
            const GeoModelMesh& mesh = geomodel_mesh( region );
            vec3 center = mesh.cells.barycenter( cell );
            const std::vector< index_t > c_in_reg =
                region.cell_nn_search().get_neighbors( center,
                    region.geomodel().epsilon() );
            ringmesh_assert( c_in_reg.size() == 1 );
            for( index_t attr_dbl_itr = 0;
                attr_dbl_itr < numeric_like_cell_attribute_names_.size();
                ++attr_dbl_itr ) {
                if( reg_cell_attr_mgr.is_defined(
                    numeric_like_cell_attribute_names_[attr_dbl_itr] ) ) {
                    const GEO::AttributeStore* attr_store =
                        reg_cell_attr_mgr.find_attribute_store(
                            numeric_like_cell_attribute_names_[attr_dbl_itr] );
                    ringmesh_assert( attr_store != nullptr );
                    ringmesh_assert(
                        attr_store->dimension()
                            == cell_attribute_dimensions_[attr_dbl_itr] );
                    ringmesh_assert(
                        GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to(
                            attr_store ) );
                    GEO::ReadOnlyScalarAttributeAdapter cur_attr( reg_cell_attr_mgr,
                        numeric_like_cell_attribute_names_[attr_dbl_itr] );
                    for( index_t dim_itr = 0;
                        dim_itr < cell_attribute_dimensions_[attr_dbl_itr];
                        ++dim_itr ) {
                        out_ << " "
                            << cur_attr[c_in_reg[0]
                                * cell_attribute_dimensions_[attr_dbl_itr] + dim_itr];
                    }
                } else {
                    write_no_data_value( cell_attribute_dimensions_[attr_dbl_itr] );
                }
            }
        }

        void export_ctetra( const Region& region, index_t c )
        {
            ringmesh_assert( out_.is_open() );
            const GeoModelMesh& mesh = geomodel_mesh( region );
            for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                out_ << " ";
                index_t polygon = NO_ID;
                bool side;
                if( mesh.cells.is_cell_facet_on_surface( c, f, polygon, side ) ) {
                    index_t surface_id = mesh.polygons.surface( polygon );
                    side ? out_ << "+" : out_ << "-";
                    out_
                        << region.geomodel().surface( surface_id ).parent( 0 ).name();
                } else {
                    out_ << "none";
                }
            }
        }

        void export_model( const GeoModel& geomodel )
        {
            ringmesh_assert( out_.is_open() );
            out_ << "MODEL" << std::endl;
            int tface_count = 1;

            const GeoModelMeshPolygons& polygons = geomodel.mesh.polygons;
            for( index_t i = 0;
                i < geomodel.nb_geological_entities( Interface::type_name_static() );
                i++ ) {
                const RINGMesh::GeoModelGeologicalEntity& interf =
                    geomodel.geological_entity( Interface::type_name_static(), i );
                out_ << "SURFACE " << interf.name() << std::endl;
                for( index_t s = 0; s < interf.nb_children(); s++ ) {
                    out_ << "TFACE " << tface_count++ << std::endl;
                    index_t surface_id = interf.child_gmme( s ).index();
                    out_ << "KEYVERTICES";
                    index_t key_polygon_id = polygons.polygon( surface_id, 0 );
                    for( index_t v = 0; v < polygons.nb_vertices( key_polygon_id );
                        v++ ) {
                        out_ << " "
                            << vertex_exported_id_[polygons.vertex( key_polygon_id,
                                v )];
                    }
                    out_ << std::endl;
                    for( index_t p = 0; p < polygons.nb_polygons( surface_id );
                        p++ ) {
                        index_t polygon_id = polygons.polygon( surface_id, p );
                        out_ << "TRGL";
                        for( index_t v = 0; v < polygons.nb_vertices( polygon_id );
                            v++ ) {
                            out_ << " "
                                << vertex_exported_id_[polygons.vertex( polygon_id,
                                    v )];
                        }
                        out_ << std::endl;
                    }
                }
            }
        }
        void export_model_region( const GeoModel& geomodel )
        {
            ringmesh_assert( out_.is_open() );
            for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                const RINGMesh::Region& region = geomodel.region( r );
                out_ << "MODEL_REGION " << region.name() << " ";
                region.side( 0 ) ? out_ << "+" : out_ << "-";
                out_ << region.boundary_gmme( 0 ).index() + 1 << std::endl;
            }
        }

        void write_no_data_value( index_t nb )
        {
            for( index_t i = 0; i < nb; ++i ) {
                out_ << " " << GEO::String::to_string( gocad_no_data_value_ );
            }
        }

    private:
        /// Attributes for save
        std::ofstream out_;
        /// numeric_like_vertex_attribute_names_ and
        /// numeric_like_cell_attribute_names_ contain all the
        /// attribute names found the regions which can be used with
        /// GEO::ReadOnlyScalarAttributeAdapter.
        /// Only vertex and cell attributes.
        /// If an attribute is defined in a region and not in another, in
        /// the region where the attribute is not defined no data values
        /// are stored.
        std::vector< std::string > numeric_like_vertex_attribute_names_;
        std::vector< index_t > vertex_attribute_dimensions_;
        std::vector< std::string > numeric_like_cell_attribute_names_;
        std::vector< index_t > cell_attribute_dimensions_;
        std::vector< bool > vertex_exported_;
        std::vector< bool > atom_exported_;
        std::vector< index_t > vertex_exported_id_;
        std::vector< index_t > atom_exported_id_;
        index_t nb_vertices_exported_ { 1 };

        /// Classical Gocad NoDataValue
        const double gocad_no_data_value_ { -99999 };
    };

}
