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

            // Print Model3d headers
            out << "GOCAD TSolid 1" << std::endl << "HEADER {" << std::endl
                << "name:" << geomodel.name() << std::endl << "}" << std::endl;

            out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl << "NAME Default"
                << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
                << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl
                << "ZPOSITIVE Elevation" << std::endl
                << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl;

            const GeoModelMesh& mesh = geomodel.mesh;
            mesh.transfert_attributes_from_gmm_to_gm() ;
            //mesh.set_duplicate_mode( GeoModelMeshCells::ALL ) ;

            std::vector< std::string > att_v_double_names;
            std::vector< index_t > vertex_attr_dims;
            fill_vertex_attribute_header( mesh, out, att_v_double_names,
                vertex_attr_dims );
            std::vector< std::string > att_c_double_names;
            std::vector< index_t > cell_attr_dims;
            fill_cell_attribute_header( mesh, out, att_c_double_names,
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
                out << "TVOLUME " << region.name() << std::endl;

                // Export not duplicated vertices
                for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
                    index_t cell = mesh.cells.cell( r, c );
                    for( index_t v = 0; v < mesh.cells.nb_vertices( cell ); v++ ) {
                        index_t atom_id;
                        if( !mesh.cells.is_corner_duplicated( cell, v, atom_id ) ) {
                            index_t vertex_id = mesh.cells.vertex( cell, v );
                            if( vertex_exported[vertex_id] ) continue;
                            vertex_exported[vertex_id] = true;
                            vertex_exported_id[vertex_id] = nb_vertices_exported;
                            // PVRTX must be used instead of VRTX because
                            // properties are not read by Gocad if it is VRTX.
                            out << "PVRTX " << nb_vertices_exported++ << " "
                                << mesh.vertices.vertex( vertex_id );
                            for( index_t attr_dbl_itr = 0;
                                attr_dbl_itr < att_v_double_names.size();
                                ++attr_dbl_itr ) {
                                GEO::Attribute< double > cur_attr(
                                    mesh.vertices.attribute_manager(),
                                    att_v_double_names[attr_dbl_itr] ) ;
                                for( index_t dim_itr = 0;
                                    dim_itr < vertex_attr_dims[attr_dbl_itr];
                                    ++dim_itr ) {
                                    out << " "
                                        << cur_attr[vertex_id
                                            * vertex_attr_dims[attr_dbl_itr]
                                            + dim_itr] ;
                                }
                            }
                            out << std::endl ;
                        }
                    }
                }

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

                /*GEO::Attribute< index_t > attribute( mesh.facet_attribute_manager(),
                 surface_att_name ) ;*/
                for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
                    out << "TETRA";
                    index_t cell = mesh.cells.cell( r, c );
                    for( index_t v = 0; v < mesh.cells.nb_vertices( cell ); v++ ) {
                        index_t atom_id;
                        if( !mesh.cells.is_corner_duplicated( cell, v, atom_id ) ) {
                            index_t vertex_id = mesh.cells.vertex( cell, v );
                            out << " " << vertex_exported_id[vertex_id];
                        } else {
                            out << " " << atom_exported_id[atom_id];
                        }
                    }
                    for( index_t attr_dbl_itr = 0;
                        attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                        GEO::Attribute< double > cur_attr(
                            mesh.cells.attribute_manager(),
                            att_c_double_names[attr_dbl_itr] ) ;
                        for( index_t dim_itr = 0;
                            dim_itr < cell_attr_dims[attr_dbl_itr]; ++dim_itr ) {
                            out << " "
                                << cur_attr[cell * cell_attr_dims[attr_dbl_itr]
                                    + dim_itr] ;
                        }
                    }
                    out << std::endl;
                    out << "# CTETRA " << region.name();
                    for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                        out << " ";
                        index_t facet = NO_ID;
                        bool side;
                        if( mesh.cells.is_cell_facet_on_surface( c, f, facet,
                            side ) ) {
                            index_t surface_id = mesh.facets.surface( facet );
                            side ? out << "+" : out << "-";
                            out << geomodel.surface( surface_id ).parent( 0 ).name();
                        } else {
                            out << "none";
                        }
                    }
                    out << std::endl;
                }
            }

            out << "MODEL" << std::endl;
            int tface_count = 1;
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
                    index_t key_facet_id = mesh.facets.facet( surface_id, 0 );
                    for( index_t v = 0; v < mesh.facets.nb_vertices( key_facet_id );
                        v++ ) {
                        out << " "
                            << vertex_exported_id[mesh.facets.vertex( key_facet_id,
                                v )];
                    }
                    out << std::endl;
                    for( index_t f = 0; f < mesh.facets.nb_facets( surface_id );
                        f++ ) {
                        index_t facet_id = mesh.facets.facet( surface_id, f );
                        out << "TRGL";
                        for( index_t v = 0; v < mesh.facets.nb_vertices( facet_id );
                            v++ ) {
                            out << " "
                                << vertex_exported_id[mesh.facets.vertex( facet_id,
                                    v )];
                        }
                        out << std::endl;
                    }
                }
            }

            for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                const RINGMesh::Region& region = geomodel.region( r );
                out << "MODEL_REGION " << region.name() << " ";
                region.side( 0 ) ? out << "+" : out << "-";
                out << region.boundary_gmme( 0 ).index() + 1 << std::endl;
            }

            out << "END" << std::endl;
        }
    private:
        void fill_vertex_attribute_header(
            const GeoModelMesh& mesh,
            std::ofstream& out,
            std::vector< std::string >& att_v_double_names,
            std::vector< index_t >& vertex_attr_dims ) const
        {
            GEO::vector< std::string > att_v_names ;
            mesh.vertices.attribute_manager().list_attribute_names( att_v_names ) ;
            for( index_t att_v = 0; att_v < mesh.vertices.attribute_manager().nb();
                att_v++ ) {

                if( att_v_names[att_v] == "point" ) {
                    continue ;
                }

                if( !GEO::Attribute< double >::is_defined(
                    mesh.vertices.attribute_manager(), att_v_names[att_v] ) ) {
                    continue ;
                }
                att_v_double_names.push_back( att_v_names[att_v] ) ;
                index_t cur_dim =
                    mesh.vertices.attribute_manager().find_attribute_store(
                        att_v_names[att_v] )->dimension() ;
                vertex_attr_dims.push_back( cur_dim ) ;
            }

            if( !att_v_double_names.empty() ) {
                out << "PROPERTIES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " " << att_v_double_names[attr_dbl_itr] ;
                }
                out << std::endl ;
                out << "PROP_LEGAL_RANGES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " **none**  **none**" ;
                }
                out << std::endl ;
                out << "NO_DATA_VALUES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " -99999" ;
                }
                out << std::endl ;
                out << "READ_ONLY" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " 1" ;
                }
                out << std::endl ;
                out << "PROPERTY_CLASSES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " " << att_v_double_names[attr_dbl_itr] ;
                }
                out << std::endl ;
                out << "PROPERTY_KINDS" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " \"Real Number\"" ;
                }
                out << std::endl ;
                out << "PROPERTY_SUBCLASSES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " QUANTITY Float" ;
                }
                out << std::endl ;
                out << "ESIZES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " "
                        << GEO::String::to_string( vertex_attr_dims[attr_dbl_itr] ) ;
                }
                out << std::endl ;
                out << "UNITS" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " unitless" ;
                }
                out << std::endl ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << "PROPERTY_CLASS_HEADER "
                        << att_v_double_names[attr_dbl_itr] << " {" << std::endl ;
                    out << "kind: Real Number" << std::endl ;
                    out << "unit: unitless" << std::endl ;
                    out << "}" << std::endl ;
                }
            }
        }

        void fill_cell_attribute_header(
            const GeoModelMesh& mesh,
            std::ofstream& out,
            std::vector< std::string >& att_c_double_names,
            std::vector< index_t >& cell_attr_dims ) const
        {
            GEO::vector< std::string > att_c_names ;
            mesh.cells.attribute_manager().list_attribute_names( att_c_names ) ;
            for( index_t att_c = 0; att_c < mesh.cells.attribute_manager().nb();
                att_c++ ) {

                if( !GEO::Attribute< double >::is_defined(
                    mesh.cells.attribute_manager(), att_c_names[att_c] ) ) {
                    continue ;
                }
                att_c_double_names.push_back( att_c_names[att_c] ) ;
                index_t cur_dim = mesh.cells.attribute_manager().find_attribute_store(
                    att_c_names[att_c] )->dimension() ;
                cell_attr_dims.push_back( cur_dim ) ;
            }

            if( !att_c_double_names.empty() ) {
                out << "TETRA_PROPERTIES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << " " << att_c_double_names[attr_dbl_itr] ;
                }
                out << std::endl ;
                out << "TETRA_PROP_LEGAL_RANGES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << " **none**  **none**" ;
                }
                out << std::endl ;
                out << "TETRA_NO_DATA_VALUES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << " -99999" ;
                }
                out << std::endl ;
                out << "READ_ONLY" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << " 1" ;
                }
                out << std::endl ;
                out << "TETRA_PROPERTY_CLASSES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << " " << att_c_double_names[attr_dbl_itr] ;
                }
                out << std::endl ;
                out << "TETRA_PROPERTY_KINDS" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << " \"Real Number\"" ;
                }
                out << std::endl ;
                out << "TETRA_PROPERTY_SUBCLASSES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << " QUANTITY Float" ;
                }
                out << std::endl ;
                out << "TETRA_ESIZES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << " "
                        << GEO::String::to_string( cell_attr_dims[attr_dbl_itr] ) ;
                }
                out << std::endl ;
                out << "TETRA_UNITS" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << " unitless" ;
                }
                out << std::endl ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_c_double_names.size(); ++attr_dbl_itr ) {
                    out << "TETRA_PROPERTY_CLASS_HEADER "
                        << att_c_double_names[attr_dbl_itr] << " {" << std::endl ;
                    out << "kind: Real Number" << std::endl ;
                    out << "unit: unitless" << std::endl ;
                    out << "}" << std::endl ;
                }
            }
        }
    };

}
