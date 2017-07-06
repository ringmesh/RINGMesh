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
    class TSolidIOHandler final: public GeoModelIOHandler {
    public:
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
            const GeoModelMeshPolygons& polygons = geomodel.mesh.polygons;
            //mesh.set_duplicate_mode( GeoModelMeshCells::ALL ) ;

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
                        index_t atom_id = mesh.cells.duplicated_corner_index( cell, v );
                        if( atom_id == NO_ID ) {
                            index_t vertex_id = mesh.cells.vertex( cell, v );
                            if( vertex_exported[vertex_id] ) continue;
                            vertex_exported[vertex_id] = true;
                            vertex_exported_id[vertex_id] = nb_vertices_exported;
                            // PVRTX must be used instead of VRTX because
                            // properties are not read by Gocad if it is VRTX.
                            out << "PVRTX " << nb_vertices_exported++ << " "
                                << mesh.vertices.vertex( vertex_id ) << std::endl;
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

                for( index_t c = 0; c < region.nb_mesh_elements(); c++ ) {
                    out << "TETRA";
                    index_t cell = mesh.cells.cell( r, c );
                    for( index_t v = 0; v < mesh.cells.nb_vertices( cell ); v++ ) {
                        index_t atom_id = mesh.cells.duplicated_corner_index( cell, v );
                        if( atom_id == NO_ID ) {
                            index_t vertex_id = mesh.cells.vertex( cell, v );
                            out << " " << vertex_exported_id[vertex_id];
                        } else {
                            out << " " << atom_exported_id[atom_id];
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
                            index_t surface_id = polygons.surface( polygon );
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

            for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                const RINGMesh::Region& region = geomodel.region( r );
                out << "MODEL_REGION " << region.name() << " ";
                region.side( 0 ) ? out << "+" : out << "-";
                out << region.boundary_gmme( 0 ).index() + 1 << std::endl;
            }

            out << "END" << std::endl;
        }
    };

}
