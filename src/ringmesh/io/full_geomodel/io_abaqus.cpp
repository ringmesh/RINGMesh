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

    struct RINGMesh2Abaqus {
        std::string entity_type;
        index_t vertices[8];
    };

    static RINGMesh2Abaqus tet_descriptor_abaqus = { "C3D4",                  // type
        { 0, 1, 2, 3 }     // vertices
    };

    static RINGMesh2Abaqus hex_descriptor_abaqus = { "C3D8",                  // type
        { 0, 4, 5, 1, 2, 6, 7, 3 }     // vertices
    };

    class AbaqusIOHandler final: public GeoModelIOHandler {
    public:
        static const index_t NB_ENTRY_PER_LINE = 16;

        virtual bool load( const std::string& filename, GeoModel< 3 >& geomodel ) final
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from abaqus not implemented yet" );
            return false;
        }
        virtual void save(
            const GeoModel< 3 >& geomodel,
            const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

            out << "*HEADING" << std::endl;
            out << "**Mesh exported from RINGMesh" << std::endl;
            out << "**https://bitbucket.org/ring_team/ringmesh" << std::endl;

            out << "*PART, name=Part-1" << std::endl;

            save_vertices( geomodel, out );
            save_nb_polygons( geomodel, out );
            save_cells( geomodel, out );

            out << "*END PART" << std::endl;
        }
    private:
        void save_vertices( const GeoModel< 3 >& geomodel, std::ofstream& out ) const
        {
            const GeoModelMeshVertices< 3 >& vertices = geomodel.mesh.vertices;
            out << "*NODE" << std::endl;
            for( index_t v = 0; v < vertices.nb(); v++ ) {
                out << v + 1;
                const vec3& vertex = vertices.vertex( v );
                for( index_t i = 0; i < 3; i++ ) {
                    out << COMMA << SPACE << vertex[i];
                }
                out << std::endl;
            }

        }
        void save_nb_polygons(
            const GeoModel< 3 >& geomodel,
            std::ofstream& out ) const
        {
            const GeologicalEntityType& type = Interface < 3 > ::type_name_static();
            index_t nb_interfaces = geomodel.nb_geological_entities( type );
            for( index_t i = 0; i < nb_interfaces; i++ ) {
                save_interface( geomodel, i, out );
            }
        }
        void save_interface(
            const GeoModel< 3 >& geomodel,
            index_t interface_id,
            std::ofstream& out ) const
        {
            const GeoModelMeshPolygons< 3 >& polygons = geomodel.mesh.polygons;
            const GeoModelGeologicalEntity< 3 >& entity = geomodel.geological_entity(
                Interface < 3 > ::type_name_static(), interface_id );
            std::string sep;
            index_t count = 0;
            std::vector< bool > vertex_exported( geomodel.mesh.vertices.nb(),
                false );
            out << "*NSET, nset=" << entity.name() << std::endl;
            for( index_t s = 0; s < entity.nb_children(); s++ ) {
                index_t surface_id = entity.child_gmme( s ).index();
                for( index_t p = 0; p < polygons.nb_polygons( surface_id ); p++ ) {
                    index_t polygon_id = polygons.polygon( surface_id, p );
                    for( index_t v = 0; v < polygons.nb_vertices( polygon_id );
                        v++ ) {
                        index_t vertex_id = polygons.vertex( polygon_id, v );
                        if( vertex_exported[vertex_id] ) continue;
                        vertex_exported[vertex_id] = true;
                        out << sep << vertex_id + 1;
                        sep = COMMA + SPACE;
                        new_line_if_needed( count, out, sep );
                    }
                }
            }
            out << std::endl;
        }

        void save_tets( const GeoModel< 3 >& geomodel, std::ofstream& out ) const
        {
            const GeoModelMeshCells< 3 >& cells = geomodel.mesh.cells;
            if( cells.nb_tet() > 0 ) {
                out << "*ELEMENT, type=" << tet_descriptor_abaqus.entity_type
                    << std::endl;
                for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                    for( index_t c = 0; c < cells.nb_tet( r ); c++ ) {
                        index_t tetra = cells.tet( r, c );
                        out << tetra + 1;
                        for( index_t v = 0; v < 4; v++ ) {
                            index_t vertex_id = tet_descriptor_abaqus.vertices[v];
                            out << COMMA << SPACE
                                << cells.vertex( tetra, vertex_id ) + 1;
                        }
                        out << std::endl;
                    }
                }
            }
        }
        void save_hex( const GeoModel< 3 >& geomodel, std::ofstream& out ) const
        {
            const GeoModelMeshCells< 3 >& cells = geomodel.mesh.cells;
            if( cells.nb_hex() > 0 ) {
                out << "*ELEMENT, type=" << hex_descriptor_abaqus.entity_type
                    << std::endl;
                for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                    for( index_t c = 0; c < cells.nb_hex( r ); c++ ) {
                        index_t hex = cells.hex( r, c );
                        out << hex + 1;
                        for( index_t v = 0; v < 8; v++ ) {
                            index_t vertex_id = hex_descriptor_abaqus.vertices[v];
                            out << COMMA << SPACE
                                << cells.vertex( hex, vertex_id ) + 1;
                        }
                        out << std::endl;
                    }
                }
            }
        }
        void save_regions( const GeoModel< 3 >& geomodel, std::ofstream& out ) const
        {
            const GeoModelMeshCells< 3 >& cells = geomodel.mesh.cells;
            for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                const std::string& name = geomodel.region( r ).name();
                out << "*ELSET, elset=" << name << std::endl;
                index_t count = 0;
                std::string sep;
                for( index_t c = 0; c < cells.nb_tet( r ); c++ ) {
                    index_t tetra = cells.tet( r, c );
                    out << sep << tetra + 1;
                    sep = COMMA + SPACE;
                    new_line_if_needed( count, out, sep );

                }
                for( index_t c = 0; c < cells.nb_hex( r ); c++ ) {
                    index_t hex = cells.hex( r, c );
                    out << sep << hex + 1;
                    sep = COMMA + SPACE;
                    new_line_if_needed( count, out, sep );
                }
                reset_line( count, out );

                out << "*NSET, nset=" << name << ", elset=" << name << std::endl;
            }
        }
        void save_cells( const GeoModel< 3 >& geomodel, std::ofstream& out ) const
        {
            save_tets( geomodel, out );
            save_hex( geomodel, out );
            save_regions( geomodel, out );
        }
        void new_line_if_needed(
            index_t& count,
            std::ofstream& out,
            std::string& sep ) const
        {
            count++;
            if( count == NB_ENTRY_PER_LINE ) {
                count = 0;
                sep = "";
                out << std::endl;
            }
        }
        void reset_line( index_t& count, std::ofstream& out ) const
        {
            if( count != 0 ) {
                count = 0;
                out << std::endl;
            }
        }
    };

}
