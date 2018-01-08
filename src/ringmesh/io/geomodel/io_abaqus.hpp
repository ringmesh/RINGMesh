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

namespace
{
    struct RINGMesh2Abaqus
    {
        std::string entity_type;
        index_t vertices[8];
    };

    static RINGMesh2Abaqus tet_descriptor_abaqus = {
        "C3D4", // type
        { 0, 1, 2, 3 } // vertices
    };

    static RINGMesh2Abaqus hex_descriptor_abaqus = {
        "C3D8", // type
        { 0, 4, 5, 1, 2, 6, 7, 3 } // vertices
    };

    class AbaqusIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        static const index_t NB_ENTRY_PER_LINE = 16;

        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

            out << "*HEADING" << EOL;
            out << "**Mesh exported from RINGMesh" << EOL;
            out << "**https://bitbucket.org/ring_team/ringmesh" << EOL;

            out << "*PART, name=Part-1" << EOL;

            save_vertices( geomodel, out );
            save_nb_polygons( geomodel, out );
            save_cells( geomodel, out );

            out << "*END PART" << EOL;
            out << std::flush;
        }

    private:
        void save_vertices(
            const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            const GeoModelMeshVertices3D& vertices = geomodel.mesh.vertices;
            out << "*NODE" << EOL;
            for( auto v : range( vertices.nb() ) )
            {
                out << v + 1;
                const vec3& vertex = vertices.vertex( v );
                for( auto i : range( 3 ) )
                {
                    out << COMMA << SPACE << vertex[i];
                }
                out << EOL;
            }
        }
        void save_nb_polygons(
            const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            const GeologicalEntityType& type = Interface3D::type_name_static();
            index_t nb_interfaces = geomodel.nb_geological_entities( type );
            for( auto i : range( nb_interfaces ) )
            {
                save_interface( geomodel, i, out );
            }
        }
        void save_interface( const GeoModel3D& geomodel,
            index_t interface_id,
            std::ofstream& out ) const
        {
            const GeoModelMeshPolygons3D& polygons = geomodel.mesh.polygons;
            const GeoModelGeologicalEntity3D& entity =
                geomodel.geological_entity(
                    Interface3D::type_name_static(), interface_id );
            std::string sep;
            index_t count = 0;
            std::vector< bool > vertex_exported(
                geomodel.mesh.vertices.nb(), false );
            out << "*NSET, nset=" << entity.name() << EOL;
            for( auto s : range( entity.nb_children() ) )
            {
                index_t surface_id = entity.child_gmme( s ).index();
                for( auto p : range( polygons.nb_polygons( surface_id ) ) )
                {
                    index_t polygon_id = polygons.polygon( surface_id, p );
                    for( auto v : range( polygons.nb_vertices( polygon_id ) ) )
                    {
                        index_t vertex_id = polygons.vertex(
                            ElementLocalVertex( polygon_id, v ) );
                        if( vertex_exported[vertex_id] )
                            continue;
                        vertex_exported[vertex_id] = true;
                        out << sep << vertex_id + 1;
                        sep = COMMA + SPACE;
                        new_line_if_needed( count, out, sep );
                    }
                }
            }
            out << EOL;
        }

        void save_tets( const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            const GeoModelMeshCells3D& cells = geomodel.mesh.cells;
            if( cells.nb_tet() > 0 )
            {
                out << "*ELEMENT, type=" << tet_descriptor_abaqus.entity_type
                    << EOL;
                for( auto r : range( geomodel.nb_regions() ) )
                {
                    for( auto c : range( cells.nb_tet( r ) ) )
                    {
                        index_t tetra = cells.tet( r, c );
                        out << tetra + 1;
                        for( auto v : range( 4 ) )
                        {
                            index_t vertex_id =
                                tet_descriptor_abaqus.vertices[v];
                            out << COMMA << SPACE
                                << cells.vertex(
                                       ElementLocalVertex( tetra, vertex_id ) )
                                       + 1;
                        }
                        out << EOL;
                    }
                }
            }
        }
        void save_hex( const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            const GeoModelMeshCells3D& cells = geomodel.mesh.cells;
            if( cells.nb_hex() > 0 )
            {
                out << "*ELEMENT, type=" << hex_descriptor_abaqus.entity_type
                    << EOL;
                for( auto r : range( geomodel.nb_regions() ) )
                {
                    for( auto c : range( cells.nb_hex( r ) ) )
                    {
                        index_t hex = cells.hex( r, c );
                        out << hex + 1;
                        for( auto v : range( 8 ) )
                        {
                            index_t vertex_id =
                                hex_descriptor_abaqus.vertices[v];
                            out << COMMA << SPACE
                                << cells.vertex(
                                       ElementLocalVertex( hex, vertex_id ) )
                                       + 1;
                        }
                        out << EOL;
                    }
                }
            }
        }
        void save_regions(
            const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            const GeoModelMeshCells3D& cells = geomodel.mesh.cells;
            for( auto r : range( geomodel.nb_regions() ) )
            {
                const std::string& name = geomodel.region( r ).name();
                out << "*ELSET, elset=" << name << EOL;
                index_t count = 0;
                std::string sep;
                for( auto c : range( cells.nb_tet( r ) ) )
                {
                    index_t tetra = cells.tet( r, c );
                    out << sep << tetra + 1;
                    sep = COMMA + SPACE;
                    new_line_if_needed( count, out, sep );
                }
                for( auto c : range( cells.nb_hex( r ) ) )
                {
                    index_t hex = cells.hex( r, c );
                    out << sep << hex + 1;
                    sep = COMMA + SPACE;
                    new_line_if_needed( count, out, sep );
                }
                reset_line( count, out );

                out << "*NSET, nset=" << name << ", elset=" << name << EOL;
            }
        }
        void save_cells( const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            save_tets( geomodel, out );
            save_hex( geomodel, out );
            save_regions( geomodel, out );
        }
        void new_line_if_needed(
            index_t& count, std::ofstream& out, std::string& sep ) const
        {
            count++;
            if( count == NB_ENTRY_PER_LINE )
            {
                count = 0;
                sep = "";
                out << EOL;
            }
        }
        void reset_line( index_t& count, std::ofstream& out ) const
        {
            if( count != 0 )
            {
                count = 0;
                out << EOL;
            }
        }
    };
}
