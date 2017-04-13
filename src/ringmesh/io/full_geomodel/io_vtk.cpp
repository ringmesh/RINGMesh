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
    struct RINGMesh2VTK {
        index_t entity_type;
        index_t vertices[8];
    };

    static RINGMesh2VTK tet_descriptor_vtk = { 10,                  // type
        { 0, 1, 2, 3 }     // vertices
    };

    static RINGMesh2VTK hex_descriptor_vtk = { 12,                         // type
        { 0, 4, 5, 1, 2, 6, 7, 3 }     // vertices
    };

    static RINGMesh2VTK prism_descriptor_vtk = { 13,                     // type
        { 0, 2, 1, 3, 5, 4 }   // vertices
    };

    static RINGMesh2VTK pyramid_descriptor_vtk = { 14,                 // type
        { 0, 1, 2, 3, 4 }  // vertices
    };

    static RINGMesh2VTK* cell_type_to_cell_descriptor_vtk[4] = {
        &tet_descriptor_vtk, &hex_descriptor_vtk, &prism_descriptor_vtk,
        &pyramid_descriptor_vtk };

    class VTKIOHandler final: public GeoModelIOHandler {
    public:
        virtual bool load( const std::string& filename, GeoModel& geomodel ) override
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from VTK not implemented yet" );
            return false;
        }
        virtual void save( const GeoModel& geomodel, const std::string& filename ) override
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

            out << "# vtk DataFile Version 2.0" << std::endl;
            out << "Unstructured Grid" << std::endl;
            out << "ASCII" << std::endl;
            out << "DATASET UNSTRUCTURED_GRID" << std::endl;

            const GeoModelMesh& mesh = geomodel.mesh;
            out << "POINTS " << mesh.vertices.nb() << " double" << std::endl;
            for( index_t v = 0; v < mesh.vertices.nb(); v++ ) {
                out << mesh.vertices.vertex( v ) << std::endl;
            }
            out << std::endl;

            index_t total_corners = ( 4 + 1 ) * mesh.cells.nb_tet()
                + ( 5 + 1 ) * mesh.cells.nb_pyramid()
                + ( 6 + 1 ) * mesh.cells.nb_prism()
                + ( 8 + 1 ) * mesh.cells.nb_hex();
            out << "CELLS " << mesh.cells.nb_cells() << SPACE << total_corners
                << std::endl;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                out << mesh.cells.nb_vertices( c );
                const RINGMesh2VTK& descriptor =
                    *cell_type_to_cell_descriptor_vtk[mesh.cells.type( c )];
                for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                    index_t vertex_id = descriptor.vertices[v];
                    out << SPACE << mesh.cells.vertex( c, vertex_id );
                }
                out << std::endl;
            }

            out << "CELL_TYPES " << mesh.cells.nb() << std::endl;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                const RINGMesh2VTK& descriptor =
                    *cell_type_to_cell_descriptor_vtk[mesh.cells.type( c )];
                out << descriptor.entity_type << std::endl;
            }
            out << std::endl;

            out << "CELL_DATA " << mesh.cells.nb() << std::endl;
            out << "SCALARS region int 1" << std::endl;
            out << "LOOKUP_TABLE default" << std::endl;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                out << mesh.cells.region( c ) << std::endl;
            }
            out << std::endl;
        }
    };

}
