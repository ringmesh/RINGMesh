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
    /// Convert the cell type of RINGMesh to the MFEM one
    /// NO_ID for pyramids and prims because there are not supported by MFEM
    static const index_t cell_type_mfem[4] = { 4, 5, NO_ID, NO_ID };

    /// Convert the polygon type of RINGMesh to the MFEM one
    /// NO_ID for polygons there are not supported by MFEM
    static const index_t polygon_type_mfem[3] = { 2, 3, NO_ID };

    /// Convert the numerotation from RINGMesh to MFEM
    /// It works for Hexaedron and also for Tetrahedron (in this
    /// case, only the first four values of this table
    /// are used while iterating on vertices)
    static const index_t cell2mfem[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };

    /// MFEM works with Surface and Region index begin with 1
    static const index_t mfem_offset = 1;

    /*!
     * Export for the MFEM format http://mfem.org/
     * Mesh file description : http://mfem.org/mesh-formats/#mfem-mesh-v10
     * "MFEM is a free, lightweight, scalable C++ library for finite element
     * methods"
     */
    class MFEMIOHandler final: public GeoModelIOHandler {
    public:
        virtual bool load( const std::string& filename, GeoModel< 3 >& geomodel ) final
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from MFEM not implemented yet" );
            return false;
        }
        virtual void save( const GeoModel< 3 >& geomodel, const std::string& filename ) final
        {
            const GeoModelMesh& geomodel_mesh = geomodel.mesh;
            index_t nb_cells = geomodel_mesh.cells.nb();
            if( geomodel_mesh.cells.nb_tet() != nb_cells
                && geomodel_mesh.cells.nb_hex() != nb_cells ) {
                throw RINGMeshException( "I/O",
                    "Export to MFEM format works only with full tet or full hex format" );
            }
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

            write_header( geomodel_mesh, out );
            write_cells( geomodel_mesh, out );
            write_polygons( geomodel_mesh, out );
            write_vertices( geomodel_mesh, out );
        }

    private:
        /*!
         * @brief Write the header for the MFEM mesh file
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_header( const GeoModelMesh& geomodel_mesh, std::ofstream& out ) const
        {
            // MFEM mesh version
            out << "MFEM mesh v1.0" << std::endl;
            out << std::endl;

            // Dimension is always 3 in our case
            out << "dimension" << std::endl;
            out << dimension << std::endl;
            out << std::endl;
        }

        /*!
         * @brief Write the cells for the MFEM mesh file (work only with
         * tetrahedra and hexaedra)
         * @details The structure of the MFEM file for cells is
         * [group_id] [cell_type] [id_vertex_0] [id_vertex_1] .....
         * cell_type is 4 for  tetrahedra and 5 for hexahedra.
         * group_id begin with 1
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_cells( const GeoModelMesh& geomodel_mesh, std::ofstream& out ) const
        {
            index_t nb_cells = geomodel_mesh.cells.nb();
            out << "elements" << std::endl;
            out << nb_cells << std::endl;
            for( index_t c = 0; c < nb_cells; c++ ) {
                out << geomodel_mesh.cells.region( c ) + mfem_offset << " ";
                out << cell_type_mfem[geomodel_mesh.cells.type( c )] << " ";
                for( index_t v = 0; v < geomodel_mesh.cells.nb_vertices( c ); v++ ) {
                    out << geomodel_mesh.cells.vertex( c, cell2mfem[v] ) << " ";
                }
                out << std::endl;
            }
            out << std::endl;
        }

        /*!
         * @brief Write the polygons for the MFEM mesh file (work only with
         * triangles and quads)
         * @details The structure of the MFEM file for polygons is
         * [group_id] [polygon_type] [id_vertex_0] [id_vertex_1] .....
         * polygon_type is 2 for triangles and 3 for the quads
         * group_id is continuous with the groupd indexes of the cells
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_polygons( const GeoModelMesh& geomodel_mesh, std::ofstream& out ) const
        {
            const GeoModelMeshPolygons& polygons = geomodel_mesh.polygons;
            out << "boundary" << std::endl;
            out << polygons.nb() << std::endl;
            for( index_t p = 0; p < polygons.nb(); p++ ) {
                index_t not_used = 0;
                out << polygons.surface( p ) + mfem_offset << " ";
                out << polygon_type_mfem[polygons.type( p, not_used )] << " ";
                for( index_t v = 0; v < polygons.nb_vertices( p ); v++ ) {
                    out << polygons.vertex( p, v ) << " ";
                }
                out << std::endl;
            }
            out << std::endl;
        }

        /*!
         * @brief Write the vertices for the MFEM mesh file
         * @details The structure of the MFEM file for vertices is
         * [x] [y] [z]
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_vertices( const GeoModelMesh& geomodel_mesh, std::ofstream& out ) const
        {
            out << "vertices" << std::endl;
            out << geomodel_mesh.vertices.nb() << std::endl;
            out << dimension << std::endl;
            for( index_t v = 0; v < geomodel_mesh.vertices.nb(); v++ ) {
                out << geomodel_mesh.vertices.vertex( v ) << std::endl;
            }
        }

    private:
        static const index_t dimension = 3;
    };

}
