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

    enum MFEM_geometry_type
    {
        POINT = 0,
        SEGMENT = 1,
        TRIANGLE = 2,
        SQUARE = 3,
        TETRAHEDRON = 4,
        CUBE = 5
    };

    /*!
     * Export for the MFEM format http://mfem.org/
     * Mesh file description : http://mfem.org/mesh-formats/#mfem-mesh-v10
     * "MFEM is a free, lightweight, scalable C++ library for finite element
     * methods"
     */
    template < index_t DIMENSION >
    class MFEMIOHandler final : public GeoModelOutputHandler< DIMENSION >
    {
    public:
        void save( const GeoModel< DIMENSION >& geomodel,
            const std::string& filename ) final
        {
            const auto& geomodel_mesh = geomodel.mesh;
            test_if_mesh_is_valid( geomodel.mesh );

            std::ofstream out( filename.c_str() );
            out.precision( 16 );

            write_header( out );
            write_elements( geomodel_mesh, out );
            write_boundaries( geomodel_mesh, out );
            write_vertices( geomodel_mesh, out );
            out << std::flush;
        }

    private:
        void test_if_mesh_is_valid(
            const GeoModelMesh< DIMENSION >& geomodel_mesh );

        /*!
         * @brief Write the header for the MFEM mesh file
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_header( std::ofstream& out ) const
        {
            out << "MFEM mesh v1.0" << EOL;
            out << EOL;
            out << "dimension" << EOL;
            out << DIMENSION << EOL;
            out << EOL;
        }

        void write_elements( const GeoModelMesh< DIMENSION >& geomodel_mesh,
            std::ofstream& out ) const;

        void write_boundaries( const GeoModelMesh< DIMENSION >& geomodel_mesh,
            std::ofstream& out ) const;

        /*!
         * @brief Write the vertices for the MFEM mesh file
         * @details The structure of the MFEM file for vertices is
         * [x] [y] [z]
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_vertices( const GeoModelMesh< DIMENSION >& geomodel_mesh,
            std::ofstream& out ) const
        {
            out << "vertices" << EOL;
            out << geomodel_mesh.vertices.nb() << EOL;
            out << DIMENSION << EOL;
            for( auto v : range( geomodel_mesh.vertices.nb() ) )
            {
                out << geomodel_mesh.vertices.vertex( v ) << EOL;
            }
        }
    };

    ALIAS_2D_AND_3D( MFEMIOHandler );

    template <>
    void MFEMIOHandler3D::test_if_mesh_is_valid(
        const GeoModelMesh3D& geomodel_mesh )
    {
        index_t nb_cells{ geomodel_mesh.cells.nb() };
        if( geomodel_mesh.cells.nb_tet() != nb_cells
            && geomodel_mesh.cells.nb_hex() != nb_cells )
        {
            throw RINGMeshException( "I/O", "Export to MFEM format works only "
                                            "with full tet or full hex "
                                            "format" );
        }
    }

    template <>
    void MFEMIOHandler2D::test_if_mesh_is_valid(
        const GeoModelMesh2D& geomodel_mesh )
    {
        if( geomodel_mesh.polygons.nb()
            != geomodel_mesh.polygons.nb_triangle() )
        {
            throw RINGMeshException( "I/O",
                "Export to MFEM format works only with triangles in 2D" );
        }
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
    template <>
    void MFEMIOHandler3D::write_elements(
        const GeoModelMesh3D& geomodel_mesh, std::ofstream& out ) const
    {
        index_t nb_cells{ geomodel_mesh.cells.nb() };
        out << "elements" << EOL;
        out << nb_cells << EOL;
        for( auto c : range( nb_cells ) )
        {
            out << geomodel_mesh.cells.region( c ) + mfem_offset << " ";
            out << cell_type_mfem[to_underlying_type(
                       geomodel_mesh.cells.type( c ) )]
                << " ";
            for( auto v : range( geomodel_mesh.cells.nb_vertices( c ) ) )
            {
                out << geomodel_mesh.cells.vertex(
                           ElementLocalVertex( c, cell2mfem[v] ) )
                    << " ";
            }
            out << EOL;
        }
        out << EOL;
    }

    template <>
    void MFEMIOHandler2D::write_elements(
        const GeoModelMesh2D& geomodel_mesh, std::ofstream& out ) const
    {
        index_t nb_triangles{ geomodel_mesh.polygons.nb_triangle() };
        out << "elements" << EOL;
        out << nb_triangles << EOL;
        for( auto c : range( nb_triangles ) )
        {
            out << geomodel_mesh.polygons.surface( c ) + mfem_offset << " ";
            out << TRIANGLE << " ";
            for( auto v : range( geomodel_mesh.polygons.nb_vertices( c ) ) )
            {
                out << geomodel_mesh.polygons.vertex(
                           ElementLocalVertex( c, v ) )
                    << " ";
            }
            out << EOL;
        }
        out << EOL;
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
    template <>
    void MFEMIOHandler3D::write_boundaries(
        const GeoModelMesh3D& geomodel_mesh, std::ofstream& out ) const
    {
        const GeoModelMeshPolygons3D& polygons = geomodel_mesh.polygons;
        out << "boundary" << EOL;
        out << polygons.nb() << EOL;
        for( auto p : range( polygons.nb() ) )
        {
            out << polygons.surface( p ) + mfem_offset << " ";
            PolygonType polygon_type;
            std::tie( polygon_type, std::ignore ) = polygons.type( p );
            out << polygon_type_mfem[to_underlying_type( polygon_type )] << " ";
            for( auto v : range( polygons.nb_vertices( p ) ) )
            {
                out << polygons.vertex( ElementLocalVertex( p, v ) ) << " ";
            }
            out << EOL;
        }
        out << EOL;
    }

    template <>
    void MFEMIOHandler2D::write_boundaries(
        const GeoModelMesh2D& geomodel_mesh, std::ofstream& out ) const
    {
        const GeoModelMeshEdges2D& edges = geomodel_mesh.edges;
        out << "boundary" << EOL;
        out << edges.nb() << EOL;
        for( auto p : range( edges.nb() ) )
        {
            out << edges.line( p ) + mfem_offset << " ";
            out << SEGMENT << " ";
            for( auto v : range( 2 ) )
            {
                out << edges.vertex( ElementLocalVertex( p, v ) ) << " ";
            }
            out << EOL;
        }
        out << EOL;
    }
}
