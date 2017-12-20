/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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
    static const std::string tet_name_in_aster_mail_file = "TETRA4";
    static const std::string hex_name_in_aster_mail_file = "HEXA10";
    static const std::string prism_name_in_aster_mail_file = "PENTA6";
    static const std::string pyr_name_in_aster_mail_file = "PYRAM5";

    static const std::string* cell_name_in_aster_mail_file[4] = {
        &tet_name_in_aster_mail_file, &hex_name_in_aster_mail_file,
        &prism_name_in_aster_mail_file, &pyr_name_in_aster_mail_file
    };

    static const std::string triangle_name_in_aster_mail_file = "TRIA3";
    static const std::string quad_name_in_aster_mail_file = "QUAD4";

    static const std::string* polygon_name_in_aster_mail_file[2] = {
        &triangle_name_in_aster_mail_file, &quad_name_in_aster_mail_file
    };
    /*!
     * @brief Export to the .mail mesh format of code aster
     * @details The descriptor of the .mail is available here:
     * http://www.code-aster.org/doc/v12/fr/man_u/u3/u3.01.00.pdf
     * Aster support multi-element mesh, so the export is region
     * based (the cells are written region by region).
     * The group of cells/polygons in aster are handle by "GROUP_MA"
     * Here, there will be one group for each Region, one group
     * for each Surface and one group for each Interfaces. It doesn't
     * matter if one polygon or cell is in several group
     *  The name of the Regions are the one given
     * by the GeoModel, the name of the Surfaces are the one given
     * by the parent Interface + the index of the child
     * @warning It supposes you have the mesh duplicate around the
     * faults if you want to use friction laws in aster
     */
    class AsterIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );
            const RINGMesh::GeoModelMesh3D& geomodel_mesh = geomodel.mesh;

            write_title( out, geomodel );
            write_vertices( out, geomodel_mesh );
            write_cells( geomodel, out );
            write_polygons( geomodel, out );
            write_regions( geomodel, out );
            write_interfaces( geomodel, out );
            out << "FIN" << EOL;
            out << std::flush;
        }

    private:
        void write_title( std::ofstream& out, const GeoModel3D& geomodel ) const
        {
            out << "TITRE" << EOL;
            out << geomodel.name() << EOL;
            out << "FINSF" << EOL;
        }
        void write_vertices( std::ofstream& out,
            const RINGMesh::GeoModelMesh3D& geomodel_mesh ) const
        {
            out << "COOR_3D" << EOL;
            for( auto v : range( geomodel_mesh.vertices.nb() ) )
            {
                out << "V" << v << " " << geomodel_mesh.vertices.vertex( v )
                    << EOL;
            }
            out << "FINSF" << EOL;
        }

        void write_cells( const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            const GeoModelMesh3D& geomodel_mesh = geomodel.mesh;
            for( auto r : range( geomodel.nb_regions() ) )
            {
                // -1 Because connectors doesn't exist in aster
                for( auto ct : range( CellType::UNCLASSIFIED ) )
                {
                    if( geomodel_mesh.cells.nb_cells( r, CellType( ct ) ) > 0 )
                    {
                        write_cells_in_region(
                            CellType( ct ), r, geomodel_mesh, out );
                    }
                }
            }
        }

        void write_polygons(
            const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            const GeoModelMesh3D& geomodel_mesh = geomodel.mesh;
            for( const auto& surface : geomodel.surfaces() )
            {
                // -1 because polygons doesn' t exist in aster
                for( auto pt :
                    range( to_underlying_type( PolygonType::UNDEFINED ) - 1 ) )
                {
                    if( geomodel_mesh.polygons.nb_polygons(
                            surface.index(), PolygonType( pt ) )
                        > 0 )
                    {
                        write_polygons_in_interface( PolygonType( pt ),
                            surface.index(), geomodel_mesh, out );
                    }
                }
            }
        }
        void write_cells_in_region( const CellType& cell_type,
            index_t region,
            const GeoModelMesh3D& geomodel_mesh,
            std::ofstream& out ) const
        {
            out << *cell_name_in_aster_mail_file[to_underlying_type(
                       cell_type )]
                << EOL;
            for( auto c :
                range( geomodel_mesh.cells.nb_cells( region, cell_type ) ) )
            {
                index_t global_id =
                    geomodel_mesh.cells.cell( region, c, cell_type );
                out << "C" << global_id << " ";
                for( auto v : range( geomodel_mesh.cells.nb_vertices( c ) ) )
                {
                    out << "V"
                        << geomodel_mesh.cells.vertex(
                               ElementLocalVertex( global_id, v ) )
                        << " ";
                }
                out << EOL;
            }
            out << "FINSF" << EOL;
        }

        void write_polygons_in_interface( const PolygonType& polygon_type,
            index_t surface,
            const RINGMesh::GeoModelMesh3D& mesh,
            std::ofstream& out ) const
        {
            out << *polygon_name_in_aster_mail_file[to_underlying_type(
                       polygon_type )]
                << EOL;
            for( auto p :
                range( mesh.polygons.nb_polygons( surface, polygon_type ) ) )
            {
                index_t global_id =
                    mesh.polygons.polygon( surface, p, polygon_type );
                out << "F" << global_id << " ";
                for( auto v : range( mesh.polygons.nb_vertices( p ) ) )
                {
                    out << "V"
                        << mesh.polygons.vertex(
                               ElementLocalVertex( global_id, v ) )
                        << " ";
                }
                out << EOL;
            }
            out << "FINSF" << EOL;
        }

        void write_regions(
            const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            for( const auto& region : geomodel.regions() )
            {
                if( region.is_meshed() )
                {
                    out << "GROUP_MA" << EOL;
                    out << region.name() << EOL;
                    for( auto c : range(
                             geomodel.mesh.cells.nb_cells( region.index() ) ) )
                    {
                        out << "C"
                            << geomodel.mesh.cells.cell( region.index(), c )
                            << EOL;
                    }
                    out << "FINSF" << EOL;
                }
            }
        }

        void write_interfaces(
            const GeoModel3D& geomodel, std::ofstream& out ) const
        {
            for( auto& cur_interface :
                geomodel.geol_entities( Interface3D::type_name_static() ) )
            {
                for( auto s : range( cur_interface.nb_children() ) )
                {
                    index_t surface_id = cur_interface.child( s ).index();
                    out << "GROUP_MA" << EOL;
                    out << cur_interface.name() << "_" << s << EOL;
                    for( auto p : range( geomodel.mesh.polygons.nb_polygons(
                             surface_id ) ) )
                    {
                        out << "F"
                            << geomodel.mesh.polygons.polygon( surface_id, p )
                            << EOL;
                    }
                    out << "FINSF" << EOL;
                }

                out << "GROUP_MA" << EOL;
                out << cur_interface.name() << EOL;
                for( auto s : range( cur_interface.nb_children() ) )
                {
                    index_t surface_id = cur_interface.child( s ).index();
                    for( auto p : range( geomodel.mesh.polygons.nb_polygons(
                             surface_id ) ) )
                    {
                        out << "F"
                            << geomodel.mesh.polygons.polygon( surface_id, p )
                            << EOL;
                    }
                }
                out << "FINSF" << EOL;
            }
        }
    };
}
