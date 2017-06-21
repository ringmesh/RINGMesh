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
    static const std::string tet_name_in_aster_mail_file = "TETRA4";
    static const std::string hex_name_in_aster_mail_file = "HEXA10";
    static const std::string prism_name_in_aster_mail_file = "PENTA6";
    static const std::string pyr_name_in_aster_mail_file = "PYRAM5";

    static const std::string* cell_name_in_aster_mail_file[4] = {
        &tet_name_in_aster_mail_file, &hex_name_in_aster_mail_file,
        &prism_name_in_aster_mail_file, &pyr_name_in_aster_mail_file };

    static const std::string triangle_name_in_aster_mail_file = "TRIA3";
    static const std::string quad_name_in_aster_mail_file = "QUAD4";

    static const std::string* polygon_name_in_aster_mail_file[2] = {
        &triangle_name_in_aster_mail_file, &quad_name_in_aster_mail_file };
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
    class AsterIOHandler final: public GeoModelIOHandler< 3 > {
    public:
        virtual bool load( const std::string& filename, GeoModel< 3 >& geomodel ) final
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from Code_Aster mesh not implemented yet" );
            return false;
        }
        virtual void save(
            const GeoModel< 3 >& geomodel,
            const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );
            const RINGMesh::GeoModelMesh< 3 >& geomodel_mesh = geomodel.mesh;

            write_title( out, geomodel );

            write_vertices( out, geomodel_mesh );

            write_cells( geomodel, out );

            write_polygons( geomodel, out );

            write_regions( geomodel, out );

            write_interfaces( geomodel, out );

            out << "FIN" << std::endl;

            out.close();
        }

    private:

        void write_title( std::ofstream& out, const GeoModel< 3 >& geomodel ) const
        {
            out << "TITRE" << std::endl;
            out << geomodel.name() << std::endl;
            out << "FINSF" << std::endl;
        }
        void write_vertices(
            std::ofstream& out,
            const RINGMesh::GeoModelMesh< 3 >& geomodel_mesh ) const
        {
            out << "COOR_3D" << std::endl;
            for( index_t v = 0; v < geomodel_mesh.vertices.nb(); v++ ) {
                out << "V" << v << " " << geomodel_mesh.vertices.vertex( v )
                    << std::endl;
            }
            out << "FINSF" << std::endl;
        }

        void write_cells( const GeoModel< 3 >& geomodel, std::ofstream& out ) const
        {
            const GeoModelMesh< 3 >& geomodel_mesh = geomodel.mesh;
            for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                // -1 Because connectors doesn't exist in aster
                for( index_t ct = 0; ct < GEO::MESH_NB_CELL_TYPES - 1; ct++ ) {
                    if( geomodel_mesh.cells.nb_cells( r, GEO::MeshCellType( ct ) )
                        > 0 ) {
                        write_cells_in_region( GEO::MeshCellType( ct ), r,
                            geomodel_mesh, out );
                    }
                }
            }
        }

        void write_polygons(
            const GeoModel< 3 >& geomodel,
            std::ofstream& out ) const
        {
            const GeoModelMesh< 3 >& geomodel_mesh = geomodel.mesh;
            for( index_t s = 0; s < geomodel.nb_surfaces(); s++ ) {
                // -1 because polygons doesn' t exist in aster
                for( index_t pt = 0; pt < GeoModelMeshPolygons < 3 > ::ALL - 1;
                    pt++ ) {
                    if( geomodel_mesh.polygons.nb_polygons( s,
                        GeoModelMeshPolygons < 3 > ::PolygonType( pt ) ) > 0 ) {
                        write_polygons_in_interface(
                            GeoModelMeshPolygons < 3 > ::PolygonType( pt ), s,
                            geomodel_mesh, out );
                    }
                }
            }
        }
        void write_cells_in_region(
            const GEO::MeshCellType& cell_type,
            index_t region,
            const GeoModelMesh< 3 >& geomodel_mesh,
            std::ofstream& out ) const
        {
            out << *cell_name_in_aster_mail_file[cell_type] << std::endl;
            for( index_t c = 0;
                c < geomodel_mesh.cells.nb_cells( region, cell_type ); c++ ) {
                index_t global_id = geomodel_mesh.cells.cell( region, c, cell_type );
                out << "C" << global_id << " ";
                for( index_t v = 0; v < geomodel_mesh.cells.nb_vertices( c ); v++ ) {
                    out << "V" << geomodel_mesh.cells.vertex( global_id, v ) << " ";
                }
                out << std::endl;
            }
            out << "FINSF" << std::endl;
        }

        void write_polygons_in_interface(
            const GeoModelMeshPolygons< 3 >::PolygonType& polygon_type,
            index_t surface,
            const RINGMesh::GeoModelMesh< 3 >& mesh,
            std::ofstream& out ) const
        {
            out << *polygon_name_in_aster_mail_file[polygon_type] << std::endl;
            for( index_t p = 0;
                p < mesh.polygons.nb_polygons( surface, polygon_type ); p++ ) {
                index_t global_id = mesh.polygons.polygon( surface, p,
                    polygon_type );
                out << "F" << global_id << " ";
                for( index_t v = 0; v < mesh.polygons.nb_vertices( p ); v++ ) {
                    out << "V" << mesh.polygons.vertex( global_id, v ) << " ";
                }
                out << std::endl;
            }
            out << "FINSF" << std::endl;
        }

        void write_regions( const GeoModel< 3 >& geomodel, std::ofstream& out ) const
        {
            for( index_t r = 0; r < geomodel.nb_regions(); r++ ) {
                if( geomodel.region( r ).is_meshed() ) {
                    out << "GROUP_MA" << std::endl;
                    out << geomodel.region( r ).name() << std::endl;
                    for( index_t c = 0; c < geomodel.mesh.cells.nb_cells( r );
                        c++ ) {
                        out << "C" << geomodel.mesh.cells.cell( r, c ) << std::endl;
                    }
                    out << "FINSF" << std::endl;
                }
            }
        }

        void write_interfaces(
            const GeoModel< 3 >& geomodel,
            std::ofstream& out ) const
        {
            for( index_t inter = 0;
                inter
                    < geomodel.nb_geological_entities(
                        Interface < 3 > ::type_name_static() ); inter++ ) {
                const GeoModelGeologicalEntity< 3 >& cur_interface =
                    geomodel.geological_entity( Interface < 3 > ::type_name_static(),
                        inter );
                for( index_t s = 0; s < cur_interface.nb_children(); s++ ) {
                    index_t surface_id = cur_interface.child( s ).index();
                    out << "GROUP_MA" << std::endl;
                    out << cur_interface.name() << "_" << s << std::endl;
                    for( index_t p = 0;
                        p < geomodel.mesh.polygons.nb_polygons( surface_id ); p++ ) {
                        out << "F" << geomodel.mesh.polygons.polygon( surface_id, p )
                            << std::endl;
                    }
                    out << "FINSF" << std::endl;
                }

                out << "GROUP_MA" << std::endl;
                out << cur_interface.name() << std::endl;
                for( index_t s = 0; s < cur_interface.nb_children(); s++ ) {
                    index_t surface_id = cur_interface.child( s ).index();
                    for( index_t p = 0;
                        p < geomodel.mesh.polygons.nb_polygons( surface_id ); p++ ) {
                        out << "F" << geomodel.mesh.polygons.polygon( surface_id, p )
                            << std::endl;
                    }
                }
                out << "FINSF" << std::endl;
            }
        }
    };

}
