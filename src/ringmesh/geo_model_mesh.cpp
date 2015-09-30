/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/geo_model_mesh.h>
#include <ringmesh/geo_model.h>

namespace {
    using namespace RINGMesh ;

    inline GeoModelMeshElement& cast_gmm_element(
        const GeoModel& M,
        GME::TYPE T,
        index_t i )
    {
        return dynamic_cast< GeoModelMeshElement& >( const_cast< GME& >( M.element(
            GME::gme_t( T, i ) ) ) ) ;
    }
}

namespace RINGMesh {

    GeoModelMeshVertices::GeoModelMeshVertices( GeoModelMesh& gmm, GEO::Mesh& mesh )
        : gmm_( gmm ), gm_( gmm.model() ), mesh_( mesh )
    {
    }

    bool GeoModelMeshVertices::is_initialized() const
    {
        return gmm_.mesh().vertices.nb() > 0 ;
    }

    void GeoModelMeshVertices::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshVertices* >( this )->initialize() ;
        }
    }

    void GeoModelMeshVertices::initialize()
    {
        mesh_.clear() ;

        // Total number of vertices in the
        // Corners, Lines, and Surfaces of the GeoModel
        index_t nb = 0 ;
        for( index_t t = GME::CORNER; t < GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            for( index_t e = 0; e < gm_.nb_elements( T ); ++e ) {
                nb += gm_.mesh_element( GME::gme_t( T, e ) ).nb_vertices() ;
            }
        }
        // Get out if no vertices
        if( nb == 0 ) {
            return ;
        }

        // Fill the vertices
        mesh_.vertices.create_vertices( nb ) ;
        gme_vertices_.resize( nb ) ;

        index_t index = 0 ;
        for( index_t t = GME::CORNER; t < GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            for( index_t e = 0; e < gm_.nb_elements( T ); ++e ) {
                GeoModelMeshElement& E = cast_gmm_element( gm_, T, e ) ;
                GEO::Memory::copy( mesh_.vertices.point_ptr( index ),
                    E.vertex( 0 ).data(), 3 * E.nb_vertices() * sizeof(double) ) ;
                for( index_t v = 0; v < E.nb_vertices(); v++ ) {
                    // Global index stored at BME level
                    E.set_model_vertex_id( v, index ) ;
                    // Index in the BME stored at global level
                    gme_vertices_[index].push_back( VertexInGME( E.gme_id(), v ) ) ;
                    // Global vertex index increment
                    index++ ;
                }
            }
        }
        // Remove colocated vertices
//        remove_colocated() ;
    }

    /*******************************************************************************/

    GeoModelMeshCells::GeoModelMeshCells( GeoModelMesh& gmm, GEO::Mesh& mesh )
        : gmm_( gmm ), gm_( gmm.model() ), mesh_( mesh )
    {
    }

    bool GeoModelMeshCells::is_initialized() const
    {
        return gmm_.mesh().cells.nb() > 0 ;
    }

    void GeoModelMeshCells::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshCells* >( this )->initialize() ;
        }
    }

    void GeoModelMeshCells::initialize()
    {

        gmm_.vertices.test_and_initialize() ;

        // Total number of  cells
        index_t nb = 0 ;
        // Table with the number of tet, hex, prism, pyramid and connector
        index_t nb_cells_by_type[GEO::MESH_NB_CELL_TYPES] ;
        for( index_t i = 0; i < GEO::MESH_NB_CELL_TYPES; ++i ) {
            nb_cells_by_type[i] = 0 ;
        }
        //
        std::vector< index_t > cells_vertices_id ;

        for( index_t r = 0; r < gm_.nb_regions(); ++r ) {
            nb += gm_.region( r ).nb_cells() ;
        }

        std::vector< GEO::MeshCellType > cells_vertices_type( nb ) ;

        index_t index = 0 ;
        for( index_t r = 0; r < gm_.nb_regions(); ++r ) {
            const GEO::Mesh& cur_region_mesh = gm_.region( r ).mesh() ;
            for( index_t c = 0; c < gm_.region( r ).nb_cells(); ++c ) {
                GEO::MeshCellType cur_cell_type = cur_region_mesh.cells.type( c ) ;
                nb_cells_by_type[cur_cell_type]++ ;
                cells_vertices_id[index] = cur_cell_type ;
            }
        }

        // Get out if no cells
        if( nb == 0 ) {
            return ;
        }

        // Fill the cells
        for( index_t i = 0; i < GEO::MESH_NB_CELL_TYPES; ++i ) {
            mesh_.cells.create_cells( nb_cells_by_type[i], GEO::MeshCellType( i ) ) ;
        }

        for( index_t r = 0; r < gm_.nb_regions(); ++r ) {
            const Region& cur_region = gm_.region( r ) ;
            const GEO::Mesh& cur_region_mesh = cur_region.mesh() ;
            for( index_t c = 0; c < gm_.region( r ).nb_cells(); ++c ) {
                GEO::MeshCellType cur_cell_type = cur_region_mesh.cells.type( c ) ;
                if( cur_cell_type == GEO::MESH_TET ) {
                    mesh_.cells.create_tet(
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 0 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 1 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 2 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 3 ) ) ) ;
                }
                if( cur_cell_type == GEO::MESH_HEX ) {
                    mesh_.cells.create_hex(
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 0 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 1 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 2 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 3 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 4 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 5 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 6 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 7 ) ) ) ;
                }
                else if( cur_cell_type == GEO::MESH_PRISM ) {
                    mesh_.cells.create_prism(
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 0 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 1 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 2 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 3 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 4 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 5 ) ) ) ;
                }
                else if( cur_cell_type == GEO::MESH_PYRAMID ) {
                    mesh_.cells.create_pyramid(
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 0 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 1 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 2 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 3 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 4 ) ) ) ;
                }
                else if( cur_cell_type == GEO::MESH_CONNECTOR ) {
                    mesh_.cells.create_connector(
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 0 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 1 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 2 ) ),
                        cur_region.model_vertex_id(
                            cur_region_mesh.cells.vertex( c, 3 ) ) ) ;
                }

                else {
                    ringmesh_assert_not_reached ;
                }

                // Retrieve the adjacencies
                mesh_.cells.connect() ;

            }
        }

    }

    /*******************************************************************************/

    GeoModelMesh::GeoModelMesh( const GeoModel& gm )
        : gm_( gm ), mesh_( new GEO::Mesh ), vertices( *this, *mesh_ ), cells(*this, *mesh_)
    {
    }

    GeoModelMesh::~GeoModelMesh()
    {
        delete mesh_ ;
    }

}
