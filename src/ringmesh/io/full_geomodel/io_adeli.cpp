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
    // The reg phys field in the GMSH format is set to 0 for each element
    static const index_t reg_phys = 0 ;

    static const index_t adeli_point_type = 15 ;
    static const index_t adeli_line_type = 1 ;
    static const index_t adeli_triangle_type = 2 ;
    static const index_t adeli_tet_type = 4 ;
    static const index_t adeli_cell_types[4] = {
        adeli_point_type, adeli_line_type, adeli_triangle_type, adeli_tet_type } ;

    // The index begins at 1.
    static const index_t id_offset_adeli = 1 ;

    /*!
     * @brief export for ADELI
     * https://sif.info-ufr.univ-montp2.fr/?q=content/adeli
     * @details This export is in fact a V1.0 of the .msh file, suitable for
     * running Finite Element Simulation with the ADELI solver.
     * The description of the output file can be found here :
     * http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-1_002e0
     * First, nodes are written, then the elements. The elements are not separated.
     * Corners are written (with vertex), then Lines (with edges), then Surfaces
     * (with surfaces, then Regions (with tetrahedron)
     */
    class AdeliIOHandler final: public GeoModelIOHandler {
    public:
        virtual bool load( const std::string& filename, GeoModel& geomodel ) override
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from Adeli .msh mesh not implemented yet" ) ;
            return false ;
        }
        virtual void save( const GeoModel& geomodel, const std::string& filename ) override
        {
            std::ofstream out( filename.c_str() ) ;
            out.precision( 16 ) ;
            const RINGMesh::GeoModelMesh& geomodel_mesh = geomodel.mesh ;
            if( geomodel_mesh.cells.nb() != geomodel_mesh.cells.nb_tet() ) {
                {
                    throw RINGMeshException( "I/O",
                        "Adeli supports only tet meshes" ) ;
                }
            }

            write_vertices( geomodel_mesh, out ) ;

            index_t elt = 1 ;
            write_corners( geomodel, out, elt ) ;
            write_mesh_elements( geomodel, out, elt ) ;
        }

    private:
        void write_vertices( const GeoModelMesh& geomodel_mesh, std::ofstream& out )
        {
            out << "$NOD" << std::endl ;
            out << geomodel_mesh.vertices.nb() << std::endl ;
            for( index_t v = 0; v < geomodel_mesh.vertices.nb(); v++ ) {
                out << v + id_offset_adeli << " "
                    << geomodel_mesh.vertices.vertex( v ) << std::endl ;
            }
            out << "$ENDNOD" << std::endl ;
        }

        void write_corners(
            const GeoModel& geomodel,
            std::ofstream& out,
            index_t& elt )
        {
            out << "$ELM" << std::endl ;
            out << nb_total_elements( geomodel ) << std::endl ;
            for( index_t corner = 0; corner < geomodel.nb_corners(); corner++ ) {
                const Corner& cur_corner = geomodel.corner( corner ) ;
                out << elt++ << " " << adeli_cell_types[0] << " " << reg_phys << " "
                    << cur_corner.index() + id_offset_adeli << " "
                    << cur_corner.nb_vertices() << " "
                    << geomodel.mesh.vertices.geomodel_vertex_id(
                        cur_corner.gmme_id(), 0 ) + id_offset_adeli << std::endl ;
            }
        }

        void write_mesh_elements(
            const GeoModel& geomodel,
            std::ofstream& out,
            index_t& elt )
        {
            // Corners are already written so we start this loop at 1
            for( index_t geomodel_mesh_entities = 1;
                geomodel_mesh_entities < MeshEntityTypeManager::nb_mesh_entity_types();
                geomodel_mesh_entities++ ) {
                for( index_t entity = 0;
                    entity
                        < geomodel.nb_mesh_entities(
                            MeshEntityTypeManager::mesh_entity_types()[geomodel_mesh_entities] );
                    entity++ ) {
                    write_mesh_elements_for_a_mesh_entity(
                        geomodel.mesh_entity(
                            MeshEntityTypeManager::mesh_entity_types()[geomodel_mesh_entities],
                            entity ), adeli_cell_types[geomodel_mesh_entities], elt,
                        out ) ;
                }
            }
            out << "$ENDELM" << std::endl ;
        }

        index_t nb_total_elements( const GeoModel& geomodel )
        {
            // Because corners does not have mesh elements, but are considered as elements
            // in adeli, we have to count the vertex of each corner in a different
            // way
            index_t nb_mesh_entities = geomodel.nb_corners() ;
            for( index_t geomodel_mesh_entities = 1;
                geomodel_mesh_entities < MeshEntityTypeManager::nb_mesh_entity_types();
                geomodel_mesh_entities++ ) {
                for( index_t entity = 0;
                    entity
                        < geomodel.nb_mesh_entities(
                            MeshEntityTypeManager::mesh_entity_types()[geomodel_mesh_entities] );
                    entity++ ) {
                    nb_mesh_entities +=
                        geomodel.mesh_entity(
                            MeshEntityTypeManager::mesh_entity_types()[geomodel_mesh_entities],
                            entity ).nb_mesh_elements() ;
                }
            }
            return nb_mesh_entities ;
        }

        void write_mesh_elements_for_a_mesh_entity(
            const GeoModelMeshEntity& geomodel_mesh_entity,
            index_t cell_descriptor,
            index_t& elt_id,
            std::ofstream& out )
        {
            for( index_t elt = 0; elt < geomodel_mesh_entity.nb_mesh_elements();
                elt++ ) {
                out << elt_id++ << " " << cell_descriptor << " " << reg_phys << " "
                    << geomodel_mesh_entity.index() + id_offset_adeli << " "
                    << geomodel_mesh_entity.nb_mesh_element_vertices( elt ) << " " ;
                for( index_t v = 0;
                    v < geomodel_mesh_entity.nb_mesh_element_vertices( elt ); v++ ) {
                    out
                        << geomodel_mesh_entity.geomodel().mesh.vertices.geomodel_vertex_id(
                            geomodel_mesh_entity.gmme_id(), elt, v ) + id_offset_adeli
                        << " " ;
                }
                out << std::endl ;
            }
        }
    } ;
}
