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





    /*!
     * @brief Export for the GMSH format 2.2 which is described here:
     * http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format
     * NB : Mesh entities are also exported
     */
    class StradivariusIOHandler final: public GeoModelIOHandler< 3 > {
    public:
        void load( const std::string& filename, GeoModel< 3 >& geomodel ) final
        {

        }
        void save( const GeoModel< 3 >& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

        }

    private:
        /*!
         * @brief Find the gmsh type on an element using
         * the number of vertices and the mesh entity index
         * in which the element belong
         */
        index_t find_gmsh_element_type(
            index_t nb_vertices,
            index_t mesh_entity_type_index )
        {
            return element_type[nb_vertices + mesh_entity_type_index];
        }

        /*!
         * @brief Find the gmsh local vertex index of an element using
         * the number of vertices and the mesh entity index
         * in which the element belong
         */
        index_t find_gmsh_element_local_vertex_id(
            index_t nb_vertices,
            index_t mesh_entity_type_index,
            index_t local_vertex_index )
        {
            return vertices_in_elements[nb_vertices + mesh_entity_type_index][local_vertex_index];
        }

        /*!
         * @brief Count all the elements in GMSH
         * an Element can be :
         * - A Corner
         * - An Edge
         * - A Triangle
         * - A Quad
         * - A Tetrahedron
         * - A Pyramid
         * - A Prism
         * - An Hexaedron
         */
        index_t count_elements( const GeoModel< 3 >& geomodel )
        {
            index_t nb_elements = 0;
            const std::vector< MeshEntityType >& gmme_types =
                geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();
            for( MeshEntityType cur_mesh_entity_type : gmme_types ) {
                for( index_t index_of_gmme_of_the_current_type : range(
                    geomodel.nb_mesh_entities( cur_mesh_entity_type ) ) ) {
                    gmme_id cur_gmme_id = gmme_id( cur_mesh_entity_type,
                        index_of_gmme_of_the_current_type );
                    const GeoModelMeshEntity< 3 >& cur_gmme = geomodel.mesh_entity(
                        cur_gmme_id );
                    nb_elements += cur_gmme.nb_mesh_elements();
                }
            }
            return nb_elements;
        }
    };
}
