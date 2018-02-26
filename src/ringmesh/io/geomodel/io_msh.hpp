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

namespace {
    // Vertices order for a vertex in GMSH
    static index_t vertices_in_vertex[1] = { 0 };
    // Vertices order for an edge in GMSH
    static index_t vertices_in_edge[2] = { 0, 1 };
    // Vertices order for a triangle in GMSH
    static index_t vertices_in_triangle[3] = { 0, 1, 2 };
    // Vertices order for a quad in GMSH
    static index_t vertices_in_quad[4] = { 0, 1, 2, 3 };
    // Vertices order for a tetrahedron in GMSH
    static index_t vertices_in_tetrahedron[4] = { 0, 1, 2, 3 };
    // Vertices order for an hexahedron in GMSH
    static index_t vertices_in_hexahedron[8] = { 4, 0, 5, 1, 7, 3, 6, 2 };
    // Vertices order for a prism in GMSH
    static index_t vertices_in_prism[6] = { 0, 1, 2, 3, 4, 5 };
    // Vertices order for a pyramid in GMSH
    static index_t vertices_in_pyramid[5] = { 0, 1, 2, 3, 4 };

    // GMSH count begin at 1
    index_t gmsh_offset = 1;

    // This is a tricky table that associate an unique id with the id of
    // elements in GMSH. The unique id is computed as the sum
    // of the MeshEntityType index and the number of vertices
    // of the elements
    // Vertex :      0 + 1 = 1
    // Edge :        1 + 2 = 3
    // Triangle :    2 + 3 = 5
    // Quad :        2 + 4 = 6
    // Tetrahedron : 3 + 4 = 7
    // Pyramid :     3 + 5 = 8
    // Prism         3 + 6 = 9
    // Hexahedron :  3 + 8 = 11
    index_t element_type[12] =
        { NO_ID, 15, NO_ID, 1, NO_ID, 2, 3, 4, 7, 6, NO_ID, 5 };

    // This is a tricky table that associate an unique id with another table
    // containing the ordering of vertices inside the elementS
    index_t* vertices_in_elements[12] = { nullptr, vertices_in_vertex, nullptr,
                                          vertices_in_edge, nullptr,
                                          vertices_in_triangle, vertices_in_quad,
                                          vertices_in_tetrahedron,
                                          vertices_in_pyramid, vertices_in_prism,
                                          nullptr, vertices_in_hexahedron };

    // The physical id is a mandatory value for GMSH which is not used
    index_t physical_id = 0;

    // in GMSH, a tag is a physical id (not used) or a geometry id
    index_t nb_of_tags = 2;

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
     * An Element can be :
     * - A Corner
     * - An Edge
     * - A Triangle
     * - A Quad
     * - A Tetrahedron
     * - A Pyramid
     * - A Prism
     * - An Hexaedron
     */
    template< index_t DIMENSION >
    index_t count_elements( const GeoModel< DIMENSION >& geomodel )
    {
        index_t nb_elements { 0 };
        const auto& gmme_types =
            geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();
        for( const auto& cur_mesh_entity_type : gmme_types ) {
            for( auto index_of_gmme_of_the_current_type : range(
                geomodel.nb_mesh_entities( cur_mesh_entity_type ) ) ) {
                gmme_id cur_gmme_id { gmme_id( cur_mesh_entity_type,
                    index_of_gmme_of_the_current_type ) };
                const GeoModelMeshEntity< DIMENSION >& cur_gmme =
                    geomodel.mesh_entity( cur_gmme_id );
                nb_elements += cur_gmme.nb_mesh_elements();
            }
        }
        return nb_elements;
    }

    vec3 vertex_in_3d( const vec2& vertex )
    {
        vec3 result;
        result.x = vertex.x;
        result.y = vertex.y;
        result.z = 0;
        return result;
    }

    vec3 vertex_in_3d( const vec3& vertex )
    {
        return vertex;
    }

    /*!
     * @brief Export GeoModel3D for the GMSH format 2.2 which is
     * described here:
     * http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format
     * NB : Mesh entities are also exported
     */
    template< index_t DIMENSION >
    class MSHIOHandler final : public GeoModelOutputHandler< DIMENSION > {
    public:
        void save(
            const GeoModel< DIMENSION >& geomodel,
            const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

            out << "$MeshFormat" << EOL;
            out << "2.2 0 8" << EOL;
            out << "$EndMeshFormat" << EOL;

            out << "$Nodes" << EOL;
            out << geomodel.mesh.vertices.nb() << EOL;
            for( auto v : range( geomodel.mesh.vertices.nb() ) ) {
                out << v + gmsh_offset << SPACE
                    << vertex_in_3d( geomodel.mesh.vertices.vertex( v ) ) << EOL;
            }
            out << "$EndNodes" << EOL;

            out << "$Elements" << EOL;
            out << count_elements( geomodel ) << EOL;
            const auto& gmme_types =
                geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types();

            index_t element_index { 1 };
            for( auto gmme_type_index : range(
                geomodel.entity_type_manager().mesh_entity_manager.nb_mesh_entity_types() ) ) {
                MeshEntityType cur_mesh_entity_type { gmme_types[gmme_type_index] };
                for( auto index_of_gmme_of_the_current_type : range(
                    geomodel.nb_mesh_entities( cur_mesh_entity_type ) ) ) {
                    gmme_id cur_gmme_id = gmme_id( cur_mesh_entity_type,
                        index_of_gmme_of_the_current_type );
                    const GeoModelMeshEntity< DIMENSION >& cur_gmme =
                        geomodel.mesh_entity( cur_gmme_id );
                    for( auto elem_in_cur_gmme : range( cur_gmme.nb_mesh_elements() ) ) {
                        index_t nb_vertices_in_cur_element =
                            cur_gmme.nb_mesh_element_vertices( elem_in_cur_gmme );
                        index_t gmsh_element_type = find_gmsh_element_type(
                            nb_vertices_in_cur_element, gmme_type_index );
                        out << element_index++ << SPACE << gmsh_element_type << SPACE
                            << nb_of_tags << SPACE << physical_id << SPACE
                            << index_of_gmme_of_the_current_type + gmsh_offset
                            << SPACE;
                        for( auto v_index_in_cur_element : range(
                            nb_vertices_in_cur_element ) ) {
                            out
                                << geomodel.mesh.vertices.geomodel_vertex_id(
                                    cur_gmme_id,
                                    cur_gmme.mesh_element_vertex_index(
                                        ElementLocalVertex( elem_in_cur_gmme,
                                            find_gmsh_element_local_vertex_id(
                                                nb_vertices_in_cur_element,
                                                gmme_type_index,
                                                v_index_in_cur_element ) ) ) )
                                    + gmsh_offset << SPACE;
                        }
                        out << EOL;
                    }
                }
            }
            out << "$EndElements" << EOL;
            out << std::flush;
        }
    };
    ALIAS_2D_AND_3D (MSHIOHandler);
}
