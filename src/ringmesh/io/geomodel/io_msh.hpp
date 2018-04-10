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
    index_t element_type[12] = { NO_ID, 15, NO_ID, 1, NO_ID, 2, 3, 4, 7, 6,
        NO_ID, 5 };

    // This is a tricky table that associate an unique id with another table
    // containing the ordering of vertices inside the elementS
    index_t* vertices_in_elements[12] = { nullptr, vertices_in_vertex, nullptr,
        vertices_in_edge, nullptr, vertices_in_triangle, vertices_in_quad,
        vertices_in_tetrahedron, vertices_in_pyramid, vertices_in_prism,
        nullptr, vertices_in_hexahedron };

    // The physical id is a mandatory value for GMSH which is not used
    index_t physical_id = 0;

    // in GMSH, a tag is a physical id (not used) or a geometry id
    index_t nb_of_tags = 2;

    enum struct GMSHElementType
    {
        EDGE = 1,
        TRIANGLE = 2,
        QUADRANGLE = 3,
        TETRAHEDRON = 4,
        HEXAHEDRON = 5,
        PRISM = 6,
        PYRAMIDE = 7,
        VERTEX = 15,
        OTHER
    };

    /*!
     * @brief Find the gmsh type on an element using
     * the number of vertices and the mesh entity index
     * in which the element belong
     */
    index_t find_gmsh_element_type(
        index_t nb_vertices, index_t mesh_entity_type_index )
    {
        return element_type[nb_vertices + mesh_entity_type_index];
    }

    /*!
     * @brief Find the gmsh local vertex index of an element using
     * the number of vertices and the mesh entity index
     * in which the element belong
     */
    index_t find_gmsh_element_local_vertex_id( index_t nb_vertices,
        index_t mesh_entity_type_index,
        index_t local_vertex_index )
    {
        return vertices_in_elements[nb_vertices + mesh_entity_type_index]
                                   [local_vertex_index];
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
    template < index_t DIMENSION >
    index_t count_elements( const GeoModel< DIMENSION >& geomodel )
    {
        index_t nb_elements{ 0 };
        const auto& gmme_types = geomodel.entity_type_manager()
                                     .mesh_entity_manager.mesh_entity_types();
        for( const auto& cur_mesh_entity_type : gmme_types )
        {
            for( auto index_of_gmme_of_the_current_type :
                range( geomodel.nb_mesh_entities( cur_mesh_entity_type ) ) )
            {
                gmme_id cur_gmme_id{ gmme_id(
                    cur_mesh_entity_type, index_of_gmme_of_the_current_type ) };
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

    bool all_third_coord_null( const std::string& filename )
    {
        GEO::LineInput file_line{ filename };
        while( !file_line.eof() && file_line.get_line() )
        {
            file_line.get_fields();
            if( file_line.nb_fields() == 1
                && file_line.field_matches( 0, "$EndNodes" ) )
            {
                // All nodes read and all the third coordinates are equal to 0
                return true;
            }
            if( file_line.nb_fields() >= 4 )
            {
                if( !file_line.field_matches( 3, "0" ) )
                {
                    return false;
                }
            }
        }
        ringmesh_assert_not_reached;
        return false;
    }

    template < index_t DIMENSION >
    std::vector< vecn< DIMENSION > > read_nodes( GEO::LineInput& line_parser )
    {
        std::vector< vecn< DIMENSION > > nodes;
        // Skip header block
        while( !line_parser.eof() && line_parser.get_line() )
        {
            line_parser.get_fields();
            if( line_parser.nb_fields() == 1
                && line_parser.field_matches( 0, "$EndMeshFormat" ) )
            {
                break;
            }
        }
        while( !line_parser.eof() && line_parser.get_line() )
        {
            line_parser.get_fields();
            if( line_parser.nb_fields() == 1
                && line_parser.field_matches( 0, "$Nodes" ) )
            {
                // Reserve capacity of vector of nodes
                line_parser.get_line();
                line_parser.get_fields();
                auto nb_nodes = line_parser.field_as_uint( 0 );
                nodes.reserve( nb_nodes );
                continue;
            }
            if( line_parser.nb_fields() == 1
                && line_parser.field_matches( 0, "$EndNodes" ) )
            {
                return nodes;
            }
            ringmesh_assert( line_parser.nb_fields() >= 4 );
            vecn< DIMENSION > cur_node;
            for( auto coord : range( DIMENSION ) )
            {
                cur_node[coord] = line_parser.field_as_double( 1 + coord );
            }
            nodes.push_back( cur_node );
        }
        ringmesh_assert_not_reached;
        return nodes;
    }

    template < index_t DIMENSION >
    void read_elements_base( GEO::LineInput& line_parser,
        GeoModel< DIMENSION >& geomodel,
        const std::vector< vecn< DIMENSION > >& nodes )
    {
        std::map< index_t, gmme_id > map_gmsh_ref_to_mesh_entities;
        std::map< index_t, std::vector< index_t > > map_gmsh_ref_to_added_nodes;
        GeoModelBuilder< DIMENSION > geomodel_builder( geomodel );
        while( !line_parser.eof() && line_parser.get_line() )
        {
            line_parser.get_fields();
            if( line_parser.nb_fields() == 1
                && line_parser.field_matches( 0, "$EndElements" ) )
            {
                return;
            }
            if( line_parser.nb_fields() < 3 )
            {
                continue;
            }
            // Read one element
            GMSHElementType cur_element_type = static_cast< GMSHElementType >(
                line_parser.field_as_uint( 1 ) );
            index_t nb_tags{ line_parser.field_as_uint( 2 ) };

            if( cur_element_type == GMSHElementType::EDGE )
            {
                index_t gmsh_entity = line_parser.field_as_uint( 4 );
                // Get or create corresponding mesh entity
                gmme_id cur_gmme_id{
                    map_gmsh_ref_to_mesh_entities[gmsh_entity]
                };
                if( !cur_gmme_id.is_defined() )
                {
                    cur_gmme_id = geomodel_builder.topology.create_mesh_entity(
                        line_type_name_static() );
                    map_gmsh_ref_to_mesh_entities[gmsh_entity] = cur_gmme_id;
                }
                ringmesh_assert( cur_gmme_id.is_defined() );

                auto line_builder =
                    geomodel_builder.geometry.create_line_builder(
                        cur_gmme_id.index() );
                auto& added_nodes = map_gmsh_ref_to_added_nodes[gmsh_entity];

                // Add corresponding nodes
                index_t edge_gmme_ids[2];
                for( auto element_vertex_id : range( 2 ) )
                {
                    index_t gmsh_node_id{ line_parser.field_as_uint(
                        3 + nb_tags + element_vertex_id ) };
                    edge_gmme_ids[element_vertex_id] =
                        find( added_nodes, gmsh_node_id );
                    if( edge_gmme_ids[element_vertex_id] == NO_ID )
                    {
                        added_nodes.push_back( gmsh_node_id );
                        edge_gmme_ids[element_vertex_id] =
                            line_builder->create_vertex(
                                nodes[gmsh_node_id - gmsh_offset] );
                    }
                    ringmesh_assert(
                        edge_gmme_ids[element_vertex_id] != NO_ID );
                }
                // Add corresponding edge
                line_builder->create_edge( edge_gmme_ids[0], edge_gmme_ids[1] );
            }
            else if( cur_element_type == GMSHElementType::TRIANGLE )
            {
                index_t gmsh_entity = line_parser.field_as_uint( 4 );
                // Get or create corresponding mesh entity
                gmme_id cur_gmme_id{
                    map_gmsh_ref_to_mesh_entities[gmsh_entity]
                };
                if( !cur_gmme_id.is_defined() )
                {
                    cur_gmme_id = geomodel_builder.topology.create_mesh_entity(
                        surface_type_name_static() );
                    map_gmsh_ref_to_mesh_entities[gmsh_entity] = cur_gmme_id;
                }
                ringmesh_assert( cur_gmme_id.is_defined() );

                auto surface_builder =
                    geomodel_builder.geometry.create_surface_builder(
                        cur_gmme_id.index() );
                auto& added_nodes = map_gmsh_ref_to_added_nodes[gmsh_entity];

                // Add corresponding nodes
                std::vector< index_t > triangle_gmme_ids( 3 );
                for( auto element_vertex_id : range( 3 ) )
                {
                    index_t gmsh_node_id{ line_parser.field_as_uint(
                        3 + nb_tags + element_vertex_id ) };
                    triangle_gmme_ids[element_vertex_id] =
                        find( added_nodes, gmsh_node_id );
                    if( triangle_gmme_ids[element_vertex_id] == NO_ID )
                    {
                        added_nodes.push_back( gmsh_node_id );
                        triangle_gmme_ids[element_vertex_id] =
                            surface_builder->create_vertex(
                                nodes[gmsh_node_id - gmsh_offset] );
                    }
                    ringmesh_assert(
                        triangle_gmme_ids[element_vertex_id] != NO_ID );
                }
                // Add corresponding edge
                surface_builder->create_polygon( triangle_gmme_ids );
            }
            else if( cur_element_type == GMSHElementType::QUADRANGLE )
            {
                throw RINGMeshException( "I/O",
                    "Import of GeoModels form GMSH with "
                    "quadrangles is not yet implemented." );
            }
        }
    }

    void read_volume_elements( const std::string& filename,
        GeoModel3D& geomodel,
        const std::vector< vec3 >& nodes )
    {
        std::map< index_t, gmme_id > map_gmsh_ref_to_mesh_entities;
        std::map< index_t, std::vector< index_t > > map_gmsh_ref_to_added_nodes;
        GeoModelBuilder3D geomodel_builder( geomodel );
        GEO::LineInput line_parser{ filename };
        // Go to Elements section
        while( !line_parser.eof() && line_parser.get_line() )
        {
            line_parser.get_fields();
            if( line_parser.nb_fields() == 1
                && line_parser.field_matches( 0, "$Elements" ) )
            {
                break;
            }
        }
        while( !line_parser.eof() && line_parser.get_line() )
        {
            line_parser.get_fields();
            if( line_parser.nb_fields() == 1
                && line_parser.field_matches( 0, "$EndElements" ) )
            {
                return;
            }
            if( line_parser.nb_fields() < 3 )
            {
                continue;
            }
            // Read one element
            GMSHElementType cur_element_type = static_cast< GMSHElementType >(
                line_parser.field_as_uint( 1 ) );
            index_t nb_tags{ line_parser.field_as_uint( 2 ) };

            if( cur_element_type == GMSHElementType::TETRAHEDRON )
            {
                index_t gmsh_entity = line_parser.field_as_uint( 4 );
                // Get or create corresponding mesh entity
                gmme_id cur_gmme_id{
                    map_gmsh_ref_to_mesh_entities[gmsh_entity]
                };
                if( !cur_gmme_id.is_defined() )
                {
                    cur_gmme_id = geomodel_builder.topology.create_mesh_entity(
                        region_type_name_static() );
                    map_gmsh_ref_to_mesh_entities[gmsh_entity] = cur_gmme_id;
                }
                ringmesh_assert( cur_gmme_id.is_defined() );

                auto region_builder =
                    geomodel_builder.geometry.create_region_builder(
                        cur_gmme_id.index() );
                auto& added_nodes = map_gmsh_ref_to_added_nodes[gmsh_entity];

                // Add corresponding edge
                index_t cell_id =
                    region_builder->create_cells( 1, CellType::TETRAHEDRON );

                // Add corresponding nodes
                std::vector< index_t > tetra_gmme_ids( 4 );
                for( auto element_vertex_id : range( 4 ) )
                {
                    index_t gmsh_node_id{ line_parser.field_as_uint(
                        3 + nb_tags + element_vertex_id ) };
                    tetra_gmme_ids[element_vertex_id] =
                        find( added_nodes, gmsh_node_id );
                    if( tetra_gmme_ids[element_vertex_id] == NO_ID )
                    {
                        added_nodes.push_back( gmsh_node_id );
                        tetra_gmme_ids[element_vertex_id] =
                            region_builder->create_vertex(
                                nodes[gmsh_node_id - gmsh_offset] );
                        region_builder->set_cell_vertex(
                            { cell_id, element_vertex_id },
                            tetra_gmme_ids[element_vertex_id] );
                    }
                    ringmesh_assert(
                        tetra_gmme_ids[element_vertex_id] != NO_ID );
                }
            }
            else if( cur_element_type == GMSHElementType::PRISM
                     || cur_element_type == GMSHElementType::PYRAMIDE
                     || cur_element_type == GMSHElementType::HEXAHEDRON )
            {
                throw RINGMeshException( "I/O",
                    "Import of GeoModels form GMSH with "
                    "non simplicial meshes is not yet implemented." );
            }
        }
    }

    void read_elements( GEO::LineInput& line_parser,
        const std::string& filename,
        GeoModel2D& geomodel,
        const std::vector< vec2 >& nodes )
    {
        ringmesh_unused( filename );
        read_elements_base( line_parser, geomodel, nodes );
    }

    void read_elements( GEO::LineInput& line_parser,
        const std::string& filename,
        GeoModel3D& geomodel,
        const std::vector< vec3 >& nodes )
    {
        read_elements_base( line_parser, geomodel, nodes );
        read_volume_elements( filename, geomodel, nodes );
    }

    void add_surface_line_relationship( GeoModel2D& geomodel,
        index_t incident_surface_id,
        index_t boundary_line_id )
    {
        GeoModelBuilder2D geomodel_builder( geomodel );
        // @todo Set the right side
        geomodel_builder.topology.add_surface_line_boundary_relation(
            incident_surface_id, boundary_line_id, false );
    }

    void add_surface_line_relationship( GeoModel3D& geomodel,
        index_t incident_surface_id,
        index_t boundary_line_id )
    {
        GeoModelBuilder3D geomodel_builder( geomodel );
        geomodel_builder.topology.add_surface_line_boundary_relation(
            incident_surface_id, boundary_line_id );
    }

    template < index_t DIMENSION >
    void complete_topology( GeoModel< DIMENSION >& geomodel )
    {
        // Remove duplicated lines
        // @todo

        // Create corners and set Corner-Line relationships
        GeoModelBuilder< DIMENSION > geomodel_builder( geomodel );
        for( auto& line : geomodel.lines() )
        {
            auto first_boundary = line.vertex( 0 );
            auto last_boundary = line.vertex( line.nb_vertices() - 1 );
            auto first_boundary_id =
                geomodel_builder.topology.find_or_create_corner(
                    first_boundary );
            auto last_boundary_id =
                geomodel_builder.topology.find_or_create_corner(
                    last_boundary );
            geomodel_builder.topology.add_line_corner_boundary_relation(
                line.index(), first_boundary_id.index() );
            geomodel_builder.topology.add_line_corner_boundary_relation(
                line.index(), last_boundary_id.index() );
        }

        // Set Line-Surface relationships
        for( auto& line : geomodel.lines() )
        {
            std::vector< index_t > added_incident_surfaces;
            added_incident_surfaces.reserve( 2 );
            if( line.nb_vertices() > 2 )
            {
                auto line_internal_vertex_global_id =
                    geomodel.mesh.vertices.geomodel_vertex_id( line.gmme(), 1 );
                auto gme_surface_vertices =
                    geomodel.mesh.vertices.gme_type_vertices(
                        surface_type_name_static(),
                        line_internal_vertex_global_id );
                for( auto gme_vertex : gme_surface_vertices )
                {
                    auto incident_surface_id = gme_vertex.gmme.index();
                    if( !contains(
                            added_incident_surfaces, incident_surface_id ) )
                    {
                        add_surface_line_relationship(
                            geomodel, incident_surface_id, line.index() );
                        added_incident_surfaces.push_back(
                            incident_surface_id );
                    }
                }
            }
            else
            {
                // Line is composed of only 2 vertices
                // Search for potential incident surfaces looking for those
                // containing the two corners
                std::set< index_t > potential_incident_surfaces;
                auto first_corner_global_id =
                    geomodel.mesh.vertices.geomodel_vertex_id( line.gmme(), 0 );
                auto first_corner_gme_surface_vertices =
                    geomodel.mesh.vertices.gme_type_vertices(
                        surface_type_name_static(), first_corner_global_id );
                auto second_corner_global_id =
                    geomodel.mesh.vertices.geomodel_vertex_id( line.gmme(), 1 );
                auto second_corner_gme_surface_vertices =
                    geomodel.mesh.vertices.gme_type_vertices(
                        surface_type_name_static(), second_corner_global_id );
                std::set< index_t > first_corner_incident_surfaces;
                std::set< index_t > second_corner_incident_surfaces;
                for( auto gme_vertex : first_corner_gme_surface_vertices )
                {
                    first_corner_incident_surfaces.insert(
                        gme_vertex.gmme.index() );
                }
                for( auto gme_vertex : second_corner_gme_surface_vertices )
                {
                    second_corner_incident_surfaces.insert(
                        gme_vertex.gmme.index() );
                }
                for( auto surface_id : first_corner_incident_surfaces )
                {
                    if( std::find( second_corner_incident_surfaces.begin(),
                            second_corner_incident_surfaces.end(), surface_id )
                        != second_corner_incident_surfaces.end() )
                    {
                        potential_incident_surfaces.insert( surface_id );
                    }
                }
                // Add relationships
                if( potential_incident_surfaces.size() <= 2 )
                {
                    for( auto incident_surface_id :
                        potential_incident_surfaces )
                    {
                        add_surface_line_relationship(
                            geomodel, incident_surface_id, line.index() );
                    }
                }
                else
                {
                    DEBUG( "problem TODO" );
                }
            }
        }

        // Connect surface polygons
        for( auto& surface : geomodel.surfaces() )
        {
            geomodel_builder.geometry.compute_surface_adjacencies(
                surface.index() );
        }
    }

    template < index_t DIMENSION >
    void complete_geometry_and_topology( GeoModel< DIMENSION >& geomodel )
    {
        complete_topology( geomodel );
    }

    template < index_t DIMENSION >
    class MSHInputHandler final : public GeoModelInputHandler< DIMENSION >
    {
    public:
        void load(
            const std::string& filename, GeoModel< DIMENSION >& geomodel ) final
        {
            GEO::LineInput file_line{ filename };
            GeoModelBuilder< DIMENSION > geomodel_builder( geomodel );
            geomodel_builder.info.set_geomodel_name( "Unnamed" );
            std::vector< vecn< DIMENSION > > nodes =
                read_nodes< DIMENSION >( file_line );
            read_elements( file_line, filename, geomodel, nodes );
            complete_geometry_and_topology( geomodel );
        }

        index_t dimension( const std::string& filename ) const final
        {
            index_t dimension{ 2 };
            bool truly_2d = all_third_coord_null( filename );
            if( !truly_2d )
            {
                dimension = 3;
            }
            DEBUG( dimension );
            return dimension;
        }
    };
    ALIAS_2D_AND_3D( MSHInputHandler );

    /*!
     * @brief Export GeoModel3D for the GMSH format 2.2 which is
     * described here:
     * http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format
     * NB : Mesh entities are also exported
     */
    template < index_t DIMENSION >
    class MSHOutputHandler final : public GeoModelOutputHandler< DIMENSION >
    {
    public:
        void save( const GeoModel< DIMENSION >& geomodel,
            const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

            out << "$MeshFormat" << EOL;
            out << "2.2 0 8" << EOL;
            out << "$EndMeshFormat" << EOL;

            out << "$Nodes" << EOL;
            out << geomodel.mesh.vertices.nb() << EOL;
            for( auto v : range( geomodel.mesh.vertices.nb() ) )
            {
                out << v + gmsh_offset << SPACE
                    << vertex_in_3d( geomodel.mesh.vertices.vertex( v ) )
                    << EOL;
            }
            out << "$EndNodes" << EOL;

            out << "$Elements" << EOL;
            out << count_elements( geomodel ) << EOL;
            const auto& gmme_types =
                geomodel.entity_type_manager()
                    .mesh_entity_manager.mesh_entity_types();

            index_t element_index{ 1 };
            for( auto gmme_type_index :
                range( geomodel.entity_type_manager()
                           .mesh_entity_manager.nb_mesh_entity_types() ) )
            {
                MeshEntityType cur_mesh_entity_type{
                    gmme_types[gmme_type_index]
                };
                for( auto index_of_gmme_of_the_current_type :
                    range( geomodel.nb_mesh_entities( cur_mesh_entity_type ) ) )
                {
                    gmme_id cur_gmme_id = gmme_id( cur_mesh_entity_type,
                        index_of_gmme_of_the_current_type );
                    const GeoModelMeshEntity< DIMENSION >& cur_gmme =
                        geomodel.mesh_entity( cur_gmme_id );
                    for( auto elem_in_cur_gmme :
                        range( cur_gmme.nb_mesh_elements() ) )
                    {
                        index_t nb_vertices_in_cur_element =
                            cur_gmme.nb_mesh_element_vertices(
                                elem_in_cur_gmme );
                        if( cur_gmme_id.type() == surface_type_name_static()
                            && nb_vertices_in_cur_element > 4 )
                        {
                            throw RINGMeshException( "I/O",
                                "Cannot export into GMSH file format a "
                                "GeoModel "
                                "with unclassified polygons (not a triangle "
                                "or a quad)." );
                        }
                        if( cur_gmme_id.type() == region_type_name_static()
                            && nb_vertices_in_cur_element > 8 )
                        {
                            throw RINGMeshException( "I/O", "Cannot export "
                                                            "into GMSH file "
                                                            "format a GeoModel "
                                                            "with polyhedra." );
                        }
                        index_t gmsh_element_type = find_gmsh_element_type(
                            nb_vertices_in_cur_element, gmme_type_index );
                        out << element_index++ << SPACE << gmsh_element_type
                            << SPACE << nb_of_tags << SPACE << physical_id
                            << SPACE
                            << index_of_gmme_of_the_current_type + gmsh_offset
                            << SPACE;
                        for( auto v_index_in_cur_element :
                            range( nb_vertices_in_cur_element ) )
                        {
                            out << geomodel.mesh.vertices.geomodel_vertex_id(
                                       cur_gmme_id,
                                       cur_gmme.mesh_element_vertex_index(
                                           ElementLocalVertex( elem_in_cur_gmme,
                                               find_gmsh_element_local_vertex_id(
                                                   nb_vertices_in_cur_element,
                                                   gmme_type_index,
                                                   v_index_in_cur_element ) ) ) )
                                       + gmsh_offset
                                << SPACE;
                        }
                        out << EOL;
                    }
                }
            }
            out << "$EndElements" << EOL;
            out << std::flush;
        }
    };
    ALIAS_2D_AND_3D( MSHOutputHandler );
}
