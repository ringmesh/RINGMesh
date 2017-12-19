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
    struct RINGMesh2CSMP
    {
        index_t entity_type;
        index_t nb_vertices;
        index_t vertices[8];
        index_t nb_polygons;
        index_t polygon[6];
    };

    static RINGMesh2CSMP tet_descriptor = {
        4, // type
        4, // nb vertices
        { 0, 1, 2, 3 }, // vertices
        4, // nb polygons
        { 0, 1, 2, 3 } // polygons
    };

    static RINGMesh2CSMP hex_descriptor = {
        6, // type
        8, // nb vertices
        { 0, 4, 5, 1, 2, 6, 7, 3 }, // vertices
        6, // nb polygons
        { 2, 0, 5, 1, 4, 3 } // polygons
    };

    static RINGMesh2CSMP prism_descriptor = {
        12, // type
        6, // nb vertices
        { 0, 1, 2, 3, 4, 5 }, // vertices
        5, // nb polygons
        { 0, 2, 4, 3, 1 } // polygons
    };

    static RINGMesh2CSMP pyramid_descriptor = {
        18, // type
        5, // nb vertices
        { 0, 1, 2, 3, 4 }, // vertices
        5, // nb polygons
        { 1, 4, 3, 2, 0 } // polygons
    };

    static RINGMesh2CSMP* cell_type_to_cell_descriptor[4] = { &tet_descriptor,
        &hex_descriptor, &prism_descriptor, &pyramid_descriptor };

    class CSMPIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        CSMPIOHandler()
        {
            clear();
        }

        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            initialize( geomodel );

            std::string directory = GEO::FileSystem::dir_name( filename );
            std::string file = GEO::FileSystem::base_name( filename );

            std::ostringstream oss_ascii;
            oss_ascii << directory << "/" << file << ".asc";
            std::ofstream ascii( oss_ascii.str().c_str() );
            ascii.precision( 16 );

            ascii << geomodel.name() << EOL;
            ascii << "Model generated from RINGMesh" << EOL;

            std::ostringstream oss_data;
            oss_data << directory << "/" << file << ".dat";
            std::ofstream data( oss_data.str().c_str() );
            data.precision( 16 );

            std::ostringstream oss_regions;
            oss_regions << directory << "/" << file << "-regions.txt";
            std::ofstream regions( oss_regions.str().c_str() );
            regions << "'" << oss_regions.str() << EOL;
            regions << "no properties" << EOL;

            const GeoModelMesh3D& mesh = geomodel.mesh;
            const GeoModelMeshPolygons3D& polygons = mesh.polygons;
            index_t count = 0;
            // Conversion from (X,Y,Z) to (X,Z,-Y)
            signed_index_t conversion_sign[3] = { 1, 1, -1 };
            index_t conversion_axis[3] = { 0, 2, 1 };
            data << mesh.vertices.nb() << " # PX, PY, PZ" << EOL;
            for( auto dim : range( 3 ) )
            {
                for( auto v : range( mesh.vertices.nb() ) )
                {
                    data << " "
                         << conversion_sign[dim]
                                * mesh.vertices.vertex(
                                      v )[conversion_axis[dim]];
                    new_line( count, 5, data );
                }
                reset_line( count, data );
            }
            reset_line( count, data );

            index_t nb_families = 0;
            index_t nb_interfaces = geomodel.nb_geological_entities(
                Interface3D::type_name_static() );
            std::vector< index_t > nb_triangle_interface( nb_interfaces, 0 );
            std::vector< index_t > nb_quad_interface( nb_interfaces, 0 );
            for( auto i : range( nb_interfaces ) )
            {
                const GeoModelGeologicalEntity3D& interf =
                    geomodel.geological_entity(
                        Interface3D::type_name_static(), i );
                for( auto s : range( interf.nb_children() ) )
                {
                    index_t s_id = interf.child_gmme( s ).index();
                    nb_triangle_interface[i] += polygons.nb_triangle( s_id );
                    nb_quad_interface[i] += polygons.nb_quad( s_id );
                }
                if( nb_triangle_interface[i] > 0 )
                    nb_families++;
                if( nb_quad_interface[i] > 0 )
                    nb_families++;
            }
            for( auto r : range( geomodel.nb_regions() ) )
            {
                if( mesh.cells.nb_tet( r ) > 0 )
                    nb_families++;
                if( mesh.cells.nb_pyramid( r ) > 0 )
                    nb_families++;
                if( mesh.cells.nb_prism( r ) > 0 )
                    nb_families++;
                if( mesh.cells.nb_hex( r ) > 0 )
                    nb_families++;
            }
            if( geomodel.wells() )
                nb_families += geomodel.wells()->nb_wells();

            ascii << nb_families << " # Number of families" << EOL;
            ascii << "# Object name" << TAB << "Entity type" << TAB
                  << "Material-ID" << TAB << "Number of entities" << EOL;
            for( const auto& region : geomodel.regions() )
            {
                regions << region.name() << EOL;
                std::string entity_type[4] = { "TETRA_4", "HEXA_8", "PENTA_6",
                    "PYRA_5" };
                for( auto type :
                    range( to_underlying_type( CellType::TETRAHEDRON ),
                        to_underlying_type( CellType::UNCLASSIFIED ) ) )
                {
                    CellType T = static_cast< CellType >( type );
                    if( mesh.cells.nb_cells( region.index(), T ) > 0 )
                    {
                        ascii << region.name() << TAB << entity_type[type]
                              << TAB << 0 << TAB
                              << mesh.cells.nb_cells( region.index(), T )
                              << EOL;
                    }
                }
            }
            for( auto i : range( nb_interfaces ) )
            {
                regions << interface_csmp_name( i, geomodel ) << EOL;
                if( nb_triangle_interface[i] > 0 )
                {
                    ascii << interface_csmp_name( i, geomodel ) << TAB
                          << "TRI_3" << TAB << 0 << TAB
                          << nb_triangle_interface[i] << EOL;
                }
                if( nb_quad_interface[i] > 0 )
                {
                    ascii << interface_csmp_name( i, geomodel ) << TAB
                          << "QUAD_4" << TAB << 0 << TAB << nb_quad_interface[i]
                          << EOL;
                }
            }
            if( geomodel.wells() )
            {
                for( auto w : range( geomodel.wells()->nb_wells() ) )
                {
                    const Well3D& well = geomodel.wells()->well( w );
                    regions << well.name() << EOL;
                    ascii << well.name() << TAB << "BAR_2" << TAB << 0 << TAB
                          << well.nb_edges() << EOL;
                }
            }

            data << "# PBFLAGS" << EOL;
            for( auto p : range( mesh.vertices.nb() ) )
            {
                data << " " << std::setw( 3 ) << point_boundary( p );
                new_line( count, 20, data );
            }
            reset_line( count, data );

            data << "# PBVALS" << EOL;
            for( auto p : range( mesh.vertices.nb() ) )
            {
                ringmesh_unused( p );
                data << " " << std::setw( 3 ) << 0;
                new_line( count, 20, data );
            }
            reset_line( count, data );

            index_t nb_total_entities = mesh.cells.nb_cells()
                                        + polygons.nb_polygons()
                                        + mesh.wells.nb_edges();
            data << nb_total_entities << " # PELEMENT" << EOL;
            for( auto r : range( geomodel.nb_regions() ) )
            {
                index_t entity_type[4] = { 4, 6, 12, 18 };
                for( auto type :
                    range( to_underlying_type( CellType::TETRAHEDRON ),
                        to_underlying_type( CellType::UNCLASSIFIED ) ) )
                {
                    CellType T = static_cast< CellType >( type );
                    for( auto el : range( mesh.cells.nb_cells( r, T ) ) )
                    {
                        ringmesh_unused( el );
                        data << " " << std::setw( 3 ) << entity_type[type];
                        new_line( count, 20, data );
                    }
                }
            }
            for( auto i : range( nb_interfaces ) )
            {
                for( auto el : range( nb_triangle_interface[i] ) )
                {
                    ringmesh_unused( el );
                    data << " " << std::setw( 3 ) << 8;
                    new_line( count, 20, data );
                }
                for( auto el : range( nb_quad_interface[i] ) )
                {
                    ringmesh_unused( el );
                    data << " " << std::setw( 3 ) << 14;
                    new_line( count, 20, data );
                }
            }
            if( geomodel.wells() )
            {
                for( auto w : range( geomodel.wells()->nb_wells() ) )
                {
                    const Well3D& well = geomodel.wells()->well( w );
                    for( auto e : range( well.nb_edges() ) )
                    {
                        ringmesh_unused( e );
                        data << " " << std::setw( 3 ) << 2;
                        new_line( count, 20, data );
                    }
                }
            }
            reset_line( count, data );

            ascii << "# now the entities which make up each object are listed "
                     "in sequence"
                  << EOL;
            index_t cur_cell = 0;
            for( auto r : range( geomodel.nb_regions() ) )
            {
                const RINGMesh::GeoModelEntity3D& region = geomodel.region( r );
                std::string entity_type[4] = { "TETRA_4", "HEXA_8", "PENTA_6",
                    "PYRA_5" };
                for( auto type :
                    range( to_underlying_type( CellType::TETRAHEDRON ),
                        to_underlying_type( CellType::UNCLASSIFIED ) ) )
                {
                    CellType T = static_cast< CellType >( type );
                    if( mesh.cells.nb_cells( r, T ) > 0 )
                    {
                        ascii << region.name() << " " << entity_type[type]
                              << " " << mesh.cells.nb_cells( r, T ) << EOL;
                        for( auto el : range( mesh.cells.nb_cells( r, T ) ) )
                        {
                            ringmesh_unused( el );
                            ascii << cur_cell++ << " ";
                            new_line( count, 10, ascii );
                        }
                        reset_line( count, ascii );
                    }
                }
            }
            for( auto i : range( nb_interfaces ) )
            {
                if( nb_triangle_interface[i] > 0 )
                {
                    ascii << interface_csmp_name( i, geomodel ) << " "
                          << "TRI_3"
                          << " " << nb_triangle_interface[i] << EOL;
                    for( auto el : range( nb_triangle_interface[i] ) )
                    {
                        ringmesh_unused( el );
                        ascii << cur_cell++ << " ";
                        new_line( count, 10, ascii );
                    }
                    reset_line( count, ascii );
                }
                if( nb_quad_interface[i] > 0 )
                {
                    ascii << interface_csmp_name( i, geomodel ) << " "
                          << "QUAD_4"
                          << " " << nb_quad_interface[i] << EOL;
                    for( auto el : range( nb_quad_interface[i] ) )
                    {
                        ringmesh_unused( el );
                        ascii << cur_cell++ << " ";
                        new_line( count, 10, ascii );
                    }
                    reset_line( count, ascii );
                }
            }
            if( geomodel.wells() )
            {
                for( auto w : range( geomodel.wells()->nb_wells() ) )
                {
                    const Well3D& well = geomodel.wells()->well( w );
                    ascii << well.name() << " "
                          << "BAR_2"
                          << " " << well.nb_edges() << EOL;
                    for( auto e : range( well.nb_edges() ) )
                    {
                        ringmesh_unused( e );
                        ascii << cur_cell++ << " ";
                        new_line( count, 10, ascii );
                    }
                    reset_line( count, ascii );
                }
            }

            index_t nb_plist =
                3 * polygons.nb_triangle() + 4 * polygons.nb_quad()
                + 4 * mesh.cells.nb_tet() + 5 * mesh.cells.nb_pyramid()
                + 6 * mesh.cells.nb_prism() + 8 * mesh.cells.nb_hex()
                + 2 * mesh.wells.nb_edges();
            data << nb_plist << " # PLIST" << EOL;
            for( auto r : range( geomodel.nb_regions() ) )
            {
                for( auto type :
                    range( to_underlying_type( CellType::TETRAHEDRON ),
                        to_underlying_type( CellType::UNCLASSIFIED ) ) )
                {
                    CellType T = static_cast< CellType >( type );
                    RINGMesh2CSMP& descriptor =
                        *cell_type_to_cell_descriptor[type];
                    for( auto el : range( mesh.cells.nb_cells( r, T ) ) )
                    {
                        index_t cell = mesh.cells.cell( r, el, T );
                        for( auto p : range( descriptor.nb_vertices ) )
                        {
                            index_t csmp_p = descriptor.vertices[p];
                            index_t vertex_id = mesh.cells.vertex(
                                ElementLocalVertex( cell, csmp_p ) );
                            data << " " << std::setw( 7 ) << vertex_id;
                            new_line( count, 10, data );
                        }
                    }
                }
            }
            for( auto i : range( nb_interfaces ) )
            {
                const GeoModelGeologicalEntity3D& interf =
                    geomodel.geological_entity(
                        Interface3D::type_name_static(), i );
                for( auto s : range( interf.nb_children() ) )
                {
                    index_t s_id = interf.child_gmme( s ).index();
                    for( auto el : range( polygons.nb_triangle( s_id ) ) )
                    {
                        index_t tri = polygons.triangle( s_id, el );
                        for( auto p : range( polygons.nb_vertices( tri ) ) )
                        {
                            index_t vertex_id =
                                polygons.vertex( ElementLocalVertex( tri, p ) );
                            data << " " << std::setw( 7 ) << vertex_id;
                            new_line( count, 10, data );
                        }
                    }
                    for( auto el : range( polygons.nb_quad( s_id ) ) )
                    {
                        index_t quad = polygons.quad( s_id, el );
                        for( auto p : range( polygons.nb_vertices( quad ) ) )
                        {
                            index_t vertex_id = polygons.vertex(
                                ElementLocalVertex( quad, p ) );
                            data << " " << std::setw( 7 ) << vertex_id;
                            new_line( count, 10, data );
                        }
                    }
                }
            }
            for( auto w : range( mesh.wells.nb_wells() ) )
            {
                for( auto e : range( mesh.wells.nb_edges( w ) ) )
                {
                    for( auto v : range( 2 ) )
                    {
                        index_t vertex_id = mesh.wells.vertex( w, e, v );
                        data << " " << std::setw( 7 ) << vertex_id;
                        new_line( count, 10, data );
                    }
                }
            }
            reset_line( count, data );

            index_t nb_polygons =
                3 * polygons.nb_triangle() + 4 * polygons.nb_quad()
                + 4 * mesh.cells.nb_tet() + 5 * mesh.cells.nb_pyramid()
                + 5 * mesh.cells.nb_prism() + 6 * mesh.cells.nb_hex()
                + 2 * mesh.wells.nb_edges();
            data << nb_polygons << " # PFVERTS" << EOL;
            for( auto r : range( geomodel.nb_regions() ) )
            {
                for( auto type :
                    range( to_underlying_type( CellType::TETRAHEDRON ),
                        to_underlying_type( CellType::UNCLASSIFIED ) ) )
                {
                    CellType T = static_cast< CellType >( type );
                    RINGMesh2CSMP& descriptor =
                        *cell_type_to_cell_descriptor[type];
                    for( auto el : range( mesh.cells.nb_cells( r, T ) ) )
                    {
                        index_t cell = mesh.cells.cell( r, el );
                        for( auto p : range( descriptor.nb_polygons ) )
                        {
                            index_t csmp_f = descriptor.polygon[p];
                            index_t adj = mesh.cells.adjacent( cell, csmp_f );
                            if( adj == NO_ID )
                            {
                                data << " " << std::setw( 7 ) << -28;
                            }
                            else
                            {
                                data << " " << std::setw( 7 ) << adj;
                            }
                            new_line( count, 10, data );
                        }
                    }
                }
            }
            for( auto i : range( nb_interfaces ) )
            {
                const GeoModelGeologicalEntity3D& interf =
                    geomodel.geological_entity(
                        Interface3D::type_name_static(), i );
                for( auto s : range( interf.nb_children() ) )
                {
                    index_t s_id = interf.child_gmme( s ).index();
                    for( auto el : range( polygons.nb_triangle( s_id ) ) )
                    {
                        index_t tri = polygons.triangle( s_id, el );
                        for( auto e : range( polygons.nb_vertices( tri ) ) )
                        {
                            index_t adj =
                                polygons.adjacent( PolygonLocalEdge( tri, e ) );
                            if( adj == NO_ID )
                            {
                                data << " " << std::setw( 7 ) << -28;
                            }
                            else
                            {
                                data << " " << std::setw( 7 ) << adj;
                            }
                            new_line( count, 10, data );
                        }
                    }
                    for( auto el : range( polygons.nb_quad( s_id ) ) )
                    {
                        index_t quad = polygons.quad( s_id, el );
                        for( auto e : range( polygons.nb_vertices( quad ) ) )
                        {
                            index_t adj = polygons.adjacent(
                                PolygonLocalEdge( quad, e ) );
                            if( adj == NO_ID )
                            {
                                data << " " << std::setw( 7 ) << -28;
                            }
                            else
                            {
                                data << " " << std::setw( 7 ) << adj;
                            }
                            new_line( count, 10, data );
                        }
                    }
                }
            }
            index_t edge_offset = polygons.nb() + mesh.cells.nb();
            index_t cur_edge = 0;
            for( auto w : range( mesh.wells.nb_wells() ) )
            {
                data << " " << std::setw( 7 ) << -28;
                new_line( count, 10, data );
                if( mesh.wells.nb_edges( w ) > 1 )
                {
                    data << " " << std::setw( 7 ) << edge_offset + cur_edge + 1;
                    cur_edge++;
                    new_line( count, 10, data );
                    for( index_t e = 1; e < mesh.wells.nb_edges( w ) - 1;
                         e++, cur_edge++ )
                    {
                        data << " " << std::setw( 7 )
                             << edge_offset + cur_edge - 1;
                        new_line( count, 10, data );
                        data << " " << std::setw( 7 )
                             << edge_offset + cur_edge + 1;
                        new_line( count, 10, data );
                    }
                    data << " " << std::setw( 7 ) << edge_offset + cur_edge - 1;
                    new_line( count, 10, data );
                }
                data << " " << std::setw( 7 ) << -28;
                cur_edge++;
                new_line( count, 10, data );
            }
            reset_line( count, data );

            data << nb_total_entities << " # PMATERIAL" << EOL;
            for( auto i : range( nb_total_entities ) )
            {
                ringmesh_unused( i );
                data << " " << std::setw( 3 ) << 0;
                new_line( count, 20, data );
            }

            ascii << std::flush;
            data << std::flush;
            regions << std::flush;
        }

    private:
        void new_line(
            index_t& count, index_t number_of_counts, std::ofstream& out ) const
        {
            count++;
            if( count == number_of_counts )
            {
                count = 0;
                out << EOL;
            }
        }
        void reset_line( index_t& count, std::ofstream& out ) const
        {
            if( count != 0 )
            {
                count = 0;
                out << EOL;
            }
        }
        void clear()
        {
            point_boundaries_.clear();
            box_model_ = false;
            back_ = NO_ID;
            top_ = NO_ID;
            front_ = NO_ID;
            bottom_ = NO_ID;
            left_ = NO_ID;
            right_ = NO_ID;
            corner_boundary_flags_.clear();
            edge_boundary_flags_.clear();
            surface_boundary_flags_.clear();
        }
        void initialize( const GeoModel3D& gm )
        {
            clear();

            const GeoModel3D& geomodel = gm;
            std::string cmsp_filename; //= GEO::CmdLine::get_arg( "out:csmp" );
            box_model_ = cmsp_filename != "";
            if( box_model_ )
            {
                GEO::LineInput parser( cmsp_filename );
                if( !parser.OK() )
                {
                    throw RINGMeshException(
                        "I/O", "Cannot open file: ", cmsp_filename );
                }
                parser.get_line();
                parser.get_fields();
                while( !parser.eof() )
                {
                    if( parser.nb_fields() == 0 )
                        continue;
                    if( parser.nb_fields() != 3 )
                        return;
                    std::string type = parser.field( 1 );
                    index_t interface_id = NO_ID;
                    if( type == "NAME" )
                    {
                        std::string name = parser.field( 2 );
                        GeologicalEntityType type =
                            Interface3D::type_name_static();
                        for( auto& cur_interface :
                            geomodel.geol_entities( type ) )
                        {
                            if( cur_interface.name() == name )
                            {
                                interface_id = cur_interface.index();
                                break;
                            }
                        }
                    }
                    else if( type == "ID" )
                    {
                        interface_id = parser.field_as_uint( 2 );
                    }
                    else
                    {
                        throw RINGMeshException(
                            "I/O", "Unknown type: ", type );
                    }

                    std::string keyword = parser.field( 0 );
                    if( keyword == "BACK" )
                    {
                        back_ = interface_id;
                    }
                    else if( keyword == "TOP" )
                    {
                        top_ = interface_id;
                    }
                    else if( keyword == "FRONT" )
                    {
                        front_ = interface_id;
                    }
                    else if( keyword == "BOTTOM" )
                    {
                        bottom_ = interface_id;
                    }
                    else if( keyword == "LEFT" )
                    {
                        left_ = interface_id;
                    }
                    else if( keyword == "RIGHT" )
                    {
                        right_ = interface_id;
                    }
                    else
                    {
                        throw RINGMeshException(
                            "I/O", "Unknown keyword: ", keyword );
                    }
                    parser.get_line();
                    parser.get_fields();
                }

                if( back_ == NO_ID || top_ == NO_ID || front_ == NO_ID
                    || bottom_ == NO_ID || left_ == NO_ID || right_ == NO_ID )
                {
                    throw RINGMeshException(
                        "I/O", "Missing box shape information" );
                }

                surface_boundary_flags_[back_] = -7;
                surface_boundary_flags_[top_] = -5;
                surface_boundary_flags_[front_] = -6;
                surface_boundary_flags_[bottom_] = -4;
                surface_boundary_flags_[left_] = -2;
                surface_boundary_flags_[right_] = -3;

                std::set< index_t > back_bottom;
                back_bottom.insert( back_ );
                back_bottom.insert( bottom_ );
                edge_boundary_flags_[back_bottom] = -16;
                std::set< index_t > back_right;
                back_right.insert( back_ );
                back_right.insert( right_ );
                edge_boundary_flags_[back_right] = -17;
                std::set< index_t > back_top;
                back_top.insert( back_ );
                back_top.insert( top_ );
                edge_boundary_flags_[back_top] = -18;
                std::set< index_t > back_left;
                back_left.insert( back_ );
                back_left.insert( left_ );
                edge_boundary_flags_[back_left] = -19;
                std::set< index_t > right_bottom;
                right_bottom.insert( right_ );
                right_bottom.insert( bottom_ );
                edge_boundary_flags_[right_bottom] = -20;
                std::set< index_t > right_top;
                right_top.insert( right_ );
                right_top.insert( top_ );
                edge_boundary_flags_[right_top] = -21;
                std::set< index_t > left_top;
                left_top.insert( left_ );
                left_top.insert( top_ );
                edge_boundary_flags_[left_top] = -22;
                std::set< index_t > left_bottom;
                left_bottom.insert( left_ );
                left_bottom.insert( bottom_ );
                edge_boundary_flags_[left_bottom] = -23;
                std::set< index_t > front_bottom;
                front_bottom.insert( front_ );
                front_bottom.insert( bottom_ );
                edge_boundary_flags_[front_bottom] = -24;
                std::set< index_t > front_right;
                front_right.insert( front_ );
                front_right.insert( right_ );
                edge_boundary_flags_[front_right] = -25;
                std::set< index_t > front_top;
                front_top.insert( front_ );
                front_top.insert( top_ );
                edge_boundary_flags_[front_top] = -26;
                std::set< index_t > front_left;
                front_left.insert( front_ );
                front_left.insert( left_ );
                edge_boundary_flags_[front_left] = -27;

                std::set< index_t > back_top_left;
                back_top_left.insert( back_ );
                back_top_left.insert( top_ );
                back_top_left.insert( left_ );
                corner_boundary_flags_[back_top_left] = -13;
                std::set< index_t > back_top_right;
                back_top_right.insert( back_ );
                back_top_right.insert( top_ );
                back_top_right.insert( right_ );
                corner_boundary_flags_[back_top_right] = -14;
                std::set< index_t > back_bottom_left;
                back_bottom_left.insert( back_ );
                back_bottom_left.insert( bottom_ );
                back_bottom_left.insert( left_ );
                corner_boundary_flags_[back_bottom_left] = -8;
                std::set< index_t > back_bottom_right;
                back_bottom_right.insert( back_ );
                back_bottom_right.insert( bottom_ );
                back_bottom_right.insert( right_ );
                corner_boundary_flags_[back_bottom_right] = -10;
                std::set< index_t > front_top_left;
                front_top_left.insert( front_ );
                front_top_left.insert( top_ );
                front_top_left.insert( left_ );
                corner_boundary_flags_[front_top_left] = -15;
                std::set< index_t > front_top_right;
                front_top_right.insert( front_ );
                front_top_right.insert( top_ );
                front_top_right.insert( right_ );
                corner_boundary_flags_[front_top_right] = -9;
                std::set< index_t > front_bottom_left;
                front_bottom_left.insert( front_ );
                front_bottom_left.insert( bottom_ );
                front_bottom_left.insert( left_ );
                corner_boundary_flags_[front_bottom_left] = -12;
                std::set< index_t > front_bottom_right;
                front_bottom_right.insert( front_ );
                front_bottom_right.insert( bottom_ );
                front_bottom_right.insert( right_ );
                corner_boundary_flags_[front_bottom_right] = -11;
            }

            point_boundaries_.resize( gm.mesh.vertices.nb() );
            for( const auto& surface : gm.surfaces() )
            {
                index_t interface_id = surface.parent_gmge( 0 ).index();
                for( auto p :
                    range( gm.mesh.polygons.nb_polygons( surface.index() ) ) )
                {
                    index_t p_id =
                        gm.mesh.polygons.polygon( surface.index(), p );
                    for( auto v :
                        range( gm.mesh.polygons.nb_vertices( p_id ) ) )
                    {
                        index_t vertex_id = gm.mesh.polygons.vertex(
                            ElementLocalVertex( p_id, v ) );
                        point_boundaries_[vertex_id].insert( interface_id );
                    }
                }
            }
        }
        std::string interface_csmp_name(
            index_t i, const GeoModel3D& geomodel ) const
        {
            if( box_model_ )
            {
                if( i == back_ )
                {
                    return "BACK";
                }
                else if( i == top_ )
                {
                    return "TOP";
                }
                else if( i == front_ )
                {
                    return "FRONT";
                }
                else if( i == bottom_ )
                {
                    return "BOTTOM";
                }
                else if( i == left_ )
                {
                    return "LEFT";
                }
                else if( i == right_ )
                {
                    return "RIGHT";
                }
            }
            return geomodel
                .geological_entity( Interface3D::type_name_static(), i )
                .name();
        }
        signed_index_t point_boundary( index_t p ) const
        {
            ringmesh_assert( p < point_boundaries_.size() );
            const std::set< unsigned int >& boundaries = point_boundaries_[p];
            if( box_model_ )
            {
                if( boundaries.size() == 1 )
                {
                    auto it =
                        surface_boundary_flags_.find( *boundaries.begin() );
                    ringmesh_assert( it != surface_boundary_flags_.end() );
                    return it->second;
                }
                else if( boundaries.size() == 2 )
                {
                    auto it = edge_boundary_flags_.find( boundaries );
                    ringmesh_assert( it != edge_boundary_flags_.end() );
                    return it->second;
                }
                else if( boundaries.size() == 3 )
                {
                    auto it = corner_boundary_flags_.find( boundaries );
                    ringmesh_assert( it != corner_boundary_flags_.end() );
                    return it->second;
                }
                else
                {
                    return 0;
                }
            }
            else
            {
                if( boundaries.empty() )
                    return 0;
                else
                    return -28;
            }
        }

    private:
        std::vector< std::set< index_t > > point_boundaries_;

        bool box_model_;
        index_t back_;
        index_t top_;
        index_t front_;
        index_t bottom_;
        index_t left_;
        index_t right_;

        std::map< std::set< index_t >, signed_index_t > corner_boundary_flags_;
        std::map< std::set< index_t >, signed_index_t > edge_boundary_flags_;
        std::map< index_t, signed_index_t > surface_boundary_flags_;
    };
}
