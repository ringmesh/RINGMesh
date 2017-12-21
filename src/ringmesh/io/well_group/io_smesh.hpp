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
    void merge_colocated_vertices( double epsilon, LineMesh3D& mesh )
    {
        std::vector< index_t > old2new;
        index_t nb_colocated = NO_ID;
        std::tie( nb_colocated, old2new ) =
            mesh.vertex_nn_search().get_colocated_index_mapping( epsilon );
        if( nb_colocated > 0 )
        {
            std::unique_ptr< LineMeshBuilder3D > builder =
                LineMeshBuilder3D::create_builder( mesh );
            for( auto e : range( mesh.nb_edges() ) )
            {
                for( auto i : range( 2 ) )
                {
                    index_t v = mesh.edge_vertex( ElementLocalVertex( e, i ) );
                    builder->set_edge_vertex(
                        EdgeLocalVertex( e, i ), old2new[v] );
                }
            }
            std::vector< bool > delete_vertices( mesh.nb_vertices(), false );
            for( auto v : range( mesh.nb_vertices() ) )
            {
                if( old2new[v] != v )
                {
                    delete_vertices[v] = true;
                }
            }
            builder->delete_vertices( delete_vertices );
        }
    }

    class SmeshIOHandler final : public WellGroupIOHandler
    {
    public:
        void load( const std::string& filename, WellGroup3D& wells ) final
        {
            GEO::LineInput in( filename );
            if( !in.OK() )
            {
                throw RINGMeshException( "I/O", "Could not open file" );
            }

            auto mesh = LineMesh3D::create_mesh();
            std::unique_ptr< LineMeshBuilder3D > builder =
                LineMeshBuilder3D::create_builder( *mesh );
            std::string name = GEO::FileSystem::base_name( filename );

            bool is_first_part = true;

            while( !in.eof() )
            {
                in.get_line();
                in.get_fields();
                if( in.nb_fields() == 0 )
                    continue;
                if( GEO::String::string_starts_with( in.field( 0 ), "#" ) )
                {
                    continue;
                }
                if( is_first_part )
                {
                    index_t nb_vertices = in.field_as_uint( 0 );
                    builder->create_vertices( nb_vertices );
                    Box3D box;

                    for( auto v : range( nb_vertices ) )
                    {
                        do
                        {
                            in.get_line();
                            in.get_fields();
                        } while( in.nb_fields() == 0 );
                        vec3 point;
                        point[0] = in.field_as_double( 1 );
                        point[1] = in.field_as_double( 2 );
                        point[2] = in.field_as_double( 3 );
                        builder->set_vertex( v, point );
                        box.add_point( point );
                    }
                    is_first_part = false;
                }
                else
                {
                    index_t nb_edges = in.field_as_uint( 0 );
                    builder->create_edges( nb_edges );
                    for( auto e : range( nb_edges ) )
                    {
                        do
                        {
                            in.get_line();
                            in.get_fields();
                        } while( in.nb_fields() == 0 );
                        builder->set_edge_vertex(
                            EdgeLocalVertex( e, 0 ), in.field_as_uint( 1 ) );
                        builder->set_edge_vertex(
                            EdgeLocalVertex( e, 1 ), in.field_as_uint( 2 ) );
                    }
                    merge_colocated_vertices(
                        wells.geomodel()->epsilon(), *mesh );
                    wells.add_well( *mesh, name );
                    break;
                }
            }
        }
        void save( const WellGroup3D& wells, const std::string& filename ) final
        {
            ringmesh_unused( wells );
            ringmesh_unused( filename );
            throw RINGMeshException(
                "I/O", "Saving of a WellGroup from Smesh not implemented yet" );
        }
    };
}
