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
    class WLIOHandler final : public WellGroupIOHandler
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
            auto builder = LineMeshBuilder3D::create_builder( *mesh );
            std::string name;
            double z_sign = 1.0;
            vec3 vertex_ref;

            while( !in.eof() )
            {
                in.get_line();
                in.get_fields();
                if( in.nb_fields() == 0 )
                    continue;
                if( in.field_matches( 0, "name:" ) )
                {
                    name = in.field( 1 );
                }
                else if( in.field_matches( 0, "ZPOSITIVE" ) )
                {
                    if( in.field_matches( 1, "Depth" ) )
                    {
                        z_sign = -1.0;
                    }
                }
                else if( in.field_matches( 0, "WREF" ) )
                {
                    vertex_ref[0] = in.field_as_double( 1 );
                    vertex_ref[1] = in.field_as_double( 2 );
                    vertex_ref[2] = z_sign * in.field_as_double( 3 );
                    builder->create_vertex( vertex_ref );
                }
                else if( in.field_matches( 0, "PATH" ) )
                {
                    if( in.field_as_double( 1 ) == 0. )
                        continue;
                    vec3 vertex;
                    vertex[2] = z_sign * in.field_as_double( 2 );
                    vertex[0] = in.field_as_double( 3 ) + vertex_ref[0];
                    vertex[1] = in.field_as_double( 4 ) + vertex_ref[1];
                    index_t id = builder->create_vertex( vertex );
                    builder->create_edge( id - 1, id );
                }
                else if( in.field_matches( 0, "END" ) )
                {
                    wells.add_well( *mesh, name );
                    mesh = LineMesh3D::create_mesh();
                    builder = LineMeshBuilder3D::create_builder( *mesh );
                }
            }
        }
        void save( const WellGroup3D& wells, const std::string& filename ) final
        {
            ringmesh_unused( wells );
            ringmesh_unused( filename );
            throw RINGMeshException(
                "I/O", "Saving of a WellGroup from Gocad not implemented yet" );
        }
    };
}
