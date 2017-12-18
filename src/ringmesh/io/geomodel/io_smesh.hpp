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
    /*!
     * @brief Save the geomodel in smesh format
     * @details No attributes and no boundary marker are transferred
     */
    class SMESHIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            if( out.bad() )
            {
                Logger::err(
                    "I/O", "Error when opening the file: ", filename.c_str() );
                return;
            }
            out.precision( 16 );

            /// 1. Write the unique vertices
            out << "# Node list" << EOL;
            out << "# node count, 3 dim, no attribute, no boundary marker"
                << EOL;
            out << geomodel.mesh.vertices.nb() << " 3 0 0" << EOL;
            out << "# node index, node coordinates " << EOL;
            for( auto p : range( geomodel.mesh.vertices.nb() ) )
            {
                const vec3& V = geomodel.mesh.vertices.vertex( p );
                out << p << " "
                    << " " << V.x << " " << V.y << " " << V.z << EOL;
            }

            /// 2. Write the triangles
            out << "# Part 2 - facet list" << EOL;
            out << "# facet count, no boundary marker" << EOL;
            out << nb_polygons( geomodel ) << "  0 " << EOL;

            for( const auto& surface : geomodel.surfaces() )
            {
                for( auto p : range( surface.nb_mesh_elements() ) )
                {
                    out << surface.nb_mesh_element_vertices( p ) << " ";
                    for( auto v :
                        range( surface.nb_mesh_element_vertices( p ) ) )
                    {
                        out << geomodel.mesh.vertices.geomodel_vertex_id(
                                   surface.gmme(), ElementLocalVertex( p, v ) )
                            << " ";
                    }
                    out << EOL;
                }
            }

            // Do not forget the stupid zeros at the end of the file
            out << EOL << "0" << EOL << "0" << EOL;
            out << std::flush;
        }
    };
}
