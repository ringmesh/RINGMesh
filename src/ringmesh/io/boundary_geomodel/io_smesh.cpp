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
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
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
     * @brief Save the geomodel in smesh format
     * @details No attributes and no boundary marker are transferred
     */
    class SMESHIOHandler final: public GeoModelIOHandler< 3 > {
    public:
        void load( const std::string& filename, GeoModel< 3 >& geomodel ) final
        {
            throw RINGMeshException( "I/O",
                "Geological model loading of a from UCD mesh not yet implemented" );
        }

        void save(
            const GeoModel< 3 >& geomodel,
            const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            if( out.bad() ) {
                Logger::err( "I/O", "Error when opening the file: ",
                    filename.c_str() );
                return;
            }
            out.precision( 16 );

            /// 1. Write the unique vertices
            out << "# Node list" << std::endl;
            out << "# node count, 3 dim, no attribute, no boundary marker"
                << std::endl;
            out << geomodel.mesh.vertices.nb() << " 3 0 0" << std::endl;
            out << "# node index, node coordinates " << std::endl;
            for( index_t p = 0; p < geomodel.mesh.vertices.nb(); p++ ) {
                const vec3& V = geomodel.mesh.vertices.vertex( p );
                out << p << " " << " " << V.x << " " << V.y << " " << V.z
                    << std::endl;
            }

            /// 2. Write the triangles
            out << "# Part 2 - facet list" << std::endl;
            out << "# facet count, no boundary marker" << std::endl;
            out << nb_polygons( geomodel ) << "  0 " << std::endl;

            for( index_t i = 0; i < geomodel.nb_surfaces(); ++i ) {
                const Surface< 3 >& S = geomodel.surface( i );
                for( index_t p = 0; p < S.nb_mesh_elements(); p++ ) {
                    out << S.nb_mesh_element_vertices( p ) << " ";
                    for( index_t v = 0; v < S.nb_mesh_element_vertices( p ); v++ ) {
                        out
                            << geomodel.mesh.vertices.geomodel_vertex_id( S.gmme(),
                                p, v ) << " ";
                    }
                    out << std::endl;
                }
            }

            // Do not forget the stupid zeros at the end of the file
            out << std::endl << "0" << std::endl << "0" << std::endl;
        }
    };

}
