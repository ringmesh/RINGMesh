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

    void save_header( const GeoModel< 3 >& geomodel, std::ostream& out )
    {
        out << "solid " << geomodel.name() << std::endl;
    }

    void save_footer( const GeoModel< 3 >& geomodel, std::ostream& out )
    {
        out << "endsolid" << geomodel.name() << std::endl;
    }

    void save_normal(
        const GeoModel< 3 >& geomodel,
        index_t triangle_id,
        std::ostream& out )
    {
        out << "facet normal " << geomodel.mesh.polygons.normal( triangle_id )
            << std::endl;
    }

    void begin_triangle( std::ostream& out )
    {
        out << "outer loop" << std::endl;
    }

    void end_triangle( std::ostream& out )
    {
        out << "endloop" << std::endl;
    }

    void save_triangle_vertex(
        const GeoModel< 3 >& geomodel,
        index_t triangle_id,
        index_t local_vertex_id,
        std::ostream& out )
    {
        out << "vertex "
            << geomodel.mesh.vertices.vertex(
                geomodel.mesh.polygons.vertex( triangle_id, local_vertex_id ) )
            << std::endl;
    }

    void save_triangle(
        const GeoModel< 3 >& geomodel,
        index_t triangle_id,
        std::ostream& out )
    {
        save_normal( geomodel, triangle_id, out );
        begin_triangle( out );
        for( index_t vertex = 0; vertex < 3; vertex++ ) {
            save_triangle_vertex( geomodel, triangle_id, vertex, out );
        }
        end_triangle( out );
    }

    void save_triangles( const GeoModel< 3 >& geomodel, std::ostream& out )
    {
        for( index_t triangle = 0; triangle < geomodel.mesh.polygons.nb_triangle();
            triangle++ ) {
            save_triangle( geomodel, triangle, out );
        }
    }

    void check_stl_validity( const GeoModel< 3 >& geomodel )
    {
        if( geomodel.mesh.polygons.nb() != geomodel.mesh.polygons.nb_triangle() ) {
            throw RINGMeshException( "I/O",
                "Geological model save in STL format support only triangles" );
        }
    }

    /*!
     * STL is an (old) file format used in CAD software
     * This is the ASCII export for this format
     * facet normal ni nj nk
     *   outer loop
     *      vertex v1x v1y v1z
     *      vertex v2x v2y v2z
     *      vertex v3x v3y v3z
     *   endloop
     * endfacet
     */
    class STLIOHandler final: public GeoModelIOHandler< 3 > {
    public:
        void load( const std::string& filename, GeoModel< 3 >& geomodel ) final
        {
            throw RINGMeshException( "I/O",
                "Geological model loading of a from STL mesh not yet implemented" );
        }

        void save(
            const GeoModel< 3 >& geomodel,
            const std::string& filename ) final
        {
            check_stl_validity( geomodel );
            std::ofstream out( filename.c_str() );
            out.precision( 17 );
            save_header( geomodel, out );
            save_triangles( geomodel, out );
            save_footer( geomodel, out );
        }
    };
}
