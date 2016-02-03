/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
 *
 *
 *
 *
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/io.h>
#include <ringmesh/well.h>

#include <geogram/basic/line_stream.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/logger.h>

namespace RINGMesh {
    static double read_double( GEO::LineInput& in, index_t field )
    {
        double result ;
        std::istringstream iss( in.field( field ) ) ;
        iss >> result >> std::ws ;
        return result ;
    }

    static std::string TAB = "\t" ;
    static std::string SPACE = " " ;

    /*!
     * This class implement the import of the well ASCII export of Gocad
     * It is a "Gocad Object", it is the file outputed by
     * File > Export > Gocad Object...
     * Please use the .wl extension to make it work well :)
     */
    class WLIOHandler: public WellGroupIOHandler {
    public:
        virtual bool load( const std::string& filename, WellGroup& wells )
        {
            GEO::LineInput in( filename ) ;
            if( !in.OK() ) {
                return false ;
            }

            GEO::Mesh mesh ;
            std::string name ;
            double z_sign = 1.0 ;
            double vertex_ref[3] ;

            while( !in.eof() ) {
                in.get_line() ;
                in.get_fields() ;
                if( in.nb_fields() == 0 ) continue ;
                if( in.field_matches( 0, "name:" ) ) {
                    name = in.field( 1 ) ;
                } else if( in.field_matches( 0, "ZPOSITIVE" ) ) {
                    if( in.field_matches( 1, "Depth" ) ) {
                        z_sign = -1.0 ;
                    }
                } else if( in.field_matches( 0, "WREF" ) ) {
                    vertex_ref[0] = read_double( in, 1 ) ;
                    vertex_ref[1] = read_double( in, 2 ) ;
                    vertex_ref[2] = z_sign * read_double( in, 3 ) ;
                    mesh.vertices.create_vertex( vertex_ref ) ;
                } else if( in.field_matches( 0, "PATH" ) ) {
                    if( read_double( in, 1 ) == 0. ) continue ;
                    double vertex[3] ;
                    vertex[2] = z_sign * read_double( in, 2 ) ;
                    vertex[0] = read_double( in, 3 ) + vertex_ref[0] ;
                    vertex[1] = read_double( in, 4 ) + vertex_ref[1] ;
                    index_t id = mesh.vertices.create_vertex( vertex ) ;
                    mesh.edges.create_edge( id - 1, id ) ;
                } else if( in.field_matches( 0, "END" ) ) {
                    wells.add_well( mesh, name ) ;
                    mesh.clear() ;
                }
            }

            return true ;
        }
        virtual bool save( const WellGroup& wells, const std::string& filename )
        {
            GEO::Logger::err( "I/O" )
                << "Saving of a WellGroup from Gocad not implemented yet"
                << std::endl ;
            return false ;
        }
    } ;

    /*!
     * This class implement the import of the well ASCII export of Gocad
     * It is NOT a "Gocad Object", it is the file outputed by
     * File > Export > Well > Well path to ASCII
     * Please use the .wg extension to make it work well :)
     */
    class WGIOHandler: public WellGroupIOHandler {
    public:
        virtual bool load( const std::string& filename, WellGroup& wells )
        {
            GEO::LineInput in( filename ) ;
            if( !in.OK() ) {
                return false ;
            }
            GEO::Mesh mesh ;
            in.get_line() ;
            in.get_fields() ;
            while( !in.eof() ) {
                in.get_line() ;
                in.get_fields() ;
                std::string well_name = in.field( 0 ) ;
                if( well_name != in.field( 0 ) ) {
                    create_well_segments(mesh) ;
                    wells.add_well( mesh, well_name ) ;
                    mesh.clear() ;
                }
                add_vertex_to_well_path(mesh,in) ;
            }

            return true ;
        }

        virtual bool save( const WellGroup& wells, const std::string& filename )
        {
            GEO::Logger::err( "I/O" )
                << "Saving of a WellGroup from Gocad not implemented yet"
                << std::endl ;
            return false ;
        }

    private:
        /*!
         * @brief create the well segments once all the vertices are set
         * @details this function has to be called only when
         * all the vertices are added to the \p mesh.
         * @param[in] mesh the mesh of the well
         */
        void create_well_segments(GEO::Mesh& mesh) {
            for( index_t v = 0; v < mesh.vertices.nb() - 1; v++ ) {
                mesh.edges.create_edge( v + 1, v ) ;
            }
        }

        /*!
         * @brief add a vertex to the well path for ASCII well export of Gocad
         * @param[in] mesh the mesh of the well
         * @param[in] in the input line matching containing the coordinate
         */
        void add_vertex_to_well_path(GEO::Mesh& mesh, GEO::LineInput& in) {
            double vertex[3] ;
            vertex[0] = read_double( in, 1 ) ;
            vertex[1] = read_double( in, 2 ) ;
            vertex[2] = read_double( in, 3 ) ;
            mesh.vertices.create_vertex( vertex ) ;
        }

    } ;

    /************************************************************************/

    /*!
     * Loads a WellGroup from a file
     * @param[in] filename the file to load
     * @param][out] wells the wells to fill
     * @return returns the success of the operation
     */
    bool well_load( const std::string& filename, WellGroup& wells )
    {
        GEO::Logger::out( "I/O" ) << "Loading file " << filename << "..."
            << std::endl ;

        WellGroupIOHandler_var handler = WellGroupIOHandler::get_handler(
            filename ) ;
        if( handler && handler->load( filename, wells ) ) {
            return true ;
        }

        GEO::Logger::err( "I/O" ) << "Could not load file: " << filename
            << std::endl ;
        return false ;
    }

    WellGroupIOHandler* WellGroupIOHandler::create( const std::string& format )
    {
        WellGroupIOHandler* handler = WellGroupIOHandlerFactory::create_object(
            format ) ;
        if( handler ) {
            return handler ;
        }

        GEO::Logger::err( "I/O" ) << "Unsupported file format: " << format
            << std::endl ;

        std::vector< std::string > names ;
        WellGroupIOHandlerFactory::list_creators( names ) ;
        GEO::Logger::out( "I/O" ) << "Currently supported file formats:" ;
        for( index_t i = 0; i < names.size(); i++ ) {
            GEO::Logger::out( "I/O" ) << " " << names[i] ;
        }
        GEO::Logger::out( "I/O" ) << std::endl ;

        return nil ;
    }

    WellGroupIOHandler* WellGroupIOHandler::get_handler(
        const std::string& filename )
    {
        std::string ext = GEO::FileSystem::extension( filename ) ;
        return create( ext ) ;
    }

    /*
     * Initializes the possible handler for IO files
     */
    void WellGroupIOHandler::initialize()
    {
        ringmesh_register_WellGroupIOHandler_creator( WLIOHandler, "wl" ) ;
        ringmesh_register_WellGroupIOHandler_creator( WGIOHandler, "wg" ) ;

    }
}
