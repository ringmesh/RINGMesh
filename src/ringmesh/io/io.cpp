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

#include <ringmesh/io/io.h>

#include <cstring>

#include <geogram/basic/file_system.h>

/*!
 * @file Please add a file description
 * @author Arnaud Botella
 * @todo Rename the file. Could'nt the 2 functions be moved somewhere else ? 
 */

namespace RINGMesh {

    /*!
     * Compares the contains of two files
     * @param[in] f1 the first filename
     * @param[in] f2 the second filename
     * @return return True if the files are identical
     */
    bool compare_files( const std::string& f1, const std::string& f2 )
    {
        const unsigned int MAX_LINE_LEN = 65535 ;

        std::ifstream lFile( f1.c_str() ) ;
        std::ifstream rFile( f2.c_str() ) ;

        char* lBuffer = new char[MAX_LINE_LEN]() ;
        char* rBuffer = new char[MAX_LINE_LEN]() ;

        do {
            lFile.read( lBuffer, MAX_LINE_LEN ) ;
            rFile.read( rBuffer, MAX_LINE_LEN ) ;
            size_t numberOfRead = static_cast< size_t >( lFile.gcount() ) ;

            if( std::memcmp( lBuffer, rBuffer, numberOfRead ) != 0 ) {
                delete[] lBuffer ;
                delete[] rBuffer ;
                return false ;
            }
        } while( lFile.good() || rFile.good() ) ;
        delete[] lBuffer ;
        delete[] rBuffer ;
        return true ;
    }

    void mesh_initialize()
    {
        GeoModelIOHandler::initialize_full_geomodel_output() ;
        GeoModelIOHandler::initialize_boundary_geomodel_output() ;
        WellGroupIOHandler::initialize() ;
    }

    /***************************************************************************/

    void zip_file( zipFile zf, const std::string& name )
    {
        std::fstream file( name.c_str(), std::ios::in | std::ios::binary ) ;
        file.seekg( 0, std::ios::end ) ;
        long size = file.tellg() ;
        file.seekg( 0, std::ios::beg ) ;
        std::vector< char > buffer( size ) ;
        file.read( &buffer[0], size ) ;
        zipOpenNewFileInZip( zf, name.c_str(), nullptr, nullptr, 0, nullptr, 0,
            nullptr, Z_DEFLATED, Z_DEFAULT_COMPRESSION ) ;
        zipWriteInFileInZip( zf, size == 0 ? "" : &buffer[0], size ) ;
        zipCloseFileInZip( zf ) ;
        file.close() ;
    }

    void unzip_file( unzFile uz, const char filename[MAX_FILENAME] )
    {
        unzLocateFile( uz, filename, 0 ) ;
        unzip_current_file( uz, filename ) ;
    }

    void unzip_current_file( unzFile uz, const char filename[MAX_FILENAME] )
    {
        char read_buffer[READ_SIZE] ;
        if( unzOpenCurrentFile( uz ) != UNZ_OK ) {
            unzClose( uz ) ;
            throw RINGMeshException( "ZLIB", "Could not open file" ) ;
        }
        FILE *out = fopen( filename, "wb" ) ;
        if( out == nullptr ) {
            unzCloseCurrentFile( uz ) ;
            unzClose( uz ) ;
            throw RINGMeshException( "ZLIB", "Could not open destination file" ) ;
        }
        int error = UNZ_OK ;
        do {
            error = unzReadCurrentFile( uz, read_buffer, READ_SIZE ) ;
            if( error < 0 ) {
                unzCloseCurrentFile( uz ) ;
                unzClose( uz ) ;
                fclose( out ) ;
                throw RINGMeshException( "ZLIB",
                    "Invalid error: " + GEO::String::to_string( error ) ) ;
            }
            if( error > 0 ) {
                fwrite( read_buffer, error, 1, out ) ;
            }
        } while( error > 0 ) ;
        fclose( out ) ;
        unzCloseCurrentFile( uz ) ;
    }

    /***************************************************************************/

    GeoModelIOHandler* GeoModelIOHandler::create( const std::string& format )
    {
        GeoModelIOHandler* handler = GeoModelIOHandlerFactory::create_object(
            format ) ;
        if( !handler ) {
            std::vector< std::string > names ;
            GeoModelIOHandlerFactory::list_creators( names ) ;
            Logger::err( "I/O", "Currently supported file formats are: " ) ;
            for( index_t i = 0; i < names.size(); i++ ) {
                Logger::err( "I/O", " ", names[i] ) ;
            }

            throw RINGMeshException( "I/O", "Unsupported file format: " + format ) ;
        }
        return handler ;
    }

    GeoModelIOHandler* GeoModelIOHandler::get_handler( const std::string& filename )
    {
        std::string ext = GEO::FileSystem::extension( filename ) ;
        return create( ext ) ;
    }

}
