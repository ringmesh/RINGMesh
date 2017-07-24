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

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_validity.h>

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
        const unsigned int MAX_LINE_LEN = std::pow( 2, 16 ) - 1;

        std::ifstream lFile( f1.c_str() );
        std::ifstream rFile( f2.c_str() );

        std::unique_ptr< char[] > lBuffer( new char[MAX_LINE_LEN]() );
        std::unique_ptr< char[] > rBuffer( new char[MAX_LINE_LEN]() );

        do {
            lFile.read( lBuffer.get(), MAX_LINE_LEN );
            rFile.read( rBuffer.get(), MAX_LINE_LEN );
            size_t numberOfRead = static_cast< size_t >( lFile.gcount() );

            if( std::memcmp( lBuffer.get(), rBuffer.get(), numberOfRead ) != 0 ) {
                return false;
            }
        } while( lFile.good() || rFile.good() );
        return true;
    }

    void mesh_initialize()
    {
        GeoModelIOHandler2D::initialize_geomodel_output();
        GeoModelIOHandler3D::initialize_geomodel_output();
        WellGroupIOHandler::initialize();
    }

    /***************************************************************************/

    void zip_file( zipFile zf, const std::string& name )
    {
        std::fstream file( name.c_str(), std::ios::in | std::ios::binary );
        file.seekg( 0, std::ios::end );
        long size = file.tellg();
        file.seekg( 0, std::ios::beg );
        std::vector< char > buffer( size );
        file.read( &buffer[0], size );
        zipOpenNewFileInZip( zf, name.c_str(), nullptr, nullptr, 0, nullptr, 0,
            nullptr, Z_DEFLATED, Z_DEFAULT_COMPRESSION );
        zipWriteInFileInZip( zf, size == 0 ? "" : &buffer[0], size );
        zipCloseFileInZip( zf );
        file.close();
    }

    void unzip_file( unzFile uz, const char filename[MAX_FILENAME] )
    {
        unzLocateFile( uz, filename, 0 );
        unzip_current_file( uz, filename );
    }

    void unzip_current_file( unzFile uz, const char filename[MAX_FILENAME] )
    {
        char read_buffer[READ_SIZE];
        if( unzOpenCurrentFile( uz ) != UNZ_OK ) {
            unzClose( uz );
            throw RINGMeshException( "ZLIB", "Could not open file" );
        }
        FILE *out = fopen( filename, "wb" );
        if( out == nullptr ) {
            unzCloseCurrentFile( uz );
            unzClose( uz );
            throw RINGMeshException( "ZLIB", "Could not open destination file" );
        }
        int error = UNZ_OK;
        do {
            error = unzReadCurrentFile( uz, read_buffer, READ_SIZE );
            if( error < 0 ) {
                unzCloseCurrentFile( uz );
                unzClose( uz );
                fclose( out );
                throw RINGMeshException( "ZLIB", "Invalid error: ", error );
            }
            if( error > 0 ) {
                fwrite( read_buffer, error, 1, out );
            }
        } while( error > 0 );
        fclose( out );
        unzCloseCurrentFile( uz );
    }

    /***************************************************************************/

    template< index_t DIMENSION >
    bool GeoModelIOHandler< DIMENSION >::load_geomodel(
        const std::string& filename,
        GeoModel< DIMENSION >& geomodel )
    {
        load( filename, geomodel );
        Logger::out( "I/O", " Loaded geomodel ", geomodel.name(), " from ",
            filename );
        return is_geomodel_valid( geomodel );
    }

    template< index_t DIMENSION >
    void GeoModelIOHandler< DIMENSION >::save_geomodel(
        const GeoModel< DIMENSION >& geomodel,
        const std::string& filename )
    {
        save( geomodel, filename );
    }

    index_t find_geomodel_dimension( const std::string& filename )
    {
        std::string ext = GEO::FileSystem::extension( filename );
        if( GeoModelIOHandlerFactory2D::has_creator( ext ) ) {
            return GeoModelIOHandler2D::get_handler( filename )->dimension(
                filename );
        } else if( GeoModelIOHandlerFactory3D::has_creator( ext ) ) {
            return GeoModelIOHandler3D::get_handler( filename )->dimension(
                filename );
        } else {
            ringmesh_assert_not_reached;
        }
        return 0;
    }

    template< index_t DIMENSION >
    bool geomodel_load(
        GeoModel< DIMENSION >& geomodel,
        const std::string& filename )
    {
        if( !GEO::FileSystem::is_file( filename ) ) {
            throw RINGMeshException( "I/O", "File does not exist: ", filename );
        }
        Logger::out( "I/O", "Loading file ", filename, "..." );

        std::unique_ptr< GeoModelIOHandler< DIMENSION > > handler(
            GeoModelIOHandler< DIMENSION >::get_handler( filename ) );
        return handler->load_geomodel( filename, geomodel );
    }

    template< index_t DIMENSION >
    void geomodel_save(
        const GeoModel< DIMENSION >& geomodel,
        const std::string& filename )
    {
        Logger::out( "I/O", "Saving file ", filename, "..." );

        std::unique_ptr< GeoModelIOHandler< DIMENSION > > handler(
            GeoModelIOHandler< DIMENSION >::get_handler( filename ) );
        handler->save_geomodel( geomodel, filename );
    }

    /************************************************************************/

    template< index_t DIMENSION >
    GeoModelIOHandler< DIMENSION >* GeoModelIOHandler< DIMENSION >::create(
        const std::string& format )
    {
        GeoModelIOHandler< DIMENSION >* handler =
            GeoModelIOHandlerFactory< DIMENSION >::create_object( format );
        if( !handler ) {
            std::vector< std::string > names;
            GeoModelIOHandlerFactory< DIMENSION >::list_creators( names );
            Logger::err( "I/O", "Currently supported file formats are: " );
            for( const std::string& name : names ) {
                Logger::err( "I/O", " ", name );
            }

            throw RINGMeshException( "I/O", "Unsupported file format: ", format );
        }
        return handler;
    }

    template< index_t DIMENSION >
    std::unique_ptr< GeoModelIOHandler< DIMENSION > > GeoModelIOHandler< DIMENSION >::get_handler(
        const std::string& filename )
    {
        std::string ext = GEO::FileSystem::extension( filename );
        return std::unique_ptr< GeoModelIOHandler< DIMENSION > >( create( ext ) );
    }

    template class RINGMESH_API GeoModelIOHandler< 2 > ;
    template class RINGMESH_API GeoModelIOHandler< 3 > ;

    template bool RINGMESH_API geomodel_load( GeoModel2D&, const std::string& );
    template void RINGMESH_API geomodel_save(
        const GeoModel2D&,
        const std::string& );

    template bool RINGMESH_API geomodel_load( GeoModel3D&, const std::string& );
    template void RINGMESH_API geomodel_save(
        const GeoModel3D&,
        const std::string& );

}
