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

#include <ringmesh/io/io.h>

#include <cstring>

#include <geogram/basic/file_system.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>

/*!
 * @file Please add a file description
 * @author Arnaud Botella
 * @todo Rename the file. Could'nt the 2 functions be moved somewhere else ?
 */

namespace RINGMesh
{
    /*!
     * Compares the contains of two files
     * @param[in] f1 the first filename
     * @param[in] f2 the second filename
     * @return return True if the files are identical
     */
    bool compare_files( const std::string& f1, const std::string& f2 )
    {
        const auto MAX_LINE_LEN =
            static_cast< unsigned int >( std::pow( 2, 16 ) - 1 );

        std::ifstream lFile( f1.c_str() );
        std::ifstream rFile( f2.c_str() );

        std::unique_ptr< char[] > lBuffer( new char[MAX_LINE_LEN]() );
        std::unique_ptr< char[] > rBuffer( new char[MAX_LINE_LEN]() );

        do
        {
            lFile.read( lBuffer.get(), MAX_LINE_LEN );
            rFile.read( rBuffer.get(), MAX_LINE_LEN );
            size_t numberOfRead = static_cast< size_t >( lFile.gcount() );

            if( std::memcmp( lBuffer.get(), rBuffer.get(), numberOfRead ) != 0 )
            {
                return false;
            }
        } while( lFile.good() || rFile.good() );
        return true;
    }

    void mesh_initialize()
    {
        GeoModelInputHandler2D::initialize();
        GeoModelInputHandler3D::initialize();
        GeoModelOutputHandler2D::initialize();
        GeoModelOutputHandler3D::initialize();
        WellGroupIOHandler::initialize();
    }

    /***************************************************************************/

    index_t find_geomodel_dimension( const std::string& filename )
    {
        auto ext = GEO::FileSystem::extension( filename );
        if( GeoModelInputHandlerFactory2D::has_creator( ext ) )
        {
            return GeoModelInputHandler2D::get_handler( filename )
                ->dimension( filename );
        }
        else if( GeoModelInputHandlerFactory3D::has_creator( ext ) )
        {
            return GeoModelInputHandler3D::get_handler( filename )
                ->dimension( filename );
        }
        else
        {
            ringmesh_assert_not_reached;
        }
        return 0;
    }

    template < index_t DIMENSION >
    bool geomodel_load(
        GeoModel< DIMENSION >& geomodel, const std::string& filename )
    {
        if( !GEO::FileSystem::is_file( filename ) )
        {
            throw RINGMeshException( "I/O", "File does not exist: ", filename );
        }
        Logger::out( "I/O", "Loading file ", filename, "..." );

        auto handler =
            GeoModelInputHandler< DIMENSION >::get_handler( filename );
        return handler->load_geomodel( filename, geomodel );
    }

    template < index_t DIMENSION >
    void geomodel_save(
        const GeoModel< DIMENSION >& geomodel, const std::string& filename )
    {
        Logger::out( "I/O", "Saving file ", filename, "..." );

        auto handler =
            GeoModelOutputHandler< DIMENSION >::get_handler( filename );
        handler->save_geomodel( geomodel, filename );
    }

    template bool io_api geomodel_load( GeoModel2D&, const std::string& );
    template void io_api geomodel_save(
        const GeoModel2D&, const std::string& );

    template bool io_api geomodel_load( GeoModel3D&, const std::string& );
    template void io_api geomodel_save(
        const GeoModel3D&, const std::string& );

} // namespace RINGMesh
