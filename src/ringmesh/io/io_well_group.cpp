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

#include <geogram/basic/file_system.h>
#include <geogram/basic/line_stream.h>

#include <ringmesh/geomodel/geomodel.h>

#include <ringmesh/mesh/geogram_mesh.h>
#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/well.h>

/*!
 * @file Implements the input - output of WellGroup
 * @author Arnaud Botella
 */

namespace {
    using namespace RINGMesh;

#include "well_group/io_smesh.cpp"
#include "well_group/io_wl.cpp"

}

namespace RINGMesh {

    void well_load( const std::string& filename, WellGroup& wells )
    {
        Logger::out( "I/O", "Loading file ", filename, "..." );

        std::unique_ptr< WellGroupIOHandler > handler =
            WellGroupIOHandler::get_handler( filename );
        handler->load( filename, wells );
    }

    WellGroupIOHandler* WellGroupIOHandler::create( const std::string& format )
    {
        WellGroupIOHandler* handler = WellGroupIOHandlerFactory::create_object(
            format );
        if( !handler ) {
            std::vector< std::string > names;
            WellGroupIOHandlerFactory::list_creators( names );
            Logger::err( "I/O", "Currently supported file formats are: " );
            for( index_t i = 0; i < names.size(); i++ ) {
                Logger::err( "I/O", " ", names[i] );
            }

            throw RINGMeshException( "I/O", "Unsupported file format: " + format );
        }
        return handler;
    }

    std::unique_ptr< WellGroupIOHandler > WellGroupIOHandler::get_handler(
        const std::string& filename )
    {
        std::string ext = GEO::FileSystem::extension( filename );
        return std::unique_ptr< WellGroupIOHandler >( create( ext ) );
    }

    /*
     * Initializes the possible handler for IO files
     */
    void WellGroupIOHandler::initialize()
    {
        ringmesh_register_WellGroupIOHandler_creator( WLIOHandler, "wl" );
        ringmesh_register_WellGroupIOHandler_creator( SmeshIOHandler, "smesh" );
    }
}
