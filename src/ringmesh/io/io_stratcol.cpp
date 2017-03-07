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
#include <ringmesh/geomodel/stratigraphic_column.h>

#include <geogram/basic/file_system.h>

namespace RINGMesh {
    class XMLStratigraphicColumnIOHandler: public StratigraphicColumnIOHandler {
    public:
        virtual void load( const std::string& filename, StratigraphicColumn& column )
        {
            throw RINGMeshException( "I/O", "This method is not implemented yet" ) ;
        }

        virtual void save(
            const StratigraphicColumn& column,
            const std::string& filename )
        {
            throw RINGMeshException( "I/O",
                "Saving of a StratigraphicColumn not implemented yet" ) ;
        }

    } ;

}

namespace RINGMesh {

    StratigraphicColumnIOHandler* StratigraphicColumnIOHandler::create(
        const std::string& format )
    {
        StratigraphicColumnIOHandler* handler =
            StratigraphicColumnIOHandlerFactory::create_object( format ) ;
        if( !handler ) {
            throw RINGMeshException( "I/O", "Unsupported file format: " + format ) ;
        }
        return handler ;
    }

    StratigraphicColumnIOHandler* StratigraphicColumnIOHandler::get_handler(
        const std::string& filename )
    {
        std::string ext = GEO::FileSystem::extension( filename ) ;
        return create( ext ) ;
    }

    /*
     * Initializes the possible handler for IO files
     */
    void StratigraphicColumnIOHandler::initialize()
    {
        ringmesh_register_StratigraphicColumnIOHandler_creator(
            XMLStratigraphicColumnIOHandler, "xml" ) ;
    }
}

