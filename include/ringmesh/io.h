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

#ifndef __RINGMESH_IO__
#define __RINGMESH_IO__

#include <ringmesh/common.h>
#include <ringmesh/geo_model_element.h>

#include <geogram/basic/string.h>


#include <third_party/zlib/unzip.h>
#include <third_party/zlib/zip.h>

#include <geogram/basic/factory.h>
#define MAX_FILENAME 512
#define READ_SIZE 8192

/*!
 * @file Global input - ouput functions of RINGMesh
 * @author Various
 */

namespace RINGMesh {
    class GeoModel ;
    class WellGroup ;
}

namespace GEO {
    class MeshSubElementsStore ;
}

namespace RINGMesh {

    bool RINGMESH_API compare_files( const std::string& f1, const std::string& f2 ) ;

    void RINGMESH_API geomodel_load(
        GeoModel& model,
        const std::string& filename ) ;

    void RINGMESH_API geomodel_save(
        const GeoModel& model,
        const std::string& filename ) ;

    void RINGMESH_API well_load( const std::string& filename, WellGroup& wells ) ;

    class RINGMESH_API GeoModelIOHandler: public GEO::Counted {
    public:
        static void initialize_full_geomodel_output() ;

        static void initialize_boundary_geomodel_output() ;

        static GeoModelIOHandler* create( const std::string& format ) ;

        static GeoModelIOHandler* get_handler( const std::string& filename ) ;

        virtual void load( const std::string& filename, GeoModel& model ) = 0 ;

        virtual void save( const GeoModel& model, const std::string& filename ) = 0 ;

    protected:
        GeoModelIOHandler()
        {
        }

        virtual ~GeoModelIOHandler()
        {
        }
    } ;

    typedef GEO::SmartPointer< GeoModelIOHandler > IOHandler_var ;
    typedef GEO::Factory0< GeoModelIOHandler > IOHandlerFactory ;

#define ringmesh_register_IOHandler_creator( type, name ) \
    geo_register_creator( IOHandlerFactory, type, name )

    /***************************************************************************/
    class RINGMESH_API WellGroupIOHandler: public GEO::Counted {
    public:
        static void initialize() ;

        static WellGroupIOHandler* create( const std::string& format ) ;

        static WellGroupIOHandler* get_handler( const std::string& filename ) ;

        virtual void load( const std::string& filename, WellGroup& mesh ) = 0 ;

        virtual void save( const WellGroup& mesh, const std::string& filename ) = 0 ;

    protected:
        WellGroupIOHandler()
        {
        }

        virtual ~WellGroupIOHandler()
        {
        }
    } ;
    typedef GEO::SmartPointer< WellGroupIOHandler > WellGroupIOHandler_var ;
    typedef GEO::Factory0< WellGroupIOHandler > WellGroupIOHandlerFactory ;

#define ringmesh_register_WellGroupIOHandler_creator( type, name ) \
    geo_register_creator( WellGroupIOHandlerFactory, type, name )

    /***************************************************************************/

    void RINGMESH_API mesh_initialize() ;

    void RINGMESH_API zip_file( zipFile zf, const std::string& name ) ;

    void RINGMESH_API unzip_file( unzFile uz, char filename[MAX_FILENAME] ) ;

}
#endif
