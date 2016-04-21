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

#include <ringmesh/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/command_line.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_builder.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/io.h>


/*!
 * @author Arnaud Botella
 */

namespace {
    using namespace RINGMesh ;

#endif
    /************************************************************************/
    /*!
     * Loads a GeoModel from a file
     * @param[in] filename the file to load
     * @param[out] model the model to fill
     */
    void geomodel_surface_load( const std::string& filename, GeoModel& model )
    {
        if( filename.empty() ) {
            throw RINGMeshException( "I/O",
                "No filename provided for structural model, use in:model" ) ;
        }
        std::ifstream input( filename.c_str() ) ;
        if( !input ) {
            throw RINGMeshException( "I/O", "Cannot open file: " + filename ) ;
        }

        GEO::Logger::out( "I/O" ) << "Loading file: " << filename << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->load( filename, model ) ;
    }

    /*!
     * Saves a GeoModel in a file
     * @param[in] model the model to save
     * @param[in] filename the filename where to save it
     */
    void geomodel_surface_save( const GeoModel& model, const std::string& filename )
    {
        GEO::Logger::out( "I/O" ) << "Saving file " << filename << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->save( model, filename ) ;
    }


    /*!
     * Loads a GeoModel from a file
     * @param[in] filename the file to load
     * @param][out] model the mesh to fill
     */
    void geomodel_volume_load( const std::string& filename, GeoModel& model )
    {
        GEO::Logger::out( "I/O" ) << "Loading file " << filename << "..."
            << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->load( filename, model ) ;
    }

    /*!
     * Saves a GeoModel in a file
     * @param[in] model the mesh to save
     * @param[in] filename the file where to save
     */
    void geomodel_volume_save( const GeoModel& model, const std::string& filename )
    {
        GEO::Logger::out( "I/O" ) << "Saving file " << filename << "..."
            << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->save( model, filename ) ;
    }


}

namespace RINGMesh {

    namespace CmdLine {
        void import_temp_in_out()
        {
            GEO::CmdLine::declare_arg( "in:model", "",
                "Filename of the input structural model" ) ;
            GEO::CmdLine::declare_arg( "in:mesh", "",
                "Filename of the input volumetric mesh" ) ;
            GEO::CmdLine::declare_arg( "out:model", "",
                "Saves the structural model" ) ;
            GEO::CmdLine::declare_arg( "out:mesh", "",
                "Saves the volumetric mesh of the structural model" ) ;
            GEO::CmdLine::declare_arg( "in:old_geomodel", "",
                "Saves the volumetric mesh of the structural model" ) ;
            ringmesh_register_IOHandler_creator( BMIOHandler, "bm" ) ;
#ifdef MINIZIP_FIXED
            ringmesh_register_IOHandler_creator( MMIOHandler, "mm" );
#endif
        }
    }
}
int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        GEO::Logger::div( "RINGMeshConvert" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshConvert !" << std::endl ;
        GEO::Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        GEO::Logger::out( "" ) << "Arnaud Botella <arnaud.botella@univ-lorraine.fr> "
            << std::endl ;

        CmdLine::import_arg_group( "in" ) ;
        CmdLine::import_arg_group( "out" ) ;
        CmdLine::import_temp_in_out() ;
        if( argc == 1 ) {
            GEO::CmdLine::show_usage() ;
            return 0 ;
        }

        std::vector< std::string > filenames ;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) ) {
            return 1 ;
        }

        GEO::Stopwatch total( "Total time" ) ;

        GeoModel model_in ;

        std::string geomodel_in_name = GEO::CmdLine::get_arg( "in:geomodel" ) ;
        if( geomodel_in_name.empty() ) {
            geomodel_in_name = GEO::CmdLine::get_arg( "in:model" ) ;
            if( geomodel_in_name.empty() ) {
                throw RINGMeshException( "I/O",
                    "Give at least a filename in in:model or in:geomodel" ) ;
            } else {
                geomodel_surface_load( geomodel_in_name, model_in ) ;
            }
        } else {
            geomodel_load( model_in, geomodel_in_name ) ;
        }

        if( GEO::CmdLine::get_arg( "in:geomodel" ).empty() ) {
            std::string mesh_in_name = GEO::CmdLine::get_arg( "in:mesh" ) ;
            if( !mesh_in_name.empty() ) {
                geomodel_volume_load( mesh_in_name, model_in ) ;
            }
        }

        std::string geomodel_out_name = GEO::CmdLine::get_arg( "out:geomodel" ) ;
        if(geomodel_out_name.empty()) {
            std::string model_out_name = GEO::CmdLine::get_arg( "out:model" ) ;
            std::string mesh_out_name = GEO::CmdLine::get_arg( "out:mesh" ) ;
            if( !model_out_name.empty() ) {
                geomodel_surface_save( model_in, model_out_name ) ;
            }
            if( !mesh_out_name.empty() ) {
                 geomodel_volume_save( model_in, mesh_out_name ) ;
             }
        }
        else {
            geomodel_save(model_in,geomodel_out_name) ;
        }

    } catch( const RINGMeshException& e ) {
        GEO::Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        GEO::Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    return 0 ;
}
