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
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geo_model_repair.h>

#include <ringmesh/io.h>

/*!
 * @author Arnaud Botella
 */

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;
        CmdLine::import_arg_group("out") ;
        CmdLine::import_arg_group("in") ;
        if( argc == 1 ) {
            GEO::CmdLine::show_usage() ;
            return 0 ;
        }

        std::vector< std::string > filenames ;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) ) {
            return 1 ;
        }

        GEO::Stopwatch total( "Total time" ) ;


        GeoModel M ;
        std::string file_name( GEO::CmdLine::get_arg( "in:geomodel" ) ) ;

        GEO::Logger::out( "RINGMesh Repair" ) << "Loading structural model:"
            << file_name << std::endl ;

        geomodel_load( M, file_name ) ;
        if( !is_geomodel_valid( M ) ) {
            // Try to repair the model if it is not valid
            geo_model_mesh_repair( M ) ;

            // Test the validity again
            if( is_geomodel_valid( M ) ) {
                std::string fixed_file_name(
                    GEO::CmdLine::get_arg( "out:geomodel" ) ) ;
                geomodel_save( M, fixed_file_name ) ;
                GEO::Logger::out( "RINGMesh Repair" ) << "Invalid geological model "
                    << M.name()
                    << " has been successfully fixed and is saved under: "
                    << fixed_file_name << std::endl ;
            } else {
                throw RINGMeshException( "RINGMesh Repair",
                    "Fixing the invalid geological model " + M.name()
                        + " failed." ) ;
            }
        } else {
            GEO::Logger::out( "RINGMesh Repair" ) << "The geological model "
                << M.name() << " is valid " << std::endl ;
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
