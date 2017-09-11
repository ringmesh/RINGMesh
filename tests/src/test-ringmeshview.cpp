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

#include <ringmesh/ringmesh_tests_config.h>

#include <chrono>
#include <future>

#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>

#include <ringmesh/basic/common.h>
#include <ringmesh/basic/logger.h>

#ifdef RINGMESH_TEST_GRAPHICS

#include <geogram/basic/process.h>

#include <ringmesh/visualization/gfx_application.h>

/*!
 * @author Pierre Anquez
 */

namespace RINGMesh
{
    void open_viewer_load_geomodel_then_close(
        const int argc, char** argv, const std::string& glup_profile )
    {
        RINGMeshApplication app( argc, argv );

        GEO::CmdLine::set_arg( "GLUP_profile", glup_profile );

        // Create the tasks for launching the app window
        // and the one for closing the window
        std::future< void > start =
            std::async( std::launch::async, &RINGMeshApplication::start, &app );
        std::this_thread::sleep_for( std::chrono::seconds( 3 ) );
        std::future< void > end =
            std::async( std::launch::async, &RINGMeshApplication::quit, &app );
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
        char ringmesh_view[] = "ringmesh-view";
        std::string input_model_file_name( ringmesh_test_data_path );
        input_model_file_name += "modelA6.ml";
        char* input_model = &input_model_file_name[0];

        char* argv[2] = { ringmesh_view, input_model };

        // Two arguments: one for 'ringmeshview' and one for the input file
        const int argc = 2;

        std::vector< std::string > GLUP_profiles( 1, "" );
        GLUP_profiles[0] = "auto";
        //        GLUP_profiles[1] = "GLUP150" ;
        //        GLUP_profiles[2] = "GLUP440" ;
        //        GLUP_profiles[3] = "VanillaGL" ;

        for( const std::string& profile : GLUP_profiles )
        {
            open_viewer_load_geomodel_then_close( argc, argv, profile );
        }
    }
    catch( const RINGMeshException& e )
    {
        Logger::err( e.category(), e.what() );
        return 1;
    }
    catch( const std::exception& e )
    {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}

#else
int main()
{
    using namespace RINGMesh;

    default_configure();
    Logger::out( "RINGMeshView",
        "To test RINGMesh viewer you need to configure ",
        "the project with the RINGMESH_TEST_GRAPHICS option ON" );
    return 0;
}
#endif
