/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/basic/common.h>
#include <ringmesh/ringmesh_tutorials_config.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @author Antoine Mazuyer
 * @author Francois Bonneau
 * @author Pierre Anquez
 */

/*!
 * Purpose of this main is to show the methods
 * to be used to build a GeoModel from scratch.
 */

int main()
{
#ifdef RINGMESH_WITH_TETGEN

    using namespace RINGMesh;

    try
    {
        // Say Hello
        print_header_information();
        Logger::div( "RINGMesh Training" );
        Logger::out( "", "Welcome to the RINGMesh training for basic "
                         "functionalities on GeoModel !" );

        // Next line is a feature of geogram which measure
        // the time of execution.
        GEO::Stopwatch total( "Total time" );

        // We instantiate the class GeoModel
        GeoModel3D geomodel;

        // load GeoModel
        // here you can load whatever the model you want in the
        // ringmesh_home/test/data directory
        std::string input_file_name( ringmesh_tutorials_data_path );
        input_file_name += "modelA1.ml";

        // We do not want to check the model validity at loading
        GEO::CmdLine::set_arg(
            "validity:do_not_check", "A" ); // "A" for all checks

        // function to load a geomodel
        geomodel_load( geomodel, input_file_name );

        // function to print the statistics of the geomodel in the command
        // terminal
        print_geomodel_mesh_stats( geomodel );

        // build volumetric mesh in regions
        GEO::CmdLine::set_arg( "algo:tet", "TetGen" );
        tetrahedralize( geomodel );

        // function to print the statistics of the geomodel in the command
        // terminal
        print_geomodel_mesh_stats( geomodel );

        // check tetrahedralized GeoModel validity
        ValidityCheckMode checks{ ValidityCheckMode::ALL };
        bool is_valid = is_geomodel_valid( geomodel, checks );
        if( !is_valid )
        {
            throw RINGMeshException( "RINGMesh Test", "Tetrahedralized model ",
                geomodel.name(), " is not valid." );
        }

        // set the name of the geomodel to output
        // you can customize the path
        std::string output_file_name = "modelA1.gm";
        geomodel_save( geomodel, output_file_name );
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

#endif

    return 0;
}
