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

#include <ringmesh/ringmesh_tests_config.h>

#include <geogram/basic/command_line.h>
#include <ringmesh/geomodel/tools/geomodel_repair.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * Load and fix a given structural model file.
 * @author Jeanne Pellerin
 */
int main()
{
    using namespace RINGMesh;

    try
    {
        std::string file_name( ringmesh_test_data_path );
        file_name += "annot.ml";

        Logger::out( "RINGMesh Test",
            "Loading and fixing structural geomodel:", file_name );

        // Check only model geometry
        GEO::CmdLine::set_arg( "validity:do_not_check", "tGI" );

        // Load the geomodel
        GeoModel3D geomodel;
        bool init_model_is_valid{ geomodel_load( geomodel, file_name ) };
        if( init_model_is_valid )
        {
            throw RINGMeshException( "RINGMesh Test", "Input test model ",
                geomodel.name(),
                " must be invalid to check the repair functionalities." );
        }

        Logger::out( "RINGMesh Test", "Repairing..." );

        // Repair the geomodel
        repair_geomodel( geomodel, RepairMode::ALL );

        // Test the validity again
        if( !is_geomodel_valid( geomodel, ValidityCheckMode::GEOMETRY ) )
        {
            throw RINGMeshException(
                "RINGMesh Test", "Fixing the invalid geological model "
                                     + geomodel.name() + " failed." );
        }

        Logger::out( "TEST", "SUCCESS" );
        return 0;
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
}
