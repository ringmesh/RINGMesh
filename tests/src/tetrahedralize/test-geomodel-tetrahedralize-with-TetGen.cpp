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

#include <geogram/basic/command_line.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_tools.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @file Test global tetrahedralization of a GeoModel
 */

int main()
{
    using namespace RINGMesh;

    try
    {
        CmdLine::import_arg_group( "global" );
        GEO::CmdLine::set_arg( "algo:tet", "TetGen" );

        std::string file_name( ringmesh_test_data_path );
        file_name += "modelA6.ml";

        // Check only model geometry
        GEO::CmdLine::set_arg( "validity:do_not_check", "tG" );

        // Loading the GeoModel
        GeoModel3D geomodel;
        bool loaded_model_is_valid = geomodel_load( geomodel, file_name );

        if( !loaded_model_is_valid )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Failed when building model ", geomodel.name(),
                ": the model geometry is not valid." );
        }

#ifdef RINGMESH_WITH_TETGEN

        // Tetrahedralize the GeoModel
        tetrahedralize( geomodel, NO_ID, false );
        for( index_t r : range( geomodel.nb_regions() ) )
        {
            if( !geomodel.region( r ).is_meshed() )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Failed when tetrahedralize model ", geomodel.name(),
                    " Region ", r, " is not meshed ",
                    "(maybe the TetGen call have failed)." );
            }
        }

        // Check validity of tetrahedralized model
        ValidityCheckMode checks{ ValidityCheckMode::GEOMETRY };

        if( !is_geomodel_valid( geomodel, checks ) )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Failed when tetrahedralize model ", geomodel.name(),
                ": the model becomes invalid." );
        }

#endif
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
