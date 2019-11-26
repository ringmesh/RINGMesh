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

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/tools/geomodel_repair.h>
#include <ringmesh/io/io.h>

/*!
 * @author Benjamin Chauvin
 */

namespace
{
    using namespace RINGMesh;

    void hello()
    {
        print_header_information();
        Logger::div( "RINGMesh-Repair" );
        Logger::out( "", "Welcome to RINGMesh-Repair !" );
    }

    void import_arg_group_repair()
    {
        GEO::CmdLine::declare_arg_group(
            "repair", "GeoModel repair processes" );
        GEO::CmdLine::declare_arg( "repair:mode", 0,
            "Repair mode: repair process to apply to the geomodel" );
    }

    void import_arg_groups()
    {
        CmdLine::import_arg_group( "in" );
        CmdLine::import_arg_group( "out" );
        import_arg_group_repair();
    }

    template < index_t DIMENSION >
    void repair_geomodel( const std::string& in_model_file_name )
    {
        GeoModel< DIMENSION > geomodel;
        geomodel_load( geomodel, in_model_file_name );

        index_t repair_mode = GEO::CmdLine::get_arg_uint( "repair:mode" );
        repair_geomodel( geomodel, static_cast< RepairMode >( repair_mode ) );

        std::string out_model_file_name =
            GEO::CmdLine::get_arg( "out:geomodel" );
        if( out_model_file_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in out:geomodel" );
        }
        geomodel_save( geomodel, out_model_file_name );
    }

    void run()
    {
        GEO::Stopwatch total( "Total time" );

        std::string in_model_file_name = GEO::CmdLine::get_arg( "in:geomodel" );
        if( in_model_file_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in in:geomodel" );
        }

        index_t dimension = find_geomodel_dimension( in_model_file_name );
        if( dimension == 2 )
        {
            repair_geomodel< 2 >( in_model_file_name );
        }
        else if( dimension == 3 )
        {
            repair_geomodel< 3 >( in_model_file_name );
        }
    }
} // namespace

int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try
    {
        hello();
        import_arg_groups();

        if( argc == 1 )
        {
            GEO::CmdLine::show_usage();
            return 0;
        }

        std::vector< std::string > filenames;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) )
        {
            return 1;
        }

        run();
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
    return 0;
}
