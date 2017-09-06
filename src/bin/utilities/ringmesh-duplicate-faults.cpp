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

#include <ringmesh/basic/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/duplicate_fntk_builder.h>
#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/io/io.h>

namespace RINGMesh
{
    void hello()
    {
        print_header_information();
        GEO::Logger::div( "RINGMesh-Duplicate-Faults" );
        GEO::Logger::out( "" ) << "Welcome to RINGMesh-Duplicate-Faults !"
                               << std::endl;
    }

    void import_arg_group_duplication_fntk()
    {
        GEO::CmdLine::declare_arg_group(
            "duplication", "Duplication of the fault network" );
        GEO::CmdLine::declare_arg( "duplication:gap", true,
            "Print statistics on the number of entities" );
    }

    void import_arg_groups()
    {
        CmdLine::import_arg_group( "in" );
        CmdLine::import_arg_group( "out" );
        import_arg_group_duplication_fntk();
    }

    void init()
    {
        default_configure();
        hello();
        import_arg_groups();
    }

    void load_input_geomodel( GeoModel3D& geomodel )
    {
        std::string input_geomodel_name =
            GEO::CmdLine::get_arg( "in:geomodel" );
        if( input_geomodel_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in in:geomodel" );
        }
        geomodel_load( geomodel, input_geomodel_name );
    }

    void duplicate_fntk( GeoModel3D& geomodel )
    {
        bool gap = GEO::CmdLine::get_arg_bool( "duplication:gap" );
        DuplicateInterfaceBuilder dib( geomodel );
        dib.duplicate_fault_network( gap );
    }

    void save_duplicated_fntk_geomodel( GeoModel3D& geomodel )
    {
        std::string output_geomodel_name =
            GEO::CmdLine::get_arg( "out:geomodel" );
        if( output_geomodel_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in out:geomodel" );
        }
        geomodel_save( geomodel, output_geomodel_name );
    }

    void run()
    {
        GEO::Stopwatch total( "Total time" );
        GeoModel3D geomodel;
        load_input_geomodel( geomodel );
        duplicate_fntk( geomodel );
        save_duplicated_fntk_geomodel( geomodel );
    }
}

/*!
 * @author Benjamin Chauvin
 */

int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try
    {
        init();

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
        GEO::Logger::err( e.category() ) << e.what() << std::endl;
        return 1;
    }
    catch( const std::exception& e )
    {
        GEO::Logger::err( "Exception" ) << e.what() << std::endl;
        return 1;
    }
    return 0;
}
