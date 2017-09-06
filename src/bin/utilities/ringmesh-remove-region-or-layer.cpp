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
#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_builder.h>
#include <ringmesh/geomodel/geomodel_builder_remove.h>
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
        Logger::div( "RINGMesh-Remove-VOI-Region-or-Layer" );
        Logger::out( "", "Welcome to RINGMesh-Remove-VOI-Region-or-Layer !" );
    }

    void import_arg_group_remove_region_or_layer()
    {
        GEO::CmdLine::declare_arg_group(
            "remove", "Remove a region or a layer (on the voi)" );
        GEO::CmdLine::declare_arg( "remove:is_region", true,
            "Type name of the GeoModelEntity to remove" );
        GEO::CmdLine::declare_arg(
            "remove:index", 0, "Index of the region or layer to remove" );
    }

    void import_arg_groups()
    {
        CmdLine::import_arg_group( "in" );
        CmdLine::import_arg_group( "out" );
        import_arg_group_remove_region_or_layer();
    }

    void get_dependent_entities_and_remove_them(
        std::set< gmme_id >& mesh_entities_to_delete,
        std::set< gmge_id >& geological_entities_to_delete,
        GeoModel3D& geomodel )
    {
        GeoModelBuilder3D builder( geomodel );
        builder.topology.get_dependent_entities(
            mesh_entities_to_delete, geological_entities_to_delete );
        builder.removal.remove_mesh_entities( mesh_entities_to_delete );
        builder.removal.remove_geological_entities(
            geological_entities_to_delete );
    }

    void save_output_geomodel( const GeoModel3D& geomodel )
    {
        const auto& out_model_file_name =
            GEO::CmdLine::get_arg( "out:geomodel" );
        if( out_model_file_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in out:geomodel" );
        }
        geomodel_save( geomodel, out_model_file_name );
    }

    void remove_a_voi_region_or_a_voi_layer(
        const std::string& in_model_file_name )
    {
        GeoModel3D geomodel;
        geomodel_load( geomodel, in_model_file_name );

        const auto is_region = GEO::CmdLine::get_arg_bool( "remove:is_region" );
        const auto gme_index = GEO::CmdLine::get_arg_uint( "remove:index" );

        std::set< gmme_id > mesh_entities_to_delete;
        std::set< gmge_id > geological_entities_to_delete;
        if( is_region )
        {
            if( gme_index
                >= geomodel.nb_mesh_entities( Region3D::type_name_static() ) )
            {
                throw RINGMeshException( "I/O", "Gme index higher than number "
                                                "of entities of the given "
                                                "type" );
            }
            // TODO check it is on the voi
            mesh_entities_to_delete.insert(
                { Region3D::type_name_static(), gme_index } );
        }
        else
        {
            if( gme_index >= geomodel.nb_geological_entities(
                                 Layer3D::type_name_static() ) )
            {
                throw RINGMeshException( "I/O", "Gme index higher than number "
                                                "of entities of the given "
                                                "type" );
            }
            // TODO check it is on the voi
            geological_entities_to_delete.insert(
                { Layer3D::type_name_static(), gme_index } );
        }
        get_dependent_entities_and_remove_them(
            mesh_entities_to_delete, geological_entities_to_delete, geomodel );

        save_output_geomodel( geomodel );
    }

    void run()
    {
        GEO::Stopwatch total( "Total time" );

        const auto& in_model_file_name = GEO::CmdLine::get_arg( "in:geomodel" );
        if( in_model_file_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in in:geomodel" );
        }
        const auto dimension = find_geomodel_dimension( in_model_file_name );
        if( dimension == 3 )
        {
            remove_a_voi_region_or_a_voi_layer( in_model_file_name );
        }
        else
        {
            throw RINGMeshException(
                "I/O", "Dimension must be 3, not ", dimension );
        }
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try
    {
        default_configure();
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
