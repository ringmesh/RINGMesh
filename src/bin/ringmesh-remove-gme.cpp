/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

namespace {
    using namespace RINGMesh;

    void hello()
    {
        print_header_information();
        Logger::div( "RINGMesh-Repair" );
        Logger::out( "", "Welcome to RINGMesh-Remove-GME !" );
    }

    void import_arg_group_remove_gme()
    {
        GEO::CmdLine::declare_arg_group( "remove", "GeoModel repair processes" );
        GEO::CmdLine::declare_arg( "remove:type", Region3D::type_name_static(),
            "Type name of the GeoModelEntity to remove" );
        GEO::CmdLine::declare_arg( "remove:index", 0,
            "Index of the entity to remove" );
    }

    void import_arg_groups()
    {
        CmdLine::import_arg_group( "in" );
        CmdLine::import_arg_group( "out" );
        import_arg_group_remove_gme();
    }

    template< index_t DIMENSION >
    void remove_gme( const std::string& in_model_file_name )
    {
        GeoModel< DIMENSION > geomodel;
        geomodel_load( geomodel, in_model_file_name );

        const std::string gme_type = GEO::CmdLine::get_arg( "remove:type" );
        index_t gme_index = GEO::CmdLine::get_arg_uint( "remove:index" );

        GeoModelBuilder< DIMENSION > builder( geomodel );
        MeshEntityType mesh_entity_type( gme_type );
        GeologicalEntityType geological_entity_type( gme_type );
        if( geomodel.entity_type_manager().mesh_entity_manager.is_valid_type(
            mesh_entity_type ) ) {
            if( gme_index >= geomodel.nb_mesh_entities( mesh_entity_type ) ) {
                throw RINGMeshException( "I/O",
                    "Gme index higher than number of entities of the given type" );
            }
            builder.removal.remove_mesh_entity_and_dependencies(
                gmme_id( mesh_entity_type, gme_index ) );
        } else if( geomodel.entity_type_manager().geological_entity_manager.is_valid_type(
            geological_entity_type ) ) {
            if( gme_index
                >= geomodel.nb_geological_entities( geological_entity_type ) ) {
                throw RINGMeshException( "I/O",
                    "Gme index higher than number of entities of the given type" );
            }
            builder.removal.remove_geological_entity_and_dependencies(
                gmge_id( geological_entity_type, gme_index ) );
        } else {
            throw RINGMeshException( "I/O", "invalid mesh type" );
        }

        std::string out_model_file_name = GEO::CmdLine::get_arg( "out:geomodel" );
        if( out_model_file_name.empty() ) {
            throw RINGMeshException( "I/O",
                "Give at least a filename in out:geomodel" );
        }
        geomodel_save( geomodel, out_model_file_name );
    }

    void run()
    {
        GEO::Stopwatch total( "Total time" );

        const std::string in_model_file_name = GEO::CmdLine::get_arg(
            "in:geomodel" );
        if( in_model_file_name.empty() ) {
            throw RINGMeshException( "I/O",
                "Give at least a filename in in:geomodel" );
        }
        index_t dimension = find_geomodel_dimension( in_model_file_name );
        if( dimension == 2 ) {
            remove_gme< 2 >( in_model_file_name );
        } else if( dimension == 3 ) {
            remove_gme< 3 >( in_model_file_name );
        } else {
            throw RINGMeshException( "I/O", "Forbidden dimension", dimension );
        }
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try {
        default_configure();
        hello();
        import_arg_groups();

        if( argc == 1 ) {
            GEO::CmdLine::show_usage();
            return 0;
        }

        std::vector< std::string > filenames;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) ) {
            return 1;
        }

        run();

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category(), e.what() );
        return 1;
    } catch( const std::exception& e ) {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    return 0;
}
