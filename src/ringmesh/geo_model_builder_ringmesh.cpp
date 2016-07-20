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

#include <ringmesh/geo_model_builder_ringmesh.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <set>
#include <stack>

#include <geogram/basic/line_stream.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/points/colocate.h>

#include <ringmesh/io.h>
#include <ringmesh/algorithm.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geometry.h>
#include <ringmesh/geogram_mesh_repair.h>
#include <ringmesh/utils.h>

/*!
 * @file ringmesh/geo_model_builder.cpp
 * @brief Implementation of the classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace {
    using namespace RINGMesh ;

    bool match_mesh_entity_type( const std::string& type )
    {
        if( type == Corner::type_name_ ) return true ;
        if( type == Line::type_name_ ) return true ;
        if( type == Surface::type_name_ ) return true ;
        if( type == Region::type_name_ ) return true ;
        return false ;
    }
}

namespace RINGMesh {

    void GeoModelBuilderGM::load_mesh_entities( GEO::LineInput& file_line )
    {
        while( !file_line.eof() && file_line.get_line() ) {

            file_line.get_fields() ;
            if( file_line.nb_fields() > 0 ) {
                // Name of the model
                if( file_line.field_matches( 0, "Name" ) ) {
                    if( file_line.nb_fields() > 1 ) {
                        set_model_name( file_line.field( 1 ) ) ;
                    }
                }
                // Number of entities of a given type
                else if( file_line.field_matches( 0, "Nb" ) ) {
                    // Allocate the space
                    create_mesh_entities( file_line.field( 1 ),
                        file_line.field_as_uint( 2 ) ) ;
                }
                // Mesh entities
                else if( match_mesh_entity_type( file_line.field( 0 ) ) ) {
                    // Read this entity
                    // First line : type - id - name - geol_feature
                    if( file_line.nb_fields() < 4 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line.line_number() )
                                + ", 4 fields are expected, the type, id, name, and geological feature" ) ;
                    }
                    const std::string type = file_line.field( 0 ) ;
                    index_t id = file_line.field_as_uint( 1 ) ;
                    GME::gme_t entity( type, id ) ;
                    set_mesh_entity_name( entity, file_line.field( 2 ) ) ;
                    set_mesh_entity_geol_feature( entity,
                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;

                    // Read second line
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    if( type == Region::type_name_ ) {
                        // Second line : signed indices of boundaries
                        for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                            bool side = false ;
                            if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                                side = true ;
                            }
                            index_t s ;
                            GEO::String::from_string( &file_line.field( c )[1], s ) ;

                            add_mesh_entity_boundary( entity,
                                GME::gme_t( Surface::type_name_, s ), side ) ;
                        }
                    } else {
                        // Second line : indices of its in boundaries
                        for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                            add_mesh_entity_in_boundary( entity,
                                GME::gme_t( model().mesh_entity( entity ).boundary_type(),
                                    file_line.field_as_uint( c ) ) ) ;
                        }
                    }
                }
                // Universe
                else if( file_line.field_matches( 0, "Universe" ) ) {
                    // Second line: signed indices of boundaries
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &file_line.field( c )[1], s ) ;

                        add_universe_boundary( GME::gme_t( Surface::type_name_, s ),
                            side ) ;
                    }
                }
            }
        }
    }
    void GeoModelBuilderGM::load_file()
    {
        unzFile uz = unzOpen( filename_.c_str() ) ;
        unz_global_info global_info ;
        if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
            unzClose( uz ) ;
            throw RINGMeshException( "ZLIB", "Could not read file global info" ) ;
        }

        const std::string mesh_entity_file( "mesh_entities.txt" ) ;
        unzip_one_file( uz, mesh_entity_file.c_str() ) ;
        GEO::LineInput line_mesh_entity( mesh_entity_file ) ;
        load_mesh_entities( line_mesh_entity ) ;
        GEO::FileSystem::delete_file( mesh_entity_file ) ;

        load_meshes( Corner::type_name_, uz ) ;
        load_meshes( Line::type_name_, uz ) ;
        load_meshes( Surface::type_name_, uz ) ;
        load_meshes( Region::type_name_, uz ) ;

        const std::string geological_entity_file( "geological_entities.txt" ) ;
        unzip_one_file( uz, geological_entity_file.c_str() ) ;
        GEO::LineInput line_geological_entity( geological_entity_file ) ;
        load_geological_entities( line_geological_entity ) ;
        GEO::FileSystem::delete_file( geological_entity_file ) ;

        unzClose( uz ) ;
    }

    void GeoModelBuilderGM::load_geological_entities( GEO::LineInput& file_line )
    {
        while( !file_line.eof() && file_line.get_line() ) {
            file_line.get_fields() ;
            if( file_line.nb_fields() > 0 ) {
                const std::string type = file_line.field( 0 ) ;
                index_t id = file_line.field_as_uint( 1 ) ;
                GME::gme_t entity( type, id ) ;
                set_geological_entity_name( entity, file_line.field( 2 ) ) ;
                set_geological_entity_geol_feature( entity,
                    GME::determine_geological_type( file_line.field( 3 ) ) ) ;
                file_line.get_line() ;
                file_line.get_fields() ;
                const GeoModelGeologicalEntity& cur_gme = model().geological_entity(
                    type, id ) ;
                const std::string& child_type = cur_gme.child_type_name() ;
                for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                    add_geological_entity_child( entity,
                        GME::gme_t( child_type, file_line.field_as_uint( in_b ) ) ) ;
                }
            }
        }
    }

    void GeoModelBuilderGM::load_meshes( const std::string& type, unzFile& uz )
    {
        for( index_t el = 0; el < model().nb_mesh_entities( type ); el++ ) {
            GME::gme_t cur_gme( type, el ) ;
            std::string file_to_extract_and_load ;
            build_string_for_geo_model_entity_export( cur_gme,
                file_to_extract_and_load ) ;
            std::string filename = file_to_extract_and_load + ".geogram" ;
            if( unzLocateFile( uz, filename.c_str(), 0 ) != UNZ_OK ) {
                if( type != Region::type_name_ ) {
                    std::string message = "Invalid format of .gm file" ;
                    message += "\n.geogram file (defining mesh) is missing." ;
                    throw RINGMeshException( "I/O", message ) ;
                }
                return ; // a region is not necessary meshed.
            }
            unzip_one_file( uz, filename.c_str() ) ;
            Mesh cur_mesh( model(), 3, false ) ;
            GEO::MeshIOFlags flags ;
            flags.set_attribute( GEO::MESH_ALL_ATTRIBUTES ) ;
            Logger::instance()->set_minimal( true ) ;
            MeshBuilder builder(cur_mesh);
            builder.load_mesh( filename, flags ) ;
            assign_mesh_to_entity( cur_mesh, cur_gme ) ;
            Logger::instance()->set_minimal( false ) ;

            GEO::FileSystem::delete_file( filename ) ;
        }

    }

    void GeoModelBuilderGM::unzip_one_file(
        unzFile& uz,
        const char filename[MAX_FILENAME] )
    {
        unzLocateFile( uz, filename, 0 ) ;
        char read_buffer[ READ_SIZE] ;

        if( unzOpenCurrentFile( uz ) != UNZ_OK ) {
            unzClose( uz ) ;
            throw RINGMeshException( "ZLIB", "Could not open file" ) ;
        }
        FILE *out = fopen( filename, "wb" ) ;
        if( out == NULL ) {
            unzCloseCurrentFile( uz ) ;
            unzClose( uz ) ;
            throw RINGMeshException( "ZLIB", "Could not open destination file" ) ;
        }
        int error = UNZ_OK ;
        do {
            error = unzReadCurrentFile( uz, read_buffer, READ_SIZE ) ;
            if( error < 0 ) {
                unzCloseCurrentFile( uz ) ;
                unzClose( uz ) ;
                fclose( out ) ;
                throw RINGMeshException( "ZLIB",
                    "Invalid error: " + GEO::String::to_string( error ) ) ;
            }
            if( error > 0 ) {
                fwrite( read_buffer, error, 1, out ) ;
            }
        } while( error > 0 ) ;
        fclose( out ) ;
        unzCloseCurrentFile( uz ) ;

    }

} // namespace
