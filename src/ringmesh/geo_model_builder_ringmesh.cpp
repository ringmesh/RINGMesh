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

namespace RINGMesh {

    void GeoModelBuilderGM::load_topology( GEO::LineInput& file_line )
    {
        while( !file_line.eof() && file_line.get_line() ) {

            file_line.get_fields() ;
            if( file_line.nb_fields() > 0 ) {
                // Name of the model
                if( file_line.field_matches( 0, "NAME" ) ) {
                    if( file_line.nb_fields() > 1 ) {
                        set_model_name( file_line.field( 1 ) ) ;
                    }
                }
                // Number of entities of a given type
                else if( match_nb_entities( file_line.field( 0 ) )
                    != GME::NO_TYPE ) {
                    // Allocate the space
                    if( file_line.nb_fields() > 1 ) {
                        create_entities( match_nb_entities( file_line.field( 0 ) ),
                            file_line.field_as_uint( 1 ) ) ;
                    }
                }

                // High-level entities
                else if( match_high_level_type( file_line.field( 0 ) ) ) {
                    // Read this entity
                    // First line : type - id - name - geol_feature
                    if( file_line.nb_fields() < 4 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line.line_number() )
                                + ", 4 fields are expected, the type, id, name, and geological feature" ) ;
                    }
                    GME::TYPE t = match_type( file_line.field( 0 ) ) ;
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t entity( t, id ) ;
                    set_entity_name( entity, file_line.field( 2 ) ) ;
                    set_entity_geol_feature( entity,
                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;
                    // Second line : indices of its children
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                        add_entity_child( entity,
                            gme_t( GME::child_type( t ),
                                file_line.field_as_uint( c ) ) ) ;
                    }
                }
                // Regions
                else if( match_type( file_line.field( 0 ) ) == GME::REGION ) {
                    // First line : type - id - name
                    if( file_line.nb_fields() < 3 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line.line_number() )
                                + ", 3 fields are expected to describe a region: REGION, id, and name" ) ;
                    }
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t entity( GME::REGION, id ) ;
                    set_entity_name( entity, file_line.field( 2 ) ) ;
                    // Second line : signed indices of boundaries
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &file_line.field( c )[1], s ) ;

                        add_entity_boundary( entity, gme_t( GME::SURFACE, s ),
                            side ) ;
                    }
                }

                // Universe
                else if( file_line.field_matches( 0, "UNIVERSE" ) ) {
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

                        add_entity_boundary( gme_t( GME::REGION, NO_ID ),
                            gme_t( GME::SURFACE, s ), side ) ;
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

        std::string topology = "topology.txt" ;
        unzip_one_file( uz, topology.c_str() ) ;

        GEO::LineInput line_topo( topology ) ;

        load_topology( line_topo ) ;
        GEO::FileSystem::delete_file( topology ) ;

        for( index_t t = GME::CORNER; t <= GME::REGION; t++ ) {
            GME::TYPE type = static_cast< GME::TYPE >( t ) ;
            load_entities( type, uz ) ;
        }

        std::string connectivity = "connectivity.txt" ;
        unzip_one_file( uz, connectivity.c_str() ) ;

        GEO::LineInput line_connectivity( connectivity ) ;
        load_connectivities( line_connectivity ) ;
        GEO::FileSystem::delete_file( connectivity ) ;

        unzClose( uz ) ;
    }

    void GeoModelBuilderGM::load_connectivities( GEO::LineInput& file_line )
    {
        while( !file_line.eof() && file_line.get_line() ) {
            file_line.get_fields() ;
            if( file_line.nb_fields() > 0 ) {
                if( file_line.field_matches( 0, "GME" ) ) {
                    GME::TYPE t = match_type( file_line.field( 1 ) ) ;
                    index_t id = file_line.field_as_uint( 2 ) ;
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    const GeoModelMeshEntity& cur_gme = model().mesh_entity( t,
                        id ) ;
                    gme_t cur_gme_type( t, id ) ;
                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                        add_entity_in_boundary( cur_gme_type,
                            gme_t( cur_gme.in_boundary_type( t ),
                                file_line.field_as_uint( in_b ) ) ) ;
                    }
                }

            }

        }

    }

    void GeoModelBuilderGM::load_entities( GME::TYPE gme_t, unzFile& uz )
    {
        for( index_t el = 0; el < model().nb_entities( gme_t ); el++ ) {
            std::string file_to_extract_and_load ;
            build_string_for_geo_model_entity_export( gme_t, el,
                file_to_extract_and_load ) ;
            std::string str_try = file_to_extract_and_load + ".geogram" ;
            if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
                str_try = file_to_extract_and_load + ".meshb" ;
                if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
                    if( gme_t != GME::REGION ) {
                        std::string message = "Invalid format of .gm file" ;
                        message += "\n.geogram file (defining mesh) is missing." ;
                        throw RINGMeshException( "I/O", message ) ;
                    }
                    return ; // a region is not necessary meshed.
                } else {
                    std::string message =
                        "Warning! you are using an old file (*.gm). \n" ;
                    message += "Please use ringmeshconvert to update this file. \n" ;
                    message +=
                        "ringmeshconvert in:geomodel=old_geomodel.gm out:geomodel=new_geomodel.gm" ;
                    GEO::Logger::warn( "I/O" ) << message << std::endl ;
                }
            }
            unzip_one_file( uz, str_try.c_str() ) ;
            Mesh cur_mesh( model(), 3, false ) ;
            GEO::MeshIOFlags flags ;
            flags.set_attribute( GEO::MESH_ALL_ATTRIBUTES ) ;
            GEO::Logger::instance()->set_minimal( true ) ;
            MeshBuilder builder(cur_mesh);
            builder.load_mesh( str_try, flags ) ;
            assign_mesh_to_entity( cur_mesh,
                model().entity( gme_t, el ).gme_id() ) ;
            GEO::Logger::instance()->set_minimal( false ) ;

            unzip_one_file( uz, str_try.c_str() ) ;

//            set_connectivities
            GEO::FileSystem::delete_file( str_try ) ;
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
