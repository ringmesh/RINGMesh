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

#include <ringmesh/geomodel/geo_model_builder_ringmesh.h>

#include <ringmesh/geomodel/geo_model_repair.h>
#include <geogram/basic/file_system.h>

/*!
 * @file ringmesh/geomodel/geo_model_builder_ringmesh.cpp
 */

#define MAX_FILENAME 512
#define READ_SIZE 8192

namespace {
    using namespace RINGMesh ;

    bool match_mesh_entity_type( const std::string& type )
    {
        if( type == Corner::type_name_static() ) return true ;
        if( type == Line::type_name_static() ) return true ;
        if( type == Surface::type_name_static() ) return true ;
        if( type == Region::type_name_static() ) return true ;
        return false ;
    }
    void build_string_for_geo_model_entity_export( gme_t id, std::string& name )
    {
        name += id.type + "_" + GEO::String::to_string( id.index ) ;
    }

    /*!
     * @brief Unzip a file in a zip file and set it to the current unZIP file
     */
    void unzip_one_file( unzFile& uz, const char filename[MAX_FILENAME] )
    {
        unzLocateFile( uz, filename, 0 ) ;
        char read_buffer[READ_SIZE] ;

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
}

namespace RINGMesh {

    void GeoModelBuilderGM::load_mesh_entities( GEO::LineInput& file_line )
    {
        while( !file_line.eof() && file_line.get_line() ) {

            file_line.get_fields() ;
            if( file_line.nb_fields() > 0 ) {
                // Name of the model
                if( file_line.field_matches( 0, "GeoModel" ) ) {
                    if( file_line.nb_fields() > 2 ) {
                        set_model_name( file_line.field( 2 ) ) ;
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
                    gme_t entity( type, id ) ;
                    set_entity_name( entity, file_line.field( 2 ) ) ;
                    set_entity_geol_feature( entity,
                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;

                    // Read second line
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    if( type == Region::type_name_static() ) {
                        // Second line : signed indices of boundaries
                        for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                            bool side = false ;
                            if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                                side = true ;
                            }
                            index_t s = NO_ID ;
                            GEO::String::from_string( &file_line.field( c )[1], s ) ;

                            add_mesh_entity_boundary( entity, s, side ) ;
                        }
                    } else {
                        // Second line : indices of its in boundaries
                        for( index_t c = 1; c < file_line.nb_fields(); c++ ) {
                            add_mesh_entity_boundary( entity,
                                file_line.field_as_uint( c ) ) ;
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
                        index_t s = NO_ID ;
                        GEO::String::from_string( &file_line.field( c )[1], s ) ;

                        add_universe_boundary( s, side ) ;
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

        load_meshes( Corner::type_name_static(), uz ) ;
        load_meshes( Line::type_name_static(), uz ) ;
        load_meshes( Surface::type_name_static(), uz ) ;
        load_meshes( Region::type_name_static(), uz ) ;

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
                // If there is no geological entity
                if( file_line.field_matches( 0, "No" )
                    && file_line.field_matches( 1, "geological" )
                    && file_line.field_matches( 2, "entity" ) ) {
                    return ;
                }
                // Number of entities of a given type
                if( file_line.field_matches( 0, "Nb" ) ) {
                    // Allocate the space
                    create_geological_entities( file_line.field( 1 ),
                        file_line.field_as_uint( 2 ) ) ;
                } else {
                    const std::string type = file_line.field( 0 ) ;
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t entity( type, id ) ;
                    set_entity_name( entity, file_line.field( 2 ) ) ;
                    set_entity_geol_feature( entity,
                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                        add_geological_entity_child( entity,
                            file_line.field_as_uint( in_b ) ) ;
                    }
                }
            }
        }
    }

    void GeoModelBuilderGM::load_meshes( const std::string& type, unzFile& uz )
    {
        for( index_t el = 0; el < model().nb_mesh_entities( type ); el++ ) {
            gme_t cur_gme( type, el ) ;
            std::string file_to_extract_and_load ;
            build_string_for_geo_model_entity_export( cur_gme,
                file_to_extract_and_load ) ;
            std::string filename = file_to_extract_and_load + ".geogram" ;
            if( unzLocateFile( uz, filename.c_str(), 0 ) != UNZ_OK ) {
                if( type != Region::type_name_static() ) {
                    std::string message = "Invalid format of .gm file" ;
                    message += "\n.geogram file (defining mesh) is missing." ;
                    throw RINGMeshException( "I/O", message ) ;
                }
                return ; // a region is not necessary meshed.
            }
            unzip_one_file( uz, filename.c_str() ) ;
            GeogramMeshAllD cur_mesh ;
            GEO::MeshIOFlags flags ;
            flags.set_attribute( GEO::MESH_ALL_ATTRIBUTES ) ;
            Logger::instance()->set_minimal( true ) ;
            GeogramMeshAllDBuilder builder( cur_mesh ) ;
            builder.load_mesh( filename, flags ) ;
            assign_mesh_to_entity( cur_mesh, cur_gme ) ;
            Logger::instance()->set_minimal( false ) ;

            GEO::FileSystem::delete_file( filename ) ;
        }

    }

    // ------------------------------------------------------------------------//

    std::string OldGeoModelBuilderGM::match_nb_entities( const char* s ) const
    {
        // Check that the first 3 characters are NB_
        if( strncmp( s, "NB_", 3 ) != 0 ) {
            return GeoModelEntity::type_name_static() ;
        } else {
            std::string old_type_name = std::string( s ).substr( 3 ) ;
            return type_name_old_to_new( old_type_name ) ;
        }
    }

    EntityType OldGeoModelBuilderGM::type_name_old_to_new(
        const std::string& old_type_name ) const
    {
        if( old_type_name == "CORNER" ) {
            return Corner::type_name_static() ;
        } else if( old_type_name == "LINE" ) {
            return Line::type_name_static() ;
        } else if( old_type_name == "SURFACE" ) {
            return Surface::type_name_static() ;
        } else if( old_type_name == "REGION" ) {
            return Region::type_name_static() ;
        } else if( old_type_name == "CONTACT" ) {
            return Contact::type_name_static() ;
        } else if( old_type_name == "INTERFACE" ) {
            return Interface::type_name_static() ;
        } else if( old_type_name == "LAYER" ) {
            return Layer::type_name_static() ;
        }
        return GeoModelEntity::type_name_static() ;
    }

    bool OldGeoModelBuilderGM::child_allowed( const char* s ) const
    {
        return entity_type_manager().is_geological_entity_type(
            type_name_old_to_new( s ) ) ;
    }

    void OldGeoModelBuilderGM::load_topology( GEO::LineInput& file_line )
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
                    != GeoModelEntity::type_name_static() ) {
                    // Allocate the space
                    if( file_line.nb_fields() > 1 ) {
                        EntityType type = match_nb_entities( file_line.field( 0 ) ) ;
                        index_t nb_entities = file_line.field_as_uint( 1 ) ;
                        if( model().is_mesh_entity_type( type ) ) {
                            create_mesh_entities( type, nb_entities ) ;
                        } else {
                            create_geological_entities( type, nb_entities ) ;
                        }

                    }
                }
                // High-level entities
                else if( child_allowed( file_line.field( 0 ) ) ) {
                    // Read this entity
                    // First line : type - id - name - geol_feature
                    if( file_line.nb_fields() < 4 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line.line_number() )
                                + ", 4 fields are expected, the type, id, name, and geological feature" ) ;
                    }
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t entity( type_name_old_to_new( file_line.field( 0 ) ),
                        id ) ;
                    set_entity_name( entity, file_line.field( 2 ) ) ;
                    set_entity_geol_feature( entity,
                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;
                    // Second line : indices of its children
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {

                        add_geological_entity_child( entity,
                            file_line.field_as_uint( c ) ) ;
                    }
                }
                // Regions
                else if( type_name_old_to_new( file_line.field( 0 ) )
                    == Region::type_name_static() ) {
                    // First line : type - id - name
                    if( file_line.nb_fields() < 3 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line.line_number() )
                                + ", 3 fields are expected to describe a region: REGION, id, and name" ) ;
                    }
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t entity( Region::type_name_static(), id ) ;
                    set_entity_name( entity, file_line.field( 2 ) ) ;
                    // Second line : signed indices of boundaries
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s = NO_ID ;
                        GEO::String::from_string( &file_line.field( c )[1], s ) ;
                        ringmesh_assert( s != NO_ID ) ;

                        add_mesh_entity_boundary( entity, s, side ) ;
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
                        index_t s = NO_ID ;
                        GEO::String::from_string( &file_line.field( c )[1], s ) ;
                        ringmesh_assert( s != NO_ID ) ;

                        add_universe_boundary( s, side ) ;
                    }
                }
            }
        }
    }

    void OldGeoModelBuilderGM::load_connectivities( GEO::LineInput& file_line )
    {
        while( !file_line.eof() && file_line.get_line() ) {
            file_line.get_fields() ;
            if( file_line.nb_fields() > 0 ) {
                if( file_line.field_matches( 0, "GME" ) ) {
                    const std::string old_name_type = file_line.field( 1 ) ;
                    index_t id = file_line.field_as_uint( 2 ) ;
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    gme_t cur_gme_type( type_name_old_to_new( old_name_type ), id ) ;
                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                        add_mesh_entity_in_boundary( cur_gme_type,
                            file_line.field_as_uint( in_b ) ) ;
                    }
                }
            }
        }
    }

    void OldGeoModelBuilderGM::load_entities(
        const std::string& old_type_name,
        unzFile& uz )
    {
        for( index_t el = 0;
            el < model().nb_mesh_entities( type_name_old_to_new( old_type_name ) );
            el++ ) {
            std::string file_to_extract_and_load = old_type_name + "_"
                + GEO::String::to_string( el ) ;
            std::string str_try = file_to_extract_and_load + ".geogram" ;
            if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
                str_try = file_to_extract_and_load + ".meshb" ;
                if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
                    if( type_name_old_to_new( old_type_name )
                        != Region::type_name_static() ) {
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
            GeogramMeshAllD cur_mesh ;
            GEO::MeshIOFlags flags ;
            flags.set_attribute( GEO::MESH_ALL_ATTRIBUTES ) ;
            GEO::Logger::instance()->set_minimal( true ) ;
            GeogramMeshAllDBuilder builder( cur_mesh ) ;
            builder.load_mesh( str_try, flags ) ;
            assign_mesh_to_entity( cur_mesh,
                model().mesh_entity( type_name_old_to_new( old_type_name ), el ).gme_id() ) ;
            GEO::Logger::instance()->set_minimal( false ) ;

            unzip_one_file( uz, str_try.c_str() ) ;

            GEO::FileSystem::delete_file( str_try ) ;
        }

    }

    void OldGeoModelBuilderGM::load_file()
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

        load_entities( "CORNER", uz ) ;
        load_entities( "LINE", uz ) ;
        load_entities( "SURFACE", uz ) ;
        load_entities( "REGION", uz ) ;

        std::string connectivity = "connectivity.txt" ;
        unzip_one_file( uz, connectivity.c_str() ) ;

        GEO::LineInput line_connectivity( connectivity ) ;
        load_connectivities( line_connectivity ) ;
        GEO::FileSystem::delete_file( connectivity ) ;

        // Repair line boundary order.
        complete_entity_connectivity() ;
        GeoModelRepair repair( model() ) ;
        repair.repair( GeoModelRepair::LINE_BOUNDARY_ORDER ) ;

        unzClose( uz ) ;
    }

} // namespace
