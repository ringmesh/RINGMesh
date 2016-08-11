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

#include <geogram/basic/file_system.h>


/*!
 * @file ringmesh/geomodel/geo_model_builder_ringmesh.cpp
 */

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

    /*!
     * @brief Load meshes of mesh entities of one type from a zip file
     * @param[in] gme_t the GeoModelMeshEntity type
     * @param[in] uz the zip file
     */
    /*void load_meshes( const std::string& type, unzFile& uz, GeoModelBuilderFile& geomodel_bf, GeoModel& geomodel )
    {
        for( index_t el = 0; el < geomodel.nb_mesh_entities( type ); el++ ) {
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
            Mesh cur_mesh( geomodel, 3, false ) ;
            GEO::MeshIOFlags flags ;
            flags.set_attribute( GEO::MESH_ALL_ATTRIBUTES ) ;
            Logger::instance()->set_minimal( true ) ;
            MeshBuilder builder( cur_mesh ) ;
            builder.load_mesh( filename, flags ) ;
            geomodel_bf.assign_mesh_to_entity( cur_mesh, cur_gme ) ;
            Logger::instance()->set_minimal( false ) ;

            GEO::FileSystem::delete_file( filename ) ;
        }

    }*/
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
                    create_entities( file_line.field( 1 ),
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

    /*void GeoModelBuilderGM::unzip_one_file(
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

    }*/

    OldGeoModelBuilderGM::TYPE OldGeoModelBuilderGM::match_nb_entities( const char* s )
    {
        // Check that the first 3 characters are NB_
        if( strncmp( s, "NB_", 3 ) != 0 ) {
            return NO_TYPE ;
        } else {
            for( index_t i = CORNER; i < NO_TYPE; i++ ) {
                TYPE type = (TYPE) i ;
                if( strstr( s, type_name( type ).data() ) != NULL ) {
                    return type ;
                }
            }
            return NO_TYPE ;
        }
    }

    std::string OldGeoModelBuilderGM::type_name( TYPE t )
    {
        switch( t ) {
            case CORNER:
                return "CORNER" ;
            case LINE:
                return "LINE" ;
            case SURFACE:
                return "SURFACE" ;
            case REGION:
                return "REGION" ;
            case CONTACT:
                return "CONTACT" ;
            case INTERFACE:
                return "INTERFACE" ;
            case LAYER:
                return "LAYER" ;
            default:
                return "NO_TYPE_NAME" ;
        }
    }

    EntityType OldGeoModelBuilderGM::old2new( TYPE type )
    {
        switch( type ) {
            case CORNER:
                return Corner::type_name_static() ;
            case LINE:
                return Line::type_name_static() ;
            case SURFACE:
                return Surface::type_name_static() ;
            case REGION:
                return Region::type_name_static() ;
            case CONTACT:
                return Contact::type_name_static() ;
            case INTERFACE:
                return Interface::type_name_static() ;
            case LAYER:
                return Layer::type_name_static() ;
            default:
                return GeoModelEntity::type_name_static() ;
        }
    }

    bool OldGeoModelBuilderGM::match_high_level_type( const char* s )
    {
        if( std::string(s) == "CONTACT" || std::string(s) == "INTERFACE" || std::string(s) == "LAYER" ) {
            return true ;
        }
        return false ;
    }

    OldGeoModelBuilderGM::TYPE OldGeoModelBuilderGM::match_type( const char* s )
    {
        for( index_t i = CORNER; i < NO_TYPE; i++ ) {
            TYPE type = (TYPE) i ;
            if( strcmp( s, type_name( type ).data() ) == 0 ) {
                return type ;
            }
        }
        return NO_TYPE ;
    }

    /*EntityType OldGeoModelBuilderGM::child_type( TYPE type ){
        switch( type ) {
            case CONTACT:
                return Line::type_name_static() ;
            case INTERFACE:
                return Surface::type_name_static() ;
            case LAYER:
                return Region::type_name_static() ;
            default:
                return GeoModelEntity::type_name_static() ;
        }
    }*/

    void OldGeoModelBuilderGM::load_topology( GEO::LineInput& file_line )
    {
        // To store the basic GeoModelGeologicalEntities in the manager.
        find_or_create_geological_entity_type( Contact::type_name_static() ) ;
        find_or_create_geological_entity_type( Interface::type_name_static() ) ;
        find_or_create_geological_entity_type( Layer::type_name_static() ) ;
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
                    != NO_TYPE ) {
                    // Allocate the space
                    if( file_line.nb_fields() > 1 ) {
//                        DEBUG(old2new(match_nb_entities( file_line.field( 0 ) ))) ;
//                        DEBUG(file_line.field_as_uint( 1 ));
//                        DEBUG(model().nb_geological_entities(Contact::type_name_static()));
                        create_entities( old2new(match_nb_entities( file_line.field( 0 ) )),
                            file_line.field_as_uint( 1 ) ) ;
//                        DEBUG(model().nb_geological_entities(Contact::type_name_static()));
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
                    TYPE t = match_type( file_line.field( 0 ) ) ;
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t entity( old2new(t), id ) ;
                    set_entity_name( entity, file_line.field( 2 ) ) ;
                    set_entity_geol_feature( entity,
                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;
                    // Second line : indices of its children
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {

//                        add_entity_child
                        add_geological_entity_child( entity,
                                file_line.field_as_uint( c )  ) ;
                    }
                }
                // Regions
                else if( match_type( file_line.field( 0 ) ) == REGION ) {
                    // First line : type - id - name
                    if( file_line.nb_fields() < 3 ) {
                        throw RINGMeshException( "I/O",
                            "Invalid line: "
                                + GEO::String::to_string( file_line.line_number() )
                                + ", 3 fields are expected to describe a region: REGION, id, and name" ) ;
                    }
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gme_t entity( old2new(REGION), id ) ;
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
                        index_t s ;
                        GEO::String::from_string( &file_line.field( c )[1], s ) ;

                        add_universe_boundary( s , side ) ;
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
                    TYPE t = match_type( file_line.field( 1 ) ) ;
                    index_t id = file_line.field_as_uint( 2 ) ;
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    /*const GeoModelMeshEntity& cur_gme = model().mesh_entity(
                        old2new( t ), id ) ;*/
                    gme_t cur_gme_type( old2new( t ), id ) ;
                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                        add_mesh_entity_in_boundary( cur_gme_type,
                            file_line.field_as_uint( in_b ) ) ;
                    }
                }

            }

        }

    }

    void OldGeoModelBuilderGM::load_entities( TYPE type, unzFile& uz )
    {
        for( index_t el = 0; el < model().nb_entities( old2new(type) ); el++ ) {
            std::string file_to_extract_and_load = type_name( type ) + "_"
                + GEO::String::to_string( el ) ;
            std::string str_try = file_to_extract_and_load + ".geogram" ;
            if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
                str_try = file_to_extract_and_load + ".meshb" ;
                if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
                    if( type != REGION ) {
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
            MeshBuilder builder( cur_mesh ) ;
            builder.load_mesh( str_try, flags ) ;
            assign_mesh_to_entity( cur_mesh,
                model().entity( old2new( type ), el ).gme_id() ) ;
            GEO::Logger::instance()->set_minimal( false ) ;

            unzip_one_file( uz, str_try.c_str() ) ;

            //            set_connectivities
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

        for( index_t t = CORNER; t <= REGION; t++ ) {
            TYPE type = static_cast< TYPE >( t ) ;
            load_entities( type, uz ) ;
        }

        std::string connectivity = "connectivity.txt" ;
        unzip_one_file( uz, connectivity.c_str() ) ;

        GEO::LineInput line_connectivity( connectivity ) ;
        load_connectivities( line_connectivity ) ;
        GEO::FileSystem::delete_file( connectivity ) ;

//        end_model() ;
        /*complete_mesh_entity_connectivity< Line >() ;
        for( index_t line_itr = 0; line_itr < model().nb_lines(); ++line_itr ) {
            const Line& cur_line = model().line( line_itr ) ;
            if( !cur_line.is_first_corner_first_vertex() ) {
                const index_t first_boundary_index = cur_line.boundary( 0 ).index() ;
                set_mesh_entity_boundary( cur_line.gme_id(), 0,
                    cur_line.boundary_gme( 1 ).index ) ;
                set_mesh_entity_boundary( cur_line.gme_id(), 1,
                    first_boundary_index ) ;
            }
        }*/

        unzClose( uz ) ;

    }

} // namespace
