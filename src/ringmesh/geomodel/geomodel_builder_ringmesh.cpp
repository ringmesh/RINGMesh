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

#include <ringmesh/geomodel/geomodel_builder_ringmesh.h>

#include <geogram/basic/file_system.h>

#include <ringmesh/io/io.h>
#include <ringmesh/geomodel/geomodel_builder_repair.h>
#include <ringmesh/mesh/geogram_mesh.h>
#include <ringmesh/mesh/geogram_mesh_builder.h>

/*!
 * @file ringmesh/geomodel/geomodel_builder_ringmesh.cpp
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
}

namespace RINGMesh {

    class GeoModelBuilderGMImpl {
    public:
        GeoModelBuilderGMImpl( GeoModelBuilderGM& builder )
            : builder_( builder )
        {
        }
        virtual ~GeoModelBuilderGMImpl()
        {
        }

        virtual void read_mesh_entity_line( GEO::LineInput& file_line ) = 0 ;

    protected:
        GeoModelBuilderGM& builder_ ;
    } ;

    class GeoModelBuilderGMImpl_0: public GeoModelBuilderGMImpl {
    public:
        GeoModelBuilderGMImpl_0( GeoModelBuilderGM& builder )
            : GeoModelBuilderGMImpl( builder )
        {
        }
        virtual ~GeoModelBuilderGMImpl_0()
        {
        }

        virtual void read_mesh_entity_line( GEO::LineInput& file_line )
        {
            // First line : type - id - name - geol_feature
            if( file_line.nb_fields() < 4 ) {
                throw RINGMeshException( "I/O",
                    "Invalid line: "
                        + GEO::String::to_string( file_line.line_number() )
                        + ", 4 fields are expected, the type, id, name, and geological feature" ) ;
            }
            gmme_t entity ;
            read_first_line( file_line, entity ) ;
            read_second_line( file_line, entity ) ;
        }

    protected:
        void read_first_line( GEO::LineInput& file_line, gmme_t& entity )
        {
            entity.type() = file_line.field( 0 ) ;
            entity.index() = file_line.field_as_uint( 1 ) ;
            builder_.info.set_mesh_entity_name( entity, file_line.field( 2 ) ) ;
            builder_.geology.set_mesh_entity_geol_feature( entity,
                GME::determine_geological_type( file_line.field( 3 ) ) ) ;
        }
        void read_second_line( GEO::LineInput& file_line, const gmme_t& entity )
        {
            file_line.get_line() ;
            file_line.get_fields() ;
            if( MeshEntityTypeManager::is_region( entity.type() ) ) {
                // Second line : signed indices of boundaries
                for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                    bool side = false ;
                    if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                        side = true ;
                    }
                    index_t s = NO_ID ;
                    GEO::String::from_string( &file_line.field( c )[1], s ) ;

                    builder_.topology.add_mesh_entity_boundary( entity, s, side ) ;
                }
            } else {
                // Second line : indices of its in boundaries
                for( index_t c = 1; c < file_line.nb_fields(); c++ ) {
                    builder_.topology.add_mesh_entity_boundary( entity,
                        file_line.field_as_uint( c ) ) ;
                }
            }
        }
    } ;

    class GeoModelBuilderGMImpl_1: public GeoModelBuilderGMImpl_0 {
    public:
        GeoModelBuilderGMImpl_1( GeoModelBuilderGM& builder )
            : GeoModelBuilderGMImpl_0( builder )
        {
        }
        virtual ~GeoModelBuilderGMImpl_1()
        {
        }

        virtual void read_mesh_entity_line( GEO::LineInput& file_line )
        {
            // Read this entity
            // First line : type - id - name - geol_feature - mesh type
            if( file_line.nb_fields() < 5 ) {
                throw RINGMeshException( "I/O",
                    "Invalid line: "
                        + GEO::String::to_string( file_line.line_number() )
                        + ", 5 fields are expected, the type, id, name, "
                        + "geological feature, and mesh type" ) ;
            }
            gmme_t entity ;
            read_first_line( file_line, entity ) ;

            const std::string mesh_type = file_line.field( 4 ) ;
            builder_.geometry.change_mesh_data_structure( entity, mesh_type ) ;

            read_second_line( file_line, entity ) ;
        }
    } ;

    GeoModelBuilderGM::GeoModelBuilderGM(
        GeoModel& geomodel,
        const std::string& filename )
        : GeoModelBuilderFile( geomodel, filename ), file_version_( 0 )
    {
        version_impl_[0] = new GeoModelBuilderGMImpl_0( *this ) ;
        version_impl_[1] = new GeoModelBuilderGMImpl_1( *this ) ;
    }

    GeoModelBuilderGM::~GeoModelBuilderGM()
    {
        for( index_t i = 0; i < NB_VERSION; i++ ) {
            delete version_impl_[i] ;
        }
    }

    void GeoModelBuilderGM::load_mesh_entities( const std::string& mesh_entity_file )
    {
        GEO::LineInput file_line( mesh_entity_file ) ;
        while( !file_line.eof() && file_line.get_line() ) {

            file_line.get_fields() ;
            if( file_line.nb_fields() > 0 ) {
                if( file_line.field_matches( 0, "Version" ) ) {
                    file_version_ = file_line.field_as_uint( 1 ) ;
                }
                // Name of the geomodel
                else if( file_line.field_matches( 0, "GeoModel" ) ) {
                    if( file_line.nb_fields() > 2 ) {
                        info.set_geomodel_name( file_line.field( 2 ) ) ;
                    }
                }
                // Number of entities of a given type
                else if( file_line.field_matches( 0, "Nb" ) ) {
                    // Allocate the space
                    topology.create_mesh_entities(
                        MeshEntityType( file_line.field( 1 ) ),
                        file_line.field_as_uint( 2 ) ) ;
                }
                // Mesh entities
                else if( match_mesh_entity_type( file_line.field( 0 ) ) ) {
                    version_impl_[file_version_]->read_mesh_entity_line(
                        file_line ) ;
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

                        topology.add_universe_boundary( s, side ) ;
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
        unzip_file( uz, mesh_entity_file.c_str() ) ;
        load_mesh_entities( mesh_entity_file ) ;
        bool ok = GEO::FileSystem::delete_file( mesh_entity_file ) ;
        ringmesh_unused( ok ) ; // avoids warning in release
        ringmesh_assert( ok ) ;
        load_meshes( uz ) ;

        const std::string geological_entity_file( "geological_entities.txt" ) ;
        unzip_file( uz, geological_entity_file.c_str() ) ;
        load_geological_entities( geological_entity_file ) ;
        ok = GEO::FileSystem::delete_file( geological_entity_file ) ;
        ringmesh_assert( ok ) ;

        unzClose( uz ) ;
    }

    void GeoModelBuilderGM::load_geological_entities(
        const std::string& geological_entity_file )
    {
        GEO::LineInput file_line( geological_entity_file ) ;
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
                    geology.create_geological_entities(
                        GeologicalEntityType( file_line.field( 1 ) ),
                        file_line.field_as_uint( 2 ) ) ;
                } else {
                    const std::string type = file_line.field( 0 ) ;
                    index_t id = file_line.field_as_uint( 1 ) ;
                    gmge_t entity( type, id ) ;
                    info.set_geological_entity_name( entity, file_line.field( 2 ) ) ;
                    geology.set_geological_entity_geol_feature( entity,
                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;
                    file_line.get_line() ;
                    file_line.get_fields() ;
                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                        geology.add_geological_entity_child( entity,
                            file_line.field_as_uint( in_b ) ) ;
                    }
                }
            }
        }
    }

    void GeoModelBuilderGM::load_meshes( unzFile& uz )
    {
        if( unzGoToFirstFile( uz ) != UNZ_OK ) {
            throw RINGMeshException( "I/O", "Unable to uncompress the first file" ) ;
        }
        std::vector< std::string > filenames ;
        do {
            char char_file_name[MAX_FILENAME] ;
            if( unzGetCurrentFileInfo64( uz, NULL, char_file_name,
            MAX_FILENAME, NULL, 0, NULL, 0 ) != UNZ_OK ) {
                throw RINGMeshException( "I/O", "Unable to get file name" ) ;
            }
            std::string file_name( char_file_name ) ;
            if( GEO::FileSystem::extension( file_name ) == "txt" ) {
                continue ;
            }

            unzip_current_file( uz, file_name.c_str() ) ;
            filenames.push_back( file_name ) ;
        } while( unzGoToNextFile( uz ) == UNZ_OK ) ;

        Logger::instance()->set_minimal( true ) ;
        RINGMESH_PARALLEL_LOOP_DYNAMIC
        for( index_t i = 0; i < filenames.size(); i++ ) {
            const std::string& file_name = filenames[i] ;
            std::string file_without_extension = GEO::FileSystem::base_name(
                file_name ) ;
            std::string entity_type, entity_id ;
            GEO::String::split_string( file_without_extension, '_', entity_type,
                entity_id ) ;
            index_t id = NO_ID ;
            GEO::String::from_string( entity_id, id ) ;
            if( MeshEntityTypeManager::is_corner( entity_type ) ) {
                Mesh0DBuilder_var builder = geometry.create_corner_builder( id ) ;
                builder->load_mesh( file_name ) ;
            } else if( MeshEntityTypeManager::is_line( entity_type ) ) {
                Mesh1DBuilder_var builder = geometry.create_line_builder( id ) ;
                builder->load_mesh( file_name ) ;
            } else if( MeshEntityTypeManager::is_surface( entity_type ) ) {
                Mesh2DBuilder_var builder = geometry.create_surface_builder( id ) ;
                builder->load_mesh( file_name ) ;
            } else if( MeshEntityTypeManager::is_region( entity_type ) ) {
                Mesh3DBuilder_var builder = geometry.create_region_builder( id ) ;
                builder->load_mesh( file_name ) ;
            }
            GEO::FileSystem::delete_file( file_name ) ;
        }
        Logger::instance()->set_minimal( false ) ;
    }

//    // ------------------------------------------------------------------------//
//
//    std::string OldGeoModelBuilderGM::match_nb_entities( const char* s ) const
//    {
//        // Check that the first 3 characters are NB_
//        if( strncmp( s, "NB_", 3 ) != 0 ) {
//
//        } else {
//            std::string old_type_name = std::string( s ).substr( 3 ) ;
//            if(MeshEntityTypeManager::is_valid_type(old_type_name)) {
//                return mesh_type_name_old_to_new( old_type_name ) ;
//            }
//            else {
//                return geological_type_name_old_to_new( old_type_name ) ;
//
//            }
//        }
//    }
//
//    MeshEntityType OldGeoModelBuilderGM::mesh_type_name_old_to_new(
//        const std::string& old_type_name ) const
//    {
//        if( old_type_name == "CORNER" ) {
//            return Corner::type_name_static() ;
//        } else if( old_type_name == "LINE" ) {
//            return Line::type_name_static() ;
//        } else if( old_type_name == "SURFACE" ) {
//            return Surface::type_name_static() ;
//        } else if( old_type_name == "REGION" ) {
//            return Region::type_name_static() ;
//        }
//        return DefaultMeshEntityType::default_entity_type() ;
//    }
//
//    GeologicalEntityType OldGeoModelBuilderGM::geological_type_name_old_to_new(
//        const std::string& old_type_name ) const
//    {
//        if( old_type_name == "CONTACT" ) {
//            return Contact::type_name_static() ;
//        } else if( old_type_name == "INTERFACE" ) {
//            return Interface::type_name_static() ;
//        } else if( old_type_name == "LAYER" ) {
//            return Layer::type_name_static() ;
//        }
//        return DefaultGeologicalEntityType::default_entity_type() ;
//    }
//
//    bool OldGeoModelBuilderGM::child_allowed( const char* s ) const
//    {
//        return geomodel_.entity_type_manager().is_geological_entity_type(
//            type_name_old_to_new( s ) ) ;
//    }
//
//    void OldGeoModelBuilderGM::load_topology( GEO::LineInput& file_line )
//    {
//        while( !file_line.eof() && file_line.get_line() ) {
//
//            file_line.get_fields() ;
//            if( file_line.nb_fields() > 0 ) {
//                // Name of the geomodel
//                if( file_line.field_matches( 0, "NAME" ) ) {
//                    if( file_line.nb_fields() > 1 ) {
//                        info.set_geomodel_name( file_line.field( 1 ) ) ;
//                    }
//                }
//                // Number of entities of a given type
//                else if( match_nb_entities( file_line.field( 0 ) )
//                    != GeoModelEntity::type_name_static() ) {
//                    // Allocate the space
//                    if( file_line.nb_fields() > 1 ) {
//                        EntityType type = match_nb_entities( file_line.field( 0 ) ) ;
//                        index_t nb_entities = file_line.field_as_uint( 1 ) ;
//                        if( geomodel_.is_mesh_entity_type( type ) ) {
//                            topology.create_mesh_entities( type, nb_entities ) ;
//                        } else {
//                            geology.create_geological_entities( type, nb_entities ) ;
//                        }
//
//                    }
//                }
//                // High-level entities
//                else if( child_allowed( file_line.field( 0 ) ) ) {
//                    // Read this entity
//                    // First line : type - id - name - geol_feature
//                    if( file_line.nb_fields() < 4 ) {
//                        throw RINGMeshException( "I/O",
//                            "Invalid line: "
//                                + GEO::String::to_string( file_line.line_number() )
//                                + ", 4 fields are expected, the type, id, name, and geological feature" ) ;
//                    }
//                    index_t id = file_line.field_as_uint( 1 ) ;
//                    gme_t entity( type_name_old_to_new( file_line.field( 0 ) ),
//                        id ) ;
//                    info.set_entity_name( entity, file_line.field( 2 ) ) ;
//                    geology.set_entity_geol_feature( entity,
//                        GME::determine_geological_type( file_line.field( 3 ) ) ) ;
//                    // Second line : indices of its children
//                    file_line.get_line() ;
//                    file_line.get_fields() ;
//                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
//
//                        geology.add_geological_entity_child( entity,
//                            file_line.field_as_uint( c ) ) ;
//                    }
//                }
//                // Regions
//                else if( type_name_old_to_new( file_line.field( 0 ) )
//                    == Region::type_name_static() ) {
//                    // First line : type - id - name
//                    if( file_line.nb_fields() < 3 ) {
//                        throw RINGMeshException( "I/O",
//                            "Invalid line: "
//                                + GEO::String::to_string( file_line.line_number() )
//                                + ", 3 fields are expected to describe a region: REGION, id, and name" ) ;
//                    }
//                    index_t id = file_line.field_as_uint( 1 ) ;
//                    gme_t entity( Region::type_name_static(), id ) ;
//                    info.set_entity_name( entity, file_line.field( 2 ) ) ;
//                    // Second line : signed indices of boundaries
//                    file_line.get_line() ;
//                    file_line.get_fields() ;
//                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
//                        bool side = false ;
//                        if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
//                            side = true ;
//                        }
//                        index_t s = NO_ID ;
//                        GEO::String::from_string( &file_line.field( c )[1], s ) ;
//                        ringmesh_assert( s != NO_ID ) ;
//
//                        topology.add_mesh_entity_boundary( entity, s, side ) ;
//                    }
//                }
//                // Universe
//                else if( file_line.field_matches( 0, "UNIVERSE" ) ) {
//                    // Second line: signed indices of boundaries
//                    file_line.get_line() ;
//                    file_line.get_fields() ;
//                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
//                        bool side = false ;
//                        if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
//                            side = true ;
//                        }
//                        index_t s = NO_ID ;
//                        GEO::String::from_string( &file_line.field( c )[1], s ) ;
//                        ringmesh_assert( s != NO_ID ) ;
//
//                        topology.add_universe_boundary( s, side ) ;
//                    }
//                }
//            }
//        }
//    }
//
//    void OldGeoModelBuilderGM::load_connectivities( GEO::LineInput& file_line )
//    {
//        while( !file_line.eof() && file_line.get_line() ) {
//            file_line.get_fields() ;
//            if( file_line.nb_fields() > 0 ) {
//                if( file_line.field_matches( 0, "GME" ) ) {
//                    const std::string old_name_type = file_line.field( 1 ) ;
//                    index_t id = file_line.field_as_uint( 2 ) ;
//                    file_line.get_line() ;
//                    file_line.get_fields() ;
//                    gme_t cur_gme_type( type_name_old_to_new( old_name_type ), id ) ;
//                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
//                        topology.add_mesh_entity_in_boundary( cur_gme_type,
//                            file_line.field_as_uint( in_b ) ) ;
//                    }
//                }
//            }
//        }
//    }
//
//    void OldGeoModelBuilderGM::load_entities(
//        const std::string& old_type_name,
//        unzFile& uz )
//    {
//        for( index_t el = 0;
//            el < geomodel_.nb_mesh_entities( type_name_old_to_new( old_type_name ) );
//            el++ ) {
//            std::string file_to_extract_and_load = old_type_name + "_"
//                + GEO::String::to_string( el ) ;
//            std::string str_try = file_to_extract_and_load + ".geogram" ;
//            if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
//                str_try = file_to_extract_and_load + ".meshb" ;
//                if( unzLocateFile( uz, str_try.c_str(), 0 ) != UNZ_OK ) {
//                    if( type_name_old_to_new( old_type_name )
//                        != Region::type_name_static() ) {
//                        std::string message = "Invalid format of .gm file" ;
//                        message += "\n.geogram file (defining mesh) is missing." ;
//                        throw RINGMeshException( "I/O", message ) ;
//                    }
//                    return ; // a region is not necessary meshed.
//                } else {
//                    std::string message =
//                        "Warning! you are using an old file (*.gm). \n" ;
//                    message += "Please use ringmeshconvert to update this file. \n" ;
//                    message +=
//                        "ringmeshconvert in:geomodel=old_geomodel.gm out:geomodel=new_geomodel.gm" ;
//                    GEO::Logger::warn( "I/O" ) << message << std::endl ;
//                }
//            }
//            unzip_file( uz, str_try.c_str() ) ;
//            GeogramMeshAllD cur_mesh ;
//            GEO::Logger::instance()->set_minimal( true ) ;
//            GeogramMeshAllDBuilder builder ;
//            builder.set_mesh( cur_mesh ) ;
//            builder.load_mesh( str_try ) ;
//            geometry.assign_mesh_to_entity( cur_mesh,
//                geomodel_.mesh_entity( type_name_old_to_new( old_type_name ), el ).gme_id() ) ;
//            GEO::Logger::instance()->set_minimal( false ) ;
//
//            unzip_file( uz, str_try.c_str() ) ;
//
//            GEO::FileSystem::delete_file( str_try ) ;
//        }
//
//    }
//
//    void OldGeoModelBuilderGM::load_file()
//    {
//        unzFile uz = unzOpen( filename_.c_str() ) ;
//        unz_global_info global_info ;
//        if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
//            unzClose( uz ) ;
//            throw RINGMeshException( "ZLIB", "Could not read file global info" ) ;
//        }
//
//        std::string topology_filename = "topology.txt" ;
//        unzip_file( uz, topology_filename.c_str() ) ;
//
//        GEO::LineInput line_topo( topology_filename ) ;
//
//        load_topology( line_topo ) ;
//        GEO::FileSystem::delete_file( topology_filename ) ;
//
//        load_entities( "CORNER", uz ) ;
//        load_entities( "LINE", uz ) ;
//        load_entities( "SURFACE", uz ) ;
//        load_entities( "REGION", uz ) ;
//
//        std::string connectivity = "connectivity.txt" ;
//        unzip_file( uz, connectivity.c_str() ) ;
//
//        GEO::LineInput line_connectivity( connectivity ) ;
//        load_connectivities( line_connectivity ) ;
//        GEO::FileSystem::delete_file( connectivity ) ;
//
//        // Repair line boundary order.
//        topology.complete_entity_connectivity() ;
//        repair.repair( GeoModelBuilderRepair::LINE_BOUNDARY_ORDER ) ;
//
//        unzClose( uz ) ;
//    }

} // namespace
