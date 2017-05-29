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
    using namespace RINGMesh;

    bool match_mesh_entity_type( const MeshEntityType& type )
    {
        if( type == Corner::type_name_static() ) return true;
        if( type == Line::type_name_static() ) return true;
        if( type == Surface::type_name_static() ) return true;
        if( type == Region::type_name_static() ) return true;
        return false;
    }
}

namespace RINGMesh {

    class GeoModelBuilderGMImpl {
    public:
        GeoModelBuilderGMImpl( GeoModelBuilderGM& builder )
            : builder_( builder )
        {
        }
        virtual ~GeoModelBuilderGMImpl() = default;

        virtual void read_mesh_entity_line( GEO::LineInput& file_line ) = 0;

    protected:
        GeoModelBuilderGM& builder_;
    };

    class GeoModelBuilderGMImpl_0: public GeoModelBuilderGMImpl {
    public:
        GeoModelBuilderGMImpl_0( GeoModelBuilderGM& builder )
            : GeoModelBuilderGMImpl( builder )
        {
        }
        virtual ~GeoModelBuilderGMImpl_0() = default;

        virtual void read_mesh_entity_line( GEO::LineInput& file_line ) override
        {
            // First line : type - id - name - geol_feature
            if( file_line.nb_fields() < 4 ) {
                throw RINGMeshException( "I/O",
                    "Invalid line: "
                        + GEO::String::to_string( file_line.line_number() )
                        + ", 4 fields are expected, the type, id, name, and geological feature" );
            }
            gmme_id entity = read_first_line( file_line );
            read_second_line( file_line, entity );
        }

    protected:
        virtual gmme_id read_first_line( GEO::LineInput& file_line )
        {

            gmme_id cur_gmme( MeshEntityType( file_line.field( 0 ) ),
                file_line.field_as_uint( 1 ) );
            builder_.info.set_mesh_entity_name( cur_gmme, file_line.field( 2 ) );
            GeoModelGeologicalEntity::GEOL_FEATURE not_used =
                GeoModelGeologicalEntity::determine_geological_type(
                    file_line.field( 3 ) );
            ringmesh_unused( not_used );
            return cur_gmme;
        }
        void read_second_line( GEO::LineInput& file_line, const gmme_id& entity )
        {
            file_line.get_line();
            file_line.get_fields();
            if( MeshEntityTypeManager::is_region( entity.type() ) ) {
                // Second line : signed indices of boundaries
                for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                    bool side = false;
                    if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                        side = true;
                    }
                    index_t s = NO_ID;
                    GEO::String::from_string( &file_line.field( c )[1], s );

                    builder_.topology.add_mesh_entity_boundary_relation( entity,
                        gmme_id( Surface::type_name_static(), s ), side );
                }
            } else {
                MeshEntityType type = MeshEntityTypeManager::boundary_type(
                    entity.type() );
                // Second line : indices of boundaries
                for( index_t c = 1; c < file_line.nb_fields(); c++ ) {
                    gmme_id boundary( type, file_line.field_as_uint( c ) );
                    builder_.topology.add_mesh_entity_boundary_relation( entity,
                        boundary );
                }
            }
        }
    };

    class GeoModelBuilderGMImpl_1: public GeoModelBuilderGMImpl_0 {
    public:
        GeoModelBuilderGMImpl_1( GeoModelBuilderGM& builder )
            : GeoModelBuilderGMImpl_0( builder )
        {
        }
        virtual ~GeoModelBuilderGMImpl_1() = default;

        virtual void read_mesh_entity_line( GEO::LineInput& file_line ) override
        {
            // Read this entity
            // First line : type - id - name - geol_feature - mesh type
            if( file_line.nb_fields() < 5 ) {
                throw RINGMeshException( "I/O",
                    "Invalid line: "
                        + GEO::String::to_string( file_line.line_number() )
                        + ", 5 fields are expected, the type, id, name, "
                        + "geological feature, and mesh type" );
            }
            gmme_id entity = read_first_line( file_line );

            const std::string mesh_type = file_line.field( 4 );
            if( GEO::String::string_starts_with( mesh_type, "Geogram" ) ) {
                builder_.geometry.change_mesh_data_structure( entity,
                    old_2_new_name( mesh_type ) );

            } else {
                builder_.geometry.change_mesh_data_structure( entity, mesh_type );

            }

            read_second_line( file_line, entity );

        }
    private:
        const std::string& old_2_new_name( const std::string& old_name )
        {

            index_t new_name_pos = GEO::String::to_int(
                GEO::String::to_string( old_name.at( old_name.length() - 2 ) ) );
            return new_names[new_name_pos];
        }
    private:
        static const std::string new_names[4];
    }
    ;
    const std::string GeoModelBuilderGMImpl_1::new_names[4] = {
        std::string( "GeogramPointSetMesh" ), std::string( "GeogramLineMesh" ),
        std::string( "GeogramSurfaceMesh" ), std::string( "GeogramVolumeMesh" ) };
    class GeoModelBuilderGMImpl_2: public GeoModelBuilderGMImpl_1 {
    public:
        GeoModelBuilderGMImpl_2( GeoModelBuilderGM& builder )
            : GeoModelBuilderGMImpl_1( builder )
        {
        }
        virtual ~GeoModelBuilderGMImpl_2() = default;

        virtual void read_mesh_entity_line( GEO::LineInput& file_line ) override
        {
            // Read this entity
            // First line : type - id - name  - mesh type
            if( file_line.nb_fields() < 4 ) {
                throw RINGMeshException( "I/O",
                    "Invalid line: "
                        + GEO::String::to_string( file_line.line_number() )
                        + ", 4 fields are expected, the type, id, name, "
                        + "geological feature, and mesh type" );
            }
            gmme_id entity = read_first_line( file_line );

            const std::string mesh_type = file_line.field( 3 );
            builder_.geometry.change_mesh_data_structure( entity, mesh_type );

            read_second_line( file_line, entity );
        }
    };

    GeoModelBuilderGM::GeoModelBuilderGM(
        GeoModel& geomodel,
        const std::string& filename )
        : GeoModelBuilderFile( geomodel, filename ), file_version_( 0 )
    {
        version_impl_[0].reset( new GeoModelBuilderGMImpl_0( *this ) );
        version_impl_[1].reset( new GeoModelBuilderGMImpl_1( *this ) );
        version_impl_[2].reset( new GeoModelBuilderGMImpl_2( *this ) );
    }

    GeoModelBuilderGM::~GeoModelBuilderGM()
    {
    }

    void GeoModelBuilderGM::load_mesh_entities( const std::string& mesh_entity_file )
    {
        GEO::LineInput file_line( mesh_entity_file );
        while( !file_line.eof() && file_line.get_line() ) {

            file_line.get_fields();
            if( file_line.nb_fields() > 0 ) {
                if( file_line.field_matches( 0, "Version" ) ) {
                    file_version_ = file_line.field_as_uint( 1 );
                }
                // Name of the geomodel
                else if( file_line.field_matches( 0, "GeoModel" ) ) {
                    if( file_line.nb_fields() > 2 ) {
                        info.set_geomodel_name( file_line.field( 2 ) );
                    }
                }
                // Number of entities of a given type
                else if( file_line.field_matches( 0, "Nb" ) ) {
                    // Allocate the space
                    topology.create_mesh_entities(
                        MeshEntityType( file_line.field( 1 ) ),
                        file_line.field_as_uint( 2 ) );
                }
                // Mesh entities
                else if( match_mesh_entity_type(
                    MeshEntityType( file_line.field( 0 ) ) ) ) {
                    version_impl_[file_version_]->read_mesh_entity_line( file_line );
                }
                // Universe
                else if( file_line.field_matches( 0, "Universe" ) ) {
                    // Second line: signed indices of boundaries
                    file_line.get_line();
                    file_line.get_fields();
                    for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                        bool side = false;
                        if( strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                            side = true;
                        }
                        index_t s = NO_ID;
                        GEO::String::from_string( &file_line.field( c )[1], s );

                        topology.add_universe_boundary( s, side );
                    }
                }
            }
        }
    }
    void GeoModelBuilderGM::load_file()
    {
        unzFile uz = unzOpen( filename_.c_str() );
        unz_global_info global_info;
        if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
            unzClose( uz );
            throw RINGMeshException( "ZLIB", "Could not read file global info" );
        }

        const std::string mesh_entity_file( "mesh_entities.txt" );
        unzip_file( uz, mesh_entity_file.c_str() );
        load_mesh_entities( mesh_entity_file );
        bool ok = GEO::FileSystem::delete_file( mesh_entity_file );
        ringmesh_unused( ok ); // avoids warning in release
        ringmesh_assert( ok );
        load_meshes( uz );

        const std::string geological_entity_file( "geological_entities.txt" );
        unzip_file( uz, geological_entity_file.c_str() );
        load_geological_entities( geological_entity_file );
        ok = GEO::FileSystem::delete_file( geological_entity_file );
        ringmesh_assert( ok );

        unzClose( uz );
    }

    void GeoModelBuilderGM::load_geological_entities(
        const std::string& geological_entity_file )
    {
        GEO::LineInput file_line( geological_entity_file );
        while( !file_line.eof() && file_line.get_line() ) {
            file_line.get_fields();
            if( file_line.nb_fields() > 0 ) {
                // If there is no geological entity
                if( file_line.field_matches( 0, "No" )
                    && file_line.field_matches( 1, "geological" )
                    && file_line.field_matches( 2, "entity" ) ) {
                    return;
                }
                // Number of entities of a given type
                if( file_line.field_matches( 0, "Nb" ) ) {
                    // Allocate the space
                    geology.create_geological_entities(
                        GeologicalEntityType( file_line.field( 1 ) ),
                        file_line.field_as_uint( 2 ) );
                } else {
                    GeologicalEntityType type( file_line.field( 0 ) );
                    index_t id = file_line.field_as_uint( 1 );
                    gmge_id entity( type, id );
                    info.set_geological_entity_name( entity, file_line.field( 2 ) );
                    geology.set_geological_entity_geol_feature( entity,
                        GeoModelGeologicalEntity::determine_geological_type(
                            file_line.field( 3 ) ) );
                    file_line.get_line();
                    file_line.get_fields();
                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                        geology.add_parent_children_relation( entity,
                            gmme_id(
                                geomodel_.entity_type_manager().relationship_manager.child_type(
                                    type ), file_line.field_as_uint( in_b ) ) );
                    }
                }
            }
        }
    }

    void GeoModelBuilderGM::load_meshes( unzFile& uz )
    {
        if( unzGoToFirstFile( uz ) != UNZ_OK ) {
            throw RINGMeshException( "I/O", "Unable to uncompress the first file" );
        }
        std::vector< std::string > filenames;
        do {
            char char_file_name[MAX_FILENAME];
            if( unzGetCurrentFileInfo64( uz, nullptr, char_file_name,
            MAX_FILENAME, nullptr, 0, nullptr, 0 ) != UNZ_OK ) {
                throw RINGMeshException( "I/O", "Unable to get file name" );
            }
            std::string file_name( char_file_name );
            if( GEO::FileSystem::extension( file_name ) == "txt" ) {
                continue;
            }

            unzip_current_file( uz, file_name.c_str() );
            filenames.push_back( file_name );
        } while( unzGoToNextFile( uz ) == UNZ_OK );

        Logger::instance()->set_minimal( true );
        RINGMESH_PARALLEL_LOOP_DYNAMIC
        for( index_t i = 0; i < filenames.size(); i++ ) {
            const std::string& file_name = filenames[i];
            std::string file_without_extension = GEO::FileSystem::base_name(
                file_name );
            std::string entity_type, entity_id;
            GEO::String::split_string( file_without_extension, '_', entity_type,
                entity_id );
            index_t id = NO_ID;
            GEO::String::from_string( entity_id, id );
            if( MeshEntityTypeManager::is_corner( entity_type ) ) {
                std::unique_ptr< PointSetMeshBuilder > builder =
                    geometry.create_corner_builder( id );
                builder->load_mesh( file_name );
            } else if( MeshEntityTypeManager::is_line( entity_type ) ) {
                std::unique_ptr< LineMeshBuilder > builder =
                    geometry.create_line_builder( id );
                builder->load_mesh( file_name );
            } else if( MeshEntityTypeManager::is_surface( entity_type ) ) {
                std::unique_ptr< SurfaceMeshBuilder > builder =
                    geometry.create_surface_builder( id );
                builder->load_mesh( file_name );
            } else if( MeshEntityTypeManager::is_region( entity_type ) ) {
                std::unique_ptr< VolumeMeshBuilder > builder =
                    geometry.create_region_builder( id );
                builder->load_mesh( file_name );
            }
            GEO::FileSystem::delete_file( file_name );
        }
        Logger::instance()->set_minimal( false );
    }

} // namespace
