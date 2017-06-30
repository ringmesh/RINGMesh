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
        if( type == Corner< 3 >::type_name_static() ) return true;
        if( type == Line< 3 >::type_name_static() ) return true;
        if( type == Surface< 3 >::type_name_static() ) return true;
        if( type == Region< 3 >::type_name_static() ) return true;
        return false;
    }
}

namespace RINGMesh {

    template< index_t DIMENSION >
    class GeoModelBuilderGMImpl {
    public:
        GeoModelBuilderGMImpl(
            GeoModelBuilderGM< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel )
            : builder_( builder ), geomodel_( geomodel )
        {
        }
        virtual ~GeoModelBuilderGMImpl() = default;

        virtual void read_mesh_entity_line( GEO::LineInput& file_line ) = 0;

    protected:
        GeoModelBuilderGM< DIMENSION >& builder_;
        GeoModel< DIMENSION >& geomodel_;
    };

    template< index_t DIMENSION >
    class GeoModelBuilderGMImpl_0: public GeoModelBuilderGMImpl< DIMENSION > {
    public:
        GeoModelBuilderGMImpl_0(
            GeoModelBuilderGM< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel )
            : GeoModelBuilderGMImpl< DIMENSION >( builder, geomodel )
        {
        }
        virtual ~GeoModelBuilderGMImpl_0() = default;

        virtual void read_mesh_entity_line( GEO::LineInput& file_line ) override
        {
            // First line : type - id - name - geol_feature
            if( file_line.nb_fields() < 4 ) {
                throw RINGMeshException( "I/O",
                    "Invalid line: " + std::to_string( file_line.line_number() )
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
            this->builder_.info.set_mesh_entity_name( cur_gmme,
                file_line.field( 2 ) );
            return cur_gmme;
        }
        void read_second_line( GEO::LineInput& file_line, const gmme_id& entity );
    private:
        template< template< index_t > class ENTITY >
        void add_relation_for_entities_with_sides(
            const gmme_id& entity,
            GEO::LineInput& file_line )
        {
            for( index_t c = 0; c < file_line.nb_fields(); c++ ) {
                bool side = false;
                if( std::strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                    side = true;
                }
                index_t s = NO_ID;
                GEO::String::from_string( &file_line.field( c )[1], s );

                this->builder_.topology.add_mesh_entity_boundary_relation( entity,
                    gmme_id( ENTITY< DIMENSION >::type_name_static(), s ), side );
            }
        }

        void add_relation_for_entities_with_no_side(
            const gmme_id& entity,
            GEO::LineInput& file_line )
        {
            const MeshEntityTypeManager< DIMENSION >& manager =
                this->geomodel_.entity_type_manager().mesh_entity_manager;
            MeshEntityType type = manager.boundary_entity_type( entity.type() );
            // Second line : indices of boundaries
            for( index_t c = 1; c < file_line.nb_fields(); c++ ) {
                gmme_id boundary( type, file_line.field_as_uint( c ) );
                this->builder_.topology.add_mesh_entity_boundary_relation( entity,
                    boundary );
            }
        }
    };

    template< >
    void GeoModelBuilderGMImpl_0< 2 >::read_second_line(
        GEO::LineInput& file_line,
        const gmme_id& entity )
    {
        file_line.get_line();
        file_line.get_fields();
        const MeshEntityTypeManager< 2 >& manager =
            this->geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_surface( entity.type() ) ) {
            add_relation_for_entities_with_sides< Line >( entity, file_line );
        } else {
            add_relation_for_entities_with_no_side( entity, file_line );
        }
    }

    template< >
    void GeoModelBuilderGMImpl_0< 3 >::read_second_line(
        GEO::LineInput& file_line,
        const gmme_id& entity )
    {
        file_line.get_line();
        file_line.get_fields();
        const MeshEntityTypeManager< 3 >& manager =
            this->geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_region( entity.type() ) ) {
            add_relation_for_entities_with_sides< Surface >( entity, file_line );
        } else {
            add_relation_for_entities_with_no_side( entity, file_line );
        }
    }

    template< index_t DIMENSION >
    class GeoModelBuilderGMImpl_1: public GeoModelBuilderGMImpl_0< DIMENSION > {
    public:
        GeoModelBuilderGMImpl_1(
            GeoModelBuilderGM< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel )
            : GeoModelBuilderGMImpl_0< DIMENSION >( builder, geomodel )
        {
        }
        virtual ~GeoModelBuilderGMImpl_1() = default;

        virtual void read_mesh_entity_line( GEO::LineInput& file_line ) override
        {
            // Read this entity
            // First line : type - id - name - geol_feature - mesh type
            if( file_line.nb_fields() < 5 ) {
                throw RINGMeshException( "I/O",
                    "Invalid line: " + std::to_string( file_line.line_number() )
                        + ", 5 fields are expected, the type, id, name, "
                        + "geological feature, and mesh type" );
            }
            gmme_id entity = this->read_first_line( file_line );

            const std::string mesh_type = file_line.field( 4 );
            if( GEO::String::string_starts_with( mesh_type, "Geogram" ) ) {
                this->builder_.geometry.change_mesh_data_structure( entity,
                    old_2_new_name( mesh_type ) );
            } else {
                this->builder_.geometry.change_mesh_data_structure( entity,
                    mesh_type );
            }

            this->read_second_line( file_line, entity );
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
    };

    template< index_t DIMENSION >
    const std::string GeoModelBuilderGMImpl_1< DIMENSION >::new_names[4] = {
        std::string( "GeogramPointSetMesh" ), std::string( "GeogramLineMesh" ),
        std::string( "GeogramSurfaceMesh" ), std::string( "GeogramVolumeMesh" ) };

    template< index_t DIMENSION >
    class GeoModelBuilderGMImpl_2: public GeoModelBuilderGMImpl_1< DIMENSION > {
    public:
        GeoModelBuilderGMImpl_2(
            GeoModelBuilderGM< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel )
            : GeoModelBuilderGMImpl_1< DIMENSION >( builder, geomodel )
        {
        }
        virtual ~GeoModelBuilderGMImpl_2() = default;

        virtual void read_mesh_entity_line( GEO::LineInput& file_line ) override
        {
            // Read this entity
            // First line : type - id - name  - mesh type
            if( file_line.nb_fields() < 4 ) {
                throw RINGMeshException( "I/O",
                    "Invalid line: " + std::to_string( file_line.line_number() )
                        + ", 4 fields are expected, the type, id, name, "
                        + "geological feature, and mesh type" );
            }
            gmme_id entity = this->read_first_line( file_line );

            const std::string mesh_type = file_line.field( 3 );
            this->builder_.geometry.change_mesh_data_structure( entity, mesh_type );

            this->read_second_line( file_line, entity );
        }
    };

    template< index_t DIMENSION >
    GeoModelBuilderGM< DIMENSION >::GeoModelBuilderGM(
        GeoModel< DIMENSION >& geomodel,
        std::string filename )
        :
            GeoModelBuilderFile< DIMENSION >( geomodel, std::move( filename ) ),
            file_version_( 0 )
    {
        version_impl_[0].reset(
            new GeoModelBuilderGMImpl_0< DIMENSION >( *this, geomodel ) );
        version_impl_[1].reset(
            new GeoModelBuilderGMImpl_1< DIMENSION >( *this, geomodel ) );
        version_impl_[2].reset(
            new GeoModelBuilderGMImpl_2< DIMENSION >( *this, geomodel ) );
    }

    template< index_t DIMENSION >
    GeoModelBuilderGM< DIMENSION >::~GeoModelBuilderGM()
    {
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGM< DIMENSION >::load_mesh_entities(
        const std::string& mesh_entity_file )
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
                        this->info.set_geomodel_name( file_line.field( 2 ) );
                    }
                }
                // Number of entities of a given type
                else if( file_line.field_matches( 0, "Nb" ) ) {
                    // Allocate the space
                    this->topology.create_mesh_entities(
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
                        if( std::strncmp( file_line.field( c ), "+", 1 ) == 0 ) {
                            side = true;
                        }
                        index_t s = NO_ID;
                        GEO::String::from_string( &file_line.field( c )[1], s );

                        this->topology.add_universe_boundary( s, side );
                    }
                }
            }
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGM< DIMENSION >::load_file()
    {
        unzFile uz = unzOpen( this->filename_.c_str() );
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

    template< index_t DIMENSION >
    void GeoModelBuilderGM< DIMENSION >::load_geological_entities(
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
                    this->geology.create_geological_entities(
                        GeologicalEntityType( file_line.field( 1 ) ),
                        file_line.field_as_uint( 2 ) );
                } else {
                    GeologicalEntityType type( file_line.field( 0 ) );
                    index_t id = file_line.field_as_uint( 1 );
                    gmge_id entity( type, id );
                    this->info.set_geological_entity_name( entity,
                        file_line.field( 2 ) );
                    this->geology.set_geological_entity_geol_feature( entity,
                        GeoModelGeologicalEntity< DIMENSION >::determine_geological_type(
                            file_line.field( 3 ) ) );
                    file_line.get_line();
                    file_line.get_fields();
                    for( index_t in_b = 0; in_b < file_line.nb_fields(); in_b++ ) {
                        this->geology.add_parent_children_relation( entity,
                            gmme_id(
                                this->geomodel_.entity_type_manager().relationship_manager.child_type(
                                    type ), file_line.field_as_uint( in_b ) ) );
                    }
                }
            }
        }
    }

    template< index_t DIMENSION >
    bool GeoModelBuilderGM< DIMENSION >::load_mesh_entity_base(
        const std::string& entity_type,
        const std::string& file_name,
        index_t id )
    {
        const MeshEntityTypeManager< DIMENSION >& manager =
            this->geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_corner( entity_type ) ) {
            std::unique_ptr< PointSetMeshBuilder< DIMENSION > > builder =
                this->geometry.create_corner_builder( id );
            builder->load_mesh( file_name );
            return true;
        } else if( manager.is_line( entity_type ) ) {
            std::unique_ptr< LineMeshBuilder< DIMENSION > > builder =
                this->geometry.create_line_builder( id );
            builder->load_mesh( file_name );
            return true;
        } else if( manager.is_surface( entity_type ) ) {
            std::unique_ptr< SurfaceMeshBuilder< DIMENSION > > builder =
                this->geometry.create_surface_builder( id );
            builder->load_mesh( file_name );
            return true;
        }
        return false;
    }

    template< >
    void GeoModelBuilderGM< 2 >::load_mesh_entity(
        const std::string& entity_type,
        const std::string& file_name,
        index_t id )
    {
        load_mesh_entity_base( entity_type, file_name, id );
    }

    template< >
    void GeoModelBuilderGM< 3 >::load_mesh_entity(
        const std::string& entity_type,
        const std::string& file_name,
        index_t id )
    {
        if( !load_mesh_entity_base( entity_type, file_name, id ) ) {
            const MeshEntityTypeManager< 3 >& manager =
                this->geomodel_.entity_type_manager().mesh_entity_manager;
            if( manager.is_region( entity_type ) ) {
                std::unique_ptr< VolumeMeshBuilder< 3 > > builder =
                    geometry.create_region_builder( id );
                builder->load_mesh( file_name );
            }
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGM< DIMENSION >::load_meshes( unzFile& uz )
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
            load_mesh_entity( entity_type, file_name, id );
            GEO::FileSystem::delete_file( file_name );
        }
        Logger::instance()->set_minimal( false );
    }

    template class RINGMESH_API GeoModelBuilderGM< 2 > ;

    template class RINGMESH_API GeoModelBuilderGM< 3 > ;

}
