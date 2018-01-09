/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

namespace
{
    using namespace RINGMesh;

    template < index_t >
    class GeoModelBuilderGM;

    bool match_mesh_entity_type( const MeshEntityType& type )
    {
        if( type == Corner3D::type_name_static() )
            return true;
        if( type == Line3D::type_name_static() )
            return true;
        if( type == Surface3D::type_name_static() )
            return true;
        if( type == Region3D::type_name_static() )
            return true;
        return false;
    }

    template < index_t DIMENSION >
    class GeoModelBuilderGMImpl
    {
    public:
        GeoModelBuilderGMImpl( GeoModelBuilderGM< DIMENSION >& builder,
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

    template < index_t DIMENSION >
    class GeoModelBuilderGMImpl_0 : public GeoModelBuilderGMImpl< DIMENSION >
    {
    public:
        GeoModelBuilderGMImpl_0( GeoModelBuilderGM< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel )
            : GeoModelBuilderGMImpl< DIMENSION >( builder, geomodel )
        {
        }
        virtual ~GeoModelBuilderGMImpl_0() = default;

        void read_mesh_entity_line( GEO::LineInput& file_line ) override
        {
            // First line : type - id - name - geol_feature
            if( file_line.nb_fields() < 4 )
            {
                throw RINGMeshException( "I/O", "Invalid line: ",
                    file_line.line_number(), ", 4 fields are expected, the "
                                             "type, id, name, and geological "
                                             "feature" );
            }
            auto entity = read_first_line( file_line );
            read_second_line( file_line, entity );
        }

    protected:
        virtual gmme_id read_first_line( GEO::LineInput& file_line )
        {
            gmme_id cur_gmme{ MeshEntityType{ file_line.field( 0 ) },
                file_line.field_as_uint( 1 ) };
            this->builder_.info.set_mesh_entity_name(
                cur_gmme, file_line.field( 2 ) );
            return cur_gmme;
        }
        void read_second_line(
            GEO::LineInput& file_line, const gmme_id& entity );
    };

    template < index_t DIMENSION >
    class GeoModelBuilderGMImpl_1 : public GeoModelBuilderGMImpl_0< DIMENSION >
    {
    public:
        GeoModelBuilderGMImpl_1( GeoModelBuilderGM< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel )
            : GeoModelBuilderGMImpl_0< DIMENSION >( builder, geomodel )
        {
        }
        virtual ~GeoModelBuilderGMImpl_1() = default;

        void read_mesh_entity_line( GEO::LineInput& file_line ) override
        {
            // Read this entity
            // First line : type - id - name - geol_feature - mesh type
            if( file_line.nb_fields() < 5 )
            {
                throw RINGMeshException( "I/O", "Invalid line: ",
                    file_line.line_number(),
                    ", 5 fields are expected, the type, id, name, ",
                    "geological feature, and mesh type" );
            }
            auto entity = this->read_first_line( file_line );

            const auto mesh_type = file_line.field( 4 );
            if( GEO::String::string_starts_with( mesh_type, "Geogram" ) )
            {
                this->builder_.geometry.change_mesh_data_structure(
                    entity, old_2_new_name( mesh_type ) );
            }
            else
            {
                this->builder_.geometry.change_mesh_data_structure(
                    entity, mesh_type );
            }

            this->read_second_line( file_line, entity );
        }

    private:
        const std::string& old_2_new_name( const std::string& old_name )
        {
            auto new_name_pos = GEO::String::to_uint( GEO::String::to_string(
                old_name.at( old_name.length() - 2 ) ) );
            return new_names[new_name_pos];
        }

    private:
        static const std::string new_names[4];
    };

    template < index_t DIMENSION >
    const std::string GeoModelBuilderGMImpl_1< DIMENSION >::new_names[4] = {
        std::string( "GeogramPointSetMesh" ), std::string( "GeogramLineMesh" ),
        std::string( "GeogramSurfaceMesh" ),
        std::string( "GeogramVolumeMesh" )
    };

    template < index_t DIMENSION >
    class GeoModelBuilderGMImpl_2 : public GeoModelBuilderGMImpl_1< DIMENSION >
    {
    public:
        GeoModelBuilderGMImpl_2( GeoModelBuilderGM< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel )
            : GeoModelBuilderGMImpl_1< DIMENSION >( builder, geomodel )
        {
        }
        virtual ~GeoModelBuilderGMImpl_2() = default;

        void read_mesh_entity_line( GEO::LineInput& file_line ) override
        {
            // Read this entity
            // First line : type - id - name  - mesh type
            if( file_line.nb_fields() < 4 )
            {
                throw RINGMeshException( "I/O", "Invalid line: ",
                    file_line.line_number(),
                    ", 4 fields are expected, the type, id, name, ",
                    "geological feature, and mesh type" );
            }
            auto entity = this->read_first_line( file_line );

            const auto mesh_type = file_line.field( 3 );
            this->builder_.geometry.change_mesh_data_structure(
                entity, mesh_type );

            this->read_second_line( file_line, entity );
        }
    };

    template < index_t DIMENSION >
    class GeoModelBuilderGM final : public GeoModelBuilderFile< DIMENSION >
    {
    public:
        static const index_t NB_VERSION = 3;
        GeoModelBuilderGM(
            GeoModel< DIMENSION >& geomodel, std::string filename )
            : GeoModelBuilderFile< DIMENSION >(
                  geomodel, std::move( filename ) )
        {
            version_impl_[0].reset(
                new GeoModelBuilderGMImpl_0< DIMENSION >( *this, geomodel ) );
            version_impl_[1].reset(
                new GeoModelBuilderGMImpl_1< DIMENSION >( *this, geomodel ) );
            version_impl_[2].reset(
                new GeoModelBuilderGMImpl_2< DIMENSION >( *this, geomodel ) );
        }
        virtual ~GeoModelBuilderGM() = default;

    private:
        void load_geological_entities(
            const std::string& geological_entity_file )
        {
            GEO::LineInput file_line{ geological_entity_file };
            while( !file_line.eof() && file_line.get_line() )
            {
                file_line.get_fields();
                if( file_line.nb_fields() > 0 )
                {
                    // If there is no geological entity
                    if( file_line.field_matches( 0, "No" )
                        && file_line.field_matches( 1, "geological" )
                        && file_line.field_matches( 2, "entity" ) )
                    {
                        return;
                    }
                    // Number of entities of a given type
                    if( file_line.field_matches( 0, "Nb" ) )
                    {
                        // Allocate the space
                        this->geology.create_geological_entities(
                            GeologicalEntityType( file_line.field( 1 ) ),
                            file_line.field_as_uint( 2 ) );
                    }
                    else
                    {
                        GeologicalEntityType type{ file_line.field( 0 ) };
                        auto id = file_line.field_as_uint( 1 );
                        gmge_id entity{ type, id };
                        this->info.set_geological_entity_name(
                            entity, file_line.field( 2 ) );
                        this->geology.set_geological_entity_geol_feature(
                            entity, GeoModelGeologicalEntity< DIMENSION >::
                                        determine_geological_type(
                                            file_line.field( 3 ) ) );
                        file_line.get_line();
                        file_line.get_fields();
                        for( auto in_b : range( file_line.nb_fields() ) )
                        {
                            this->geology.add_parent_children_relation( entity,
                                { this->geomodel_.entity_type_manager()
                                        .relationship_manager.child_type(
                                            type ),
                                    file_line.field_as_uint( in_b ) } );
                        }
                    }
                }
            }
        }

        /*!
         * @brief Load meshes of all the mesh entities from a zip file
         * @param[in] uz the zip file
         */
        void load_meshes( UnZipFile& uz )
        {
            uz.start_extract();

            Logger::instance()->set_minimal( true );
            std::vector< std::future< void > > files;
            do
            {
                if( GEO::FileSystem::extension( uz.get_current_filename() )
                    == "txt" )
                {
                    continue;
                }

                const auto file_name = uz.get_current_file();
                files.push_back(
                    std::async( std::launch::deferred, [file_name, this] {
                        auto file_without_extension =
                            GEO::FileSystem::base_name( file_name );
                        std::string entity_type, entity_id;
                        GEO::String::split_string( file_without_extension, '_',
                            entity_type, entity_id );
                        index_t id{ NO_ID };
                        GEO::String::from_string( entity_id, id );
                        load_mesh_entity(
                            MeshEntityType{ entity_type }, file_name, id );
                        GEO::FileSystem::delete_file( file_name );
                    } ) );
            } while( uz.next_file() );

            for( auto& file : files )
            {
                file.wait();
            }
            Logger::instance()->set_minimal( false );
        }

        void load_mesh_entity( const MeshEntityType& entity_type,
            const std::string& file_name,
            index_t id );

        bool load_mesh_entity_base( const MeshEntityType& entity_type,
            const std::string& file_name,
            index_t id )
        {
            const auto& manager =
                this->geomodel_.entity_type_manager().mesh_entity_manager;
            if( manager.is_corner( entity_type ) )
            {
                auto builder = this->geometry.create_corner_builder( id );
                builder->load_mesh( file_name );
                return true;
            }
            else if( manager.is_line( entity_type ) )
            {
                auto builder = this->geometry.create_line_builder( id );
                builder->load_mesh( file_name );
                return true;
            }
            else if( manager.is_surface( entity_type ) )
            {
                auto builder = this->geometry.create_surface_builder( id );
                builder->load_mesh( file_name );
                return true;
            }
            return false;
        }

        void load_file() final
        {
            const std::string directory_to_unzip{
                GEO::FileSystem::normalized_path(
                    GEO::FileSystem::dir_name( this->filename() ) )
                + "/" + std::to_string(
                            std::hash< std::string >()( this->filename() ) )
            };
            UnZipFile uz{ this->filename(), directory_to_unzip };

            const auto mesh_entity_file = uz.get_file( "mesh_entities.txt" );
            load_mesh_entities( mesh_entity_file );
            auto ok = GEO::FileSystem::delete_file( mesh_entity_file );
            ringmesh_unused( ok ); // avoids warning in release
            ringmesh_assert( ok );
            load_meshes( uz );

            const auto geological_entity_file =
                uz.get_file( "geological_entities.txt" );
            load_geological_entities( geological_entity_file );
            ok = GEO::FileSystem::delete_file( geological_entity_file );
            ringmesh_assert( ok );

            GEO::FileSystem::delete_directory( directory_to_unzip );
        }

        void load_mesh_entities( const std::string& mesh_entity_file )
        {
            GEO::LineInput file_line{ mesh_entity_file };
            while( !file_line.eof() && file_line.get_line() )
            {
                file_line.get_fields();
                if( file_line.nb_fields() > 0 )
                {
                    if( file_line.field_matches( 0, "Version" ) )
                    {
                        file_version_ = file_line.field_as_uint( 1 );
                    }
                    // Name of the geomodel
                    else if( file_line.field_matches( 0, "GeoModel" ) )
                    {
                        if( file_line.nb_fields() > 2 )
                        {
                            this->info.set_geomodel_name(
                                file_line.field( 2 ) );
                        }
                    }
                    // Number of entities of a given type
                    else if( file_line.field_matches( 0, "Nb" ) )
                    {
                        // Allocate the space
                        this->topology.create_mesh_entities(
                            MeshEntityType( file_line.field( 1 ) ),
                            file_line.field_as_uint( 2 ) );
                    }
                    // Mesh entities
                    else if( match_mesh_entity_type(
                                 MeshEntityType( file_line.field( 0 ) ) ) )
                    {
                        version_impl_[file_version_]->read_mesh_entity_line(
                            file_line );
                    }
                }
            }
        }

        void load_region_if_entity_is_region( const std::string& entity_type,
            index_t id,
            const std::string& file_name,
            const MeshEntityTypeManager< DIMENSION >& manager );

    private:
        index_t file_version_{ 0 };
        std::unique_ptr< GeoModelBuilderGMImpl< DIMENSION > >
            version_impl_[NB_VERSION];
    };

    template <>
    void GeoModelBuilderGM< 2 >::load_mesh_entity(
        const MeshEntityType& entity_type,
        const std::string& file_name,
        index_t id )
    {
        load_mesh_entity_base( entity_type, file_name, id );
    }

    template <>
    void GeoModelBuilderGM< 3 >::load_mesh_entity(
        const MeshEntityType& entity_type,
        const std::string& file_name,
        index_t id )
    {
        if( !load_mesh_entity_base( entity_type, file_name, id ) )
        {
            const auto& manager =
                this->geomodel_.entity_type_manager().mesh_entity_manager;
            if( manager.is_region( entity_type ) )
            {
                auto builder = geometry.create_region_builder( id );
                builder->load_mesh( file_name );
            }
        }
    }

    /*!
     * @brief Write in the out stream things to save for CONTACT, INTERFACE and
     * LAYERS
     */
    template < index_t DIMENSION >
    void save_geological_entity(
        std::ofstream& out, const GeoModelGeologicalEntity< DIMENSION >& E )
    {
        /// First line:  TYPE - ID - NAME - GEOL
        out << E.gmge() << " " << E.name() << " ";
        out << GeoModelGeologicalEntity< DIMENSION >::geol_name(
                   E.geological_feature() )
            << EOL;

        /// Second line:  IDS of children
        for( auto j : range( E.nb_children() ) )
        {
            out << E.child_gmme( j ).index() << " ";
        }
        out << EOL;
    }

    /*!
     * @brief Save the connectivity of a GeoModel in a file
     * @param[in] M the GeoModel
     * @param[in] file_name path to the input file
     */
    template < index_t DIMENSION >
    void save_geological_entities(
        const GeoModel< DIMENSION >& geomodel, const std::string& file_name )
    {
        std::ofstream out{ file_name.c_str() };
        out.precision( 16 );
        if( out.bad() )
        {
            throw RINGMeshException(
                "I/O", "Error when opening the file: ", file_name );
        }

        if( geomodel.nb_geological_entity_types() == 0 )
        {
            // Compression of an empty files crashes ? (in debug on windows at
            // least)
            out << "No geological entity in the geomodel" << EOL;
            return;
        }
        for( auto i : range( geomodel.nb_geological_entity_types() ) )
        {
            const auto& type = geomodel.geological_entity_type( i );
            auto nb = geomodel.nb_geological_entities( type );
            out << "Nb " << type << " " << nb << EOL;
        }
        for( auto i : range( geomodel.nb_geological_entity_types() ) )
        {
            const auto& type = geomodel.geological_entity_type( i );
            auto nb = geomodel.nb_geological_entities( type );
            for( auto j : range( nb ) )
            {
                save_geological_entity(
                    out, geomodel.geological_entity( type, j ) );
            }
        }
        out << std::flush;
    }

    template < typename ENTITY, index_t DIMENSION >
    void save_mesh_entities_of_type(
        const GeoModel< DIMENSION >& geomodel, std::ofstream& out )
    {
        const MeshEntityType& type = ENTITY::type_name_static();
        for( auto e : range( geomodel.nb_mesh_entities( type ) ) )
        {
            const ENTITY& cur_mesh_entity = dynamic_cast< const ENTITY& >(
                geomodel.mesh_entity( type, e ) );
            out << type << " " << e << " " << cur_mesh_entity.name() << " "
                << cur_mesh_entity.mesh().type_name() << EOL;
            out << "boundary ";
            for( auto b : range( cur_mesh_entity.nb_boundaries() ) )
            {
                out << cur_mesh_entity.boundary_gmme( b ).index() << " ";
            }
            out << EOL;
        }
    }

    void save_dimension( index_t dimension, std::ofstream& out )
    {
        out << "Dimension " << dimension << EOL;
    }

    template < index_t DIMENSION >
    void save_version_and_name(
        const GeoModel< DIMENSION >& geomodel, std::ofstream& out )
    {
        out << "Version 2" << EOL;
        out << "GeoModel name " << geomodel.name() << EOL;
    }

    template < index_t DIMENSION >
    void save_number_of_mesh_entities_base(
        const GeoModel< DIMENSION >& geomodel, std::ofstream& out )
    {
        // Numbers of the different types of mesh entities
        out << "Nb " << Corner< DIMENSION >::type_name_static() << " "
            << geomodel.nb_corners() << EOL;
        out << "Nb " << Line< DIMENSION >::type_name_static() << " "
            << geomodel.nb_lines() << EOL;
        out << "Nb " << Surface< DIMENSION >::type_name_static() << " "
            << geomodel.nb_surfaces() << EOL;
    }

    template < index_t DIMENSION >
    void save_number_of_mesh_entities(
        const GeoModel< DIMENSION >& M, std::ofstream& out );

    template <>
    void save_number_of_mesh_entities(
        const GeoModel2D& geomodel, std::ofstream& out )
    {
        save_number_of_mesh_entities_base( geomodel, out );
    }

    template <>
    void save_number_of_mesh_entities(
        const GeoModel3D& geomodel, std::ofstream& out )
    {
        save_number_of_mesh_entities_base( geomodel, out );
        out << "Nb " << Region< 3 >::type_name_static() << " "
            << geomodel.nb_regions() << EOL;
    }

    template < index_t DIMENSION >
    void save_mesh_entities_topology_and_sides(
        const GeoModel< DIMENSION >& geomodel, std::ofstream& out );

    template < template < index_t > class ENTITY, index_t DIMENSION >
    void save_mesh_entities_topology_and_sides_impl(
        const GeoModel< DIMENSION >& geomodel, std::ofstream& out )
    {
        for( auto i : range( geomodel.nb_mesh_entities(
                 ENTITY< DIMENSION >::type_name_static() ) ) )
        {
            const ENTITY< DIMENSION >& E =
                static_cast< const ENTITY< DIMENSION >& >( geomodel.mesh_entity(
                    ENTITY< DIMENSION >::type_name_static(), i ) );
            // Save ID - NAME
            out << E.gmme() << " " << E.name() << " " << E.mesh().type_name()
                << EOL;
            // Second line Signed ids of boundary surfaces
            for( auto j : range( E.nb_boundaries() ) )
            {
                if( E.side( j ) )
                {
                    out << "+";
                }
                else
                {
                    out << "-";
                }
                out << E.boundary_gmme( j ).index() << " ";
            }
            out << EOL;
        }
    }

    template <>
    void save_mesh_entities_topology_and_sides(
        const GeoModel2D& geomodel, std::ofstream& out )
    {
        save_mesh_entities_topology_and_sides_impl< Surface >( geomodel, out );
    }

    template <>
    void save_mesh_entities_topology_and_sides(
        const GeoModel3D& geomodel, std::ofstream& out )
    {
        save_mesh_entities_topology_and_sides_impl< Region >( geomodel, out );
    }

    template < index_t DIMENSION >
    void save_mesh_entities_topology_base(
        const GeoModel< DIMENSION >& geomodel, std::ofstream& out )
    {
        save_mesh_entities_of_type< Corner< DIMENSION > >( geomodel, out );
        save_mesh_entities_of_type< Line< DIMENSION > >( geomodel, out );
    }

    template < index_t DIMENSION >
    void save_mesh_entities_topology(
        const GeoModel< DIMENSION >& geomodel, std::ofstream& out );

    template <>
    void save_mesh_entities_topology(
        const GeoModel2D& geomodel, std::ofstream& out )
    {
        save_mesh_entities_topology_base( geomodel, out );
    }
    template <>
    void save_mesh_entities_topology(
        const GeoModel3D& geomodel, std::ofstream& out )
    {
        save_mesh_entities_topology_base( geomodel, out );
        save_mesh_entities_of_type< Surface3D >( geomodel, out );
    }

    /*!
     * @brief Save the topology of a GeoModelin a file
     * @param[in] geomodel the GeoModel
     * @param[in] file_name the output file name
     */
    template < index_t DIMENSION >
    void save_mesh_entities(
        const GeoModel< DIMENSION >& geomodel, const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() );
        out.precision( 16 );
        if( out.bad() )
        {
            throw RINGMeshException(
                "I/O", "Error when opening the file: ", file_name );
        }
        save_dimension( DIMENSION, out );
        save_version_and_name( geomodel, out );
        save_number_of_mesh_entities( geomodel, out );
        save_mesh_entities_topology( geomodel, out );
        save_mesh_entities_topology_and_sides( geomodel, out );
        out << std::flush;
    }

    template < index_t DIMENSION >
    bool save_mesh( const GeoModelMeshEntity< DIMENSION >& geomodel_entity_mesh,
        const std::string& name );

    template <>
    bool save_mesh( const GeoModelMeshEntity3D& geomodel_entity_mesh,
        const std::string& name )
    {
        if( geomodel_entity_mesh.type_name() == Region3D::type_name_static() )
        {
            const auto& region = geomodel_entity_mesh.geomodel().region(
                geomodel_entity_mesh.index() );
            if( !region.is_meshed() )
            {
                // a region is not necessary meshed.
                return false;
            }
        }
        geomodel_entity_mesh.save( name );
        return true;
    }

    template <>
    bool save_mesh( const GeoModelMeshEntity2D& geomodel_entity_mesh,
        const std::string& name )
    {
        if( geomodel_entity_mesh.type_name() == Surface2D::type_name_static() )
        {
            const auto& surface = geomodel_entity_mesh.geomodel().surface(
                geomodel_entity_mesh.index() );
            if( !surface.is_meshed() )
            {
                return false;
            }
        }
        geomodel_entity_mesh.save( name );
        return true;
    }

    template < typename ENTITY >
    std::string build_string_for_geomodel_entity_export( const ENTITY& entity )
    {
        const auto& id = entity.gmme();
        std::string base_name =
            id.type().string() + "_" + std::to_string( id.index() );
        return base_name + "." + entity.mesh().default_extension();
    }

    /*!
     * @brief Save the GeoModelMeshEntity in a meshb file
     * @param[in] geomodel_entity_mesh the GeoModelMeshEntity you want to save
     * @param[in] zf the zip file
     */
    template < typename ENTITY >
    void save_geomodel_mesh_entity( const ENTITY& geomodel_entity_mesh,
        std::vector< std::string >& filenames )
    {
        static std::mutex lock;
        auto name =
            build_string_for_geomodel_entity_export( geomodel_entity_mesh );
        if( save_mesh( geomodel_entity_mesh, name ) )
        {
            std::lock_guard< std::mutex > locking( lock );
            filenames.push_back( name );
        }
    }

    void zip_files( const std::vector< std::string >& filenames, ZipFile& zf )
    {
        for( const std::string& name : filenames )
        {
            zf.add_file( name );
            GEO::FileSystem::delete_file( name );
        }
    }

    template < template < index_t > class ENTITY, index_t DIMENSION >
    void save_geomodel_mesh_entities( const GeoModel< DIMENSION >& geomodel,
        std::vector< std::string >& filenames )
    {
        const auto& type = ENTITY< DIMENSION >::type_name_static();
        auto* logger = Logger::instance();
        auto logger_status = logger->is_quiet();
        logger->set_quiet( true );
        parallel_for( geomodel.nb_mesh_entities( type ),
            [&geomodel, &type, &filenames]( index_t i ) {
                const auto& entity = dynamic_cast< const ENTITY< DIMENSION >& >(
                    geomodel.mesh_entity( type, i ) );
                save_geomodel_mesh_entity< ENTITY< DIMENSION > >(
                    entity, filenames );
            } );
        logger->set_quiet( logger_status );
    }

    template < index_t DIMENSION >
    void save_all_geomodel_mesh_entities_base(
        const GeoModel< DIMENSION >& geomodel,
        std::vector< std::string >& filenames )
    {
        save_geomodel_mesh_entities< Corner, DIMENSION >( geomodel, filenames );
        save_geomodel_mesh_entities< Line, DIMENSION >( geomodel, filenames );
        save_geomodel_mesh_entities< Surface, DIMENSION >(
            geomodel, filenames );
    }

    template < index_t DIMENSION >
    void save_all_geomodel_mesh_entities( const GeoModel< DIMENSION >& geomodel,
        std::vector< std::string >& filenames );

    template <>
    void save_all_geomodel_mesh_entities(
        const GeoModel2D& geomodel, std::vector< std::string >& filenames )
    {
        save_all_geomodel_mesh_entities_base( geomodel, filenames );
    }
    template <>
    void save_all_geomodel_mesh_entities(
        const GeoModel3D& geomodel, std::vector< std::string >& filenames )
    {
        save_all_geomodel_mesh_entities_base( geomodel, filenames );
        save_geomodel_mesh_entities< Region, 3 >( geomodel, filenames );
    }

    template < index_t DIMENSION >
    index_t nb_mesh_entities( const GeoModel< DIMENSION >& geomodel );

    template <>
    index_t nb_mesh_entities( const GeoModel2D& geomodel )
    {
        return geomodel.nb_corners() + geomodel.nb_lines()
               + geomodel.nb_surfaces();
    }

    template <>
    index_t nb_mesh_entities( const GeoModel3D& geomodel )
    {
        return geomodel.nb_corners() + geomodel.nb_lines()
               + geomodel.nb_surfaces() + geomodel.nb_regions();
    }

    index_t find_dimension( const std::string& mesh_entities_filename )
    {
        GEO::LineInput file_line{ mesh_entities_filename };
        while( !file_line.eof() && file_line.get_line() )
        {
            file_line.get_fields();
            if( file_line.nb_fields() == 2 )
            {
                if( file_line.field_matches( 0, "Dimension" ) )
                {
                    return file_line.field_as_uint( 1 );
                }
            }
        }
        return 3;
    }

    template < index_t DIMENSION >
    class GeoModelHandlerGM final : public GeoModelOutputHandler< DIMENSION >,
                                    public GeoModelInputHandler< DIMENSION >
    {
    public:
        void load(
            const std::string& filename, GeoModel< DIMENSION >& geomodel ) final
        {
            GeoModelBuilderGM< DIMENSION > builder{ geomodel, filename };
            builder.build_geomodel();
        }

        void save( const GeoModel< DIMENSION >& geomodel,
            const std::string& filename ) final
        {
            ZipFile zf{ filename };

            const std::string mesh_entity_file{ "mesh_entities.txt" };
            save_mesh_entities( geomodel, mesh_entity_file );
            zf.add_file( mesh_entity_file );
            GEO::FileSystem::delete_file( mesh_entity_file );

            const std::string geological_entity_file{
                "geological_entities.txt"
            };
            save_geological_entities( geomodel, geological_entity_file );
            zf.add_file( geological_entity_file );
            GEO::FileSystem::delete_file( geological_entity_file );

            auto nb_mesh_entites = nb_mesh_entities( geomodel );
            std::vector< std::string > filenames;
            filenames.reserve( nb_mesh_entites );
            Logger::instance()->set_quiet( true );
            save_all_geomodel_mesh_entities( geomodel, filenames );
            Logger::instance()->set_quiet( false );
            std::sort( filenames.begin(), filenames.end() );
            zip_files( filenames, zf );
        }

        index_t dimension( const std::string& filename ) const final
        {
            UnZipFile uz{ filename };
            const std::string mesh_entity_file( "mesh_entities.txt" );
            uz.get_file( mesh_entity_file );
            auto dimension = find_dimension( mesh_entity_file );
            auto ok = GEO::FileSystem::delete_file( mesh_entity_file );
            ringmesh_unused( ok );
            return dimension;
        }

        virtual ~GeoModelHandlerGM() = default;

    private:
        void save_geomodel_regions( const GeoModel< DIMENSION >& geomodel,
            std::vector< std::string >& filenames );
    };

    ALIAS_2D_AND_3D( GeoModelHandlerGM );

    template <>
    void GeoModelBuilderGMImpl_0< 2 >::read_second_line(
        GEO::LineInput& file_line, const gmme_id& entity )
    {
        file_line.get_line();
        file_line.get_fields();
        const auto& manager =
            this->geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_surface( entity.type() ) )
        {
            for( auto c : range( file_line.nb_fields() ) )
            {
                bool side{ false };
                if( std::strncmp( file_line.field( c ), "+", 1 ) == 0 )
                {
                    side = true;
                }
                index_t s{ NO_ID };
                GEO::String::from_string( &file_line.field( c )[1], s );

                this->builder_.topology.add_surface_line_boundary_relation(
                    entity.index(), s, side );
            }
        }
        else
        {
            // Second line : indices of boundaries
            for( auto c : range( 1, file_line.nb_fields() ) )
            {
                ringmesh_assert( entity.type() == Line2D::type_name_static() );
                index_t boundary_id{ file_line.field_as_uint( c ) };
                this->builder_.topology.add_line_corner_boundary_relation(
                    entity.index(), boundary_id );
            }
        }
    }

    template <>
    void GeoModelBuilderGMImpl_0< 3 >::read_second_line(
        GEO::LineInput& file_line, const gmme_id& entity )
    {
        file_line.get_line();
        file_line.get_fields();
        const auto& manager =
            this->geomodel_.entity_type_manager().mesh_entity_manager;
        if( manager.is_region( entity.type() ) )
        {
            for( auto c : range( file_line.nb_fields() ) )
            {
                bool side{ false };
                if( std::strncmp( file_line.field( c ), "+", 1 ) == 0 )
                {
                    side = true;
                }
                index_t s{ NO_ID };
                GEO::String::from_string( &file_line.field( c )[1], s );

                this->builder_.topology.add_region_surface_boundary_relation(
                    entity.index(), s, side );
            }
        }
        else
        {
            // Second line : indices of boundaries
            // Corners are skipped
            for( auto c : range( 1, file_line.nb_fields() ) )
            {
                ringmesh_assert( manager.is_line( entity.type() )
                                 || manager.is_surface( entity.type() ) );
                index_t boundary_id{ file_line.field_as_uint( c ) };
                if( manager.is_line( entity.type() ) )
                {
                    this->builder_.topology.add_line_corner_boundary_relation(
                        entity.index(), boundary_id );
                }
                else
                {
                    ringmesh_assert( manager.is_surface( entity.type() ) );
                    this->builder_.topology.add_surface_line_boundary_relation(
                        entity.index(), boundary_id );
                }
            }
        }
    }
}
