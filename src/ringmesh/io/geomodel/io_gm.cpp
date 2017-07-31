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
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
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

namespace {
    /*!
     * @brief Write in the out stream things to save for CONTACT, INTERFACE and LAYERS
     */
    template< index_t DIMENSION >
    void save_geological_entity(
        std::ofstream& out,
        const GeoModelGeologicalEntity< DIMENSION >& E )
    {
        /// First line:  TYPE - ID - NAME - GEOL
        out << E.gmge() << " " << E.name() << " ";
        out << GeoModelGeologicalEntity < DIMENSION
            > ::geol_name( E.geological_feature() ) << EOL;

        /// Second line:  IDS of children
        for( index_t j : range( E.nb_children() ) ) {
            out << E.child_gmme( j ).index() << " ";
        }
        out << EOL;
    }

    /*!
     * @brief Save the connectivity of a GeoModel in a file
     * @param[in] M the GeoModel
     * @param[in] file_name path to the input file
     */
    template< index_t DIMENSION >
    void save_geological_entities(
        const GeoModel< DIMENSION >& geomodel,
        const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() );
        out.precision( 16 );
        if( out.bad() ) {
            throw RINGMeshException( "I/O", "Error when opening the file: ",
                file_name );
        }

        if( geomodel.nb_geological_entity_types() == 0 ) {
            // Compression of an empty files crashes ? (in debug on windows at least)
            out << "No geological entity in the geomodel" << EOL;
            return;
        }
        for( index_t i : range( geomodel.nb_geological_entity_types() ) ) {
            const GeologicalEntityType& type = geomodel.geological_entity_type( i );
            index_t nb = geomodel.nb_geological_entities( type );
            out << "Nb " << type << " " << nb << EOL;
        }
        for( index_t i : range( geomodel.nb_geological_entity_types() ) ) {
            const GeologicalEntityType& type = geomodel.geological_entity_type( i );
            index_t nb = geomodel.nb_geological_entities( type );
            for( index_t j : range( nb ) ) {
                save_geological_entity( out, geomodel.geological_entity( type, j ) );
            }
        }
        out << std::flush;
    }

    template< typename ENTITY, index_t DIMENSION >
    void save_mesh_entities_of_type(
        const GeoModel< DIMENSION >& geomodel,
        std::ofstream& out )
    {
        const MeshEntityType& type = ENTITY::type_name_static();
        for( index_t e : range( geomodel.nb_mesh_entities( type ) ) ) {
            const ENTITY& cur_mesh_entity =
                dynamic_cast< const ENTITY& >( geomodel.mesh_entity( type, e ) );
            out << type << " " << e << " " << cur_mesh_entity.name() << " "
                << cur_mesh_entity.low_level_mesh_storage().type_name() << EOL;
            out << "boundary ";
            for( index_t b : range( cur_mesh_entity.nb_boundaries() ) ) {
                out << cur_mesh_entity.boundary_gmme( b ).index() << " ";
            }
            out << EOL;
        }
    }

    template< index_t DIMENSION >
    void save_dimension( const GeoModel< DIMENSION >& geomodel, std::ofstream& out )
    {
        out << "Dimension " << DIMENSION << EOL;
    }

    template< index_t DIMENSION >
    void save_version_and_name(
        const GeoModel< DIMENSION >& geomodel,
        std::ofstream& out )
    {
        out << "Version 2" << EOL;
        out << "GeoModel name " << geomodel.name() << EOL;
    }

    template< index_t DIMENSION >
    void save_number_of_mesh_entities_base(
        const GeoModel< DIMENSION >& geomodel,
        std::ofstream& out )
    {
        // Numbers of the different types of mesh entities
        out << "Nb " << Corner < DIMENSION
            > ::type_name_static() << " " << geomodel.nb_corners() << EOL;
        out << "Nb " << Line < DIMENSION
            > ::type_name_static() << " " << geomodel.nb_lines() << EOL;
        out << "Nb " << Surface < DIMENSION
            > ::type_name_static() << " " << geomodel.nb_surfaces() << EOL;
    }

    template< index_t DIMENSION >
    void save_number_of_mesh_entities(
        const GeoModel< DIMENSION >& M,
        std::ofstream& out );

    template< >
    void save_number_of_mesh_entities(
        const GeoModel2D& geomodel,
        std::ofstream& out )
    {
        save_number_of_mesh_entities_base( geomodel, out );
    }

    template< >
    void save_number_of_mesh_entities(
        const GeoModel3D& geomodel,
        std::ofstream& out )
    {
        save_number_of_mesh_entities_base( geomodel, out );
        out << "Nb " << Region < 3
            > ::type_name_static() << " " << geomodel.nb_regions() << EOL;
    }

    template< index_t DIMENSION >
    void save_mesh_entities_topology_and_sides(
        const GeoModel< DIMENSION >& geomodel,
        std::ofstream& out );

    template< template< index_t > class ENTITY, index_t DIMENSION >
    void save_mesh_entities_topology_and_sides_impl(
        const GeoModel< DIMENSION >& geomodel,
        std::ofstream& out )
    {
        for( index_t i : range(
            geomodel.nb_mesh_entities( ENTITY< DIMENSION >::type_name_static() ) ) ) {
            const ENTITY< DIMENSION >& E =
                static_cast< const ENTITY< DIMENSION >& >( geomodel.mesh_entity(
                    ENTITY< DIMENSION >::type_name_static(), i ) );
            // Save ID - NAME
            out << E.gmme() << " " << E.name() << " "
                << E.low_level_mesh_storage().type_name() << EOL;
            // Second line Signed ids of boundary surfaces
            for( index_t j : range( E.nb_boundaries() ) ) {
                if( E.side( j ) ) {
                    out << "+";
                } else {
                    out << "-";
                }
                out << E.boundary_gmme( j ).index() << " ";
            }
            out << EOL;
        }
    }

    template< >
    void save_mesh_entities_topology_and_sides(
        const GeoModel2D& geomodel,
        std::ofstream& out )
    {
        save_mesh_entities_topology_and_sides_impl< Surface >( geomodel, out );
    }

    template< >
    void save_mesh_entities_topology_and_sides(
        const GeoModel3D& geomodel,
        std::ofstream& out )
    {
        save_mesh_entities_topology_and_sides_impl< Region >( geomodel, out );

    }

    template< index_t DIMENSION >
    void save_mesh_entities_topology_base(
        const GeoModel< DIMENSION >& geomodel,
        std::ofstream& out )
    {
        save_mesh_entities_of_type< Corner< DIMENSION > >( geomodel, out );
        save_mesh_entities_of_type< Line< DIMENSION > >( geomodel, out );
    }

    template< index_t DIMENSION >
    void save_mesh_entities_topology(
        const GeoModel< DIMENSION >& geomodel,
        std::ofstream& out );

    template< >
    void save_mesh_entities_topology(
        const GeoModel2D& geomodel,
        std::ofstream& out )
    {
        save_mesh_entities_topology_base( geomodel, out );

    }
    template< >
    void save_mesh_entities_topology(
        const GeoModel3D& geomodel,
        std::ofstream& out )
    {
        save_mesh_entities_topology_base( geomodel, out );
        save_mesh_entities_of_type< Surface3D >( geomodel, out );
    }

    template< index_t DIMENSION >
    void save_universe( const GeoModel< DIMENSION >& M, std::ofstream& out )
    {
        out << "Universe " << EOL;
        for( index_t j : range( M.universe().nb_boundaries() ) ) {
            if( M.universe().side( j ) ) {
                out << "+";
            } else {
                out << "-";
            }
            out << M.universe().boundary_gmme( j ).index() << " ";
        }
        out << EOL;
    }

    /*!
     * @brief Save the topology of a GeoModelin a file
     * @param[in] geomodel the GeoModel
     * @param[in] file_name the output file name
     */
    template< index_t DIMENSION >
    void save_mesh_entities(
        const GeoModel< DIMENSION >& geomodel,
        const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() );
        out.precision( 16 );
        if( out.bad() ) {
            throw RINGMeshException( "I/O", "Error when opening the file: ",
                file_name );
        }
        save_dimension( geomodel, out );
        save_version_and_name( geomodel, out );
        save_number_of_mesh_entities( geomodel, out );
        save_mesh_entities_topology( geomodel, out );
        save_mesh_entities_topology_and_sides( geomodel, out );
        save_universe( geomodel, out );
        out << std::flush;
    }

    template< index_t DIMENSION >
    bool save_mesh(
        const GeoModelMeshEntity< DIMENSION >& geomodel_entity_mesh,
        const std::string& name );

    template< >
    bool save_mesh(
        const GeoModelMeshEntity3D& geomodel_entity_mesh,
        const std::string& name )
    {
        if( geomodel_entity_mesh.type_name() == Region3D::type_name_static() ) {
            const Region3D& region = geomodel_entity_mesh.geomodel().region(
                geomodel_entity_mesh.index() );
            if( !region.is_meshed() ) {
                // a region is not necessary meshed.
                return false;
            }
        }
        geomodel_entity_mesh.save( name );
        return true;
    }

    template< >
    bool save_mesh(
        const GeoModelMeshEntity2D& geomodel_entity_mesh,
        const std::string& name )
    {
        geomodel_entity_mesh.save( name );
        return true;
    }

    template< typename ENTITY >
    std::string build_string_for_geomodel_entity_export( const ENTITY& entity )
    {
        const gmme_id& id = entity.gmme();
        std::string base_name = id.type().to_string() + "_"
            + std::to_string( id.index() );
        return base_name + "." + entity.low_level_mesh_storage().default_extension();
    }

    /*!
     * @brief Save the GeoModelMeshEntity in a meshb file
     * @param[in] geomodel_entity_mesh the GeoModelMeshEntity you want to save
     * @param[in] zf the zip file
     */
    template< typename ENTITY >
    void save_geomodel_mesh_entity(
        const ENTITY& geomodel_entity_mesh,
        std::vector< std::string >& filenames )
    {
        std::string name = build_string_for_geomodel_entity_export(
            geomodel_entity_mesh );
        if( save_mesh( geomodel_entity_mesh, name ) ) {
#pragma omp critical
            {
                filenames.push_back( name );
            }
        }
    }

    void zip_files( const std::vector< std::string >& filenames, zipFile& zf )
    {
        for( const std::string& name : filenames ) {
            zip_file( zf, name );
            GEO::FileSystem::delete_file( name );
        }
    }

    template< template< index_t > class ENTITY, index_t DIMENSION >
    void save_geomodel_mesh_entities(
        const GeoModel< DIMENSION >& geomodel,
        std::vector< std::string >& filenames )
    {
        const MeshEntityType& type = ENTITY< DIMENSION >::type_name_static();
        GEO::Logger* logger = Logger::instance();
        bool logger_status = logger->is_quiet();
        logger->set_quiet( true );
        RINGMESH_PARALLEL_LOOP_DYNAMIC
        for( index_t e = 0; e < geomodel.nb_mesh_entities( type ); e++ ) {
            const ENTITY< DIMENSION >& entity =
                dynamic_cast< const ENTITY< DIMENSION >& >( geomodel.mesh_entity(
                    type, e ) );
            save_geomodel_mesh_entity< ENTITY< DIMENSION > >( entity, filenames );
        }
        logger->set_quiet( logger_status );
    }

    template< index_t DIMENSION >
    void save_all_geomodel_mesh_entities_base(
        const GeoModel< DIMENSION >& geomodel,
        std::vector< std::string >& filenames )
    {
        save_geomodel_mesh_entities< Corner, DIMENSION >( geomodel, filenames );
        save_geomodel_mesh_entities< Line, DIMENSION >( geomodel, filenames );
        save_geomodel_mesh_entities< Surface, DIMENSION >( geomodel, filenames );
    }

    template< index_t DIMENSION >
    void save_all_geomodel_mesh_entities(
        const GeoModel< DIMENSION >& geomodel,
        std::vector< std::string >& filenames );

    template< >
    void save_all_geomodel_mesh_entities(
        const GeoModel2D& geomodel,
        std::vector< std::string >& filenames )
    {
        save_all_geomodel_mesh_entities_base( geomodel, filenames );
    }
    template< >
    void save_all_geomodel_mesh_entities(
        const GeoModel3D& geomodel,
        std::vector< std::string >& filenames )
    {
        save_all_geomodel_mesh_entities_base( geomodel, filenames );
        save_geomodel_mesh_entities< Region, 3 >( geomodel, filenames );
    }

    template< index_t DIMENSION >
    index_t nb_mesh_entities( const GeoModel< DIMENSION >& geomodel );

    template< >
    index_t nb_mesh_entities( const GeoModel2D& geomodel )
    {
        return geomodel.nb_corners() + geomodel.nb_lines() + geomodel.nb_surfaces();
    }

    template< >
    index_t nb_mesh_entities( const GeoModel3D& geomodel )
    {
        return geomodel.nb_corners() + geomodel.nb_lines() + geomodel.nb_surfaces()
            + geomodel.nb_regions();
    }

    index_t find_dimension( const std::string& mesh_entities_filename )
    {
        GEO::LineInput file_line( mesh_entities_filename );
        while( !file_line.eof() && file_line.get_line() ) {
            file_line.get_fields();
            if( file_line.nb_fields() == 2 ) {
                if( file_line.field_matches( 0, "Dimension" ) ) {
                    return file_line.field_as_int( 1 );
                }
            }
        }
        return 3;
    }

    template< index_t DIMENSION >
    class GeoModelHandlerGM final : public GeoModelIOHandler< DIMENSION > {
    public:
        void load( const std::string& filename, GeoModel< DIMENSION >& geomodel ) final
        {
            std::string pwd = GEO::FileSystem::get_current_working_directory();
            GEO::FileSystem::set_current_working_directory(
                GEO::FileSystem::dir_name( filename ) );
            GeoModelBuilderGM < DIMENSION
                > builder( geomodel, GEO::FileSystem::base_name( filename, false ) );
            builder.build_geomodel();
            GEO::FileSystem::set_current_working_directory( pwd );
        }

        void save(
            const GeoModel< DIMENSION >& geomodel,
            const std::string& filename ) final
        {
            const std::string pwd = GEO::FileSystem::get_current_working_directory();
            bool valid_new_working_directory =
                GEO::FileSystem::set_current_working_directory(
                    GEO::FileSystem::dir_name( filename ) );
            if( !valid_new_working_directory ) {
                throw RINGMeshException( "I/O", "Output directory does not exist" );
            }

            zipFile zf = zipOpen(
                GEO::FileSystem::base_name( filename, false ).c_str(),
                APPEND_STATUS_CREATE );
            ringmesh_assert( zf != nil );

            const std::string mesh_entity_file( "mesh_entities.txt" );
            save_mesh_entities( geomodel, mesh_entity_file );
            zip_file( zf, mesh_entity_file );
            GEO::FileSystem::delete_file( mesh_entity_file );

            const std::string geological_entity_file( "geological_entities.txt" );
            save_geological_entities( geomodel, geological_entity_file );
            zip_file( zf, geological_entity_file );
            GEO::FileSystem::delete_file( geological_entity_file );

            index_t nb_mesh_entites = nb_mesh_entities( geomodel );
            std::vector< std::string > filenames;
            filenames.reserve( nb_mesh_entites );
            save_all_geomodel_mesh_entities( geomodel, filenames );
            std::sort( filenames.begin(), filenames.end() );
            zip_files( filenames, zf );

            zipClose( zf, NULL );
            GEO::FileSystem::set_current_working_directory( pwd );
        }

        index_t dimension( const std::string& filename ) const final
        {
            unzFile uz = unzOpen( filename.c_str() );
            unz_global_info global_info;
            if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
                unzClose( uz );
                throw RINGMeshException( "ZLIB", "Could not read file global info" );
            }
            const std::string mesh_entity_file( "mesh_entities.txt" );
            unzip_file( uz, mesh_entity_file.c_str() );
            index_t dimension = find_dimension( mesh_entity_file );
            bool ok = GEO::FileSystem::delete_file( mesh_entity_file );
            ringmesh_unused( ok );
            return dimension;
        }

        virtual ~GeoModelHandlerGM() = default;
    private:
        void save_geomodel_regions(
            const GeoModel< DIMENSION >& geomodel,
            std::vector< std::string >& filenames );
    };

    CLASS_DIMENSION_ALIASES (GeoModelHandlerGM);
}
