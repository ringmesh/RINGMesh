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
    void save_geological_entity(
        std::ofstream& out,
        const GeoModelGeologicalEntity& E )
    {
        /// First line:  TYPE - ID - NAME - GEOL
        out << E.gmge() << " " << E.name() << " ";
        out << GeoModelGeologicalEntity::geol_name( E.geological_feature() )
            << std::endl;

        /// Second line:  IDS of children
        for( index_t j = 0; j < E.nb_children(); ++j ) {
            out << E.child_gmme( j ).index() << " ";
        }
        out << std::endl;
    }

    /*!
     * @brief Save the connectivity of a GeoModel in a file
     * @param[in] M the GeoModel
     * @param[in] file_name path to the input file
     */
    void save_geological_entities( const GeoModel& M, const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() );
        out.precision( PRECISION );
        if( out.bad() ) {
            throw RINGMeshException( "I/O",
                "Error when opening the file: " + file_name );
        }

        if( M.nb_geological_entity_types() == 0 ) {
            // Compression of an empty files crashes ? (in debug on windows at least)
            out << "No geological entity in the geomodel" << std::endl;
            return;
        }
        for( index_t i = 0; i < M.nb_geological_entity_types(); i++ ) {
            const std::string& type = M.geological_entity_type( i );
            index_t nb = M.nb_geological_entities( type );
            out << "Nb " << type << " " << nb << std::endl;
        }
        for( index_t i = 0; i < M.nb_geological_entity_types(); i++ ) {
            const std::string& type = M.geological_entity_type( i );
            index_t nb = M.nb_geological_entities( type );
            for( index_t j = 0; j < nb; ++j ) {
                save_geological_entity( out, M.geological_entity( type, j ) );
            }
        }
    }

    template< typename ENTITY >
    void save_mesh_entities_of_type( const GeoModel& M, std::ofstream& out )
    {
        const std::string& type = ENTITY::type_name_static();
        for( index_t e = 0; e < M.nb_mesh_entities( type ); e++ ) {
            const ENTITY& cur_mesh_entity =
                dynamic_cast< const ENTITY& >( M.mesh_entity( type, e ) );
            out << type << " " << e << " " << cur_mesh_entity.name() << " "
                << cur_mesh_entity.low_level_mesh_storage().type_name() << std::endl;
            out << "boundary ";
            for( index_t b = 0; b < cur_mesh_entity.nb_boundaries(); b++ ) {
                out << cur_mesh_entity.boundary_gmme( b ).index() << " ";
            }
            out << std::endl;
        }
    }

    /*!
     * @brief Save the topology of a GeoModelin a file
     * @param[in] M the GeoModel
     * @param[in] file_name the output file name
     */
    void save_mesh_entities( const GeoModel& M, const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() );
        out.precision( PRECISION );
        if( out.bad() ) {
            throw RINGMeshException( "I/O",
                "Error when opening the file: " + file_name );
        }

        out << "Version 2" << std::endl;
        out << "GeoModel name " << M.name() << std::endl;

        // Numbers of the different types of mesh entities
        out << "Nb " << Corner::type_name_static() << " " << M.nb_corners()
            << std::endl;
        out << "Nb " << Line::type_name_static() << " " << M.nb_lines() << std::endl;
        out << "Nb " << Surface::type_name_static() << " " << M.nb_surfaces()
            << std::endl;
        out << "Nb " << Region::type_name_static() << " " << M.nb_regions()
            << std::endl;

        save_mesh_entities_of_type< Corner >( M, out );
        save_mesh_entities_of_type< Line >( M, out );
        save_mesh_entities_of_type< Surface >( M, out );

        // Regions
        for( index_t i = 0; i < M.nb_regions(); ++i ) {
            const Region& E = M.region( i );
            // Save ID - NAME
            out << Region::type_name_static() << " " << i << " " << E.name() << " "
                << E.low_level_mesh_storage().type_name() << std::endl;
            // Second line Signed ids of boundary surfaces
            for( index_t j = 0; j < E.nb_boundaries(); ++j ) {
                if( E.side( j ) ) {
                    out << "+";
                } else {
                    out << "-";
                }
                out << E.boundary_gmme( j ).index() << " ";
            }
            out << std::endl;
        }

        // Universe
        out << "Universe " << std::endl;
        for( index_t j = 0; j < M.universe().nb_boundaries(); ++j ) {
            if( M.universe().side( j ) ) {
                out << "+";
            } else {
                out << "-";
            }
            out << M.universe().boundary_gmme( j ).index() << " ";
        }
        out << std::endl;
    }

    bool save_mesh(
        const GeoModelMeshEntity& geomodel_entity_mesh,
        const std::string& name )
    {
        if( geomodel_entity_mesh.type_name() == Region::type_name_static() ) {
            const Region& region = geomodel_entity_mesh.geomodel().region(
                geomodel_entity_mesh.index() );
            if( !region.is_meshed() ) {
                // a region is not necessary meshed.
                return false;
            }
        }
        geomodel_entity_mesh.save( name );
        return true;
    }

    template< typename ENTITY >
    std::string build_string_for_geomodel_entity_export( const ENTITY& entity )
    {
        const gmme_id& id = entity.gmme();
        std::string base_name = static_cast< std::string >( id.type() ) + "_"
            + GEO::String::to_string( id.index() );
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
        for( index_t i = 0; i < filenames.size(); i++ ) {
            const std::string& name = filenames[i];
            zip_file( zf, name );
            GEO::FileSystem::delete_file( name );
        }
    }

    template< typename ENTITY >
    void save_geomodel_mesh_entities(
        const GeoModel& geomodel,
        std::vector< std::string >& filenames )
    {
        const std::string& type = ENTITY::type_name_static();
        GEO::Logger* logger = Logger::instance();
        bool logger_status = logger->is_quiet();
        logger->set_quiet( true );
        RINGMESH_PARALLEL_LOOP_DYNAMIC
        for( index_t e = 0; e < geomodel.nb_mesh_entities( type ); e++ ) {
            const ENTITY& entity =
                dynamic_cast< const ENTITY& >( geomodel.mesh_entity( type, e ) );
            save_geomodel_mesh_entity< ENTITY >( entity, filenames );
        }
        logger->set_quiet( logger_status );
    }

    class GeoModelHandlerGM final: public GeoModelIOHandler {
    public:
        virtual bool load( const std::string& filename, GeoModel& geomodel ) final
        {
            std::string pwd = GEO::FileSystem::get_current_working_directory();
            GEO::FileSystem::set_current_working_directory(
                GEO::FileSystem::dir_name( filename ) );
            GeoModelBuilderGM builder( geomodel,
                GEO::FileSystem::base_name( filename, false ) );
            builder.build_geomodel();
            Logger::out( "I/O", " Loaded geomodel ", geomodel.name(), " from ",
                filename );
            print_geomodel( geomodel );
            bool is_valid = is_geomodel_valid( geomodel );
            GEO::FileSystem::set_current_working_directory( pwd );
            return is_valid;

        }
        virtual void save( const GeoModel& geomodel, const std::string& filename ) final
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

            index_t nb_mesh_entites = geomodel.nb_corners() + geomodel.nb_lines()
                + geomodel.nb_surfaces() + geomodel.nb_regions();
            std::vector< std::string > filenames;
            filenames.reserve( nb_mesh_entites );
            save_geomodel_mesh_entities< Corner >( geomodel, filenames );
            save_geomodel_mesh_entities< Line >( geomodel, filenames );
            save_geomodel_mesh_entities< Surface >( geomodel, filenames );
            save_geomodel_mesh_entities< Region >( geomodel, filenames );
            std::sort( filenames.begin(), filenames.end() );
            zip_files( filenames, zf );

            zipClose( zf, NULL );
            GEO::FileSystem::set_current_working_directory( pwd );
        }
    };

}
