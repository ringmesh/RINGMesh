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

#include <ringmesh/io.h>

#include <geogram/basic/file_system.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/string.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/geo_model_builder.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geo_model_api.h>

#include <cstring>
#include <fstream>

/*!
 * @file Please add a file description
 * @author Arnaud Botella
 * @todo Rename the file. Could'nt the 2 functions be moved somewhere else ? 
 */

namespace {
    using namespace RINGMesh ;

    /*!
     * @brief Write in the out stream things to save for CONTACT, INTERFACE and LAYERS
     */
    void save_high_level_bme( std::ofstream& out, const GeoModelEntity& E )
    {
        /// First line:  TYPE - ID - NAME - GEOL
        out << E.gme_id() << " " ;
        if( E.has_name() ) {
            out << E.name() << " " ;
        } else {
            out << "no_name " ;
        }
        out << GeoModelEntity::geol_name( E.geological_feature() ) << std::endl ;

        /// Second line:  IDS of children
        for( index_t j = 0; j < E.nb_children(); ++j ) {
            out << " " << E.child_id( j ).index ;
        }
        out << std::endl ;
    }

    void save_topology( const GeoModel& M, const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() ) ;
        out.precision( 16 ) ;
        if( out.bad() ) {
            throw RINGMeshException( "I/O",
                "Error when opening the file: " + file_name ) ;
        }

        out << "RINGMESH BOUNDARY MODEL" << std::endl ;
        out << "NAME " << M.name() << std::endl ;

        // Numbers of the different types of entities
        for( index_t i = GME::CORNER; i < GME::NO_TYPE; i++ ) {
            GME::TYPE type = static_cast< GME::TYPE >( i ) ;
            out << "NB_" << GME::type_name( type ) << " " << M.nb_entities( type )
                << std::endl ;
        }
        // Write high-level entities
        for( index_t i = GME::CONTACT; i < GME::NO_TYPE; i++ ) {
            GME::TYPE type = static_cast< GME::TYPE >( i ) ;
            index_t nb = M.nb_entities( type ) ;
            for( index_t j = 0; j < nb; ++j ) {
                save_high_level_bme( out, M.entity( GME::gme_t( type, j ) ) ) ;
            }
        }
        // Regions
        for( index_t i = 0; i < M.nb_regions(); ++i ) {
            const Region& E = M.region( i ) ;
            // Save ID - NAME
            out << E.gme_id() << " " ;
            if( E.has_name() ) {
                out << E.name() ;
            } else {
                out << "no_name" ;
            }
            out << std::endl ;
            // Second line Signed ids of boundary surfaces
            for( index_t j = 0; j < E.nb_boundaries(); ++j ) {
                if( E.side( j ) ) {
                    out << "+" ;
                } else {
                    out << "-" ;
                }
                out << E.boundary_gme( j ).index << " " ;
            }
            out << std::endl ;
        }

        // Universe
        out << "UNIVERSE " << std::endl ;
        for( index_t j = 0; j < M.universe().nb_boundaries(); ++j ) {
            if( M.universe().side( j ) ) {
                out << "+" ;
            } else {
                out << "-" ;
            }
            out << M.universe().boundary_gme( j ).index << " " ;
        }
        out << std::endl ;
    }

    /*!
     * @brief Save the GeoModel into a dedicated format bm
     * @todo Write the description of the BM format
     * @todo We need a generic read/write for the attributes !!
     */
    void save_bm_file( const GeoModel& M, const std::string& file_name )
    {
        save_topology( M, file_name ) ;
        std::ofstream out( file_name.c_str(), std::ios::out | std::ios::app ) ;

        out.precision( 16 ) ;

        // Corners
        for( index_t i = 0; i < M.nb_corners(); ++i ) {
            out << M.corner( i ).gme_id() << " " << M.corner( i ).vertex(0)
                << std::endl ;
        }
        // Lines
        for( index_t i = 0; i < M.nb_lines(); ++i ) {
            const Line& L = M.line( i ) ;
            out << L.gme_id() << std::endl ;
            out << "LINE_VERTICES " << L.nb_vertices() << std::endl ;
            for( index_t j = 0; j < L.nb_vertices(); ++j ) {
                out << L.vertex( j ) << std::endl ;
            }
            out << "IN_BOUNDARY " ;
            for( index_t j = 0; j < L.nb_in_boundary(); ++j ) {
                out << L.in_boundary_gme( j ).index << " " ;
            }
            out << std::endl ;
        }

        // Surfaces
        for( index_t i = 0; i < M.nb_surfaces(); ++i ) {
            const Surface& S = M.surface( i ) ;
            out << S.gme_id() << std::endl ;
            out << "SURFACE_VERTICES " << S.nb_vertices() << std::endl ;
            for( index_t j = 0; j < S.nb_vertices(); ++j ) {
                out << S.vertex( j ) << std::endl ;
            }

            //to remove or porte
        //    out << "SURFACE_CORNERS " << S.nb_facet_corners() << std::endl ;
        /*    out << "SURFACE_FACETS " << S.nb_mesh_element() << std::endl ;
            for( index_t j = 0; j < S.nb_mesh_element(); ++j ) {
                out << S.nb_mesh_element_vertices( j ) << " " ;
                for( index_t v = 0; v < S.nb_mesh_element_vertices( j ); ++v ) {
                    out << S.vertex_indexsurf_vertex_id( j, v ) << " " ;
                }
                out << std::endl ;
            }*/
        }
    }
}

namespace RINGMesh {

    /*!
    * Compares the contains of two files
    * @param[in] f1 the first filename
    * @param[in] f2 the second filename
    * @return return True if the files are identical
    */
    bool compare_files( const std::string& f1, const std::string& f2 )
    {
        const unsigned int MAX_LINE_LEN = 65535 ;

        std::ifstream lFile( f1.c_str() ) ;
        std::ifstream rFile( f2.c_str() ) ;

        char* lBuffer = new char[ MAX_LINE_LEN ]() ;
        char* rBuffer = new char[ MAX_LINE_LEN ]() ;

        do {
            lFile.read( lBuffer, MAX_LINE_LEN ) ;
            rFile.read( rBuffer, MAX_LINE_LEN ) ;
            size_t numberOfRead = static_cast< size_t >( lFile.gcount() ) ;

            if( std::memcmp( lBuffer, rBuffer, numberOfRead ) != 0 ) {
                delete[] lBuffer ;
                delete[] rBuffer ;
                return false ;
            }
        } while( lFile.good() || rFile.good() ) ;
        delete[] lBuffer ;
        delete[] rBuffer ;
        return true ;
    }

    void mesh_initialize()
    {
        GeoModelIOHandler::initialize_full_geomodel_output() ;
        GeoModelIOHandler::initialize_boundary_geomodel_output() ;
        WellGroupIOHandler::initialize() ;
    }


    /***************************************************************************/



    void zip_file( zipFile zf, const std::string& name )
    {
        zip_fileinfo zfi = {} ;
        std::fstream file( name.c_str(), std::ios::in | std::ios::binary ) ;
        file.seekg( 0, std::ios::end ) ;
        long size = file.tellg() ;
        file.seekg( 0, std::ios::beg ) ;
        std::vector< char > buffer( size ) ;
        file.read( &buffer[0], size ) ;
        zipOpenNewFileInZip( zf, name.c_str(), &zfi,
        NULL, 0, NULL, 0, NULL, Z_DEFLATED, Z_DEFAULT_COMPRESSION ) ;
        zipWriteInFileInZip( zf, size == 0 ? "" : &buffer[0], size ) ;
        zipCloseFileInZip( zf ) ;
        file.close() ;
    }

    void unzip_file( unzFile uz, char filename[MAX_FILENAME] )
    {
        char read_buffer[ READ_SIZE] ;
        unz_file_info file_info ;
        if( unzGetCurrentFileInfo( uz, &file_info, filename,
        MAX_FILENAME,
        NULL, 0, NULL, 0 ) != UNZ_OK ) {
            unzClose( uz ) ;
            throw RINGMeshException( "ZLIB", "Could not read file global info" ) ;
        }
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

    /***************************************************************************/



    GeoModelIOHandler* GeoModelIOHandler::create( const std::string& format )
    {
        GeoModelIOHandler* handler = IOHandlerFactory::create_object(
            format ) ;
        if( !handler ) {
            std::vector< std::string > names ;
            IOHandlerFactory::list_creators( names ) ;
            GEO::Logger::err( "I/O" ) << "Currently supported file formats are: " ;
            for( index_t i = 0; i < names.size(); i++ ) {
                GEO::Logger::err( "I/O" ) << " " << names[i] ;
            }
            GEO::Logger::err( "I/O" ) << std::endl ;

            throw RINGMeshException( "I/O", "Unsupported file format: " + format ) ;
        }
        return handler ;
    }

    GeoModelIOHandler* GeoModelIOHandler::get_handler(
        const std::string& filename )
    {
        std::string ext = GEO::FileSystem::extension( filename ) ;
        return create( ext ) ;
    }


    void BMIOHandler::load( const std::string& filename, GeoModel& model )
        {
            std::ifstream input( filename.c_str() ) ;
            if( !input ) {
                throw RINGMeshException( "I/O", "Failed to open file " + filename ) ;
            }
            GeoModelBuilderBM builder( model, filename ) ;
            builder.build_model() ;
            GEO::Logger::out( "I/O" ) << " Loaded model " << model.name() << " from "
                << filename << std::endl ;
            print_geomodel( model ) ;
            is_geomodel_valid( model ) ;
        }

    void BMIOHandler::save( const GeoModel& model, const std::string& filename )
    {
        save_bm_file( model, filename ) ;
    } ;


    void MMIOHandler::load( const std::string& filename, GeoModel& gm )
        {/*
            unzFile uz = unzOpen( filename.c_str() ) ;
            unz_global_info global_info ;
            if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
                unzClose( uz ) ;
                throw RINGMeshException( "ZLIB",
                    "Could not read file global info" ) ;
            }
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                char filename[MAX_FILENAME] ;
                unzip_file( uz, filename ) ;
                GEO::MeshIOFlags flags ;
                flags.set_entity( GEO::MESH_FACETS ) ;
                flags.set_entity( GEO::MESH_CELLS ) ;
                flags.set_entity( GEO::MESH_EDGES ) ;
                flags.set_attribute( GEO::MESH_FACET_REGION ) ;
                GEO::Mesh& m = gm.region( r ).mesh() ;
                std::string ext = GEO::FileSystem::extension( filename ) ;
                if( ext == "meshb" ) {
                    GEO::Logger::instance()->set_minimal( true ) ;
                    GEO::mesh_load( GEO::String::to_string( filename ), m, flags ) ;
                    GEO::Logger::instance()->set_minimal( false ) ;
                } else {
                    ringmesh_assert_not_reached;
                }
                GEO::FileSystem::delete_file( filename ) ;  // WHY ?? [Jeanne]

                if( ( r + 1 ) < global_info.number_entry ) {
                    if( unzGoToNextFile( uz ) != UNZ_OK ) {
                        unzClose( uz ) ;
                        throw RINGMeshException( "ZLIB",
                            "Could not read next file" ) ;
                    }
                }
            }
            unzClose( uz ) ;*/
        }

    void MMIOHandler::save( const GeoModel& gm, const std::string& filename )
        {/*
            std::string pwd = GEO::FileSystem::get_current_working_directory() ;
            GEO::FileSystem::set_current_working_directory(
                GEO::FileSystem::dir_name( filename ) ) ;
            zipFile zf = zipOpen( filename.c_str(), APPEND_STATUS_CREATE ) ;
            for( index_t m = 0; m < gm.nb_regions(); m++ ) {
                GEO::MeshIOFlags flags ;
                flags.set_entity( GEO::MESH_FACETS ) ;
                flags.set_entity( GEO::MESH_CELLS ) ;
                flags.set_entity( GEO::MESH_EDGES ) ;
                flags.set_attribute( GEO::MESH_FACET_REGION ) ;

                const GEO::Mesh& cur_mesh = gm.region( m ).mesh() ;
                std::string name_mesh_file = "region_" + GEO::String::to_string( m )
                    + ".meshb" ;

                GEO::Logger::instance()->set_quiet( true ) ;
                GEO::mesh_save( cur_mesh, name_mesh_file, flags ) ;
                GEO::Logger::instance()->set_quiet( false ) ;

                zip_file( zf, name_mesh_file ) ;

                GEO::FileSystem::delete_file( name_mesh_file ) ;

            }
            zipClose( zf, NULL ) ;
            GEO::FileSystem::set_current_working_directory( pwd ) ;*/
        }
}
