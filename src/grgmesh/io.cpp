/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#include <grgmesh/io.h>
#include <grgmesh/boundary_model.h>

#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/string.h>
#include <geogram/basic/file_system.h>

#include <third_party/zlib/zip.h>
#include <third_party/zlib/unzip.h>

#include <iostream>
#include <fstream>
#include <string>

namespace GRGMesh {
    namespace GRGMeshIO {

        bool load_BoundaryModel_from_Model3D(
            const std::string& filename,
            BoundaryModel& model )
        {
            std::ifstream input( filename.c_str() ) ;
            if( !input ) {
                std::cout << "cannot open file:" << filename << std::endl ;
                return false ;
            }

            BoundaryModelBuilder builder( model ) ;
            builder.load_file( input ) ;
            return true ;
        }
        /// Save a \param[in] macro mesh in a .zip file which contains all the mesh file. Type of the export is
        /// determined by the extension given in \param[in] filename
        bool save_macro_mesh( const MacroMesh& mm, const std::string& filename )
        {
            std::string pwd = GEO::FileSystem::get_current_working_directory() ;
            GEO::FileSystem::set_current_working_directory(
                GEO::FileSystem::dir_name( filename ) ) ;
            zipFile zf = zipOpen( filename.c_str(), APPEND_STATUS_CREATE ) ;
            zip_fileinfo zfi = { 0 } ;
            for( index_t i = 0; i < mm.nb_meshes(); i++ ) {
                GEO::MeshIOFlags flags ;
                flags.set_element( GEO::MESH_CELLS ) ;
                const GEO::Mesh& m = mm.mesh( i ) ;
                std::string name = GEO::String::to_string( i ) + ".meshb" ;
                GEO::mesh_save( m, name, flags ) ;
                std::fstream file( name.c_str(), std::ios::in ) ;
                file.seekg( 0, std::ios::end ) ;
                long size = file.tellg() ;
                file.seekg( 0, std::ios::beg ) ;
                std::vector< char > buffer( size ) ;
                file.read( &buffer[0], size ) ;
                zipOpenNewFileInZip( zf, name.c_str(), &zfi,
                NULL, 0, NULL, 0, NULL, Z_DEFLATED, Z_DEFAULT_COMPRESSION ) ;
                zipWriteInFileInZip( zf, size == 0 ? "" : &buffer[0], size ) ;
                zipCloseFileInZip( file ) ;
                file.close() ;
                GEO::FileSystem::delete_file( name ) ;
            }
            zipClose( zf, NULL ) ;
            GEO::FileSystem::set_current_working_directory( pwd ) ;
            return true ;
        }

        bool load_macro_mesh( MacroMesh& mm, const std::string& mesh_file )
        {
            unzFile uz = unzOpen( mesh_file.c_str() ) ;
            unz_global_info global_info ;
            if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
                GEO::Logger::err( "could not read file global info" ) ;
                unzClose( uz ) ;
                return false ;
            }
            char read_buffer[ READ_SIZE] ;
            for( index_t i = 0; i < mm.model()->nb_regions(); i++ ) {
                unz_file_info file_info ;
                char filename[MAX_FILENAME] ;
                if( unzGetCurrentFileInfo( uz, &file_info, filename,
                MAX_FILENAME,
                NULL, 0, NULL, 0 ) != UNZ_OK ) {
                    GEO::Logger::err( "could not read file global info" ) ;
                    unzClose( uz ) ;
                    return false ;
                }
                if( unzOpenCurrentFile( uz ) != UNZ_OK ) {
                    GEO::Logger::err( "could not open file" ) ;
                    unzClose( uz ) ;
                    return false ;
                }
                FILE *out = fopen( filename, "wb" ) ;
                if( out == NULL ) {
                    GEO::Logger::err( "could not open destination file" ) ;
                    unzCloseCurrentFile( uz ) ;
                    unzClose( uz ) ;
                    return false ;
                }
                int error = UNZ_OK ;
                do {
                    error = unzReadCurrentFile( uz, read_buffer, READ_SIZE ) ;
                    if( error < 0 ) {
                        GEO::Logger::err(
                            "Invalid error: " + GEO::String::to_string( error ) ) ;
                        unzCloseCurrentFile( uz ) ;
                        unzClose( uz ) ;
                        return false ;
                    }
                    if( error > 0 ) {
                        fwrite( read_buffer, error, 1, out ) ;
                    }
                } while( error > 0 ) ;
                fclose( out ) ;
                unzCloseCurrentFile( uz ) ;
                GEO::MeshIOFlags flags ;
                flags.set_element( GEO::MESH_CELLS ) ;
                GEO::Mesh& m = mm.mesh( i ) ;
                GEO::mesh_load( GEO::String::to_string( filename ), m, flags ) ;
                GEO::FileSystem::delete_file( filename ) ;

                if( ( i + 1 ) < global_info.number_entry ) {
                    if( unzGoToNextFile( uz ) != UNZ_OK ) {
                        GEO::Logger::err( "Could not read next file" ) ;
                        unzClose( uz ) ;
                        return false ;
                    }
                }
            }
            unzClose( uz ) ;

        }
    }
}
