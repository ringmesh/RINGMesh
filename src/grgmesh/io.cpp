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
#include <geogram/basic/line_stream.h>

#include <geogram/mesh/mesh_private.h>


#include <third_party/zlib/zip.h>
#include <third_party/zlib/unzip.h>

#include <iostream>
#include <fstream>
#include <string>

namespace GRGMesh {
    namespace GRGMeshIO {
        void zip_file( zipFile zf, const std::string& name )
        {
            zip_fileinfo zfi = { 0 } ;
            std::fstream file( name.c_str(), std::ios::in ) ;
            file.seekg( 0, std::ios::end ) ;
            long size = file.tellg() ;
            std::cerr << "SIZE " << size << " - " << name << std::endl ;
            file.seekg( 0, std::ios::beg ) ;
            std::vector< char > buffer( size ) ;
            file.read( &buffer[0], size ) ;
            zipOpenNewFileInZip( zf, name.c_str(), &zfi,
            NULL, 0, NULL, 0, NULL, Z_DEFLATED, Z_DEFAULT_COMPRESSION ) ;
            zipWriteInFileInZip( zf, size == 0 ? "" : &buffer[0], size ) ;
            zipCloseFileInZip( file ) ;
            file.close() ;
        }

        bool unzip_file( unzFile uz, char filename[MAX_FILENAME] )
        {
            char read_buffer[ READ_SIZE] ;
            unz_file_info file_info ;
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
        }
        bool load_BoundaryModel_from_Model3D(
            const std::string& filename,
            BoundaryModel& model )
        {
            std::ifstream input( filename.c_str() ) ;
            if( !input ) {
                GEO::Logger::err( "I/O" ) << "Cannot open file : " << filename
                    << std::endl ;
                return false ;
            }

            BoundaryModelBuilder builder( model ) ;
            builder.load_file( input ) ;
            return true ;
        }


        bool save_BoundaryModel(
            BoundaryModel& model,
            const std::string& filename )
        {
            std::ofstream out( filename.c_str() ) ;
           return model.save_gocad_model3d( out ) ;
        }

        /// Save a \param[in] macro mesh in a .zip file which contains all the mesh file. Type of the export is
        /// determined by the extension given in \param[in] filename
        bool save_macro_mesh( const MacroMesh& mm, const std::string& filename )
        {
            std::string pwd = GEO::FileSystem::get_current_working_directory() ;
            GEO::FileSystem::set_current_working_directory(
                GEO::FileSystem::dir_name( filename ) ) ;
            zipFile zf = zipOpen( filename.c_str(), APPEND_STATUS_CREATE ) ;
            for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                GEO::MeshIOFlags flags ;
                flags.set_element( GEO::MESH_CELLS ) ;
                const GEO::Mesh& cur_mesh = mm.mesh( m ) ;
                std::string name_mesh_file = GEO::String::to_string( m ) + ".meshb" ;
                std::string name_facet_file = GEO::String::to_string( m )
                    + ".facets" ;

                GEO::mesh_save( cur_mesh, name_mesh_file, flags ) ;
                std::ofstream out( name_facet_file.c_str(), std::ios::out ) ;

                out << cur_mesh.nb_facets() << std::endl ;
                for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                    out << cur_mesh.facet_region( f ) << std::endl ;
                }
                out.close() ;

                zip_file( zf, name_mesh_file ) ;
                zip_file( zf, name_facet_file ) ;

                GEO::FileSystem::delete_file( name_mesh_file ) ;
                GEO::FileSystem::delete_file( name_facet_file ) ;

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
            for( index_t r = 0; r < mm.model()->nb_regions(); r++ ) {
                char filename[MAX_FILENAME] ;
                unzip_file( uz, filename ) ;
                GEO::MeshIOFlags flags ;
                flags.set_element( GEO::MESH_CELLS ) ;
                GEO::Mesh& m = mm.mesh( r ) ;
                if( !GEO::mesh_load( GEO::String::to_string( filename ), m,
                    flags ) ) {
                    GEO::LineInput line( GEO::String::to_string( filename ) ) ;
                    line.get_line() ; line.get_fields() ;
                    index_t nb_facets = line.field_as_uint( 0 ) ;
                    GEO::vector< signed_index_t >& facet_regions =
                        GEO::MeshMutator::facet_regions( m ) ;
                    facet_regions.resize( nb_facets ) ;
                    for( index_t f = 0; f < nb_facets; f++ ) {
                        line.get_line() ; line.get_fields() ;
                        facet_regions[f] = line.field_as_int( 0 ) ;
                    }
                    GEO::MeshMutator::set_attributes( m, GEO::MESH_FACET_REGION ) ;
                }
                GEO::FileSystem::delete_file( filename ) ;

                if( ( r + 1 ) < global_info.number_entry ) {
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
