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

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/string.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_builder.h>

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

        //    ___                   _               __  __         _     _
        //   | _ ) ___ _  _ _ _  __| |__ _ _ _ _  _|  \/  |___  __| |___| |
        //   | _ \/ _ \ || | ' \/ _` / _` | '_| || | |\/| / _ \/ _` / -_) |
        //   |___/\___/\_,_|_||_\__,_\__,_|_|  \_, |_|  |_\___/\__,_\___|_|
        //                                     |__/

        bool load( const std::string& filename, BoundaryModel& model )
        {
            GEO::Logger::out( "I/O" ) << "Loading file " << filename << "..."
                << std::endl ;

            BoundaryModelIOHandler_var handler = BoundaryModelIOHandler::get_handler(
                filename ) ;
            if( handler && handler->load( filename, model ) ) {
                return true ;
            }

            GEO::Logger::err( "I/O" ) << "Could not load file: " << filename
                << std::endl ;
            return false ;
        }

        bool save( BoundaryModel& model, const std::string& filename )
        {
            GEO::Logger::out( "I/O" ) << "Saving file " << filename << "..."
                << std::endl ;

            BoundaryModelIOHandler_var handler = BoundaryModelIOHandler::get_handler(
                filename ) ;
            if( handler && handler->save( model, filename ) ) {
                return true ;
            }

            GEO::Logger::err( "I/O" ) << "Could not save file: " << filename
                << std::endl ;
            return false ;
        }

        /************************************************************************/

        class MLIOHandler: public BoundaryModelIOHandler {
        public:
            virtual bool load( const std::string& filename, BoundaryModel& model )
            {
                if( filename.empty() ) {
                    GEO::Logger::err( "I/O" )
                        << "No filename provided for structural model, use in:model"
                        << std::endl ;
                    return false ;
                }

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

            virtual bool save( BoundaryModel& model, const std::string& filename )
            {
                std::ofstream out( filename.c_str() ) ;
                return model.save_gocad_model3d( out ) ;
            }
        } ;

        /************************************************************************/

        BoundaryModelIOHandler* BoundaryModelIOHandler::create(
            const std::string& format )
        {
            grgmesh_register_BoundaryModelIOHandler_creator( MLIOHandler, "ml" ) ;

            BoundaryModelIOHandler* handler = BoundaryModelIOHandlerFactory::create_object(format) ;
            if( handler ) {
                return handler ;
            }

            GEO::Logger::err("I/O")
                << "Unsupported file format: " << format
                << std::endl ;
            return nil ;
        }

        BoundaryModelIOHandler* BoundaryModelIOHandler::get_handler(
            const std::string& filename )
        {
            std::string ext = GEO::FileSystem::extension( filename ) ;
            return create( ext ) ;
        }

        //    __  __                 __  __        _
        //   |  \/  |__ _ __ _ _ ___|  \/  |___ __| |_
        //   | |\/| / _` / _| '_/ _ \ |\/| / -_|_-< ' \
        //   |_|  |_\__,_\__|_| \___/_|  |_\___/__/_||_|
        //

        bool load( const std::string& filename, MacroMesh& model )
        {
            GEO::Logger::out( "I/O" ) << "Loading file " << filename << "..."
                << std::endl ;

            MacroMeshIOHandler_var handler = MacroMeshIOHandler::get_handler(
                filename ) ;
            if( handler && handler->load( filename, model ) ) {
                return true ;
            }

            GEO::Logger::err( "I/O" ) << "Could not load file: " << filename
                << std::endl ;
            return false ;
        }

        bool save( const MacroMesh& model, const std::string& filename )
        {
            GEO::Logger::out( "I/O" ) << "Saving file " << filename << "..."
                << std::endl ;

            MacroMeshIOHandler_var handler = MacroMeshIOHandler::get_handler(
                filename ) ;
            if( handler && handler->save( model, filename ) ) {
                return true ;
            }

            GEO::Logger::err( "I/O" ) << "Could not save file: " << filename
                << std::endl ;
            return false ;
        }


        /************************************************************************/

        class MMIOHandler: public MacroMeshIOHandler {
        public:
            virtual bool load( const std::string& filename, MacroMesh& mm )
            {
                unzFile uz = unzOpen( filename.c_str() ) ;
                unz_global_info global_info ;
                if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
                    GEO::Logger::err( "could not read file global info" ) ;
                    unzClose( uz ) ;
                    return false ;
                }
                for( index_t r = 0; r < mm.model().nb_regions(); r++ ) {
                    char filename[MAX_FILENAME] ;
                    unzip_file( uz, filename ) ;
                    GEO::MeshIOFlags flags ;
                    flags.set_element( GEO::MESH_FACETS ) ;
                    flags.set_element( GEO::MESH_CELLS ) ;
                    flags.set_attribute( GEO::MESH_FACET_REGION ) ;
                    GEO::Mesh& m = mm.mesh( r ) ;
                    std::string ext = GEO::FileSystem::extension( filename ) ;
                    if( ext == "meshb" ) {
                        GEO::Logger::instance()->set_quiet( true ) ;
                        GEO::mesh_load( GEO::String::to_string( filename ), m,
                            flags ) ;
                        GEO::Logger::instance()->set_quiet( false ) ;
                    } else {
                        grgmesh_assert_not_reached;
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
                return true ;
            }

            /// Save a \param[in] macro mesh in a .zip file which contains all the mesh file. Type of the export is
            /// determined by the extension given in \param[in] filename
            virtual bool save( const MacroMesh& mm, const std::string& filename )
            {
                std::string pwd = GEO::FileSystem::get_current_working_directory() ;
                GEO::FileSystem::set_current_working_directory(
                    GEO::FileSystem::dir_name( filename ) ) ;
                zipFile zf = zipOpen( filename.c_str(), APPEND_STATUS_CREATE ) ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    GEO::MeshIOFlags flags ;
                    flags.set_element( GEO::MESH_FACETS ) ;
                    flags.set_element( GEO::MESH_CELLS ) ;
                    flags.set_attribute( GEO::MESH_FACET_REGION ) ;

                    const GEO::Mesh& cur_mesh = mm.mesh( m ) ;
                    std::string name_mesh_file = "region_"
                        + GEO::String::to_string( m ) + ".meshb" ;

                    GEO::Logger::instance()->set_quiet( true ) ;
                    GEO::mesh_save( cur_mesh, name_mesh_file, flags ) ;
                    GEO::Logger::instance()->set_quiet( false ) ;

                    zip_file( zf, name_mesh_file ) ;

                    GEO::FileSystem::delete_file( name_mesh_file ) ;

                }
                zipClose( zf, NULL ) ;
                GEO::FileSystem::set_current_working_directory( pwd ) ;
                return true ;

            }
        } ;

        /************************************************************************/

        class MESHBIOHandler: public MacroMeshIOHandler {
        public:
            virtual bool load( const std::string& filename, MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                    << "Loading of a MacroMesh from a meshb not implemented yet"
                    << std::endl ;
                return false ;
            }
            virtual bool save( const MacroMesh& mm, const std::string& filename )
            {
                GEO::Mesh mesh( 3 ) ;
                GEO::MeshBuilder builder( &mesh ) ;
                builder.begin_mesh() ;
                GEO::MeshMutator::set_attributes( mesh, GEO::MESH_FACET_REGION ) ;

                std::vector< vec3 > unique_vertices ;
                std::vector< index_t > indices ;
                mm.unique_points( unique_vertices, indices ) ;
                for( index_t p = 0; p < unique_vertices.size(); p++ ) {
                    builder.add_vertex( unique_vertices[p] ) ;
                }

                index_t vertex_offset = 0 ;
                index_t cell_offset = 0 ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& cur_mesh = mm.mesh( m ) ;
                    for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                        builder.begin_facet() ;
                        for( index_t v = cur_mesh.facet_begin( f );
                            v < cur_mesh.facet_end( f ); v++ ) {
                            builder.add_vertex_to_facet(
                                indices[vertex_offset + cur_mesh.corner_vertex_index( v )] ) ;
                        }
                        builder.end_facet( cur_mesh.facet_region( f ) ) ;
                    }

                    for( index_t c = 0; c < cur_mesh.nb_cells(); c++ ) {
                        index_t vertex_indices[8] ;
                        for( unsigned int v = 0; v < cur_mesh.cell_nb_vertices( c );
                            v++ ) {
                            vertex_indices[v] = indices[vertex_offset
                                +cur_mesh.cell_vertex_index( c, v )] ;
                        }
                        signed_index_t adj_indices[6] ;
                        for( index_t f = 0; f < cur_mesh.cell_nb_facets( c ); f++ ) {
                            signed_index_t adj = cur_mesh.cell_adjacent( c, f )  ;
                            adj_indices[f] = adj == -1 ? adj : cell_offset + adj ;
                        }
                        mesh.add_cell( cur_mesh.cell_type( c ), vertex_indices, adj_indices ) ;
                    }
                    vertex_offset += cur_mesh.nb_vertices() ;
                    cell_offset += cur_mesh.nb_cells() ;
                }

                builder.end_mesh( false ) ;

                GEO::MeshIOFlags flags ;
                flags.set_element( GEO::MESH_FACETS ) ;
                flags.set_element( GEO::MESH_CELLS ) ;
                flags.set_attribute( GEO::MESH_FACET_REGION ) ;
                GEO::Logger::instance()->set_quiet( true ) ;
                GEO::mesh_save( mesh, filename, flags ) ;
                GEO::Logger::instance()->set_quiet( false ) ;

                return true ;
            }
        } ;

        /************************************************************************/

        MacroMeshIOHandler* MacroMeshIOHandler::create(
            const std::string& format )
        {
            grgmesh_register_MacroMeshIOHandler_creator( MMIOHandler, "mm" ) ;
            grgmesh_register_MacroMeshIOHandler_creator( MESHBIOHandler, "meshb" ) ;

            MacroMeshIOHandler* handler = MacroMeshIOHandlerFactory::create_object(format) ;
            if( handler ) {
                return handler ;
            }

            GEO::Logger::err("I/O")
                << "Unsupported file format: " << format
                << std::endl ;
            return nil ;
        }

        MacroMeshIOHandler* MacroMeshIOHandler::get_handler(
            const std::string& filename )
        {
            std::string ext = GEO::FileSystem::extension( filename ) ;
            return create( ext ) ;
        }

        MacroMeshExport::MacroMeshExport( const MacroMesh& mm )
            : mm_( mm )
        {
        }

    }

}
