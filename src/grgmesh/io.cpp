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
#include <geogram/mesh/mesh_geometry.h>

#include <third_party/zlib/zip.h>
#include <third_party/zlib/unzip.h>

#include <iostream>
#include <fstream>
#include <string>
#include <stack>

namespace GRGMesh {
    namespace GRGMeshIO {

        static std::string TAB = "\t" ;
        static std::string SPACE = " " ;

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
            return true ;
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

                BoundaryModelBuilderGocad builder( model ) ;
                builder.load_ml_file( input ) ;
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
                MacroMesh& MM = const_cast <MacroMesh&> (mm) ;
                GEO::MeshMutator::set_attributes( mesh, GEO::MESH_FACET_REGION ) ;
                for( index_t p = 0; p < MM.nb_vertices(); p++ ) {
                    builder.add_vertex( MM.global_vertex(p) ) ;
                }

                index_t cell_offset = 0 ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& cur_mesh = mm.mesh( m ) ;
                    for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                        builder.begin_facet() ;
                        for( index_t v = cur_mesh.facet_begin( f );
                            v < cur_mesh.facet_end( f ); v++ ) {
                            builder.add_vertex_to_facet(MM.global_vertex_id(m,cur_mesh.corner_vertex_index( v )));
                        }
                        builder.end_facet( cur_mesh.facet_region( f ) ) ;
                    }

                    for( index_t c = 0; c < cur_mesh.nb_cells(); c++ ) {
                        index_t vertex_indices[8] ;
                        for( unsigned int v = 0; v < cur_mesh.cell_nb_vertices( c );
                            v++ ) {
                            vertex_indices[v] = MM.global_vertex_id(m,
                                cur_mesh.cell_vertex_index( c, v )) ;
                        }
                        signed_index_t adj_indices[6] ;
                        for( index_t f = 0; f < cur_mesh.cell_nb_facets( c ); f++ ) {
                            signed_index_t adj = cur_mesh.cell_adjacent( c, f )  ;
                            adj_indices[f] = adj == -1 ? adj : cell_offset + adj ;
                        }
                        mesh.add_cell( cur_mesh.cell_type( c ), vertex_indices, adj_indices ) ;
                    }
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

        class TetGenIOHandler: public MacroMeshIOHandler {
        public:
            virtual bool load( const std::string& filename, MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                    << "Loading of a MacroMesh from TetGen not implemented yet"
                    << std::endl ;
                return false ;
            }
            virtual bool save( const MacroMesh& mm, const std::string& filename )
            {
                MacroMeshExport db( mm ) ;
                db.compute_database() ;

                std::string directory = GEO::FileSystem::dir_name( filename ) ;
                std::string file = GEO::FileSystem::base_name( filename ) ;

                std::ostringstream oss_node ;
                oss_node << directory << "/" << file << ".node" ;
                std::ofstream node( oss_node.str().c_str() ) ;
                node.precision( 16 ) ;
                node << mm.nb_vertices() << " 3 0 0" << std::endl ;
                for( index_t v = 0; v < mm.nb_vertices(); v++ ) {
                    const vec3& point = mm.global_vertex( v ) ;
                    node << v << SPACE << point.x << SPACE << point.y << SPACE
                        << point.z << std::endl ;
                }

                std::ostringstream oss_ele ;
                oss_ele << directory << "/" << file << ".ele" ;
                std::ofstream ele( oss_ele.str().c_str() ) ;
                std::ostringstream oss_neigh ;
                oss_neigh << directory << "/" << file << ".neigh" ;
                std::ofstream neigh( oss_neigh.str().c_str() ) ;

                ele << db.nb_cells() << " 4 1" << std::endl ;
                neigh << db.nb_cells() << " 4" << std::endl ;
                grgmesh_debug_assert( db.nb_cells() == db.nb_tet() ) ;
                index_t nb_tet_exported = 0 ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& mesh = mm.mesh( m ) ;
                    for( index_t t = 0; t < db.nb_tet( m ); t++ ) {
                        index_t tet = db.local_tet_id( m, t ) ;
                        ele << nb_tet_exported + tet << SPACE
                            << mm.global_vertex_id( m,
                                mesh.cell_vertex_index( tet, 0 ) ) << SPACE
                            << mm.global_vertex_id( m,
                                mesh.cell_vertex_index( tet, 1 ) ) << SPACE
                            << mm.global_vertex_id( m,
                                mesh.cell_vertex_index( tet, 2 ) ) << SPACE
                            << mm.global_vertex_id( m,
                                mesh.cell_vertex_index( tet, 3 ) ) << SPACE << m + 1
                            << std::endl ;
                        neigh << nb_tet_exported + tet ;
                        for( index_t f = 0; f < mesh.cell_nb_facets( tet ); f++ ) {
                            neigh << SPACE << mm.global_cell_adjacent( m, t, f ) ;
                        }
                        neigh << std::endl ;
                    }
                }
                return true ;
            }
        } ;

        /************************************************************************/

        MacroMeshIOHandler* MacroMeshIOHandler::create(
            const std::string& format )
        {
            grgmesh_register_MacroMeshIOHandler_creator( MMIOHandler, "mm" ) ;
            grgmesh_register_MacroMeshIOHandler_creator( MESHBIOHandler, "meshb" ) ;
            grgmesh_register_MacroMeshIOHandler_creator( TetGenIOHandler, "tetgen" ) ;

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

        /************************************************************************/

        MacroMeshExport::MacroMeshExport( const MacroMesh& mm )
            :
                mm_( mm ),
                facet_ptr_( ( NB_FACET_TYPES + 1 ) * mm.model().nb_surfaces(), 0 ),
                mesh_facet_ptr_( mm.model().nb_surfaces() + 1, 0 ),
                cell_ptr_( ( NB_CELL_TYPES + 1 ) * mm.model().nb_regions(), 0 ),
                mesh_cell_ptr_( mm.model().nb_regions() + 1, 0 ),
                mesh_corner_ptr_( mm.model().nb_regions() + 1, 0 ),
                surface2mesh_( mm.model().nb_surfaces(), Surface::NO_ID ),
                first_duplicated_vertex_id_( 0 ),
                nb_triangle_( 0 ),
                nb_quad_( 0 ),
                nb_tet_( 0 ),
                nb_pyramid_( 0 ),
                nb_prism_( 0 ),
                nb_hex_( 0 )
        {
        }

        void MacroMeshExport::compute_database( const DuplicateMode& mode )
        {
            fill_with_geometry() ;
            if( mode != NONE ) {
                duplicate_vertices( mode ) ;
            }
        }

        void MacroMeshExport::fill_with_geometry()
        {
            /// 1 - fill facet information
            index_t facet_access[5] = { -1, -1, -1, 0, 1 } ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& cur_mesh = mm_.mesh( m ) ;
                std::vector< signed_index_t > surface_proccessed ;
                for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                    signed_index_t surface_id = cur_mesh.facet_region( f ) ;
                    if( surface2mesh_[surface_id] != Surface::NO_ID ) continue ;
                    if( !Utils::contains( surface_proccessed, surface_id ) ) {
                        surface_proccessed.push_back( surface_id ) ;
                    }
                    facet_ptr_[(NB_FACET_TYPES+1)*surface_id
                        + facet_access[cur_mesh.facet_size( f )] + 1]++ ;
                }
                for( index_t s = 0; s < surface_proccessed.size(); s++ ) {
                    surface2mesh_[surface_proccessed[s]] = m ;
                }
            }

            for( index_t s = 0; s < mm_.model().nb_surfaces(); s++ ) {
                for( index_t t = 1; t < NB_FACET_TYPES; t++ ) {
                    facet_ptr_[( NB_FACET_TYPES + 1 ) * s + t + 1] +=
                        facet_ptr_[( NB_FACET_TYPES + 1 ) * s + t] ;
                }
                mesh_facet_ptr_[s + 1] += mesh_facet_ptr_[s] ;
                index_t nb_facets_in_surface = facet_ptr_[( NB_FACET_TYPES + 1 ) * s
                    + NB_FACET_TYPES] - facet_ptr_[( NB_FACET_TYPES + 1 ) * s] ;
                mesh_facet_ptr_[s + 1] += nb_facets_in_surface ;
            }
            facets_.resize( mesh_facet_ptr_.back() ) ;

            std::vector< index_t > cur_facet_index_type( NB_FACET_TYPES*mm_.model().nb_surfaces(), 0 ) ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& cur_mesh = mm_.mesh( m ) ;
                for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                    signed_index_t surface_id = cur_mesh.facet_region( f ) ;
                    if( surface2mesh_[surface_id] != m ) continue ;
                    index_t type_access = facet_access[cur_mesh.facet_size( f )] ;
                    facets_[mesh_facet_ptr_[surface_id] + facet_ptr_[(NB_FACET_TYPES+1)*surface_id + type_access]
                        + cur_facet_index_type[NB_FACET_TYPES*surface_id+type_access]++ ] = f ;
                }
            }


            /// 2 - fill cell/corner information
            index_t cell_access[4] = { 0, 3, 2, 1 } ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& mesh = mm_.mesh( m ) ;
                for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                    cell_ptr_[(NB_CELL_TYPES+1) * m + cell_access[mesh.cell_type( c )] + 1]++ ;
                }
                for( index_t type = 1; type < NB_CELL_TYPES; type++ ) {
                    cell_ptr_[(NB_CELL_TYPES+1)  * m + type + 1] += cell_ptr_[(NB_CELL_TYPES+1)  * m + type] ;
                }
            }

            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& mesh = mm_.mesh( m ) ;
                mesh_cell_ptr_[m + 1] += mesh_cell_ptr_[m] + mesh.nb_cells() ;
                mesh_corner_ptr_[m + 1] += mesh_corner_ptr_[m] + mesh.nb_corners() ;
            }
            cells_.resize( mesh_cell_ptr_.back() ) ;
            corners_.resize( mesh_corner_ptr_.back() ) ;

            std::vector< index_t >cur_cell_index_type( NB_CELL_TYPES*mm_.nb_meshes(), 0 ) ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                GEO::Mesh& mesh = const_cast< GEO::Mesh& >( mm_.mesh( m ) ) ;
                for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                    index_t type_access = cell_access[mesh.cell_type( c )] ;
//                    std::cout << "mesh_cell_ptr_[m] " << mesh_cell_ptr_[m] << std::endl ;
//                    std::cout << "cell_ptr_[(NB_CELL_TYPES+1) * m + type_access " << cell_ptr_[(NB_CELL_TYPES+1) * m + type_access] <<std::endl ;
//                    std::cout << "cur_index_type[NB_CELL_TYPES*m+type_access] " << cur_index_type[NB_CELL_TYPES*m+type_access] << std::endl ;
//                    std::cout << "total " << mesh_cell_ptr_[m] + cell_ptr_[(NB_CELL_TYPES+1) * m + type_access]
//                                                                           + cur_index_type[NB_CELL_TYPES*m+type_access] << std::endl ;
//                    std::cout << "size " << cells_.size() << std::endl ;
                    cells_[mesh_cell_ptr_[m] + cell_ptr_[(NB_CELL_TYPES+1) * m + type_access]
                        + cur_cell_index_type[NB_CELL_TYPES*m+type_access]++ ] = c ;
                }
                std::copy( mesh.corner_vertex_index_ptr( 0 ),
                    mesh.corner_vertex_index_ptr( mesh.nb_corners() - 1 ),
                    &corners_[mesh_corner_ptr_[m]] ) ;
            }

            /// 3 - fill vertex information
            first_duplicated_vertex_id_ = mm_.nb_vertex_indices() ;

            /// 4 - update cached values
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                nb_triangle_ += nb_triangle( m ) ;
                nb_quad_ += nb_quad( m ) ;
                nb_tet_ += nb_tet( m ) ;
                nb_pyramid_ += nb_pyramid( m ) ;
                nb_prism_ += nb_prism( m ) ;
                nb_hex_ += nb_hex( m ) ;
            }

        }

        void MacroMeshExport::duplicate_vertices( const DuplicateMode& mode )
        {
            /// 1 - Get all the corner vertices (redondant information)
            std::vector< vec3 > corner_vertices( corners_.size() ) ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& mesh = mm_.mesh( m ) ;
                index_t mesh_start = mesh_corner_ptr_[m] ;
                for( index_t c = 0; c < mesh.nb_corners(); c++ ) {
                    corner_vertices[mesh_start + c] = GEO::Geom::mesh_corner_vertex(
                        mesh, c ) ;
                }
            }

            /// 2 - Tag all corners to duplicate
            std::vector< bool > corner_to_duplicate( corner_vertices.size(), false ) ;
            {
                ColocaterANN ann( corner_vertices ) ;
                const BoundaryModel& model = mm_.model() ;
                for( index_t s = 0; s < model.nb_surfaces(); s++ ) {
                    if( !is_surface_to_duplicate( s, mode ) ) continue ;
                    const Surface& surface = model.surface( s ) ;
                    for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
                        std::vector< index_t > colocated_corners ;
                        ann.get_colocated( surface.vertex( v ), colocated_corners ) ;
                        for( index_t co = 0; co < colocated_corners.size(); co++ ) {
                            corner_to_duplicate[colocated_corners[co]] = true ;
                        }
                    }
                }
            }
            corner_vertices.clear() ;

            /// 3 - Duplicate the corners
//            typedef std::pair< index_t, bool > surface_side ;
            std::vector< bool > is_vertex_to_duplicate( mm_.nb_vertices(), false ) ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& mesh = mm_.mesh( m ) ;
                ColocaterANN ann( mesh, ColocaterANN::FACETS ) ;
                for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                    for( index_t co = mesh.cell_vertices_begin( c );
                        co < mesh.cell_vertices_begin( c + 1 ); co++ ) {
                        if( !corner_to_duplicate[co] ) continue ;

                        index_t vertex_id = mesh.corner_vertex_index( co ) ;
                        std::vector< index_t > corner_used ;
                        std::vector< index_t > cell_added ;
//                        std::set< surface_side > surfaces ;
                        std::stack< index_t > S ;
                        S.push( c ) ;
                        cell_added.push_back( c ) ;
                        do {
                            index_t cur_c = S.top() ;
                            S.pop() ;
                            index_t cur_co = mesh.cell_vertices_begin( cur_c ) ;
                            for( ; cur_co < mesh.cell_vertices_begin( cur_c + 1 );
                                cur_co++ ) {
                                if( mesh.corner_vertex_index( cur_co )
                                    == vertex_id ) {
                                    break ;
                                }
                            }
                            corner_to_duplicate[cur_co] = false ;
                            corner_used.push_back( cur_co ) ;
                            for( index_t cur_f = 0;
                                cur_f < mesh.cell_nb_facets( cur_c ); cur_f++ ) {
                                for( index_t cur_v = 0;
                                    cur_v
                                        < mesh.cell_facet_nb_vertices( cur_c, cur_f );
                                    cur_v++ ) {
                                    if( mesh.cell_facet_vertex_index( cur_c, cur_f,
                                        cur_v ) == vertex_id ) {
                                        std::vector< index_t > result ;
                                        if( !ann.get_colocated(
                                            Utils::mesh_cell_facet_center( mesh,
                                                cur_c, cur_f ), result ) ) {
//                                            index_t surface_id = mesh.facet_region( result[0] ) ;
//                                            vec3 facet_normal =
//                                                GEO::Geom::mesh_facet_normal( mesh,
//                                                    cur_f ) ;
//                                            vec3 cell_facet_normal =
//                                                Utils::mesh_cell_facet_normal( mesh,
//                                                    cur_c, cur_f ) ;
//                                            bool side = dot( facet_normal,
//                                                cell_facet_normal ) > 0 ;
//                                            surfaces.insert(
//                                                surface_side( surface_id, side ) ) ;
                                            signed_index_t cur_adj =
                                                mesh.cell_adjacent( cur_c, cur_f ) ;
                                            if( cur_adj != -1
                                                && !Utils::contains( cell_added,
                                                    index_t( cur_adj ) ) ) {
                                                cell_added.push_back( cur_adj ) ;
                                            }
                                        }
                                        break ;
                                    }
                                }
                            }
                        } while( !S.empty() ) ;

                        index_t global_vertex_id = mm_.global_vertex_id( m, vertex_id ) ;
                        if( is_vertex_to_duplicate[global_vertex_id] ) {
                            index_t duplicated_vertex_id =
                                first_duplicated_vertex_id_
                                    + duplicated_vertex_indices_.size() ;
                            duplicated_vertex_indices_.push_back( global_vertex_id ) ;
                            for( index_t cur_co = 0; cur_co < corner_used.size();
                                cur_co++ ) {
                                corners_[mesh_corner_ptr_[m] + corner_used[cur_co]] =
                                    duplicated_vertex_id ;
                            }
                        } else {
                            is_vertex_to_duplicate[global_vertex_id] = true ;
                        }
                    }
                }
            }

        }

        bool MacroMeshExport::is_surface_to_duplicate(
            index_t s,
            const DuplicateMode& mode ) const
        {
            BoundaryModelElement::GEOL_FEATURE feature = mm_.model().surface( s ).geological_feature() ;
            if( mode == FAULT && feature == BoundaryModelElement::FAULT ) return true ;

            return false ;
        }

    }

}
