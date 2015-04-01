/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr 
 *     Antoine.Mazuyer@univ-lorraine.fr 
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Sup�rieure de G�ologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
 */

#include <ringmesh/io.h>
#include <ringmesh/boundary_model.h>
#include <ringmesh/boundary_model_builder.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/string.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>

#include <third_party/zlib/zip.h>
#include <third_party/zlib/unzip.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stack>

namespace RINGMesh {
    namespace RINGMeshIO {

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
                    fclose( out ) ;
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
                builder.load_ml_file( filename ) ;
                return true ;
            }

            virtual bool save( BoundaryModel& model, const std::string& filename )
            {
                std::ofstream out( filename.c_str() ) ;
                return model.save_gocad_model3d( out ) ;
            }
        } ;

        class BMIOHandler: public BoundaryModelIOHandler {
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

                BoundaryModelBuilderBM builder( model ) ;
                builder.load_file( filename ) ;
                return true ;
            }

            virtual bool save( BoundaryModel& model, const std::string& filename )
            {
                model.save_bm_file( filename ) ;
                return true ;
            }
        } ;

        /************************************************************************/

        BoundaryModelIOHandler* BoundaryModelIOHandler::create(
            const std::string& format )
        {
            ringmesh_register_BoundaryModelIOHandler_creator( MLIOHandler, "ml" ) ;
            ringmesh_register_BoundaryModelIOHandler_creator( BMIOHandler, "bm" );

            BoundaryModelIOHandler* handler =
                BoundaryModelIOHandlerFactory::create_object( format ) ;
            if( handler ) {
                return handler ;
            }

            GEO::Logger::err( "I/O" ) << "Unsupported file format: " << format
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
        class AsterIOHandler: public MacroMeshIOHandler {
        public:
            virtual bool load( const std::string& filename, MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                    << "Loading of a MacroMesh from Code_Aster mesh not implemented yet"
                    << std::endl ;
                return false ;
            }
            virtual bool save( const MacroMesh& mm, const std::string& filename )
            {
                const BoundaryModel& model = mm.model() ;
                std::vector< index_t > vertex_exported_id( mm.vertices.nb_vertices(),
                    NO_ID ) ;
                std::vector< index_t > atom_exported_id( mm.vertices.nb_duplicated_vertices(),
                    NO_ID ) ;
                std::ofstream out( filename.c_str() ) ;
                out.precision( 16 ) ;
                std::vector< bool > vertex_exported( mm.vertices.nb_vertices(), false ) ;
                std::vector< bool > atom_exported( mm.vertices.nb_duplicated_vertices(),
                    false ) ;


                std::cout << "nb vertices " << mm.vertices.nb_vertices() << std::endl ;
                std::cout << "nb total vertices " << mm.vertices.nb_total_vertices()
                    << std::endl ;

                index_t nb_vertices_exported = 0 ;
                index_t cur_cell = 0 ;
                index_t cur_facet = 0 ;

                /// 1. Write the vertices coordinates (with the duplicate ones)
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const GEO::Mesh& mesh = mm.mesh( r ) ;
                    out << "COOR_3D" << std::endl ;
                    for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                        for( index_t co = mesh.cells.corners_begin( c );
                            co < mesh.cells.corners_end( c ); co++ ) {
                            index_t vertex_id, atom_id ;
                            if( mm.vertices.vertex_id( r, co, vertex_id, atom_id ) ) {
                                if( vertex_exported[vertex_id] ) continue ;
                                vertex_exported[vertex_id] = true ;
                                vertex_exported_id[vertex_id] =
                                    nb_vertices_exported ;
                                out << "V" << nb_vertices_exported++ << " "
                                    << mm.vertices.vertex( vertex_id ) << std::endl ;
                            }
                        }
                    }
                    out << "FINSF" << std::endl ;

                    out << "TETRA4" << std::endl ;

                    for( index_t c = 0; c < mm.cells.nb_tet( r ); c++ ) {
                        index_t cur_tet = mm.cells.tet_id( r, c ) ;
                        out << "C" << cur_cell++ << " " ;

                        for( index_t co = mesh.cells.corners_begin( cur_tet );
                            co < mesh.cells.corners_end( cur_tet ); co++ ) {
                            index_t vertex_id, atom_id ;
                            if( mm.vertices.vertex_id( r, co, vertex_id, atom_id ) ) {
                                out << "V" << vertex_exported_id[vertex_id] << " " ;
                            } else {
                                out << "V" << atom_exported_id[atom_id] << " " ;
                            }
                        }
                        out << std::endl ;

                    }
                    out << "FINSF" << std::endl ;
                }

                cur_cell = 0 ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    out << "GROUP_MA" << std::endl ;
                    out << model.region(r).name() ;
                    for(index_t c = 0 ; c < mm.mesh(r).cells.nb() ; c++) {
                        out << "C" << cur_cell++ << std::endl ;
                    }
                }
                out << "FINSF" << std::endl ;

                out << "TRIA3" << std::endl ;
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    const RINGMesh::BoundaryModelElement& interf =
                        model.one_interface( i ) ;
                    for( index_t s = 0; s < interf.nb_children(); s++ ) {
                        index_t surface_id = interf.child_id( s ) ;
                        index_t mesh_id = mm.facets.mesh( surface_id ) ;
                        const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                        for( index_t f = 0; f < mm.facets.nb_facets( surface_id );
                            f++ ) {
                            index_t facet_id = mm.facets.facet( surface_id, f ) ;
                            out << "F" << facet_id ;
                            for( index_t v = mesh.facets.corners_begin( facet_id );
                                v < mesh.facets.corners_end( facet_id ); v++ ) {
                                out << " V"
                                    << vertex_exported_id[mm.vertices.vertex_id(
                                        mesh_id, mesh.facet_corners.vertex( v ) )] ;
                            }
                            out << std::endl ;
                        }
                    }
                }
                out << "FINSF" << std::endl ;

                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    const RINGMesh::BoundaryModelElement& interf =
                        model.one_interface( i ) ;
                    for( index_t s = 0; s < interf.nb_children(); s++ ) {
                        index_t surface_id = interf.child_id( s ) ;
                        out << "GROUP_MA" <<  std::endl ;
                        out << "S" << surface_id << std::endl ;

                        index_t mesh_id = mm.facets.mesh( surface_id ) ;
                        const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                        for( index_t f = 0; f < mm.facets.nb_facets( surface_id );
                            f++ ) {
                            index_t facet_id = mm.facets.facet( surface_id, f ) ;
                            out << "F" << facet_id ;
                            out << std::endl ;
                        }
                        out << "FINSF" << std::endl ;
                    }
                }



                out << "FIN" << std::endl ;

                out.close() ;

                return true ;

            }
        } ;
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
                        ringmesh_assert_not_reached;
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
                mesh.vertices.create_vertices( mm.vertices.nb_vertices() ) ;
                for( index_t p = 0; p < mm.vertices.nb_vertices(); p++ ) {
                    mesh.vertices.point( p ) = mm.vertices.vertex( p ) ;
                }

                GEO::Attribute< index_t > attribute( mesh.facets.attributes(),
                    surface_att_name ) ;
                index_t cell_offset = 0 ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& cur_mesh = mm.mesh( m ) ;
                    GEO::Attribute< index_t > cur_attribute(
                        cur_mesh.facets.attributes(), surface_att_name ) ;
                    for( index_t f = 0; f < cur_mesh.facets.nb(); f++ ) {
                        GEO::vector< index_t > vertices ;
                        vertices.reserve( cur_mesh.facets.nb_vertices( f ) ) ;
                        for( index_t v = cur_mesh.facets.corners_begin( f );
                            v < cur_mesh.facets.corners_end( f ); v++ ) {
                            vertices.push_back(
                                mm.vertices.vertex_id( m,
                                    cur_mesh.facet_corners.vertex( v ) ) ) ;
                        }
                        index_t f_id = mesh.facets.create_polygon( vertices ) ;
                        attribute[f_id] = cur_attribute[f] ;
                    }

                    for( index_t c = 0; c < cur_mesh.cells.nb(); c++ ) {
                        std::vector< index_t > vertex_indices ;
                        vertex_indices.reserve( 8 ) ;
                        for( index_t v = 0; v < cur_mesh.cells.nb_vertices( c );
                            v++ ) {
                            vertex_indices.push_back(
                                mm.vertices.vertex_id( m,
                                    cur_mesh.cells.vertex( c, v ) ) ) ;
                        }
                        std::vector< index_t > adj_indices ;
                        adj_indices.reserve( 6 ) ;
                        for( index_t f = 0; f < cur_mesh.cells.nb_facets( c );
                            f++ ) {
                            index_t adj = cur_mesh.cells.adjacent( c, f ) ;
                            adj_indices.push_back(
                                adj == GEO::NO_CELL ? adj : cell_offset + adj ) ;
                        }
                        if( cur_mesh.cells.type( c ) == GEO::MESH_TET ) {
                            mesh.cells.create_tet( vertex_indices[0],
                                vertex_indices[1], vertex_indices[2],
                                vertex_indices[3], adj_indices[0], adj_indices[1],
                                adj_indices[2], adj_indices[3] ) ;
                        } else if( cur_mesh.cells.type( c ) == GEO::MESH_PRISM ) {
                            mesh.cells.create_prism( vertex_indices[0],
                                vertex_indices[1], vertex_indices[2],
                                vertex_indices[3], vertex_indices[4],
                                vertex_indices[5], adj_indices[0], adj_indices[1],
                                adj_indices[2], adj_indices[3], adj_indices[4] ) ;
                        } else if( cur_mesh.cells.type( c ) == GEO::MESH_PYRAMID ) {
                            mesh.cells.create_pyramid( vertex_indices[0],
                                vertex_indices[1], vertex_indices[2],
                                vertex_indices[3], vertex_indices[4], adj_indices[0],
                                adj_indices[1], adj_indices[2], adj_indices[3],
                                adj_indices[4] ) ;
                        } else if( cur_mesh.cells.type( c ) == GEO::MESH_HEX ) {
                            mesh.cells.create_hex( vertex_indices[0],
                                vertex_indices[1], vertex_indices[2],
                                vertex_indices[3], vertex_indices[4],
                                vertex_indices[5], vertex_indices[6],
                                vertex_indices[7], adj_indices[0], adj_indices[1],
                                adj_indices[2], adj_indices[3], adj_indices[4],
                                adj_indices[5] ) ;
                        } else {
                            ringmesh_assert_not_reached;
                        }
                    }
                    cell_offset += cur_mesh.cells.nb() ;
                }

                mesh.facets.connect() ;
                mesh.cells.connect() ;

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
                std::string directory = GEO::FileSystem::dir_name( filename ) ;
                std::string file = GEO::FileSystem::base_name( filename ) ;

                std::ostringstream oss_node ;
                oss_node << directory << "/" << file << ".node" ;
                std::ofstream node( oss_node.str().c_str() ) ;
                node.precision( 16 ) ;
                node << mm.vertices.nb_vertices() << " 3 0 0" << std::endl ;
                for( index_t v = 0; v < mm.vertices.nb_vertices(); v++ ) {
                    const vec3& point = mm.vertices.vertex( v ) ;
                    node << v << SPACE << point.x << SPACE << point.y << SPACE
                        << point.z << std::endl ;
                }

                std::ostringstream oss_ele ;
                oss_ele << directory << "/" << file << ".ele" ;
                std::ofstream ele( oss_ele.str().c_str() ) ;
                std::ostringstream oss_neigh ;
                oss_neigh << directory << "/" << file << ".neigh" ;
                std::ofstream neigh( oss_neigh.str().c_str() ) ;

                ele << mm.cells.nb_cells() << " 4 1" << std::endl ;
                neigh << mm.cells.nb_cells() << " 4" << std::endl ;
                index_t nb_tet_exported = 0 ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& mesh = mm.mesh( m ) ;
                    for( index_t tet = 0; tet < mesh.cells.nb(); tet++ ) {
                        ele << nb_tet_exported + tet << SPACE
                            << mm.vertices.vertex_id( m, mesh.cells.vertex( tet, 0 ) )
                            << SPACE
                            << mm.vertices.vertex_id( m, mesh.cells.vertex( tet, 1 ) )
                            << SPACE
                            << mm.vertices.vertex_id( m, mesh.cells.vertex( tet, 2 ) )
                            << SPACE
                            << mm.vertices.vertex_id( m, mesh.cells.vertex( tet, 3 ) )
                            << SPACE << m + 1 << std::endl ;
                        neigh << nb_tet_exported + tet ;
                        for( index_t f = 0; f < mesh.cells.nb_facets( tet ); f++ ) {
                            neigh << SPACE << mm.cells.cell_adjacent( m, tet, f ) ;
                        }
                        neigh << std::endl ;
                    }
                    nb_tet_exported += mesh.cells.nb() ;
                }
                return true ;
            }
        } ;

        /************************************************************************/

        class VTKIOHandler: public MacroMeshIOHandler {
        public:
            virtual bool load( const std::string& filename, MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                    << "Loading of a MacroMesh from VTK not implemented yet"
                    << std::endl ;
                return false ;
            }
            virtual bool save( const MacroMesh& mm, const std::string& filename )
            {
                std::ofstream out( filename.c_str() ) ;
                out.precision( 16 ) ;

                out << "# vtk DataFile Version 2.0" << std::endl ;
                out << "Unstructured Grid" << std::endl ;
                out << "ASCII" << std::endl ;
                out << "DATASET UNSTRUCTURED_GRID" << std::endl ;

                out << "POINTS " << mm.vertices.nb_vertices() << " double" << std::endl ;
                for( index_t v = 0; v < mm.vertices.nb_vertices(); v++ ) {
                    const vec3& point = mm.vertices.vertex( v ) ;
                    out << point.x << SPACE << point.y << SPACE << point.z << std::endl ;
                }
                out << std::endl ;

                index_t total_corners = ( 4 + 1 ) * mm.cells.nb_tet()
                    + ( 5 + 1 ) * mm.cells.nb_pyramid() + ( 6 + 1 ) * mm.cells.nb_prism()
                    + ( 8 + 1 ) * mm.cells.nb_hex() ;
                out << "CELLS " << mm.cells.nb_cells() << SPACE << total_corners
                    << std::endl ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& mesh = mm.mesh( m ) ;
                    for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                        out << mesh.cells.nb_vertices( c ) ;
                        for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                            out << SPACE << mm.vertices.vertex_id( m, mesh.cells.vertex( c, v ) ) ;
                        }
                        out << std::endl ;
                    }
                }

                out << "CELL_TYPES " << mm.cells.nb_cells() << std::endl ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& mesh = mm.mesh( m ) ;
                    for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                        out << cell_type( mesh.cells.type( c ) ) << std::endl ;
                    }
                }
                out << std::endl ;

                out << "CELL_DATA " << mm.cells.nb_cells() << std::endl ;
                out << "SCALARS region int 1" << std::endl ;
                out << "LOOKUP_TABLE default" << std::endl ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& mesh = mm.mesh( m ) ;
                    for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                        out << m << std::endl ;
                    }
                }
                out << std::endl ;
                return true ;
            }

        private:
            index_t cell_type( GEO::MeshCellType t ) const {
                switch( t ) {
                    case GEO::MESH_TET:
                        return 10 ;
                    case GEO::MESH_PYRAMID:
                        return 14 ;
                    case GEO::MESH_PRISM:
                        return 13 ;
                    case GEO::MESH_HEX:
                        return 12 ;
                    default:
                        ringmesh_assert_not_reached ;
                        return NO_ID ;
                }
            }
        } ;

        /************************************************************************/

        class TSolidIOHandler: public MacroMeshIOHandler {
        public:
            virtual bool load( const std::string& filename, MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                    << "Loading of a MacroMesh from TSolid not implemented yet"
                    << std::endl ;
                return false ;
            }
            virtual bool save( const MacroMesh& mm, const std::string& filename )
            {
                std::ofstream out( filename.c_str() ) ;
                out.precision( 16 ) ;

                const BoundaryModel& model = mm.model() ;
                // Print Model3d headers
                out << "GOCAD TSolid 1" << std::endl << "HEADER {" << std::endl
                    << "name:" << model.name() << std::endl << "}" << std::endl ;

                out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl
                    << "NAME Default" << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\""
                    << std::endl << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl
                    << "ZPOSITIVE Elevation" << std::endl
                    << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl ;

                mm.set_duplicate_mode( ALL ) ;

                std::vector< bool > vertex_exported( mm.vertices.nb_vertices(), false ) ;
                std::vector< bool > atom_exported( mm.vertices.nb_duplicated_vertices(),
                    false ) ;
                std::vector< index_t > vertex_exported_id( mm.vertices.nb_vertices(),
                    NO_ID ) ;
                std::vector< index_t > atom_exported_id( mm.vertices.nb_duplicated_vertices(),
                    NO_ID ) ;
                index_t nb_vertices_exported = 1 ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const RINGMesh::BoundaryModelElement& region = model.region(
                        r ) ;
                    out << "TVOLUME " << region.name() << std::endl ;

                    const GEO::Mesh& mesh = mm.mesh( r ) ;
                    for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                        for( index_t co = mesh.cells.corners_begin( c );
                            co < mesh.cells.corners_end( c ); co++ ) {
                            index_t vertex_id, atom_id ;
                            if( mm.vertices.vertex_id( r, co, vertex_id, atom_id ) ) {
                                if( vertex_exported[vertex_id] ) continue ;
                                vertex_exported[vertex_id] = true ;
                                vertex_exported_id[vertex_id] =
                                    nb_vertices_exported ;
                                out << "VRTX " << nb_vertices_exported++ << " "
                                    << mm.vertices.vertex( vertex_id ) << std::endl ;
                            }
                        }
                    }

                    for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                        for( index_t co = mesh.cells.corners_begin( c );
                            co < mesh.cells.corners_end( c ); co++ ) {
                            index_t vertex_id, atom_id ;
                            if( !mm.vertices.vertex_id( r, co, vertex_id, atom_id ) ) {
                                if( atom_exported[atom_id] ) continue ;
                                atom_exported[atom_id] = true ;
                                atom_exported_id[atom_id] = nb_vertices_exported ;
                                out << "ATOM " << nb_vertices_exported++ << " "
                                    << vertex_exported_id[vertex_id] << std::endl ;
                            }
                        }
                    }

                    std::map< index_t, index_t > sides ;
                    for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                        if( sides.count( region.boundary_id( s ) ) > 0 )
                            sides[region.boundary_id( s )] = 2 ;
                        else
                            sides[region.boundary_id( s )] = region.side( s ) ;
                    }

                    GEO::Attribute< index_t > attribute( mesh.facets.attributes(),
                        surface_att_name ) ;
                    ColocaterANN ann( mesh, ColocaterANN::FACETS ) ;
                    for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                        out << "TETRA" ;
                        for( index_t co = mesh.cells.corners_begin( c );
                            co < mesh.cells.corners_end( c ); co++ ) {
                            index_t vertex_id, atom_id ;
                            if( mm.vertices.vertex_id( r, co, vertex_id, atom_id ) )
                                out << " " << vertex_exported_id[vertex_id] ;
                            else
                                out << " " << atom_exported_id[atom_id] ;
                        }
                        out << std::endl ;
                        out << "# CTETRA " << region.name() ;
                        for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                            out << " " ;
                            vec3 facet_center = Utils::mesh_cell_facet_center( mesh,
                                c, f ) ;
                            std::vector< index_t > result ;
                            if( ann.get_colocated( facet_center, result ) ) {
                                ringmesh_debug_assert( result.size() == 1 ) ;
                                index_t surface_id = attribute[result[0]] ;
                                index_t side = sides[surface_id] ;
                                if( side < 2 )
                                    side ? out << "+" : out << "-" ;
                                else {
                                    vec3 cell_facet_normal =
                                        Utils::mesh_cell_facet_normal( mesh, c, f ) ;
                                    vec3 facet_normal = GEO::Geom::mesh_facet_normal(
                                        mesh, result[0] ) ;
                                    double d = dot( cell_facet_normal,
                                        facet_normal ) ;
                                    d > 0 ? out << "+" : out << "-" ;
                                }
                                out << model.surface( surface_id ).parent().name() ;
                            } else {
                                out << "none" ;
                            }
                        }
                        out << std::endl ;
                    }
                }

                out << "MODEL" << std::endl ;
                int tface_count = 1 ;
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    const RINGMesh::BoundaryModelElement& interf =
                        model.one_interface( i ) ;
                    out << "SURFACE " << interf.name() << std::endl ;
                    for( index_t s = 0; s < interf.nb_children(); s++ ) {
                        out << "TFACE " << tface_count++ << std::endl ;
                        index_t surface_id = interf.child_id( s ) ;
                        index_t mesh_id = mm.facets.mesh( surface_id ) ;
                        const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                        out << "KEYVERTICES" ;
                        index_t key_facet_id = mm.facets.facet( surface_id, 0 ) ;
                        for( index_t v = mesh.facets.corners_begin( key_facet_id );
                            v < mesh.facets.corners_end( key_facet_id ); v++ ) {
                            out << " "
                                << vertex_exported_id[mm.vertices.vertex_id( mesh_id,
                                    mesh.facet_corners.vertex( v ) )] ;
                        }
                        out << std::endl ;
                        for( index_t f = 0; f < mm.facets.nb_facets( surface_id );
                            f++ ) {
                            index_t facet_id = mm.facets.facet( surface_id, f ) ;
                            out << "TRGL" ;
                            for( index_t v = mesh.facets.corners_begin( facet_id );
                                v < mesh.facets.corners_end( facet_id ); v++ ) {
                                out << " "
                                    << vertex_exported_id[mm.vertices.vertex_id(
                                        mesh_id, mesh.facet_corners.vertex( v ) )] ;
                            }
                            out << std::endl ;
                        }
                    }
                }

                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const RINGMesh::BoundaryModelElement& region = model.region(
                        r ) ;
                    out << "MODEL_REGION " << region.name() << " " ;
                    region.side( 0 ) ? out << "+" : out << "-" ;
                    out << region.boundary_id( 0 ) + 1 << std::endl ;
                }

                out << "END" << std::endl ;
                return true ;
            }
        } ;

        /************************************************************************/
        struct RINGMesh2CSMP {
            index_t element_type ;
            index_t nb_vertices ;
            index_t vertices[8] ;
            index_t nb_facets ;
            index_t nb_vertices_in_facet[6] ;
            index_t facet[6] ;
            index_t vertices_in_facet[6][4] ;
        } ;

        static RINGMesh2CSMP tet_descriptor = { 4,                  // type
            4,                  // nb vertices
            { 0, 1, 2, 3 },     // vertices
            4,                  // nb facets
            { 3, 3, 3, 3 },     // nb vertices in facet
            { 0, 1, 2, 3 },     // facets
            { { 1, 3, 2 }, { 0, 2, 3 }, { 3, 1, 0 }, { 0, 1, 2 } } } ;

        static RINGMesh2CSMP hex_descriptor = { 6,                         // type
            8,                              // nb vertices
            { 0, 1, 3, 2, 4, 5, 7, 6 },     // vertices
            6,                              // nb facets
            { 4, 4, 4, 4, 4, 4 },           // nb vertices in facet
            { 4, 2, 1, 3, 0, 5 },           // facets
            { { 0, 3, 7, 4 }, { 2, 1, 5, 6 }, { 1, 0, 4, 5 }, { 3, 2, 6, 7 }, {
                1, 2, 3, 0 }, { 4, 7, 6, 5 } } } ;

        static RINGMesh2CSMP prism_descriptor = { 12,                     // type
            6,                      // nb vertices
            { 0, 1, 2, 3, 4, 5 },   // vertices
            5,                      // nb facets
            { 3, 4, 4, 4, 3 },      // nb vertices in facet
            { 0, 2, 4, 3, 1 },      // facets
            {
                { 0, 1, 2 }, { 3, 5, 4 }, { 0, 3, 4, 1 }, { 0, 2, 5, 3 }, {
                    1, 4, 5, 2 } } } ;

        static RINGMesh2CSMP pyramid_descriptor = { 18,                 // type
            5,                  // nb vertices
            { 0, 1, 2, 3, 4 },  // vertices
            5,                  // nb facets
            { 3, 3, 3, 3, 4 },  // nb vertices in facet
            { 1, 3, 4, 2, 0 },  // facets
            { { 0, 1, 2, 3 }, { 0, 4, 1 }, { 0, 3, 4 }, { 2, 4, 3 }, { 2, 1, 4 } } } ;

        class CSMPIOHandler: public MacroMeshIOHandler {
        public:
            virtual bool load( const std::string& filename, MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                    << "Loading of a MacroMesh from CSMP not implemented yet"
                    << std::endl ;
                return false ;
            }
            virtual bool save( const MacroMesh& mm, const std::string& filename )
            {
                if( !initialize( mm ) ) {
                    return false ;
                }

                std::string directory = GEO::FileSystem::dir_name( filename ) ;
                std::string file = GEO::FileSystem::base_name( filename ) ;

                std::ostringstream oss_ascii ;
                oss_ascii << directory << "/" << file << ".asc" ;
                std::ofstream ascii( oss_ascii.str().c_str() ) ;
                ascii.precision( 16 ) ;

                const BoundaryModel& model = mm.model() ;
                ascii << model.name() << std::endl ;
                ascii << "Model generated from RINGMesh" << std::endl ;

                std::ostringstream oss_data ;
                oss_data << directory << "/" << file << ".dat" ;
                std::ofstream data( oss_data.str().c_str() ) ;
                data.precision( 16 ) ;

                std::ostringstream oss_regions ;
                oss_regions << directory << "/" << file << "-regions.txt" ;
                std::ofstream regions( oss_regions.str().c_str() ) ;
                regions << "'" << oss_regions.str() <<std::endl ;
                regions << "no properties" << std::endl ;

                index_t count = 0 ;
                // Conversion from (X,Y,Z) to (X,Z,-Y)
                signed_index_t conversion_sign[3] = { 1, 1, -1 } ;
                index_t conversion_axis[3] = { 0, 2, 1 } ;
                data << mm.vertices.nb_vertices() << " # PX, PY, PZ" << std::endl ;
                for( index_t dim = 0; dim < 3; dim++ ) {
                    for( index_t v = 0; v < mm.vertices.nb_vertices(); v++ ) {
                        data << " "
                            << conversion_sign[dim]
                                * mm.vertices.vertex( v )[conversion_axis[dim]] ;
                        count++ ;
                        if( count == 5 ) {
                            count = 0 ;
                            data << std::endl ;
                        }
                    }
                    if( count != 0 ) {
                        count = 0 ;
                        data << std::endl ;
                    }
                }
                if( count != 0 ) {
                    count = 0 ;
                    data << std::endl ;
                }

                index_t nb_families = 0 ;
                std::vector< index_t > nb_triangle_interface( model.nb_interfaces(),
                    0 ) ;
                std::vector< index_t > nb_quad_interface( model.nb_interfaces(),
                    0 ) ;
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    const BoundaryModelElement& interf = model.one_interface( i ) ;
                    for( index_t s = 0; s < interf.nb_children(); s++ ) {
                        index_t s_id = interf.child_id( s ) ;
                        nb_triangle_interface[i] += mm.facets.nb_triangle( s_id ) ;
                        nb_quad_interface[i] += mm.facets.nb_quad( s_id ) ;
                    }
                    if( nb_triangle_interface[i] > 0 ) nb_families++ ;
                    if( nb_quad_interface[i] > 0 ) nb_families++ ;
                }
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    if( mm.cells.nb_tet( r ) > 0 ) nb_families++ ;
                    if( mm.cells.nb_pyramid( r ) > 0 ) nb_families++ ;
                    if( mm.cells.nb_prism( r ) > 0 ) nb_families++ ;
                    if( mm.cells.nb_hex( r ) > 0 ) nb_families++ ;
                }

                ascii << nb_families << " # Number of families" << std::endl ;
                ascii << "# Object name" << TAB << "Element type" << TAB
                    << "Material-ID" << TAB << "Number of elements" << std::endl ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const RINGMesh::BoundaryModelElement& region = model.region(
                        r ) ;
                    regions << region.name() << std::endl ;
                    if( mm.cells.nb_tet( r ) > 0 ) {
                        ascii << region.name() << TAB << "TETRA_4" << TAB << 0 << TAB
                            << mm.cells.nb_tet( r ) << std::endl ;
                    }
                    if( mm.cells.nb_pyramid( r ) > 0 ) {
                        ascii << region.name() << TAB << "PYRA_5" << TAB << 0 << TAB
                            << mm.cells.nb_pyramid( r ) << std::endl ;
                    }
                    if( mm.cells.nb_prism( r ) > 0 ) {
                        ascii << region.name() << TAB << "PENTA_6" << TAB << 0 << TAB
                            << mm.cells.nb_prism( r ) << std::endl ;
                    }
                    if( mm.cells.nb_hex( r ) > 0 ) {
                        ascii << region.name() << TAB << "HEXA_8" << TAB << 0 << TAB
                            << mm.cells.nb_hex( r ) << std::endl ;
                    }
                }
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    regions << interface_name( i, mm ) << std::endl ;
                    if( nb_triangle_interface[i] > 0 ) {
                        ascii << interface_name( i, mm ) << TAB << "TRI_3" << TAB
                            << 0 << TAB << nb_triangle_interface[i] << std::endl ;
                    }
                    if( nb_quad_interface[i] > 0 ) {
                        ascii << interface_name( i, mm ) << TAB << "QUAD_4" << TAB
                            << 0 << TAB << nb_quad_interface[i] << std::endl ;
                    }
                }

                data << "# PBFLAGS" << std::endl ;
                for( index_t p = 0; p < mm.vertices.nb_vertices(); p++ ) {
                    data << " " << std::setw( 3 ) << point_boundary( p ) ;
                    count++ ;
                    if( count == 20 ) {
                        data << std::endl ;
                        count = 0 ;
                    }
                }
                if( count != 0 ) {
                    count = 0 ;
                    data << std::endl ;
                }

                data << "# PBVALS" << std::endl ;
                for( index_t p = 0; p < mm.vertices.nb_vertices(); p++ ) {
                    data << " " << std::setw( 3 ) << 0 ;
                    count++ ;
                    if( count == 20 ) {
                        data << std::endl ;
                        count = 0 ;
                    }
                }
                if( count != 0 ) {
                    count = 0 ;
                    data << std::endl ;
                }

                index_t nb_total_elements = mm.cells.nb_cells() + mm.facets.nb_facets() ;
                data << nb_total_elements << " # PELEMENT" << std::endl ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    for( index_t el = 0; el < mm.cells.nb_tet( r ); el++ ) {
                        data << " " << std::setw( 3 ) << 4 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                    for( index_t el = 0; el < mm.cells.nb_pyramid( r ); el++ ) {
                        data << " " << std::setw( 3 ) << 18 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                    for( index_t el = 0; el < mm.cells.nb_prism( r ); el++ ) {
                        data << " " << std::setw( 3 ) << 12 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                    for( index_t el = 0; el < mm.cells.nb_hex( r ); el++ ) {
                        data << " " << std::setw( 3 ) << 6 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                }
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    for( index_t el = 0; el < nb_triangle_interface[i]; el++ ) {
                        data << " " << std::setw( 3 ) << 8 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                    for( index_t el = 0; el < nb_quad_interface[i]; el++ ) {
                        data << " " << std::setw( 3 ) << 14 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                }
                if( count != 0 ) {
                    count = 0 ;
                    data << std::endl ;
                }

                ascii
                    << "# now the elements which make up each object are listed in sequence"
                    << std::endl ;
                index_t cur_cell = 0 ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const RINGMesh::BoundaryModelElement& region = model.region(
                        r ) ;
                    if( mm.cells.nb_tet( r ) > 0 ) {
                        ascii << region.name() << " " << "TETRA_4" << " "
                            << mm.cells.nb_tet( r ) << std::endl ;
                        for( index_t el = 0; el < mm.cells.nb_tet( r ); el++ ) {
                            ascii << cur_cell++ << " " ;
                            count++ ;
                            if( count == 10 ) {
                                ascii << std::endl ;
                                count = 0 ;
                            }
                        }
                        if( count != 0 ) {
                            count = 0 ;
                            ascii << std::endl ;
                        }
                    }
                    if( mm.cells.nb_pyramid( r ) > 0 ) {
                        ascii << region.name() << " " << "PYRA_5" << " "
                            << mm.cells.nb_pyramid( r ) << std::endl ;
                        for( index_t el = 0; el < mm.cells.nb_pyramid( r ); el++ ) {
                            ascii << cur_cell++ << " " ;
                            count++ ;
                            if( count == 10 ) {
                                ascii << std::endl ;
                                count = 0 ;
                            }
                        }
                        if( count != 0 ) {
                            count = 0 ;
                            ascii << std::endl ;
                        }
                    }
                    if( mm.cells.nb_prism( r ) > 0 ) {
                        ascii << region.name() << " " << "PENTA_6" << " "
                            << mm.cells.nb_prism( r ) << std::endl ;
                        for( index_t el = 0; el < mm.cells.nb_prism( r ); el++ ) {
                            ascii << cur_cell++ << " " ;
                            count++ ;
                            if( count == 10 ) {
                                ascii << std::endl ;
                                count = 0 ;
                            }
                        }
                        if( count != 0 ) {
                            count = 0 ;
                            ascii << std::endl ;
                        }
                    }
                    if( mm.cells.nb_hex( r ) > 0 ) {
                        ascii << region.name() << " " << "HEXA_8" << " "
                            << mm.cells.nb_hex( r ) << std::endl ;
                        for( index_t el = 0; el < mm.cells.nb_hex( r ); el++ ) {
                            ascii << cur_cell++ << " " ;
                            count++ ;
                            if( count == 10 ) {
                                ascii << std::endl ;
                                count = 0 ;
                            }
                        }
                        if( count != 0 ) {
                            count = 0 ;
                            ascii << std::endl ;
                        }
                    }
                }
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    if( nb_triangle_interface[i] > 0 ) {
                        ascii << interface_name( i, mm ) << " " << "TRI_3" << " "
                            << nb_triangle_interface[i] << std::endl ;
                        for( index_t el = 0; el < nb_triangle_interface[i]; el++ ) {
                            ascii << cur_cell++ << " " ;
                            count++ ;
                            if( count == 10 ) {
                                ascii << std::endl ;
                                count = 0 ;
                            }
                        }
                        if( count != 0 ) {
                            count = 0 ;
                            ascii << std::endl ;
                        }
                    }
                    if( nb_quad_interface[i] > 0 ) {
                        ascii << interface_name( i, mm ) << " " << "QUAD_4" << " "
                            << nb_quad_interface[i] << std::endl ;
                        for( index_t el = 0; el < nb_quad_interface[i]; el++ ) {
                            ascii << cur_cell++ << " " ;
                            count++ ;
                            if( count == 10 ) {
                                ascii << std::endl ;
                                count = 0 ;
                            }
                        }
                        if( count != 0 ) {
                            count = 0 ;
                            ascii << std::endl ;
                        }
                    }
                }

                index_t nb_plist = 3 * mm.facets.nb_triangle() + 4 * mm.facets.nb_quad()
                    + 4 * mm.cells.nb_tet() + 5 * mm.cells.nb_pyramid() + 6 * mm.cells.nb_prism()
                    + 8 * mm.cells.nb_hex() ; //+ 2 * nb_edges ;
                data << nb_plist << " # PLIST" << std::endl ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const GEO::Mesh& mesh = mm.mesh( r ) ;
                    for( index_t el = 0; el < mm.cells.nb_tet( r ); el++ ) {
                        index_t tet = mm.cells.tet_id( r, el ) ;
                        for( index_t p = 0; p < 4; p++ ) {
                            index_t vertex_id ;
                            index_t csmp_p = tet_descriptor.vertices[p] ;
                            mm.vertices.vertex_id( r,
                                mesh.cells.corners_begin( tet ) + csmp_p,
                                vertex_id ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                    for( index_t el = 0; el < mm.cells.nb_pyramid( r ); el++ ) {
                        index_t py = mm.cells.pyramid_id( r, el ) ;
                        for( index_t p = 0; p < 5; p++ ) {
                            index_t vertex_id ;
                            index_t csmp_p = pyramid_descriptor.vertices[p] ;
                            mm.vertices.vertex_id( r, mesh.cells.corners_begin( py ) + csmp_p,
                                vertex_id ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                    for( index_t el = 0; el < mm.cells.nb_prism( r ); el++ ) {
                        index_t prism = mm.cells.prism_id( r, el ) ;
                        for( index_t p = 0; p < 6; p++ ) {
                            index_t vertex_id ;
                            index_t csmp_p = prism_descriptor.vertices[p] ;
                            mm.vertices.vertex_id( r,
                                mesh.cells.corners_begin( prism ) + csmp_p,
                                vertex_id ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                    for( index_t el = 0; el < mm.cells.nb_hex( r ); el++ ) {
                        index_t hex = mm.cells.prism_id( r, el ) ;
                        for( index_t p = 0; p < 8; p++ ) {
                            index_t vertex_id ;
                            index_t csmp_p = hex_descriptor.vertices[p] ;
                            mm.vertices.vertex_id( r,
                                mesh.cells.corners_begin( hex ) + csmp_p,
                                vertex_id ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                }
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    const BoundaryModelElement& interf = model.one_interface( i ) ;
                    for( index_t s = 0; s < interf.nb_children(); s++ ) {
                        index_t s_id = interf.child_id( s ) ;
                        index_t mesh_id = mm.facets.mesh( s_id ) ;
                        const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                        for( index_t el = 0; el < mm.facets.nb_triangle( s_id ); el++ ) {
                            index_t tri = mm.facets.triangle_id( mesh_id, el ) ;
                            for( index_t p = mesh.facets.corners_begin( tri );
                                p < mesh.facets.corners_end( tri ); p++ ) {
                                index_t vertex_id = mm.vertices.vertex_id( mesh_id,
                                    mesh.facet_corners.vertex( p ) ) ;
                                data << " " << std::setw( 7 ) << vertex_id ;
                                count++ ;
                                if( count == 10 ) {
                                    data << std::endl ;
                                    count = 0 ;
                                }
                            }
                        }
                        for( index_t el = 0; el < mm.facets.nb_quad( s_id ); el++ ) {
                            index_t quad = mm.facets.quad_id( mesh_id, el ) ;
                            for( index_t p = mesh.facets.corners_begin( quad );
                                p < mesh.facets.corners_end( quad ); p++ ) {
                                index_t vertex_id = mm.vertices.vertex_id( mesh_id,
                                    mesh.facet_corners.vertex( p ) ) ;
                                data << " " << std::setw( 7 ) << vertex_id ;
                                count++ ;
                                if( count == 10 ) {
                                    data << std::endl ;
                                    count = 0 ;
                                }
                            }
                        }
                    }
                }
                if( count != 0 ) {
                    count = 0 ;
                    data << std::endl ;
                }

                int nb_facets = 3 * mm.facets.nb_triangle() + 4 * mm.facets.nb_quad()
                    + 4 * mm.cells.nb_tet() + 5 * mm.cells.nb_pyramid() + 5 * mm.cells.nb_prism()
                    + 6 * mm.cells.nb_hex() ; //+ 2 * nb_edges ;
                data << nb_facets << " # PFVERTS" << std::endl ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    for( index_t el = 0; el < mm.cells.nb_tet( r ); el++ ) {
                        index_t tet = mm.cells.tet_id( r, el ) ;
                        for( index_t f = 0; f < 4; f++ ) {
                            index_t csmp_f = tet_descriptor.facet[f] ;
                            signed_index_t adj = mm.cells.cell_adjacent( r, tet,
                                csmp_f ) ;
                            if( adj == -1 ) {
                                data << " " << std::setw( 7 ) << -28 ;
                            } else {
                                data << " " << std::setw( 7 ) << adj ;
                            }
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                    for( index_t el = 0; el < mm.cells.nb_pyramid( r ); el++ ) {
                        index_t py = mm.cells.pyramid_id( r, el ) ;
                        for( index_t f = 0; f < 5; f++ ) {
                            index_t csmp_f = pyramid_descriptor.facet[f] ;
                            signed_index_t adj = mm.cells.cell_adjacent( r, py,
                                csmp_f ) ;
                            if( adj == -1 ) {
                                data << " " << std::setw( 7 ) << -28 ;
                            } else {
                                data << " " << std::setw( 7 ) << adj ;
                            }
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                    for( index_t el = 0; el < mm.cells.nb_prism( r ); el++ ) {
                        index_t prism = mm.cells.prism_id( r, el ) ;
                        for( index_t f = 0; f < 5; f++ ) {
                            index_t csmp_f = prism_descriptor.facet[f] ;
                            signed_index_t adj = mm.cells.cell_adjacent( r, prism,
                                csmp_f ) ;
                            if( adj == -1 ) {
                                data << " " << std::setw( 7 ) << -28 ;
                            } else {
                                data << " " << std::setw( 7 ) << adj ;
                            }
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                    for( index_t el = 0; el < mm.cells.nb_hex( r ); el++ ) {
                        index_t hex = mm.cells.hex_id( r, el ) ;
                        for( index_t f = 0; f < 6; f++ ) {
                            index_t csmp_f = hex_descriptor.facet[f] ;
                            signed_index_t adj = mm.cells.cell_adjacent( r, hex,
                                csmp_f ) ;
                            if( adj == -1 ) {
                                data << " " << std::setw( 7 ) << -28 ;
                            } else {
                                data << " " << std::setw( 7 ) << adj ;
                            }
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                }
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    const BoundaryModelElement& interf = model.one_interface( i ) ;
                    for( index_t s = 0; s < interf.nb_children(); s++ ) {
                        index_t s_id = interf.child_id( s ) ;
                        index_t mesh_id = mm.facets.mesh( s_id ) ;
                        const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                        for( index_t el = 0; el < mm.facets.nb_triangle( s_id ); el++ ) {
                            index_t tri = mm.facets.triangle_id( mesh_id, el ) ;
                            for( index_t f = mesh.facets.corners_begin( tri );
                                f < mesh.facets.corners_end( tri ); f++ ) {
                                index_t adj = mesh.facet_corners.adjacent_facet(
                                    f ) ;
                                if( adj == GEO::NO_FACET ) {
                                    data << " " << std::setw( 7 ) << -28 ;
                                } else {
                                    data << " " << std::setw( 7 ) << adj ;
                                }
                                count++ ;
                                if( count == 10 ) {
                                    data << std::endl ;
                                    count = 0 ;
                                }
                            }
                        }
                        for( index_t el = 0; el < mm.facets.nb_quad( s_id ); el++ ) {
                            index_t quad = mm.facets.quad_id( mesh_id, el ) ;
                            for( index_t f = mesh.facets.corners_begin( quad );
                                f < mesh.facets.corners_end( quad ); f++ ) {
                                index_t adj = mesh.facet_corners.adjacent_facet(
                                    f ) ;
                                if( adj == GEO::NO_FACET ) {
                                    data << " " << std::setw( 7 ) << -28 ;
                                } else {
                                    data << " " << std::setw( 7 ) << adj ;
                                }
                                count++ ;
                                if( count == 10 ) {
                                    data << std::endl ;
                                    count = 0 ;
                                }
                            }
                        }
                    }
                }
                if( count != 0 ) {
                    data << std::endl ;
                    count = 0 ;
                }

                data << nb_total_elements << " # PMATERIAL" << std::endl ;
                for( unsigned int i = 0; i < nb_total_elements; i++ ) {
                    data << " " << std::setw( 3 ) << 0 ;
                    count++ ;
                    if( count == 20 ) {
                        data << std::endl ;
                        count = 0 ;
                    }
                }

                return true ;
            }

        private:

            bool initialize( const MacroMesh& mm )
            {
                const BoundaryModel& model = mm.model() ;
                std::string cmsp_filename = GEO::CmdLine::get_arg( "out:csmp" ) ;
                box_model_ = cmsp_filename != "" ;
                if( box_model_ ) {
                    GEO::LineInput parser( cmsp_filename ) ;
                    if( !parser.OK() ) {
                        GEO::Logger::err( "I/O" ) << "Cannot open file: "
                            << cmsp_filename << std::endl ;
                        return false ;
                    }
                    parser.get_line() ;
                    parser.get_fields() ;
                    while( !parser.eof() ) {
                        if( parser.nb_fields() == 0 ) continue ;
                        if( parser.nb_fields() != 3 ) return false ;
                        std::string type = parser.field( 1 ) ;
                        index_t interface_id = NO_ID ;
                        if( type == "NAME" ) {
                            std::string name = parser.field( 2 ) ;
                            for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                                if( model.one_interface( i ).name() == name ) {
                                    interface_id = i ;
                                    break ;
                                }
                            }
                        } else if( type == "ID" ) {
                            interface_id = parser.field_as_uint( 2 ) ;
                        } else {
                            GEO::Logger::err( "I/O" ) << "Unknown type: " << type
                                << std::endl ;
                            return false ;
                        }

                        std::string keyword = parser.field( 0 ) ;
                        if( keyword == "BACK" ) {
                            back_ = interface_id ;
                        } else if( keyword == "TOP" ) {
                            top_ = interface_id ;
                        } else if( keyword == "FRONT" ) {
                            front_ = interface_id ;
                        } else if( keyword == "BOTTOM" ) {
                            bottom_ = interface_id ;
                        } else if( keyword == "LEFT" ) {
                            left_ = interface_id ;
                        } else if( keyword == "RIGHT" ) {
                            right_ = interface_id ;
                        } else {
                            GEO::Logger::err( "I/O" ) << "Unknown keyword: "
                                << keyword << std::endl ;
                            return false ;
                        }
                        parser.get_line() ;
                        parser.get_fields() ;
                    }

                    if( back_ == NO_ID || top_ == NO_ID || front_ == NO_ID
                        || bottom_ == NO_ID || left_ == NO_ID || right_ == NO_ID ) {
                        GEO::Logger::err( "I/O" ) << "Missing box shape information"
                            << std::endl ;
                        return false ;
                    }

                    surface_boundary_flags_[back_] = -7 ;
                    surface_boundary_flags_[top_] = -5 ;
                    surface_boundary_flags_[front_] = -6 ;
                    surface_boundary_flags_[bottom_] = -4 ;
                    surface_boundary_flags_[left_] = -2 ;
                    surface_boundary_flags_[right_] = -3 ;

                    std::set< index_t > back_bottom ;
                    back_bottom.insert( back_ ) ;
                    back_bottom.insert( bottom_ ) ;
                    edge_boundary_flags_[back_bottom] = -16 ;
                    std::set< index_t > back_right ;
                    back_right.insert( back_ ) ;
                    back_right.insert( right_ ) ;
                    edge_boundary_flags_[back_right] = -17 ;
                    std::set< index_t > back_top ;
                    back_top.insert( back_ ) ;
                    back_top.insert( top_ ) ;
                    edge_boundary_flags_[back_top] = -18 ;
                    std::set< index_t > back_left ;
                    back_left.insert( back_ ) ;
                    back_left.insert( left_ ) ;
                    edge_boundary_flags_[back_left] = -19 ;
                    std::set< index_t > right_bottom ;
                    right_bottom.insert( right_ ) ;
                    right_bottom.insert( bottom_ ) ;
                    edge_boundary_flags_[right_bottom] = -20 ;
                    std::set< index_t > right_top ;
                    right_top.insert( right_ ) ;
                    right_top.insert( top_ ) ;
                    edge_boundary_flags_[right_top] = -21 ;
                    std::set< index_t > left_top ;
                    left_top.insert( left_ ) ;
                    left_top.insert( top_ ) ;
                    edge_boundary_flags_[left_top] = -22 ;
                    std::set< index_t > left_bottom ;
                    left_bottom.insert( left_ ) ;
                    left_bottom.insert( bottom_ ) ;
                    edge_boundary_flags_[left_bottom] = -23 ;
                    std::set< index_t > front_bottom ;
                    front_bottom.insert( front_ ) ;
                    front_bottom.insert( bottom_ ) ;
                    edge_boundary_flags_[front_bottom] = -24 ;
                    std::set< index_t > front_right ;
                    front_right.insert( front_ ) ;
                    front_right.insert( right_ ) ;
                    edge_boundary_flags_[front_right] = -25 ;
                    std::set< index_t > front_top ;
                    front_top.insert( front_ ) ;
                    front_top.insert( top_ ) ;
                    edge_boundary_flags_[front_top] = -26 ;
                    std::set< index_t > front_left ;
                    front_left.insert( front_ ) ;
                    front_left.insert( left_ ) ;
                    edge_boundary_flags_[front_left] = -27 ;

                    std::set< index_t > back_top_left ;
                    back_top_left.insert( back_ ) ;
                    back_top_left.insert( top_ ) ;
                    back_top_left.insert( left_ ) ;
                    corner_boundary_flags_[back_top_left] = -13 ;
                    std::set< index_t > back_top_right ;
                    back_top_right.insert( back_ ) ;
                    back_top_right.insert( top_ ) ;
                    back_top_right.insert( right_ ) ;
                    corner_boundary_flags_[back_top_right] = -14 ;
                    std::set< index_t > back_bottom_left ;
                    back_bottom_left.insert( back_ ) ;
                    back_bottom_left.insert( bottom_ ) ;
                    back_bottom_left.insert( left_ ) ;
                    corner_boundary_flags_[back_bottom_left] = -8 ;
                    std::set< index_t > back_bottom_right ;
                    back_bottom_right.insert( back_ ) ;
                    back_bottom_right.insert( bottom_ ) ;
                    back_bottom_right.insert( right_ ) ;
                    corner_boundary_flags_[back_bottom_right] = -10 ;
                    std::set< index_t > front_top_left ;
                    front_top_left.insert( front_ ) ;
                    front_top_left.insert( top_ ) ;
                    front_top_left.insert( left_ ) ;
                    corner_boundary_flags_[front_top_left] = -15 ;
                    std::set< index_t > front_top_right ;
                    front_top_right.insert( front_ ) ;
                    front_top_right.insert( top_ ) ;
                    front_top_right.insert( right_ ) ;
                    corner_boundary_flags_[front_top_right] = -9 ;
                    std::set< index_t > front_bottom_left ;
                    front_bottom_left.insert( front_ ) ;
                    front_bottom_left.insert( bottom_ ) ;
                    front_bottom_left.insert( left_ ) ;
                    corner_boundary_flags_[front_bottom_left] = -12 ;
                    std::set< index_t > front_bottom_right ;
                    front_bottom_right.insert( front_ ) ;
                    front_bottom_right.insert( bottom_ ) ;
                    front_bottom_right.insert( right_ ) ;
                    corner_boundary_flags_[front_bottom_right] = -11 ;
                }

                point_boundaries_.resize( mm.vertices.nb_vertices() ) ;
                for( index_t s = 0; s < model.nb_surfaces(); s++ ) {
                    index_t interface_id = model.surface( s ).parent_id() ;
                    index_t mesh_id = mm.facets.mesh( s ) ;
                    const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                    for( index_t f = 0; f < mm.facets.nb_facets( s ); f++ ) {
                        index_t f_id = mm.facets.facet( s, f ) ;
                        for( index_t c = mesh.facets.corners_begin( f_id );
                            c < mesh.facets.corners_end( f_id ); c++ ) {
                            index_t vertex_id = mm.vertices.vertex_id( mesh_id,
                                mesh.facet_corners.vertex( c ) ) ;
                            point_boundaries_[vertex_id].insert( interface_id ) ;
                        }
                    }
                }

                return true ;
            }
            std::string interface_name( index_t i, const MacroMesh& mm )
            {
                if( box_model_ ) {
                    if( i == back_ ) {
                        return "BACK" ;
                    } else if( i == top_ ) {
                        return "TOP" ;
                    } else if( i == front_ ) {
                        return "FRONT" ;
                    } else if( i == bottom_ ) {
                        return "BOTTOM" ;
                    } else if( i == left_ ) {
                        return "LEFT" ;
                    } else if( i == right_ ) {
                        return "RIGHT" ;
                    }
                }
                return mm.model().one_interface( i ).name() ;
            }
            signed_index_t point_boundary( index_t p )
            {
                ringmesh_debug_assert( p < point_boundaries_.size() ) ;
                const std::set< unsigned int >& boundaries = point_boundaries_[p] ;
                if( box_model_ ) {
                    if( boundaries.size() == 1 ) {
                        std::map< unsigned int, int >::iterator it =
                            surface_boundary_flags_.find( *boundaries.begin() ) ;
                        ringmesh_debug_assert( it != surface_boundary_flags_.end() ) ;
                        return it->second ;
                    } else if( boundaries.size() == 2 ) {
                        std::map< std::set< unsigned int >, int >::iterator it =
                            edge_boundary_flags_.find( boundaries ) ;
                        ringmesh_debug_assert( it != edge_boundary_flags_.end() ) ;
                        return it->second ;
                    } else if( boundaries.size() == 3 ) {
                        std::map< std::set< unsigned int >, int >::iterator it =
                            corner_boundary_flags_.find( boundaries ) ;
                        ringmesh_debug_assert( it != corner_boundary_flags_.end() ) ;
                        return it->second ;
                    } else {
                        return 0 ;
                    }
                } else {
                    if( boundaries.empty() )
                        return 0 ;
                    else
                        return -28 ;
                }
            }
        private:
            std::vector< std::set< index_t > > point_boundaries_ ;

            bool box_model_ ;
            index_t back_ ;
            index_t top_ ;
            index_t front_ ;
            index_t bottom_ ;
            index_t left_ ;
            index_t right_ ;

            std::map< std::set< index_t >, signed_index_t > corner_boundary_flags_ ;
            std::map< std::set< index_t >, signed_index_t > edge_boundary_flags_ ;
            std::map< index_t, signed_index_t > surface_boundary_flags_ ;
        } ;

        /************************************************************************/

        class GPRSIOHandler: public MacroMeshIOHandler {
        public:
            struct Pipe {
                Pipe( index_t v0_in, index_t v1_in )
                    : v0( v0_in ), v1( v1_in )
                {
                }
                index_t v0 ;
                index_t v1 ;
            } ;
            virtual bool load( const std::string& filename, MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                    << "Loading of a MacroMesh from GPRS not implemented yet"
                    << std::endl ;
                return false ;
            }
            virtual bool save( const MacroMesh& mm, const std::string& filename )
            {
                std::string path = GEO::FileSystem::dir_name( filename ) ;
                std::string directory = GEO::FileSystem::base_name( filename ) ;
                if( path == "." ) {
                    path = GEO::FileSystem::get_current_working_directory() ;
                }
                std::ostringstream oss ;
                oss << path << "/" << directory ;
                std::string full_path = oss.str() ;
                GEO::FileSystem::create_directory( full_path ) ;

                std::ostringstream oss_pipes ;
                oss_pipes << full_path << "/pipes.in" ;
                std::ofstream out_pipes( oss_pipes.str().c_str() ) ;

                std::ostringstream oss_xyz ;
                oss_xyz << full_path << "/gprs.xyz" ;
                std::ofstream out_xyz( oss_xyz.str().c_str() ) ;
                out_xyz.precision( 16 ) ;

                const BoundaryModel& model = mm.model() ;
                std::vector< index_t > region_offsets( mm.nb_meshes(), 0 ) ;
                std::vector< index_t > surface_offsets( model.nb_surfaces(), 0 ) ;
                for( index_t r = 0; r < mm.nb_meshes()-1; r++ ) {
                    region_offsets[r+1] += region_offsets[r] + mm.mesh( r ).cells.nb() ;
                }
                index_t last_region = mm.nb_meshes()-1 ;
                surface_offsets[0] += region_offsets[last_region] + mm.mesh( last_region ).cells.nb() ;
                for( index_t s = 0; s < model.nb_surfaces()-1; s++ ) {
                    surface_offsets[s+1] += surface_offsets[s] + model.surface( s ).nb_cells() ;
                }

                std::deque< Pipe > pipes ;
                for( index_t r = 0; r < mm.nb_meshes(); r++ ) {
                    const GEO::Mesh& mesh = mm.mesh( r ) ;
                    const BoundaryModelElement& region = model.region( r ) ;
                    std::vector< index_t > boundary_ids( region.nb_boundaries() ) ;
                    for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                        boundary_ids[s] = region.boundary_id( s ) ;
                    }
                    index_t cell_offset = region_offsets[r] ;
                    std::stack< index_t > S ;
                    std::vector< bool > visited( mesh.cells.nb(), false ) ;
                    S.push( 0 ) ;
                    visited[0] = true ;
                    do {
                        index_t c = S.top() ;
                        S.pop() ;
                        for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                            index_t adj = mesh.cells.adjacent( c, f ) ;
                            if( adj != GEO::NO_CELL ) {
                                if( visited[adj] ) continue ;
                                pipes.push_back(
                                    Pipe( c + cell_offset, adj + cell_offset ) ) ;
                                S.push( adj ) ;
                                visited[adj] = true ;
                            } else {
                                vec3 query = Utils::mesh_cell_facet_center( mesh, c, f ) ;
                                for( index_t s = 0; s < boundary_ids.size(); s++ ) {
                                    index_t s_id = boundary_ids[s] ;
                                    index_t surface_offset = surface_offsets[s_id] ;
                                    std::vector< index_t > results ;
                                    if( model.surface( s_id ).tools.ann().get_colocated(
                                        query, results ) ) {
                                        pipes.push_back( Pipe( c + cell_offset, results[0] + surface_offset ) ) ;
                                    }
                                }
                            }
                        }
                    } while( !S.empty() ) ;
                }
                for( index_t s = 0; s < model.nb_surfaces(); s++ ) {
                    const Surface& surface = model.surface( s ) ;
                    index_t surface_offset = surface_offsets[s] ;
                    std::stack< index_t > S ;
                    std::vector< bool > visited( surface.nb_cells(), false ) ;
                    S.push( 0 ) ;
                    visited[0] = true ;
                    do {
                        index_t f = S.top() ;
                        S.pop() ;
                        for( index_t e = 0; e < surface.nb_vertices_in_facet( f ); e++ ) {
                            index_t adj = surface.adjacent( f, e ) ;
                            if( adj != GEO::NO_CELL ) {
                                if( visited[adj] ) continue ;
                                pipes.push_back(
                                    Pipe( f + surface_offset, adj + surface_offset ) ) ;
                                S.push( adj ) ;
                                visited[adj] = true ;
                            }
                        }
                    } while( !S.empty() ) ;
                }

                out_pipes << pipes.size() << std::endl ;
                for( index_t p = 0; p < pipes.size(); p++ ) {
                    const Pipe& pipe = pipes[p] ;
                    out_pipes << pipe.v0 << SPACE << pipe.v1 << std::endl ;
                }


                out_xyz
                    << "Node geometry, not used by GPRS but useful to reconstruct a pipe-network"
                    << std::endl ;
                for( index_t r = 0; r < mm.nb_meshes(); r++ ) {
                    const GEO::Mesh& mesh = mm.mesh( r ) ;
                    for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                        out_xyz << Utils::mesh_cell_center( mesh, c ) << std::endl ;
                    }
                }


                return true ;
            }

        } ;

        /************************************************************************/

        MacroMeshIOHandler* MacroMeshIOHandler::create( const std::string& format )
        {
            ringmesh_register_MacroMeshIOHandler_creator( MMIOHandler, "mm" ) ;
            ringmesh_register_MacroMeshIOHandler_creator( MESHBIOHandler, "meshb" );
            ringmesh_register_MacroMeshIOHandler_creator( TetGenIOHandler, "tetgen" );
            ringmesh_register_MacroMeshIOHandler_creator( TSolidIOHandler, "so" );
            ringmesh_register_MacroMeshIOHandler_creator( CSMPIOHandler, "csmp" );
            ringmesh_register_MacroMeshIOHandler_creator( AsterIOHandler, "mail" );
            ringmesh_register_MacroMeshIOHandler_creator( VTKIOHandler, "vtk" );
            ringmesh_register_MacroMeshIOHandler_creator( GPRSIOHandler, "gprs" );

            MacroMeshIOHandler* handler = MacroMeshIOHandlerFactory::create_object(
                format ) ;
            if( handler ) {
                return handler ;
            }

            GEO::Logger::err( "I/O" ) << "Unsupported file format: " << format
                << std::endl ;
            return nil ;
        }

        MacroMeshIOHandler* MacroMeshIOHandler::get_handler(
            const std::string& filename )
        {
            std::string ext = GEO::FileSystem::extension( filename ) ;
            return create( ext ) ;
        }
    }
}
