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
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_builder.h>
#include <geogram/mesh/mesh_private.h>
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

        void zip_file(
            zipFile zf,
            const std::string& name )
        {
            zip_fileinfo zfi = {
                0
            } ;
            std::fstream file( name.c_str(), std::ios::in ) ;
            file.seekg( 0, std::ios::end ) ;
            long size = static_cast<long>( file.tellg() ) ;
            file.seekg( 0, std::ios::beg ) ;
            std::vector< char > buffer( size ) ;
            file.read( &buffer[ 0 ], size ) ;
            zipOpenNewFileInZip( zf, name.c_str(), &zfi,
                NULL, 0, NULL, 0, NULL, Z_DEFLATED, Z_DEFAULT_COMPRESSION ) ;
            zipWriteInFileInZip( zf, size == 0 ? "" : &buffer[ 0 ], size ) ;
            zipCloseFileInZip( file ) ;
            file.close() ;
        }


        bool unzip_file(
            unzFile uz,
            char filename[ MAX_FILENAME ] )
        {
            char read_buffer[ READ_SIZE ] ;
            unz_file_info file_info ;
            if( unzGetCurrentFileInfo( uz, &file_info, filename,
                    MAX_FILENAME,
                    NULL, 0, NULL, 0 ) != UNZ_OK )
            {
                GEO::Logger::err( "could not read file global info" ) ;
                unzClose( uz ) ;
                return false ;
            }
            if( unzOpenCurrentFile( uz ) != UNZ_OK ) {
                GEO::Logger::err( "could not open file" ) ;
                unzClose( uz ) ;
                return false ;
            }
            FILE* out = fopen( filename, "wb" ) ;
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

        bool load(
            const std::string& filename,
            BoundaryModel& model )
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


        bool save(
            BoundaryModel& model,
            const std::string& filename )
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

        class MLIOHandler : public BoundaryModelIOHandler {
    public:
            virtual bool load(
                const std::string& filename,
                BoundaryModel& model )
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

            virtual bool save(
                BoundaryModel& model,
                const std::string& filename )
            {
                std::ofstream out( filename.c_str() ) ;
                return model.save_gocad_model3d( out ) ;
            }
        } ;

        class BMIOHandler : public BoundaryModelIOHandler {
    public:
            virtual bool load(
                const std::string& filename,
                BoundaryModel& model )
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

            virtual bool save(
                BoundaryModel& model,
                const std::string& filename )
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
            ringmesh_register_BoundaryModelIOHandler_creator( BMIOHandler, "bm" ) ;

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

        bool load(
            const std::string& filename,
            MacroMesh& model )
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


        bool save(
            const MacroMesh& model,
            const std::string& filename )
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
        class AsterIOHandler : public MacroMeshIOHandler {
    public:
            virtual bool load(
                const std::string& filename,
                MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                << "Loading of a MacroMesh from Code_Aster mesh not implemented yet"
                << std::endl ;
                return false ;
            }

            virtual bool save(
                const MacroMesh& mm,
                const std::string& filename )
            {
                MacroMeshExport db( mm ) ;
                db.compute_database( MacroMeshExport::FAULT ) ;
                const BoundaryModel& model = mm.model() ;
                std::vector< index_t > vertex_exported_id( mm.nb_vertices(),
                                                           NO_ID ) ;
                std::vector< index_t > atom_exported_id( db.nb_duplicated_vertices(),
                                                         NO_ID ) ;
                std::ofstream out( filename.c_str() ) ;
                out.precision( 16 ) ;
                std::vector< bool > vertex_exported( mm.nb_vertices(), false ) ;
                std::vector< bool > atom_exported( db.nb_duplicated_vertices(),
                                                   false ) ;

                // 1. Write the vertices coordinates (with the duplicate ones)
//                for( index_t v = 0; v < db.nb_total_vertices(); v++ ) {
//                    const vec3& cur_v = db. ;
//                    out << "V" << v << SPACE << cur_v.x << SPACE << cur_v.y << SPACE
//                        << cur_v.z << std::endl ;
//
//                }
                std::cout << "nb vertices " << mm.nb_vertices() << std::endl ;
                std::cout << "nb total vertices " << db.nb_total_vertices()
                          << std::endl ;

                index_t nb_vertices_exported = 0 ;
                index_t cur_cell = 0 ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const GEO::Mesh& mesh = mm.mesh( r ) ;
                    out << "COOR_3D" << std::endl ;
                    for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                        for( index_t co = mesh.cell_vertices_begin( c );
                             co
                             < mesh.cell_vertices_begin( c )
                             + mesh.cell_nb_vertices( c ); co++ )
                        {
                            index_t vertex_id, atom_id ;
                            if( db.vertex_id( r, co, vertex_id, atom_id ) ) {
                                if( vertex_exported[ vertex_id ] ) {continue ;}
                                vertex_exported[ vertex_id ] = true ;
                                vertex_exported_id[ vertex_id ] =
                                    nb_vertices_exported ;
                                out << "V" << nb_vertices_exported++ << " "
                                    << mm.global_vertex( vertex_id ) << std::endl ;
                            }
                        }
                    }
                    for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                        for( index_t co = mesh.cell_vertices_begin( c );
                             co
                             < mesh.cell_vertices_begin( c )
                             + mesh.cell_nb_vertices( c ); co++ )
                        {
                            index_t vertex_id, atom_id ;
                            if( !db.vertex_id( r, co, vertex_id, atom_id ) ) {
                                if( atom_exported[ atom_id ] ) {continue ;}
                                atom_exported[ atom_id ] = true ;
                                atom_exported_id[ atom_id ] = nb_vertices_exported ;
                                out << "V" << nb_vertices_exported++ << " "
                                    << mm.global_vertex( vertex_id ) << std::endl ;
                            }
                        }
                    }
                    out << "FINSF" << std::endl ;

                    out << "TETRA4" << std::endl ;

                    ColocaterANN ann( mesh, ColocaterANN::FACETS ) ;
                    for( index_t c = 0; c < db.nb_tet( r ); c++ ) {
                        index_t cur_tet = db.local_tet_id( r, c ) ;
                        out << "C" << cur_cell++ << " " ;

                        for( index_t co = mesh.cell_vertices_begin( cur_tet );
                             co
                             < mesh.cell_vertices_begin( cur_tet )
                             + mesh.cell_nb_vertices( cur_tet ); co++ )
                        {
                            index_t vertex_id, atom_id ;
                            if( db.vertex_id( r, co, vertex_id, atom_id ) ) {
                                out << "V" << vertex_exported_id[ vertex_id ] << " " ;
                            } else {
                                out << "V" << atom_exported_id[ atom_id ] << " " ;
                            }
                        }
                        out << std::endl ;
                    }
                    out << "FINSF" << std::endl ;

//                    int triangle_count = 0 ;
//                    out << "TRI3" << std::endl ;
//
//                    for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
//                        const RINGMesh::BoundaryModelElement& interf =
//                            model.one_interface( i ) ;
//                        for( index_t s = 0; s < interf.nb_children(); s++ ) {
//
//                            index_t surface_id = interf.child_id( s ) ;
//                            index_t mesh_id = mm.surface_mesh( surface_id ) ;
//                            const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
//                            out << std::endl ;
//                            for( index_t f = 0;
//                                f < mm.nb_surface_facets( surface_id ); f++ ) {
//                                index_t facet_id = mm.surface_facet( surface_id,
//                                    f ) ;
//                                out << "F" << triangle_count++ << " " ;
//                                for( index_t v = mesh.facet_begin( facet_id );
//                                    v < mesh.facet_end( facet_id ); v++ ) {
//                                    out
//                                        << vertex_exported_id[mm.global_vertex_id(
//                                            mesh_id, mesh.corner_vertex_index( v ) )]
//                                        << " " ;
//                                }
//                                out << std::endl ;
//                            }
//                        }
//                    }
//                    out << "FINSF" << std::endl ;
//
//                    triangle_count = 0 ;

//                    for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
//                        const RINGMesh::BoundaryModelElement& interf =
//                            model.one_interface( i ) ;
//                        for( index_t s = 0; s < interf.nb_children(); s++ ) {
//
//                            index_t surface_id = interf.child_id( s ) ;
//                            index_t mesh_id = mm.surface_mesh( surface_id ) ;
//                            out << "GROUP_MA nom = s" << surface_id << std::endl ;
//                            const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
//                            out << std::endl ;
//                            for( index_t f = 0;
//                                f < mm.nb_surface_facets( surface_id ); f++ ) {
//                                index_t facet_id = mm.surface_facet( surface_id,
//                                    f ) ;
//                                out << "F" << triangle_count++ << " " ;
//                                out << std::endl ;
//                            }
//                            out << "FINSF" << std::endl ;
//
//                        }
//                    }
                }
                out << "FIN" << std::endl ;

//                    //2. Write the cells per type
//                    //2.1 Tetraedrons
//                    out << "TETRA4" << std::endl ;

//                    index_t cell_id = 0 ;
//                    for( index_t r = 0; r < mm.nb_meshes(); r++ ) {
//                        const GEO::Mesh& cur_r = mm.mesh( r ) ;
//                        for( index_t tet = 0; tet < db.nb_tet( r ); tet++ ) {
//                            index_t local_tet = db.local_tet_id( r, tet ) ;
//                            out << "C" << cell_id << SPACE ;
//                            for( index_t co = cur_r.cell_vertices_begin( local_tet );
//                                co < cur_r.cell_vertices_begin( local_tet ) + 4;
//                                co++ ) {
//                                index_t global_co ;
//                                index_t orig_co ;
//                                if( db.vertex_id( r, co, orig_co, global_co ) ) {
//                                    out << "V" << orig_co << SPACE ;
//                                } else {
//                                    out << "V" << global_co << SPACE ;
//
//                                }
//                            }
//                            out << std::endl ;
//                            cell_id++ ;
//                        }
//                    }
//                    out << "FINSF" << std::endl ;

                out.close() ;

                return true ;
            }
        } ;

        /************************************************************************/

        class MMIOHandler : public MacroMeshIOHandler {
    public:
            virtual bool load(
                const std::string& filename,
                MacroMesh& mm )
            {
                unzFile uz = unzOpen( filename.c_str() ) ;
                unz_global_info global_info ;
                if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
                    GEO::Logger::err( "could not read file global info" ) ;
                    unzClose( uz ) ;
                    return false ;
                }
                for( index_t r = 0; r < mm.model().nb_regions(); r++ ) {
                    char filename[ MAX_FILENAME ] ;
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
                        ringmesh_assert_not_reached ;
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

            // / Save a \param[in] macro mesh in a .zip file which contains all the mesh file. Type of the export is
            // / determined by the extension given in \param[in] filename
            virtual bool save(
                const MacroMesh& mm,
                const std::string& filename )
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
                                                 + GEO::String::to_string( m ) +
                                                 ".meshb" ;

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

        class MESHBIOHandler : public MacroMeshIOHandler {
    public:
            virtual bool load(
                const std::string& filename,
                MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                << "Loading of a MacroMesh from a meshb not implemented yet"
                << std::endl ;
                return false ;
            }

            virtual bool save(
                const MacroMesh& mm,
                const std::string& filename )
            {
                GEO::Mesh mesh( 3 ) ;
                GEO::MeshBuilder builder( &mesh ) ;
                builder.begin_mesh() ;
                MacroMesh& MM = const_cast< MacroMesh& >( mm ) ;
                GEO::MeshMutator::set_attributes( mesh, GEO::MESH_FACET_REGION ) ;
                for( index_t p = 0; p < MM.nb_vertices(); p++ ) {
                    builder.add_vertex( MM.global_vertex( p ) ) ;
                }

                index_t cell_offset = 0 ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& cur_mesh = mm.mesh( m ) ;
                    for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                        builder.begin_facet() ;
                        for( index_t v = cur_mesh.facet_begin( f );
                             v < cur_mesh.facet_end( f ); v++ )
                        {
                            builder.add_vertex_to_facet(
                                MM.global_vertex_id( m,
                                    cur_mesh.corner_vertex_index( v ) ) ) ;
                        }
                        builder.end_facet( cur_mesh.facet_region( f ) ) ;
                    }

                    for( index_t c = 0; c < cur_mesh.nb_cells(); c++ ) {
                        index_t vertex_indices[ 8 ] ;
                        for( unsigned int v = 0; v < cur_mesh.cell_nb_vertices( c );
                             v++ )
                        {
                            vertex_indices[ v ] = MM.global_vertex_id( m,
                                cur_mesh.cell_vertex_index( c, v ) ) ;
                        }
                        signed_index_t adj_indices[ 6 ] ;
                        for( index_t f = 0; f < cur_mesh.cell_nb_facets( c );
                             f++ )
                        {
                            signed_index_t adj = cur_mesh.cell_adjacent( c, f ) ;
                            adj_indices[ f ] = adj == - 1 ? adj : cell_offset + adj ;
                        }
                        mesh.add_cell( cur_mesh.cell_type( c ), vertex_indices,
                            adj_indices ) ;
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

        class TetGenIOHandler : public MacroMeshIOHandler {
    public:
            virtual bool load(
                const std::string& filename,
                MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                << "Loading of a MacroMesh from TetGen not implemented yet"
                << std::endl ;
                return false ;
            }

            virtual bool save(
                const MacroMesh& mm,
                const std::string& filename )
            {
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

                ele << mm.nb_cells() << " 4 1" << std::endl ;
                neigh << mm.nb_cells() << " 4" << std::endl ;
                index_t nb_tet_exported = 0 ;
                for( index_t m = 0; m < mm.nb_meshes(); m++ ) {
                    const GEO::Mesh& mesh = mm.mesh( m ) ;
                    for( index_t tet = 0; tet < mesh.nb_cells(); tet++ ) {
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
                            neigh << SPACE << mm.global_cell_adjacent( m, tet, f ) ;
                        }
                        neigh << std::endl ;
                    }
                    nb_tet_exported += mesh.nb_cells() ;
                }
                return true ;
            }
        } ;

        /************************************************************************/

        class TSolidIOHandler : public MacroMeshIOHandler {
    public:
            virtual bool load(
                const std::string& filename,
                MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                << "Loading of a MacroMesh from TSolid not implemented yet"
                << std::endl ;
                return false ;
            }

            virtual bool save(
                const MacroMesh& mm,
                const std::string& filename )
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

                MacroMeshExport db( mm ) ;
                db.compute_database( MacroMeshExport::ALL ) ;

                std::vector< bool > vertex_exported( mm.nb_vertices(), false ) ;
                std::vector< bool > atom_exported( db.nb_duplicated_vertices(),
                                                   false ) ;
                std::vector< index_t > vertex_exported_id( mm.nb_vertices(),
                                                           NO_ID ) ;
                std::vector< index_t > atom_exported_id( db.nb_duplicated_vertices(),
                                                         NO_ID ) ;
                index_t nb_vertices_exported = 1 ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const RINGMesh::BoundaryModelElement& region = model.region( r ) ;
                    out << "TVOLUME " << region.name() << std::endl ;

                    const GEO::Mesh& mesh = mm.mesh( r ) ;
                    for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                        for( index_t co = mesh.cell_vertices_begin( c );
                             co
                             < mesh.cell_vertices_begin( c )
                             + mesh.cell_nb_vertices( c ); co++ )
                        {
                            index_t vertex_id, atom_id ;
                            if( db.vertex_id( r, co, vertex_id, atom_id ) ) {
                                if( vertex_exported[ vertex_id ] ) {continue ;}
                                vertex_exported[ vertex_id ] = true ;
                                vertex_exported_id[ vertex_id ] =
                                    nb_vertices_exported ;
                                out << "VRTX " << nb_vertices_exported++ << " "
                                    << mm.global_vertex( vertex_id ) << std::endl ;
                            }
                        }
                    }

                    for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                        for( index_t co = mesh.cell_vertices_begin( c );
                             co
                             < mesh.cell_vertices_begin( c )
                             + mesh.cell_nb_vertices( c ); co++ )
                        {
                            index_t vertex_id, atom_id ;
                            if( !db.vertex_id( r, co, vertex_id, atom_id ) ) {
                                if( atom_exported[ atom_id ] ) {continue ;}
                                atom_exported[ atom_id ] = true ;
                                atom_exported_id[ atom_id ] = nb_vertices_exported ;
                                out << "ATOM " << nb_vertices_exported++ << " "
                                    << vertex_exported_id[ vertex_id ] << std::endl ;
                            }
                        }
                    }

                    std::map< index_t, index_t > sides ;
                    for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                        if( sides.count( region.boundary_id( s ) ) > 0 ) {
                            sides[ region.boundary_id( s ) ] = 2 ;
                        } else {
                            sides[ region.boundary_id( s ) ] = region.side( s ) ;
                        }
                    }

                    ColocaterANN ann( mesh, ColocaterANN::FACETS ) ;
                    for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                        out << "TETRA" ;
                        for( index_t co = mesh.cell_vertices_begin( c );
                             co
                             < mesh.cell_vertices_begin( c )
                             + mesh.cell_nb_vertices( c ); co++ )
                        {
                            index_t vertex_id, atom_id ;
                            if( db.vertex_id( r, co, vertex_id, atom_id ) ) {
                                out << " " << vertex_exported_id[ vertex_id ] ;
                            } else {
                                out << " " << atom_exported_id[ atom_id ] ;
                            }
                        }
                        out << std::endl ;
                        out << "# CTETRA " << region.name() ;
                        for( index_t f = 0; f < mesh.cell_nb_facets( c ); f++ ) {
                            out << " " ;
                            vec3 facet_center = Utils::mesh_cell_facet_center( mesh,
                                c, f ) ;
                            std::vector< index_t > result ;
                            if( ann.get_colocated( facet_center, result ) ) {
                                ringmesh_debug_assert( result.size() == 1 ) ;
                                index_t surface_id = mesh.facet_region( result[ 0 ] ) ;
                                index_t side = sides[ surface_id ] ;
                                if( side < 2 ) {
                                    side ? out << "+" : out << "-" ;
                                } else {
                                    vec3 cell_facet_normal =
                                        Utils::mesh_cell_facet_normal( mesh, c, f ) ;
                                    vec3 facet_normal = GEO::Geom::mesh_facet_normal(
                                        mesh, result[ 0 ] ) ;
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
                        index_t mesh_id = mm.surface_mesh( surface_id ) ;
                        const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                        out << "KEYVERTICES" ;
                        index_t key_facet_id = mm.surface_facet( surface_id, 0 ) ;
                        for( index_t v = mesh.facet_begin( key_facet_id );
                             v < mesh.facet_end( key_facet_id ); v++ )
                        {
                            out << " "
                                << vertex_exported_id[ mm.global_vertex_id( mesh_id,
                                                       mesh.corner_vertex_index( v ) )
                            ] ;
                        }
                        out << std::endl ;
                        for( index_t f = 0; f < mm.nb_surface_facets( surface_id );
                             f++ )
                        {
                            index_t facet_id = mm.surface_facet( surface_id, f ) ;
                            out << "TRGL" ;
                            for( index_t v = mesh.facet_begin( facet_id );
                                 v < mesh.facet_end( facet_id ); v++ )
                            {
                                out << " "
                                    << vertex_exported_id[ mm.global_vertex_id(
                                                           mesh_id,
                                                           mesh.corner_vertex_index(
                                                               v ) ) ] ;
                            }
                            out << std::endl ;
                        }
                    }
                }

                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const RINGMesh::BoundaryModelElement& region = model.region( r ) ;
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
            index_t vertices[ 8 ] ;
            index_t nb_facets ;
            index_t nb_vertices_in_facet[ 6 ] ;
            index_t facet[ 6 ] ;
            index_t vertices_in_facet[ 6 ][ 4 ] ;
        } ;

        static RINGMesh2CSMP tet_descriptor = {
            4,                                                      // type
            4,                  // nb vertices
            { 0, 1, 2, 3 },     // vertices
            4,                  // nb facets
            { 3, 3, 3, 3 },     // nb vertices in facet
            { 0, 1, 2, 3 },     // facets
            {
                { 1, 3, 2 }, { 0, 2, 3 }, { 3, 1, 0 }, { 0, 1, 2 }
            }
        } ;

        static RINGMesh2CSMP hex_descriptor = {
            6,                                                             // type
            8,                              // nb vertices
            { 0, 1, 3, 2, 4, 5, 7, 6 },     // vertices
            6,                              // nb facets
            { 4, 4, 4, 4, 4, 4 },           // nb vertices in facet
            { 4, 2, 1, 3, 0, 5 },           // facets
            {
                { 0, 3, 7, 4 }, { 2, 1, 5, 6 }, { 1, 0, 4, 5 }, { 3, 2, 6, 7 }, {
                    1, 2, 3, 0
                }, { 4, 7, 6, 5 }
            }
        } ;

        static RINGMesh2CSMP prism_descriptor = {
            12,                                                           // type
            6,                      // nb vertices
            { 0, 1, 2, 3, 4, 5 },   // vertices
            5,                      // nb facets
            { 3, 4, 4, 4, 3 },      // nb vertices in facet
            { 0, 2, 4, 3, 1 },      // facets
            {
                { 0, 1, 2 }, { 3, 5, 4 }, { 0, 3, 4, 1 }, { 0, 2, 5, 3 }, {
                    1, 4, 5, 2
                }
            }
        } ;

        static RINGMesh2CSMP pyramid_descriptor = {
            18,                                                         // type
            5,                  // nb vertices
            { 0, 1, 2, 3, 4 },  // vertices
            5,                  // nb facets
            { 3, 3, 3, 3, 4 },  // nb vertices in facet
            { 1, 3, 4, 2, 0 },  // facets
            {
                { 0, 1, 2, 3 }, { 0, 4, 1 }, { 0, 3, 4 }, { 2, 4, 3 }, { 2, 1, 4 }
            }
        } ;

        class CSMPIOHandler : public MacroMeshIOHandler {
    public:
            virtual bool load(
                const std::string& filename,
                MacroMesh& mesh )
            {
                GEO::Logger::err( "I/O" )
                << "Loading of a MacroMesh from CSMP not implemented yet"
                << std::endl ;
                return false ;
            }

            virtual bool save(
                const MacroMesh& mm,
                const std::string& filename )
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

                index_t count = 0 ;

                // Conversion from (X,Y,Z) to (X,Z,-Y)
                signed_index_t conversion_sign[ 3 ] = {
                    1, 1, - 1
                } ;
                index_t conversion_axis[ 3 ] = {
                    0, 2, 1
                } ;
                data << mm.nb_vertices() << " # PX, PY, PZ" << std::endl ;
                for( index_t dim = 0; dim < 3; dim++ ) {
                    for( index_t v = 0; v < mm.nb_vertices(); v++ ) {
                        data << " "
                             << conversion_sign[ dim ]
                        * mm.global_vertex( v )[ conversion_axis[ dim ] ] ;
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

                MacroMeshExport db( mm ) ;
                db.compute_database() ;

                index_t nb_families = 0 ;
                std::vector< index_t > nb_triangle_interface( model.nb_interfaces(),
                                                              0 ) ;
                std::vector< index_t > nb_quad_interface( model.nb_interfaces(),
                                                          0 ) ;
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    const BoundaryModelElement& interf = model.one_interface( i ) ;
                    for( index_t s = 0; s < interf.nb_children(); s++ ) {
                        index_t s_id = interf.child_id( s ) ;
                        nb_triangle_interface[ i ] += db.nb_triangle( s_id ) ;
                        nb_quad_interface[ i ] += db.nb_quad( s_id ) ;
                    }
                    if( nb_triangle_interface[ i ] > 0 ) {nb_families++ ;}
                    if( nb_quad_interface[ i ] > 0 ) {nb_families++ ;}
                }
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    if( db.nb_tet( r ) > 0 ) {nb_families++ ;}
                    if( db.nb_pyramid( r ) > 0 ) {nb_families++ ;}
                    if( db.nb_prism( r ) > 0 ) {nb_families++ ;}
                    if( db.nb_hex( r ) > 0 ) {nb_families++ ;}
                }

                ascii << nb_families << " # Number of families" << std::endl ;
                ascii << "# Object name" << TAB << "Element type" << TAB
                      << "Material-ID" << TAB << "Number of elements" << std::endl ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const RINGMesh::BoundaryModelElement& region = model.region( r ) ;
                    if( db.nb_tet( r ) > 0 ) {
                        ascii << region.name() << TAB << "TETRA_4" << TAB << 0 <<
                        TAB
                              << db.nb_tet( r ) << std::endl ;
                    }
                    if( db.nb_pyramid( r ) > 0 ) {
                        ascii << region.name() << TAB << "PYRA_5" << TAB << 0 << TAB
                              << db.nb_pyramid( r ) << std::endl ;
                    }
                    if( db.nb_prism( r ) > 0 ) {
                        ascii << region.name() << TAB << "PENTA_6" << TAB << 0 <<
                        TAB
                              << db.nb_prism( r ) << std::endl ;
                    }
                    if( db.nb_hex( r ) > 0 ) {
                        ascii << region.name() << TAB << "HEXA_8" << TAB << 0 << TAB
                              << db.nb_hex( r ) << std::endl ;
                    }
                }
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    if( nb_triangle_interface[ i ] > 0 ) {
                        ascii << interface_name( i, mm ) << TAB << "TRI_3" << TAB
                              << 0 << TAB << nb_triangle_interface[ i ] << std::endl ;
                    }
                    if( nb_quad_interface[ i ] > 0 ) {
                        ascii << interface_name( i, mm ) << TAB << "QUAD_4" << TAB
                              << 0 << TAB << nb_quad_interface[ i ] << std::endl ;
                    }
                }

                data << "# PBFLAGS" << std::endl ;
                for( index_t p = 0; p < mm.nb_vertices(); p++ ) {
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
                for( index_t p = 0; p < mm.nb_vertices(); p++ ) {
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

                index_t nb_total_elements = db.nb_cells() + db.nb_facets() ;
                data << nb_total_elements << " # PELEMENT" << std::endl ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    for( index_t el = 0; el < db.nb_tet( r ); el++ ) {
                        data << " " << std::setw( 3 ) << 4 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                    for( index_t el = 0; el < db.nb_pyramid( r ); el++ ) {
                        data << " " << std::setw( 3 ) << 18 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                    for( index_t el = 0; el < db.nb_prism( r ); el++ ) {
                        data << " " << std::setw( 3 ) << 12 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                    for( index_t el = 0; el < db.nb_hex( r ); el++ ) {
                        data << " " << std::setw( 3 ) << 6 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                }
                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                    for( index_t el = 0; el < nb_triangle_interface[ i ]; el++ ) {
                        data << " " << std::setw( 3 ) << 8 ;
                        count++ ;
                        if( count == 20 ) {
                            data << std::endl ;
                            count = 0 ;
                        }
                    }
                    for( index_t el = 0; el < nb_quad_interface[ i ]; el++ ) {
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
                <<
                "# now the elements which make up each object are listed in sequence"
                << std::endl ;
                index_t cur_cell = 0 ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const RINGMesh::BoundaryModelElement& region = model.region( r ) ;
                    if( db.nb_tet( r ) > 0 ) {
                        ascii << region.name() << " " << "TETRA_4" << " "
                              << db.nb_tet( r ) << std::endl ;
                        for( index_t el = 0; el < db.nb_tet( r ); el++ ) {
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
                    if( db.nb_pyramid( r ) > 0 ) {
                        ascii << region.name() << " " << "PYRA_5" << " "
                              << db.nb_pyramid( r ) << std::endl ;
                        for( index_t el = 0; el < db.nb_pyramid( r ); el++ ) {
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
                    if( db.nb_prism( r ) > 0 ) {
                        ascii << region.name() << " " << "PENTA_6" << " "
                              << db.nb_prism( r ) << std::endl ;
                        for( index_t el = 0; el < db.nb_prism( r ); el++ ) {
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
                    if( db.nb_hex( r ) > 0 ) {
                        ascii << region.name() << " " << "HEXA_8" << " "
                              << db.nb_hex( r ) << std::endl ;
                        for( index_t el = 0; el < db.nb_hex( r ); el++ ) {
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
                    if( nb_triangle_interface[ i ] > 0 ) {
                        ascii << interface_name( i, mm ) << " " << "TRI_3" << " "
                              << nb_triangle_interface[ i ] << std::endl ;
                        for( index_t el = 0; el < nb_triangle_interface[ i ];
                             el++ )
                        {
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
                    if( nb_quad_interface[ i ] > 0 ) {
                        ascii << interface_name( i, mm ) << " " << "QUAD_4" << " "
                              << nb_quad_interface[ i ] << std::endl ;
                        for( index_t el = 0; el < nb_quad_interface[ i ]; el++ ) {
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

                index_t nb_plist = 3 * db.nb_triangle() + 4 * db.nb_quad()
                                   + 4 * db.nb_tet() + 5 * db.nb_pyramid() + 6 *
                                   db.nb_prism()
                                   + 8 * db.nb_hex() ; // + 2 * nb_edges ;
                data << nb_plist << " # PLIST" << std::endl ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    const GEO::Mesh& mesh = mm.mesh( r ) ;
                    for( index_t el = 0; el < db.nb_tet( r ); el++ ) {
                        index_t tet = db.local_tet_id( r, el ) ;
                        for( index_t p = 0; p < 4; p++ ) {
                            index_t vertex_id ;
                            index_t csmp_p = tet_descriptor.vertices[ p ] ;
                            db.vertex_id( r,
                                mesh.cell_vertices_begin( tet ) + csmp_p,
                                vertex_id ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                    for( index_t el = 0; el < db.nb_pyramid( r ); el++ ) {
                        index_t py = db.local_pyramid_id( r, el ) ;
                        for( index_t p = 0; p < 5; p++ ) {
                            index_t vertex_id ;
                            index_t csmp_p = pyramid_descriptor.vertices[ p ] ;
                            db.vertex_id( r, mesh.cell_vertices_begin( py ) + csmp_p,
                                vertex_id ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                    for( index_t el = 0; el < db.nb_prism( r ); el++ ) {
                        index_t prism = db.local_prism_id( r, el ) ;
                        for( index_t p = 0; p < 6; p++ ) {
                            index_t vertex_id ;
                            index_t csmp_p = prism_descriptor.vertices[ p ] ;
                            db.vertex_id( r,
                                mesh.cell_vertices_begin( prism ) + csmp_p,
                                vertex_id ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            count++ ;
                            if( count == 10 ) {
                                data << std::endl ;
                                count = 0 ;
                            }
                        }
                    }
                    for( index_t el = 0; el < db.nb_hex( r ); el++ ) {
                        index_t hex = db.local_prism_id( r, el ) ;
                        for( index_t p = 0; p < 8; p++ ) {
                            index_t vertex_id ;
                            index_t csmp_p = hex_descriptor.vertices[ p ] ;
                            db.vertex_id( r,
                                mesh.cell_vertices_begin( hex ) + csmp_p,
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
                        index_t mesh_id = mm.surface_mesh( s_id ) ;
                        const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                        for( index_t el = 0; el < db.nb_triangle( s_id ); el++ ) {
                            index_t tri = db.local_triangle_id( mesh_id, el ) ;
                            for( index_t p = mesh.facet_begin( tri );
                                 p < mesh.facet_end( tri ); p++ )
                            {
                                index_t vertex_id = mm.global_vertex_id( mesh_id,
                                    mesh.corner_vertex_index( p ) ) ;
                                data << " " << std::setw( 7 ) << vertex_id ;
                                count++ ;
                                if( count == 10 ) {
                                    data << std::endl ;
                                    count = 0 ;
                                }
                            }
                        }
                        for( index_t el = 0; el < db.nb_quad( s_id ); el++ ) {
                            index_t quad = db.local_quad_id( mesh_id, el ) ;
                            for( index_t p = mesh.facet_begin( quad );
                                 p < mesh.facet_end( quad ); p++ )
                            {
                                index_t vertex_id = mm.global_vertex_id( mesh_id,
                                    mesh.corner_vertex_index( p ) ) ;
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

                int nb_facets = 3 * db.nb_triangle() + 4 * db.nb_quad()
                                + 4 * db.nb_tet() + 5 * db.nb_pyramid() + 5 *
                                db.nb_prism()
                                + 6 * db.nb_hex() ; // + 2 * nb_edges ;
                data << nb_facets << " # PFVERTS" << std::endl ;
                for( index_t r = 0; r < model.nb_regions(); r++ ) {
                    for( index_t el = 0; el < db.nb_tet( r ); el++ ) {
                        index_t tet = db.local_tet_id( r, el ) ;
                        for( index_t f = 0; f < 4; f++ ) {
                            index_t csmp_f = tet_descriptor.facet[ f ] ;
                            signed_index_t adj = mm.global_cell_adjacent( r, tet,
                                csmp_f ) ;
                            if( adj == - 1 ) {
                                data << " " << std::setw( 7 ) << - 28 ;
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
                    for( index_t el = 0; el < db.nb_pyramid( r ); el++ ) {
                        index_t py = db.local_pyramid_id( r, el ) ;
                        for( index_t f = 0; f < 5; f++ ) {
                            index_t csmp_f = pyramid_descriptor.facet[ f ] ;
                            signed_index_t adj = mm.global_cell_adjacent( r, py,
                                csmp_f ) ;
                            if( adj == - 1 ) {
                                data << " " << std::setw( 7 ) << - 28 ;
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
                    for( index_t el = 0; el < db.nb_prism( r ); el++ ) {
                        index_t prism = db.local_prism_id( r, el ) ;
                        for( index_t f = 0; f < 5; f++ ) {
                            index_t csmp_f = prism_descriptor.facet[ f ] ;
                            signed_index_t adj = mm.global_cell_adjacent( r, prism,
                                csmp_f ) ;
                            if( adj == - 1 ) {
                                data << " " << std::setw( 7 ) << - 28 ;
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
                    for( index_t el = 0; el < db.nb_hex( r ); el++ ) {
                        index_t hex = db.local_hex_id( r, el ) ;
                        for( index_t f = 0; f < 6; f++ ) {
                            index_t csmp_f = hex_descriptor.facet[ f ] ;
                            signed_index_t adj = mm.global_cell_adjacent( r, hex,
                                csmp_f ) ;
                            if( adj == - 1 ) {
                                data << " " << std::setw( 7 ) << - 28 ;
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
                        index_t mesh_id = mm.surface_mesh( s_id ) ;
                        const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                        for( index_t el = 0; el < db.nb_triangle( s_id ); el++ ) {
                            index_t tri = db.local_triangle_id( mesh_id, el ) ;
                            for( index_t f = mesh.facet_begin( tri );
                                 f < mesh.facet_end( tri ); f++ )
                            {
                                signed_index_t adj = mesh.corner_adjacent_facet(
                                    f ) ;
                                if( adj == - 1 ) {
                                    data << " " << std::setw( 7 ) << - 28 ;
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
                        for( index_t el = 0; el < db.nb_quad( s_id ); el++ ) {
                            index_t quad = db.local_quad_id( mesh_id, el ) ;
                            for( index_t f = mesh.facet_begin( quad );
                                 f < mesh.facet_end( quad ); f++ )
                            {
                                signed_index_t adj = mesh.corner_adjacent_facet(
                                    f ) ;
                                if( adj == - 1 ) {
                                    data << " " << std::setw( 7 ) << - 28 ;
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
                        if( parser.nb_fields() == 0 ) {continue ;}
                        if( parser.nb_fields() != 3 ) {return false ;}
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
                        || bottom_ == NO_ID || left_ == NO_ID || right_ == NO_ID )
                    {
                        GEO::Logger::err( "I/O" ) << "Missing box shape information"
                                                  << std::endl ;
                        return false ;
                    }

                    surface_boundary_flags_[ back_ ] = - 7 ;
                    surface_boundary_flags_[ top_ ] = - 5 ;
                    surface_boundary_flags_[ front_ ] = - 6 ;
                    surface_boundary_flags_[ bottom_ ] = - 4 ;
                    surface_boundary_flags_[ left_ ] = - 2 ;
                    surface_boundary_flags_[ right_ ] = - 3 ;

                    std::set< index_t > back_bottom ;
                    back_bottom.insert( back_ ) ;
                    back_bottom.insert( bottom_ ) ;
                    edge_boundary_flags_[ back_bottom ] = - 16 ;
                    std::set< index_t > back_right ;
                    back_right.insert( back_ ) ;
                    back_right.insert( right_ ) ;
                    edge_boundary_flags_[ back_right ] = - 17 ;
                    std::set< index_t > back_top ;
                    back_top.insert( back_ ) ;
                    back_top.insert( top_ ) ;
                    edge_boundary_flags_[ back_top ] = - 18 ;
                    std::set< index_t > back_left ;
                    back_left.insert( back_ ) ;
                    back_left.insert( left_ ) ;
                    edge_boundary_flags_[ back_left ] = - 19 ;
                    std::set< index_t > right_bottom ;
                    right_bottom.insert( right_ ) ;
                    right_bottom.insert( bottom_ ) ;
                    edge_boundary_flags_[ right_bottom ] = - 20 ;
                    std::set< index_t > right_top ;
                    right_top.insert( right_ ) ;
                    right_top.insert( top_ ) ;
                    edge_boundary_flags_[ right_top ] = - 21 ;
                    std::set< index_t > left_top ;
                    left_top.insert( left_ ) ;
                    left_top.insert( top_ ) ;
                    edge_boundary_flags_[ left_top ] = - 22 ;
                    std::set< index_t > left_bottom ;
                    left_bottom.insert( left_ ) ;
                    left_bottom.insert( bottom_ ) ;
                    edge_boundary_flags_[ left_bottom ] = - 23 ;
                    std::set< index_t > front_bottom ;
                    front_bottom.insert( front_ ) ;
                    front_bottom.insert( bottom_ ) ;
                    edge_boundary_flags_[ front_bottom ] = - 24 ;
                    std::set< index_t > front_right ;
                    front_right.insert( front_ ) ;
                    front_right.insert( right_ ) ;
                    edge_boundary_flags_[ front_right ] = - 25 ;
                    std::set< index_t > front_top ;
                    front_top.insert( front_ ) ;
                    front_top.insert( top_ ) ;
                    edge_boundary_flags_[ front_top ] = - 26 ;
                    std::set< index_t > front_left ;
                    front_left.insert( front_ ) ;
                    front_left.insert( left_ ) ;
                    edge_boundary_flags_[ front_left ] = - 27 ;

                    std::set< index_t > back_top_left ;
                    back_top_left.insert( back_ ) ;
                    back_top_left.insert( top_ ) ;
                    back_top_left.insert( left_ ) ;
                    corner_boundary_flags_[ back_top_left ] = - 13 ;
                    std::set< index_t > back_top_right ;
                    back_top_right.insert( back_ ) ;
                    back_top_right.insert( top_ ) ;
                    back_top_right.insert( right_ ) ;
                    corner_boundary_flags_[ back_top_right ] = - 14 ;
                    std::set< index_t > back_bottom_left ;
                    back_bottom_left.insert( back_ ) ;
                    back_bottom_left.insert( bottom_ ) ;
                    back_bottom_left.insert( left_ ) ;
                    corner_boundary_flags_[ back_bottom_left ] = - 8 ;
                    std::set< index_t > back_bottom_right ;
                    back_bottom_right.insert( back_ ) ;
                    back_bottom_right.insert( bottom_ ) ;
                    back_bottom_right.insert( right_ ) ;
                    corner_boundary_flags_[ back_bottom_right ] = - 10 ;
                    std::set< index_t > front_top_left ;
                    front_top_left.insert( front_ ) ;
                    front_top_left.insert( top_ ) ;
                    front_top_left.insert( left_ ) ;
                    corner_boundary_flags_[ front_top_left ] = - 15 ;
                    std::set< index_t > front_top_right ;
                    front_top_right.insert( front_ ) ;
                    front_top_right.insert( top_ ) ;
                    front_top_right.insert( right_ ) ;
                    corner_boundary_flags_[ front_top_right ] = - 9 ;
                    std::set< index_t > front_bottom_left ;
                    front_bottom_left.insert( front_ ) ;
                    front_bottom_left.insert( bottom_ ) ;
                    front_bottom_left.insert( left_ ) ;
                    corner_boundary_flags_[ front_bottom_left ] = - 12 ;
                    std::set< index_t > front_bottom_right ;
                    front_bottom_right.insert( front_ ) ;
                    front_bottom_right.insert( bottom_ ) ;
                    front_bottom_right.insert( right_ ) ;
                    corner_boundary_flags_[ front_bottom_right ] = - 11 ;
                }

                point_boundaries_.resize( mm.nb_vertices() ) ;
                for( index_t s = 0; s < model.nb_surfaces(); s++ ) {
                    index_t interface_id = model.surface( s ).parent_id() ;
                    index_t mesh_id = mm.surface_mesh( s ) ;
                    const GEO::Mesh& mesh = mm.mesh( mesh_id ) ;
                    for( index_t f = 0; f < mm.nb_surface_facets( s ); f++ ) {
                        index_t f_id = mm.surface_facet( s, f ) ;
                        for( index_t c = mesh.facet_begin( f_id );
                             c < mesh.facet_end( f_id ); c++ )
                        {
                            index_t vertex_id = mm.global_vertex_id( mesh_id,
                                mesh.corner_vertex_index( c ) ) ;
                            point_boundaries_[ vertex_id ].insert( interface_id ) ;
                        }
                    }
                }

                return true ;
            }

            std::string interface_name(
                index_t i,
                const MacroMesh& mm )
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
                const std::set< unsigned int >& boundaries = point_boundaries_[ p ] ;
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
                    if( boundaries.empty() ) {
                        return 0 ;
                    } else {
                        return - 28 ;
                    }
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

        MacroMeshIOHandler* MacroMeshIOHandler::create( const std::string& format )
        {
            ringmesh_register_MacroMeshIOHandler_creator( MMIOHandler, "mm" ) ;
            ringmesh_register_MacroMeshIOHandler_creator( MESHBIOHandler, "meshb" ) ;
            ringmesh_register_MacroMeshIOHandler_creator( TetGenIOHandler, "tetgen" ) ;
            ringmesh_register_MacroMeshIOHandler_creator( TSolidIOHandler, "so" ) ;
            ringmesh_register_MacroMeshIOHandler_creator( CSMPIOHandler, "csmp" ) ;
            ringmesh_register_MacroMeshIOHandler_creator( AsterIOHandler, "mail" ) ;

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


        /************************************************************************/

        MacroMeshExport::MacroMeshExport( const MacroMesh& mm )
              :
                mm_( mm ),
                facet_ptr_( ( NB_FACET_TYPES + 1 ) * mm.model().nb_surfaces(), 0 ),
                mesh_facet_ptr_( mm.model().nb_surfaces() + 1, 0 ),
                surface2mesh_( mm.model().nb_surfaces(), Surface::NO_ID ),
                cell_ptr_( ( NB_CELL_TYPES + 1 ) * mm.model().nb_regions(), 0 ),
                mesh_cell_ptr_( mm.model().nb_regions() + 1, 0 ),
                mesh_corner_ptr_( mm.model().nb_regions() + 1, 0 ),
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
            // / 1 - fill facet information
            index_t facet_access[ 5 ] = {
                - 1, - 1, - 1, 0, 1
            } ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& cur_mesh = mm_.mesh( m ) ;
                std::vector< signed_index_t > surface_proccessed ;
                for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                    signed_index_t surface_id = cur_mesh.facet_region( f ) ;
                    if( surface2mesh_[ surface_id ] != Surface::NO_ID ) {continue ;}
                    if( !Utils::contains( surface_proccessed, surface_id ) ) {
                        surface_proccessed.push_back( surface_id ) ;
                    }
                    facet_ptr_[ ( NB_FACET_TYPES + 1 ) * surface_id
                                + facet_access[ cur_mesh.facet_size( f ) ] + 1 ]++ ;
                }
                for( index_t s = 0; s < surface_proccessed.size(); s++ ) {
                    surface2mesh_[ surface_proccessed[ s ] ] = m ;
                }
            }

            for( index_t s = 0; s < mm_.model().nb_surfaces(); s++ ) {
                for( index_t t = 1; t < NB_FACET_TYPES; t++ ) {
                    facet_ptr_[ ( NB_FACET_TYPES + 1 ) * s + t + 1 ] +=
                        facet_ptr_[ ( NB_FACET_TYPES + 1 ) * s + t ] ;
                }
                mesh_facet_ptr_[ s + 1 ] += mesh_facet_ptr_[ s ] ;
                index_t nb_facets_in_surface =
                    facet_ptr_[ ( NB_FACET_TYPES + 1 ) * s
                                + NB_FACET_TYPES ] -
                    facet_ptr_[ ( NB_FACET_TYPES + 1 ) * s ] ;
                mesh_facet_ptr_[ s + 1 ] += nb_facets_in_surface ;
            }
            facets_.resize( mesh_facet_ptr_.back() ) ;

            std::vector< index_t > cur_facet_index_type(
                NB_FACET_TYPES * mm_.model().nb_surfaces(), 0 ) ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& cur_mesh = mm_.mesh( m ) ;
                for( index_t f = 0; f < cur_mesh.nb_facets(); f++ ) {
                    signed_index_t surface_id = cur_mesh.facet_region( f ) ;
                    if( surface2mesh_[ surface_id ] != m ) {continue ;}
                    index_t type_access = facet_access[ cur_mesh.facet_size( f ) ] ;
                    facets_[ mesh_facet_ptr_[ surface_id ]
                             + facet_ptr_[ ( NB_FACET_TYPES + 1 ) * surface_id
                                           + type_access ]
                             + cur_facet_index_type[ NB_FACET_TYPES * surface_id
                                                     + type_access ]++ ] = f ;
                }
            }

            // / 2 - fill cell/corner information
            index_t cell_access[ 4 ] = {
                0, 3, 2, 1
            } ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& mesh = mm_.mesh( m ) ;
                for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                    cell_ptr_[ ( NB_CELL_TYPES + 1 ) * m
                               + cell_access[ mesh.cell_type( c ) ] + 1 ]++ ;
                }
                for( index_t type = 1; type < NB_CELL_TYPES; type++ ) {
                    cell_ptr_[ ( NB_CELL_TYPES + 1 ) * m + type + 1 ] +=
                        cell_ptr_[ ( NB_CELL_TYPES + 1 ) * m + type ] ;
                }
            }

            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& mesh = mm_.mesh( m ) ;
                for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                    mesh_corner_ptr_[ m + 1 ] += mesh.cell_nb_vertices( c ) ;
                }
                mesh_cell_ptr_[ m + 1 ] += mesh_cell_ptr_[ m ] + mesh.nb_cells() ;
                mesh_corner_ptr_[ m + 1 ] += mesh_corner_ptr_[ m ] ;
            }
            cells_.resize( mesh_cell_ptr_.back() ) ;
            corners_.resize( mesh_corner_ptr_.back() ) ;

            std::vector< index_t > cur_cell_index_type(
                NB_CELL_TYPES * mm_.nb_meshes(), 0 ) ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                GEO::Mesh& mesh = const_cast< GEO::Mesh& >( mm_.mesh( m ) ) ;
                for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                    index_t type_access = cell_access[ mesh.cell_type( c ) ] ;
                    cells_[ mesh_cell_ptr_[ m ]
                            + cell_ptr_[ ( NB_CELL_TYPES + 1 ) * m + type_access ]
                            + cur_cell_index_type[ NB_CELL_TYPES * m +
                                                   type_access ]++ ] =
                        c ;
                    for( index_t v = 0; v < mesh.cell_nb_vertices( c ); v++ ) {
                        corners_[ mesh_corner_ptr_[ m ] +
                                  mesh.cell_vertices_begin( c )
                                  + v ] = mesh.cell_vertex_index( c, v ) ;
                    }
                }
            }

            // / 3 - fill vertex information
            first_duplicated_vertex_id_ = mm_.nb_vertices() ;

            // / 4 - update cached values
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
            // / 1 - Get all the corner vertices (redondant information)
            std::vector< vec3 > corner_vertices( corners_.size() ) ;
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& mesh = mm_.mesh( m ) ;
                index_t mesh_start = mesh_corner_ptr_[ m ] ;
                for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                    for( index_t v = 0; v < mesh.cell_nb_vertices( c ); v++ ) {
                        corner_vertices[ mesh_start + mesh.cell_vertices_begin( c )
                                         + v ] = GEO::Geom::mesh_vertex( mesh,
                            mesh.cell_vertex_index( c, v ) ) ;
                    }
                }
            }

            // / 2 - Tag all corners to duplicate
            const BoundaryModel& model = mm_.model() ;
            std::vector< SurfaceAction > surface_actions( model.nb_surfaces(),
                                                          SKIP ) ;
            std::vector< bool > corner_to_duplicate( corner_vertices.size(),
                                                     false ) ;
            {
                ColocaterANN ann( corner_vertices ) ;
                for( index_t s = 0; s < model.nb_surfaces(); s++ ) {
                    if( !is_surface_to_duplicate( s, mode ) ) {continue ;}
                    surface_actions[ s ] = TO_PROCESS ;
                    const Surface& surface = model.surface( s ) ;
                    for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
                        std::vector< index_t > colocated_corners ;
                        ann.get_colocated( surface.vertex( v ), colocated_corners ) ;
                        for( index_t co = 0; co < colocated_corners.size(); co++ ) {
                            corner_to_duplicate[ colocated_corners[ co ] ] = true ;
                        }
                    }
                }
            }
            corner_vertices.clear() ;

            // / 3 - Duplicate the corners
            for( index_t m = 0; m < mm_.nb_meshes(); m++ ) {
                const GEO::Mesh& mesh = mm_.mesh( m ) ;
                ColocaterANN ann( mesh, ColocaterANN::FACETS ) ;
                for( index_t c = 0; c < mesh.nb_cells(); c++ ) {
                    for( index_t v = 0; v < mesh.cell_nb_vertices( c ); v++ ) {
                        index_t co = mesh_corner_ptr_[ m ]
                                     + mesh.cell_vertices_begin( c ) + v ;
                        if( !corner_to_duplicate[ co ] ) {continue ;}

                        index_t vertex_id = mesh.cell_vertex_index( c, v ) ;
                        std::vector< index_t > corner_used ;
                        std::vector< index_t > cell_added ;
                        std::set< surface_side > surfaces ;
                        std::stack< index_t > S ;
                        S.push( c ) ;
                        cell_added.push_back( c ) ;
                        do {
                            index_t cur_c = S.top() ;
                            S.pop() ;
                            for( index_t cur_v = 0;
                                 cur_v < mesh.cell_nb_vertices( cur_c ); cur_v++ )
                            {
                                if( mesh.cell_vertex_index( cur_c, cur_v )
                                    == vertex_id )
                                {
                                    index_t cur_co = mesh_corner_ptr_[ m ]
                                                     + mesh.cell_vertices_begin(
                                        cur_c ) + cur_v ;
                                    corner_to_duplicate[ cur_co ] = false ;
                                    corner_used.push_back( cur_co ) ;
                                    break ;
                                }
                            }
                            for( index_t cur_f = 0;
                                 cur_f < mesh.cell_nb_facets( cur_c ); cur_f++ )
                            {
                                for( index_t cur_v = 0;
                                     cur_v
                                     < mesh.cell_facet_nb_vertices( cur_c, cur_f );
                                     cur_v++ )
                                {
                                    if( mesh.cell_facet_vertex_index( cur_c, cur_f,
                                            cur_v ) == vertex_id )
                                    {
                                        std::vector< index_t > result ;
                                        if( ann.get_colocated(
                                                Utils::mesh_cell_facet_center( mesh,
                                                    cur_c, cur_f ), result ) )
                                        {
                                            index_t surface_id = mesh.facet_region(
                                                result[ 0 ] ) ;
                                            vec3 facet_normal =
                                                GEO::Geom::mesh_facet_normal( mesh,
                                                    result[ 0 ] ) ;
                                            vec3 cell_facet_normal =
                                                Utils::mesh_cell_facet_normal( mesh,
                                                    cur_c, cur_f ) ;
                                            SurfaceAction side = SurfaceAction(
                                                dot( facet_normal,
                                                    cell_facet_normal ) > 0 ) ;
                                            surfaces.insert(
                                                surface_side( surface_id, side ) ) ;
                                        } else {
                                            signed_index_t cur_adj =
                                                mesh.cell_adjacent( cur_c, cur_f ) ;
                                            if( cur_adj != - 1
                                                && !Utils::contains( cell_added,
                                                    index_t( cur_adj ) ) )
                                            {
                                                cell_added.push_back( cur_adj ) ;
                                                S.push( cur_adj ) ;
                                            }
                                        }
                                        break ;
                                    }
                                }
                            }
                        } while( !S.empty() ) ;

                        if( duplicate_corner( surfaces, surface_actions ) ) {
                            index_t duplicated_vertex_id =
                                first_duplicated_vertex_id_
                                + duplicated_vertex_indices_.size() ;
                            index_t global_vertex_id = mm_.global_vertex_id( m,
                                vertex_id ) ;
                            duplicated_vertex_indices_.push_back(
                                global_vertex_id ) ;
                            for( index_t cur_co = 0; cur_co < corner_used.size();
                                 cur_co++ )
                            {
                                corners_[ corner_used[ cur_co ] ] =
                                    duplicated_vertex_id ;
                            }
                        }
                    }
                }
            }
        }


        bool MacroMeshExport::duplicate_corner(
            const std::set< surface_side >& surfaces,
            std::vector< SurfaceAction >& info )
        {
            std::vector< SurfaceAction > temp_info( info.size(), TO_PROCESS ) ;
            for( std::set< surface_side >::const_iterator it( surfaces.begin() );
                 it != surfaces.end(); it++ )
            {
                index_t surface_id = it->first ;
                if( info[ surface_id ] == SKIP || temp_info[ surface_id ] ==
                    SKIP )
                {
                    continue ;
                }
                if( temp_info[ surface_id ] == TO_PROCESS ) {
                    temp_info[ surface_id ] = it->second ;
                } else {
                    if( temp_info[ surface_id ] != it->second ) {
                        // Free border
                        temp_info[ surface_id ] = SKIP ;
                    }
                }
            }

            bool result = false ;
            for( index_t s = 0; s < info.size(); s++ ) {
                if( temp_info[ s ] < 0 ) {continue ;}
                ringmesh_debug_assert( info[ s ] != SKIP ) ;
                if( info[ s ] == TO_PROCESS ) {
                    info[ s ] = SurfaceAction( !temp_info[ s ] ) ;
                } else {
                    if( info[ s ] == temp_info[ s ] ) {
                        result = true ;
                    }
                }
            }

            return result ;
        }


        bool MacroMeshExport::is_surface_to_duplicate(
            index_t s,
            const DuplicateMode& mode ) const
        {
            BoundaryModelElement::GEOL_FEATURE feature =
                mm_.model().surface( s ).geological_feature() ;
            if( mode == ALL && feature != BoundaryModelElement::VOI ) {return true ;}
            if( mode == FAULT && feature == BoundaryModelElement::FAULT ) {
                return true ;
            }

            return false ;
        }


        bool MacroMeshExport::vertex_id(
            index_t r,
            index_t corner,
            index_t& vertex_id,
            index_t& duplicated_vertex_id ) const
        {
            index_t corner_value = corners_[ mesh_corner_ptr_[ r ] + corner ] ;
            if( corner_value < first_duplicated_vertex_id_ ) {
                vertex_id = mm_.global_vertex_id( r, corner_value ) ;
                return true ;
            } else {
                duplicated_vertex_id = corner_value - first_duplicated_vertex_id_ ;
                vertex_id = duplicated_vertex_indices_[ duplicated_vertex_id ] ;
                return false ;
            }
        }
    }
}
