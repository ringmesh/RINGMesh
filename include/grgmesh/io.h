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
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
*/


#ifndef __GRGMESH_IO__
#define __GRGMESH_IO__

#include <grgmesh/common.h>
#include <grgmesh/macro_mesh.h>
#include <grgmesh/boundary_model.h>


#define MAX_FILENAME 512
#define READ_SIZE 8192
class string ;

namespace GRGMesh {
    class BoundaryModel ;
}

namespace GRGMesh {
    namespace GRGMeshIO {

        //    ___                   _               __  __         _     _
        //   | _ ) ___ _  _ _ _  __| |__ _ _ _ _  _|  \/  |___  __| |___| |
        //   | _ \/ _ \ || | ' \/ _` / _` | '_| || | |\/| / _ \/ _` / -_) |
        //   |___/\___/\_,_|_||_\__,_\__,_|_|  \_, |_|  |_\___/\__,_\___|_|
        //                                     |__/
        bool GRGMESH_API load( const std::string& filename, BoundaryModel& model ) ;
        bool GRGMESH_API save( BoundaryModel& model, const std::string& filename ) ;

        //    __  __                 __  __        _
        //   |  \/  |__ _ __ _ _ ___|  \/  |___ __| |_
        //   | |\/| / _` / _| '_/ _ \ |\/| / -_|_-< ' \
        //   |_|  |_\__,_\__|_| \___/_|  |_\___/__/_||_|
        //
        bool GRGMESH_API load( const std::string& mesh_file, MacroMesh& mm ) ;
        bool GRGMESH_API save( const MacroMesh& mm, const std::string& filename ) ;


        class GRGMESH_API BoundaryModelIOHandler: public GEO::Counted {
        public:
            static BoundaryModelIOHandler* create( const std::string& format ) ;
            static BoundaryModelIOHandler* get_handler( const std::string& filename ) ;

            virtual bool load(
                const std::string& filename,
                BoundaryModel& model ) = 0 ;
            virtual bool save(
                BoundaryModel& model,
                const std::string& filename ) = 0 ;

        protected:
            BoundaryModelIOHandler()
            {
            }
            virtual ~BoundaryModelIOHandler()
            {
            }
        } ;
        typedef GEO::SmartPointer< BoundaryModelIOHandler > BoundaryModelIOHandler_var ;
        typedef GEO::Factory0< BoundaryModelIOHandler > BoundaryModelIOHandlerFactory;
#define grgmesh_register_BoundaryModelIOHandler_creator(type, name) \
    geo_register_creator(BoundaryModelIOHandlerFactory, type, name)


        class GRGMESH_API MacroMeshIOHandler: public GEO::Counted {
        public:
            static MacroMeshIOHandler* create( const std::string& format ) ;
            static MacroMeshIOHandler* get_handler( const std::string& filename ) ;

            virtual bool load(
                const std::string& filename,
                MacroMesh& mesh ) = 0 ;
            virtual bool save(
                const MacroMesh& mesh,
                const std::string& filename ) = 0 ;

        protected:
            MacroMeshIOHandler()
            {
            }
            virtual ~MacroMeshIOHandler()
            {
            }
        } ;
        typedef GEO::SmartPointer< MacroMeshIOHandler > MacroMeshIOHandler_var ;
        typedef GEO::Factory0< MacroMeshIOHandler > MacroMeshIOHandlerFactory ;
#define grgmesh_register_MacroMeshIOHandler_creator(type, name) \
    geo_register_creator(MacroMeshIOHandlerFactory, type, name)

        class GRGMESH_API MacroMeshExport {
        private:
            const static index_t NB_FACET_TYPES = 2 ;
            const static index_t NB_CELL_TYPES = 4 ;

        public:
            enum DuplicateMode {
                NONE = 0, FAULT = 1, HORIZON = 2
            } ;

            MacroMeshExport( const MacroMesh& mm ) ;

            void compute_database( const DuplicateMode& mode = NONE ) ;

            index_t nb_triangle() const { return nb_triangle_ ; }
            index_t nb_triangle( index_t r ) const {
                return facet_ptr_[(NB_FACET_TYPES+1) * r + 1] - facet_ptr_[(NB_FACET_TYPES+1) * r] ;
            }
            index_t nb_quad() const { return nb_quad_ ; }
            index_t nb_quad( index_t r ) const {
                return facet_ptr_[(NB_FACET_TYPES+1) * r + 2] - facet_ptr_[(NB_FACET_TYPES+1) * r + 1] ;
            }
            index_t nb_facets() const { return facets_.size() ; }
            index_t nb_facets( index_t r ) const {
                return mesh_facet_end( r ) -  mesh_facet_begin( r ) ;
            }

            index_t nb_tet() const { return nb_tet_ ; }
            index_t nb_tet( index_t r ) const {
                return cell_ptr_[(NB_CELL_TYPES+1) * r + 1] - cell_ptr_[(NB_CELL_TYPES+1) * r] ;
            }
            index_t nb_pyramid() const { return nb_pyramid_ ; }
            index_t nb_pyramid( index_t r ) const {
                return cell_ptr_[(NB_CELL_TYPES+1) * r + 2] - cell_ptr_[(NB_CELL_TYPES+1) * r + 1] ;
            }
            index_t nb_prism() const { return nb_prism_ ; }
            index_t nb_prism( index_t r ) const {
                return cell_ptr_[(NB_CELL_TYPES+1) * r + 3] - cell_ptr_[(NB_CELL_TYPES+1) * r + 2] ;
            }
            index_t nb_hex() const { return nb_hex_ ; }
            index_t nb_hex( index_t r ) const {
                return cell_ptr_[(NB_CELL_TYPES+1) * r + 4] - cell_ptr_[(NB_CELL_TYPES+1) * r + 3] ;
            }
            index_t nb_cells() const { return cells_.size() ; }
            index_t nb_cells( index_t r ) const {
                return mesh_cell_end( r ) -  mesh_cell_begin( r ) ;
            }


            index_t local_triangle_id( index_t r, index_t t ) const {
                return facet( mesh_facet_begin( r ) + facet_ptr_[(NB_FACET_TYPES+1)* r] + t ) ;
            }
            index_t local_quad_id( index_t r, index_t q ) const {
                return facet( mesh_facet_begin( r ) + facet_ptr_[(NB_FACET_TYPES+1) * r + 1] + q ) ;
            }

            index_t local_tet_id( index_t r, index_t t ) const {
                return cell( mesh_facet_begin( r ) + cell_ptr_[(NB_CELL_TYPES+1) * r] + t ) ;
            }
            index_t local_pyramid_id( index_t r, index_t p ) const {
                return cell( mesh_facet_begin( r ) + cell_ptr_[(NB_CELL_TYPES+1) * r + 1] + p ) ;
            }
            index_t local_prism_id( index_t r, index_t p ) const {
                return cell( mesh_facet_begin( r ) + cell_ptr_[(NB_CELL_TYPES+1) * r + 2] + p ) ;
            }
            index_t local_hex_id( index_t r, index_t h ) const {
                return cell( mesh_facet_begin( r ) + cell_ptr_[(NB_CELL_TYPES+1) * r + 3] + h ) ;
            }

        private:
            void fill_with_geometry() ;
            void duplicate_vertices( const DuplicateMode& mode ) ;
            bool is_surface_to_duplicate(
                index_t s,
                const DuplicateMode& mode ) const ;

            index_t mesh_facet_begin( index_t r ) const { return mesh_facet_ptr_[r] ; }
            index_t mesh_facet_end( index_t r ) const { return mesh_facet_ptr_[r+1] ; }
            index_t facet( index_t f ) const { return facets_[f] ; }

            index_t mesh_cell_begin( index_t r ) const { return mesh_cell_ptr_[r] ; }
            index_t mesh_cell_end( index_t r ) const { return mesh_cell_ptr_[r+1] ; }
            index_t cell( index_t c ) const { return cells_[c] ; }

        private:
            const MacroMesh& mm_ ;

            std::vector< index_t > facets_ ;
            // [TRGL/QUAD]
            std::vector< index_t > facet_ptr_ ;
            std::vector< index_t > mesh_facet_ptr_ ;
            std::vector< index_t > surface2mesh_ ;

            std::vector< index_t > cells_ ;
            // [TET/PY/PRISM/HEX]
            std::vector< index_t > cell_ptr_ ;
            std::vector< index_t > mesh_cell_ptr_ ;

            std::vector< index_t > corners_ ;
            std::vector< index_t > mesh_corner_ptr_ ;

            std::vector< index_t > duplicated_vertex_indices_ ;
            index_t first_duplicated_vertex_id_ ;

            index_t nb_triangle_ ;
            index_t nb_quad_ ;
            index_t nb_tet_ ;
            index_t nb_pyramid_ ;
            index_t nb_prism_ ;
            index_t nb_hex_ ;
        } ;
    }
}
#endif
