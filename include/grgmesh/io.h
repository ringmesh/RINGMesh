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





        private:
            const MacroMesh& mm_ ;

            std::vector< index_t > facets_ ;
            std::vector< index_t > facet_ptr_ ;
            std::vector< index_t > mesh_facet_ptr_ ;

            std::vector< index_t > cells_ ;
            std::vector< index_t > cell_ptr_ ;
            std::vector< index_t > mesh_cell_ptr_ ;

            std::vector< index_t > cell_corners_ ;
        } ;
    }
}
#endif
