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
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_IO__
#define __RINGMESH_IO__


#include <ringmesh/common.h>

#include <geogram/basic/factory.h>

EXPIMP_TEMPLATE template class RINGMESH_API std::basic_string<char>;

namespace RINGMesh {
    class BoundaryModel ;
    class MacroMesh ;
    class WellGroup ;
}

namespace RINGMesh {
    namespace RINGMeshIO {
        //    ___                   _               __  __         _     _
        //   | _ ) ___ _  _ _ _  __| |__ _ _ _ _  _|  \/  |___  __| |___| |
        //   | _ \/ _ \ || | ' \/ _` / _` | '_| || | |\/| / _ \/ _` / -_) |
        //   |___/\___/\_,_|_||_\__,_\__,_|_|  \_, |_|  |_\___/\__,_\___|_|
        //                                     |__/

        bool RINGMESH_API load(
            const std::string& filename,
            BoundaryModel& model ) ;

        bool RINGMESH_API save(
            BoundaryModel& model,
            const std::string& filename ) ;

        //    __  __                 __  __        _
        //   |  \/  |__ _ __ _ _ ___|  \/  |___ __| |_
        //   | |\/| / _` / _| '_/ _ \ |\/| / -_|_-< ' \
        //   |_|  |_\__,_\__|_| \___/_|  |_\___/__/_||_|
        //

        bool RINGMESH_API load(
            const std::string& mesh_file,
            MacroMesh& mm ) ;

        bool RINGMESH_API save(
            const MacroMesh& mm,
            const std::string& filename ) ;

        //   __      __   _ _  ___
        //   \ \    / /__| | |/ __|_ _ ___ _  _ _ __
        //    \ \/\/ / -_) | | (_ | '_/ _ \ || | '_ \
        //     \_/\_/\___|_|_|\___|_| \___/\_,_| .__/
        //                                     |_|

        bool RINGMESH_API load(
            const std::string& mesh_file,
            WellGroup& wells ) ;


        class RINGMESH_API BoundaryModelIOHandler: public GEO::Counted {
        public:
            static BoundaryModelIOHandler* create( const std::string& format ) ;

            static BoundaryModelIOHandler* get_handler(
                const std::string& filename ) ;

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
        typedef GEO::Factory0< BoundaryModelIOHandler > BoundaryModelIOHandlerFactory ;
#define ringmesh_register_BoundaryModelIOHandler_creator( type, name ) \
    geo_register_creator( BoundaryModelIOHandlerFactory, type, name )

        class RINGMESH_API MacroMeshIOHandler: public GEO::Counted {
        public:
            static MacroMeshIOHandler* create( const std::string& format ) ;

            static MacroMeshIOHandler* get_handler( const std::string& filename ) ;

            virtual bool load( const std::string& filename, MacroMesh& mesh ) = 0 ;

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
#define ringmesh_register_MacroMeshIOHandler_creator( type, name ) \
    geo_register_creator( MacroMeshIOHandlerFactory, type, name )

        class RINGMESH_API WellGroupIOHandler: public GEO::Counted {
        public:
            static WellGroupIOHandler* create( const std::string& format ) ;

            static WellGroupIOHandler* get_handler( const std::string& filename ) ;

            virtual bool load( const std::string& filename, WellGroup& mesh ) = 0 ;

            virtual bool save(
                const WellGroup& mesh,
                const std::string& filename ) = 0 ;

        protected:
            WellGroupIOHandler()
            {
            }

            virtual ~WellGroupIOHandler()
            {
            }
        } ;

        typedef GEO::SmartPointer< WellGroupIOHandler > WellGroupIOHandler_var ;
        typedef GEO::Factory0< WellGroupIOHandler > WellGroupIOHandlerFactory ;
#define ringmesh_register_WellGroupIOHandler_creator( type, name ) \
    geo_register_creator( WellGroupIOHandlerFactory, type, name )

    }
}
#endif
