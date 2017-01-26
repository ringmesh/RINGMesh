/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

/*! \author Francois Bonneau */

#include <ringmesh/mesh/mesh.h>

#include <ringmesh/mesh/geogram_mesh.h>

namespace RINGMesh {
    MeshBase::~MeshBase()
    {
        if( vertices_nn_search_ ) delete vertices_nn_search_ ;
    }

    Mesh0D* Mesh0D::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMesh0D::type_name_static() ;
        }
        Mesh0D* mesh = Mesh0DFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "Mesh0D" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "Mesh0D" ) << "Falling back to GeogramMesh0D data structure"
                << std::endl ;

            mesh = new GeogramMesh0D ;
        }
        return mesh ;
    }

    Mesh1D* Mesh1D::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMesh1D::type_name_static() ;
        }
        Mesh1D* mesh = Mesh1DFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "Mesh1D" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "Mesh1D" ) << "Falling back to GeogramMesh1D data structure"
                << std::endl ;

            mesh = new GeogramMesh1D ;
        }
        return mesh ;
    }

    Mesh2D* Mesh2D::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMesh2D::type_name_static() ;
        }
        Mesh2D* mesh = Mesh2DFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "Mesh2D" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "Mesh2D" ) << "Falling back to GeogramMesh2D data structure"
                << std::endl ;

            mesh = new GeogramMesh2D ;
        }
        return mesh ;
    }

    Mesh3D* Mesh3D::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMesh3D::type_name_static() ;
        }
        Mesh3D* mesh = Mesh3DFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "Mesh3D" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "Mesh3D" ) << "Falling back to GeogramMesh3D data structure"
                << std::endl ;

            mesh = new GeogramMesh3D ;
        }
        return mesh ;
    }

    MeshAllD* MeshAllD::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMeshAllD::type_name_static() ;
        }
        MeshAllD* mesh = MeshAllDFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "MeshAllD" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "MeshAllD" ) << "Falling back to GeogramMeshAllD data structure"
                << std::endl ;

            mesh = new GeogramMeshAllD ;
        }
        return mesh ;
    }


} // namespace
