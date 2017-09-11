/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include <ringmesh/mesh/mesh_builder.h>

#include <ringmesh/basic/nn_search.h>

#include <ringmesh/mesh/aabb.h>
#include <ringmesh/mesh/geogram_mesh_builder.h>

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    std::unique_ptr< PointSetMeshBuilder< DIMENSION > >
        create_point_mesh_builder( PointSetMesh< DIMENSION >& mesh )
    {
        return PointSetMeshBuilderFactory< DIMENSION >::create(
            mesh.type_name(), mesh );
    }

    template < index_t DIMENSION >
    std::unique_ptr< LineMeshBuilder< DIMENSION > > create_line_mesh_builder(
        LineMesh< DIMENSION >& mesh )
    {
        return LineMeshBuilderFactory< DIMENSION >::create(
            mesh.type_name(), mesh );
    }

    template < index_t DIMENSION >
    std::unique_ptr< SurfaceMeshBuilder< DIMENSION > >
        create_surface_mesh_builder( SurfaceMesh< DIMENSION >& mesh )
    {
        return SurfaceMeshBuilderFactory< DIMENSION >::create(
            mesh.type_name(), mesh );
    }

    template < index_t DIMENSION >
    std::unique_ptr< VolumeMeshBuilder< DIMENSION > >
        create_volume_mesh_builder( VolumeMesh< DIMENSION >& mesh )
    {
        return VolumeMeshBuilderFactory< DIMENSION >::create(
            mesh.type_name(), mesh );
    }

    template < index_t DIMENSION >
    std::unique_ptr< MeshBaseBuilder< DIMENSION > > create_pointset_builder(
        MeshBase< DIMENSION >& mesh )
    {
        auto point_set = dynamic_cast< PointSetMesh< DIMENSION >* >( &mesh );
        if( point_set )
        {
            return create_point_mesh_builder( *point_set );
        }
        else
        {
            auto line = dynamic_cast< LineMesh< DIMENSION >* >( &mesh );
            if( line )
            {
                return create_line_mesh_builder( *line );
            }
            else
            {
                auto surface =
                    dynamic_cast< SurfaceMesh< DIMENSION >* >( &mesh );
                if( surface )
                {
                    return create_surface_mesh_builder( *surface );
                }
            }
        }
        return {};
    }
} // namespace

namespace RINGMesh
{
    template <>
    std::unique_ptr< MeshBaseBuilder< 2 > >
        RINGMESH_API MeshBaseBuilder< 2 >::create_builder( MeshBase< 2 >& mesh )
    {
        auto builder = create_pointset_builder( mesh );
        if( !builder )
        {
            throw RINGMeshException( "MeshBaseBuilder",
                "Could not create mesh builder of data structure: ",
                mesh.type_name() );
        }
        return builder;
    }

    template <>
    std::unique_ptr< MeshBaseBuilder< 3 > >
        RINGMESH_API MeshBaseBuilder< 3 >::create_builder( MeshBase< 3 >& mesh )
    {
        auto builder = create_pointset_builder( mesh );
        if( !builder )
        {
            auto volume = dynamic_cast< VolumeMesh< 3 >* >( &mesh );
            if( volume )
            {
                builder = create_volume_mesh_builder( *volume );
            }
        }
        if( !builder )
        {
            throw RINGMeshException( "MeshBaseBuilder",
                "Could not create mesh builder of data structure: ",
                mesh.type_name() );
        }
        return builder;
    }

    template < index_t DIMENSION >
    std::unique_ptr< PointSetMeshBuilder< DIMENSION > >
        PointSetMeshBuilder< DIMENSION >::create_builder(
            PointSetMesh< DIMENSION >& mesh )
    {
        auto builder = create_point_mesh_builder( mesh );
        if( !builder )
        {
            throw RINGMeshException( "PointSet",
                "Could not create mesh builder of data structure: ",
                mesh.type_name() );
        }
        return builder;
    }

    template < index_t DIMENSION >
    std::unique_ptr< LineMeshBuilder< DIMENSION > >
        LineMeshBuilder< DIMENSION >::create_builder(
            LineMesh< DIMENSION >& mesh )
    {
        auto builder = create_line_mesh_builder( mesh );
        if( !builder )
        {
            Logger::warn( "LineMeshBuilder",
                "Could not create mesh builder of data structure: ",
                mesh.type_name() );
        }
        return builder;
    }

    template < index_t DIMENSION >
    std::unique_ptr< SurfaceMeshBuilder< DIMENSION > >
        SurfaceMeshBuilder< DIMENSION >::create_builder(
            SurfaceMesh< DIMENSION >& mesh )
    {
        auto builder = create_surface_mesh_builder( mesh );
        if( !builder )
        {
            Logger::warn( "SurfaceMeshBuilder",
                "Could not create mesh builder of data structure: ",
                mesh.type_name() );
        }
        return builder;
    }

    template < index_t DIMENSION >
    std::unique_ptr< VolumeMeshBuilder< DIMENSION > >
        VolumeMeshBuilder< DIMENSION >::create_builder(
            VolumeMesh< DIMENSION >& mesh )
    {
        auto builder = create_volume_mesh_builder( mesh );
        if( !builder )
        {
            Logger::warn( "VolumeMeshBuilder",
                "Could not create mesh builder of data structure: ",
                mesh.type_name() );
        }
        return builder;
    }

    template std::unique_ptr< PointSetMeshBuilder< 2 > > RINGMESH_API
        PointSetMeshBuilder< 2 >::create_builder( PointSetMesh< 2 >& );
    template std::unique_ptr< LineMeshBuilder< 2 > >
        RINGMESH_API LineMeshBuilder< 2 >::create_builder( LineMesh< 2 >& );
    template std::unique_ptr< SurfaceMeshBuilder< 2 > > RINGMESH_API
        SurfaceMeshBuilder< 2 >::create_builder( SurfaceMesh< 2 >& );

    template std::unique_ptr< PointSetMeshBuilder< 3 > > RINGMESH_API
        PointSetMeshBuilder< 3 >::create_builder( PointSetMesh< 3 >& );
    template std::unique_ptr< LineMeshBuilder< 3 > >
        RINGMESH_API LineMeshBuilder< 3 >::create_builder( LineMesh< 3 >& );
    template std::unique_ptr< SurfaceMeshBuilder< 3 > > RINGMESH_API
        SurfaceMeshBuilder< 3 >::create_builder( SurfaceMesh< 3 >& );
    template std::unique_ptr< VolumeMeshBuilder< 3 > >
        RINGMESH_API VolumeMeshBuilder< 3 >::create_builder( VolumeMesh< 3 >& );
} // namespace RINGMesh
