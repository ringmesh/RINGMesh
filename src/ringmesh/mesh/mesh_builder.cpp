/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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
#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_index.h>

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
        auto line = dynamic_cast< LineMesh< DIMENSION >* >( &mesh );
        if( line )
        {
            return create_line_mesh_builder( *line );
        }
        auto surface = dynamic_cast< SurfaceMesh< DIMENSION >* >( &mesh );
        if( surface )
        {
            return create_surface_mesh_builder( *surface );
        }
        return {};
    }
} // namespace

namespace RINGMesh
{
    template <>
    std::unique_ptr< MeshBaseBuilder< 2 > >
        mesh_api MeshBaseBuilder< 2 >::create_builder( MeshBase< 2 >& mesh )
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
        mesh_api MeshBaseBuilder< 3 >::create_builder( MeshBase< 3 >& mesh )
    {
        auto builder = create_pointset_builder( mesh );
        if( !builder )
        {
            auto volume = dynamic_cast< VolumeMesh< 3 >* >( &mesh );
            if( volume != nullptr )
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

    template < index_t DIMENSION >
    void MeshBaseBuilder< DIMENSION >::delete_vertex_nn_search()
    {
        mesh_base_.vertex_nn_search_.reset();
    }

    template < index_t DIMENSION >
    void MeshBaseBuilder< DIMENSION >::copy( const MeshBase< DIMENSION >& rhs, bool copy_attributes )
    {
        do_copy( rhs, copy_attributes );
        clear_vertex_linked_objects();
    }

    template < index_t DIMENSION >
    void MeshBaseBuilder< DIMENSION >::clear( bool keep_attributes, bool keep_memory )
    {
        do_clear( keep_attributes, keep_memory );
        clear_vertex_linked_objects();
    }

    template < index_t DIMENSION >
    void MeshBaseBuilder< DIMENSION >::set_vertex( index_t v_id, const vecn< DIMENSION >& vertex )
    {
        do_set_vertex( v_id, vertex );
        clear_vertex_linked_objects();
    }

    template < index_t DIMENSION >
    index_t MeshBaseBuilder< DIMENSION >::create_vertex()
    {
        index_t index = do_create_vertex();
        clear_vertex_linked_objects();
        return index;
    }

    template < index_t DIMENSION >
    index_t MeshBaseBuilder< DIMENSION >::create_vertex( const vecn< DIMENSION >& vertex )
    {
        index_t index = create_vertex();
        set_vertex( index, vertex );
        return index;
    }

    template < index_t DIMENSION >
    index_t MeshBaseBuilder< DIMENSION >::create_vertices( index_t nb )
    {
        index_t index = do_create_vertices( nb );
        clear_vertex_linked_objects();
        return index;
    }

    template < index_t DIMENSION >
    void MeshBaseBuilder< DIMENSION >::assign_vertices( const std::vector< double >& point_coordinates )
    {
        do_assign_vertices( point_coordinates );
        clear_vertex_linked_objects();
    }

    template < index_t DIMENSION >
    void MeshBaseBuilder< DIMENSION >::delete_vertices( const std::vector< bool >& to_delete )
    {
        do_delete_vertices( to_delete );
        clear_vertex_linked_objects();
    }

    template < index_t DIMENSION >
    void MeshBaseBuilder< DIMENSION >::clear_vertices( bool keep_attributes, bool keep_memory )
    {
        do_clear_vertices( keep_attributes, keep_memory );
        clear_vertex_linked_objects();
    }

    template < index_t DIMENSION >
    void MeshBaseBuilder< DIMENSION >::permute_vertices( const std::vector< index_t >& permutation )
    {
        do_permute_vertices( permutation );
        clear_vertex_linked_objects();
    }

    template < index_t DIMENSION >
    void LineMeshBuilder< DIMENSION >::remove_isolated_vertices()
    {
        std::vector< bool > to_delete( line_mesh_.nb_vertices(), true );
        for( auto e : range( line_mesh_.nb_edges() ) )
        {
            for( auto v : range( 2 ) )
            {
                auto vertex_id = line_mesh_.edge_vertex( { e, v } );
                to_delete[vertex_id] = false;
            }
        }
        this->delete_vertices( to_delete );
    }

    template < index_t DIMENSION >
    void LineMeshBuilder< DIMENSION >::create_edge( index_t v1_id, index_t v2_id )
    {
        do_create_edge( v1_id, v2_id );
        clear_edge_linked_objects();
    }

    template < index_t DIMENSION >
    index_t LineMeshBuilder< DIMENSION >::create_edges( index_t nb_edges )
    {
        index_t index = do_create_edges( nb_edges );
        clear_edge_linked_objects();
        return index;
    }

    template < index_t DIMENSION >
    void LineMeshBuilder< DIMENSION >::set_edge_vertex(
        const EdgeLocalVertex& edge_local_vertex, index_t vertex_id )
    {
        do_set_edge_vertex( edge_local_vertex, vertex_id );
        clear_edge_linked_objects();
    }

    template < index_t DIMENSION >
    void LineMeshBuilder< DIMENSION >::delete_edges( const std::vector< bool >& to_delete,
        bool remove_isolated_vertices )
    {
        do_delete_edges( to_delete );
        if( remove_isolated_vertices )
        {
            this->remove_isolated_vertices();
        }
        clear_edge_linked_objects();
    }

    template < index_t DIMENSION >
    void LineMeshBuilder< DIMENSION >::clear_edges( bool keep_attributes, bool keep_memory )
    {
        do_clear_edges( keep_attributes, keep_memory );
        clear_edge_linked_objects();
    }

    template < index_t DIMENSION >
    void LineMeshBuilder< DIMENSION >::permute_edges( const std::vector< index_t >& permutation )
    {
        do_permute_edges( permutation );
        clear_edge_linked_objects();
    }

    template < index_t DIMENSION >
    void SurfaceMeshBuilder< DIMENSION >::remove_isolated_vertices()
    {
        std::vector< bool > to_delete( surface_mesh_.nb_vertices(), true );
        for( auto p : range( surface_mesh_.nb_polygons() ) )
        {
            for( auto v : range( surface_mesh_.nb_polygon_vertices( p ) ) )
            {
                auto vertex_id = surface_mesh_.polygon_vertex( { p, v } );
                to_delete[vertex_id] = false;
            }
        }
        this->delete_vertices( to_delete );
    }

    template < index_t DIMENSION >
    void VolumeMeshBuilder< DIMENSION >::remove_isolated_vertices()
    {
        std::vector< bool > to_delete( volume_mesh_.nb_vertices(), true );
        for( auto c : range( volume_mesh_.nb_cells() ) )
        {
            for( auto v : range( volume_mesh_.nb_cell_vertices( c ) ) )
            {
                auto vertex_id = volume_mesh_.cell_vertex( { c, v } );
                to_delete[vertex_id] = false;
            }
        }
        this->delete_vertices( to_delete );
    }

    template < index_t DIMENSION >
    void VolumeMeshBuilder< DIMENSION >::delete_cell_nn_search()
    {
        volume_mesh_.cell_nn_search_.reset();
        volume_mesh_.cell_facet_nn_search_.reset();
    }

    template < index_t DIMENSION >
    void VolumeMeshBuilder< DIMENSION >::delete_cell_aabb()
    {
        volume_mesh_.cell_aabb_.reset();
    }

    template class mesh_api MeshBaseBuilder< 2 >;
    template class mesh_api PointSetMeshBuilder< 2 >;
    template class mesh_api LineMeshBuilder< 2 >;
    template class mesh_api SurfaceMeshBuilder< 2 >;

    template class mesh_api MeshBaseBuilder< 3 >;
    template class mesh_api PointSetMeshBuilder< 3 >;
    template class mesh_api LineMeshBuilder< 3 >;
    template class mesh_api SurfaceMeshBuilder< 3 >;
    template class mesh_api VolumeMeshBuilder< 3 >;
} // namespace RINGMesh
