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

#include <ringmesh/geomodel/geomodel_mesh.h>

#include <numeric>
#include <stack>

#include <geogram/basic/algorithm.h>

#include <geogram/mesh/mesh_geometry.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_builder.h>

#include <ringmesh/mesh/well.h>
#include <ringmesh/mesh/geogram_mesh.h>

/*!
 * @author Arnaud Botella - Jeanne Pellerin - Antoine Mazuyer
 */

namespace {
    using namespace RINGMesh;

    template< index_t DIMENSION >
    class GeoModelMeshPolygonsBaseSort {
    public:
        GeoModelMeshPolygonsBaseSort(
            const SurfaceMesh< DIMENSION >& mesh,
            const GEO::Attribute< index_t >& surface_id )
            : mesh_( mesh ), surface_id_( surface_id )
        {
        }

        bool operator()( index_t i, index_t j ) const
        {
            if( surface_id_[i] != surface_id_[j] ) {
                return surface_id_[i] < surface_id_[j];
            } else {
                return mesh_.nb_polygon_vertices( i )
                    < mesh_.nb_polygon_vertices( j );
            }
        }
    private:
        const SurfaceMesh< DIMENSION >& mesh_;
        const GEO::Attribute< index_t >& surface_id_;
    };

    template< index_t DIMENSION >
    class GeoModelMeshCellsSort {
    public:
        GeoModelMeshCellsSort(
            const VolumeMesh< DIMENSION >& mesh,
            const GEO::Attribute< index_t >& region_id )
            : mesh_( mesh ), region_id_( region_id )
        {
        }

        bool operator()( index_t i, index_t j ) const
        {
            if( region_id_[i] != region_id_[j] ) {
                return region_id_[i] < region_id_[j];
            } else {
                return mesh_.cell_type( i ) < mesh_.cell_type( j );
            }
        }
    private:
        const VolumeMesh< DIMENSION >& mesh_;
        const GEO::Attribute< index_t >& region_id_;
    };

    std::string vertex_map_name()
    {
        return "model_vertex_map";
    }

    template< index_t DIMENSION >
    std::vector< index_t > cell_facets_around_vertex(
        const VolumeMesh< DIMENSION >& mesh,
        index_t cell,
        index_t vertex_id )
    {
        std::vector< index_t > facets;
        facets.reserve( mesh.nb_cell_facets( cell ) );
        for( index_t f : range( mesh.nb_cell_facets( cell ) ) ) {
            for( index_t v : range( mesh.nb_cell_facet_vertices( cell, f ) ) ) {
                if( mesh.cell_facet_vertex( cell, f, v ) == vertex_id ) {
                    facets.push_back( f );
                    break;
                }
            }
        }
        return facets;
    }

    template< index_t DIMENSION >
    void copy_vertices(
        MeshBaseBuilder< DIMENSION >* builder,
        const MeshBase< DIMENSION >& mesh )
    {
        builder->clear_vertices( true, false );
        builder->create_vertices( mesh.nb_vertices() );
        for( index_t v : range( mesh.nb_vertices() ) ) {
            builder->set_vertex( v, mesh.vertex( v ) );
        }
    }
}

namespace RINGMesh {

    template< index_t DIMENSION >
    GeoModelMeshCommon< DIMENSION >::GeoModelMeshCommon(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm )
        : gmm_( gmm ), gm_( gm ), mesh_base_( nullptr )
    {
    }

    template< index_t DIMENSION >
    GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::GeoModelVertexMapper(
        GeoModelMeshVerticesBase& geomodel_vertices,
        const GeoModel< DIMENSION >& geomodel )
        : geomodel_vertices_( geomodel_vertices ), geomodel_( geomodel )
    {
        vertex_maps_[Corner< DIMENSION >::type_name_static()] = &corner_vertex_maps_;
        vertex_maps_[Line< DIMENSION >::type_name_static()] = &line_vertex_maps_;
        vertex_maps_[Surface< DIMENSION >::type_name_static()] =
            &surface_vertex_maps_;
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::clear()
    {
        gme_vertices_.clear();
        clear_all_mesh_entity_vertex_map();
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::geomodel_vertex_index(
        const gmme_id& mesh_entity_id,
        index_t mesh_entity_vertex_index ) const
    {
        ringmesh_assert(
            mesh_entity_vertex_index
                < geomodel_.mesh_entity( mesh_entity_id ).nb_vertices() );

        return vertex_map( mesh_entity_id )[mesh_entity_vertex_index];
    }

    template< index_t DIMENSION >
    const std::vector< GMEVertex >& GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::mesh_entity_vertex_indices(
        index_t v ) const
    {
        ringmesh_assert( v < gme_vertices_.size() );
        return gme_vertices_[v];
    }

    template< index_t DIMENSION >
    std::vector< GMEVertex > GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::mesh_entity_vertex_indices(
        index_t v,
        const MeshEntityType& mesh_entity_type ) const
    {
        const std::vector< GMEVertex >& all_gmes = mesh_entity_vertex_indices( v );
        std::vector< GMEVertex > result;
        result.reserve( all_gmes.size() );
        for( const GMEVertex& vertex : all_gmes ) {
            if( vertex.gmme.type() == mesh_entity_type ) {
                result.push_back( vertex );
            }
        }
        return result;
    }

    template< index_t DIMENSION >
    std::vector< index_t > GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::mesh_entity_vertex_indices(
        index_t v,
        const gmme_id& mesh_entity_id ) const
    {
        std::vector< index_t > result;
        std::vector< GMEVertex > all_gmes = mesh_entity_vertex_indices( v );
        for( const GMEVertex& vertex : all_gmes ) {
            if( vertex.gmme == mesh_entity_id ) {
                result.push_back( vertex.v_index );
            }
        }
        return result;
    }

    template< index_t DIMENSION >
    const GEO::Attribute< index_t >&
    GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::vertex_map(
        const gmme_id& mesh_entity_id ) const
    {
        return ( *vertex_maps_.at( mesh_entity_id.type() ) )[mesh_entity_id.index()];
    }

    template< index_t DIMENSION >
    GEO::Attribute< index_t >& GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::vertex_map(
        const gmme_id& mesh_entity_id )
    {
        return ( *vertex_maps_[mesh_entity_id.type()] )[mesh_entity_id.index()];
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::set_vertex_map_value(
        const gmme_id& mesh_entity_id,
        index_t mesh_entity_vertex_index,
        index_t geomodel_entity_vertex_index )
    {
        test_and_initialize_mesh_entity_vertex_map( mesh_entity_id );
        vertex_map( mesh_entity_id )[mesh_entity_vertex_index] =
            geomodel_entity_vertex_index;
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::add_to_gme_vertices(
        const GMEVertex& gme_vertex,
        index_t geomodel_vertex_index )
    {
        gme_vertices_[geomodel_vertex_index].push_back( gme_vertex );
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::bind_all_mesh_entity_vertex_maps()
    {
        const std::vector< MeshEntityType >& all_mesh_entity_types =
            geomodel_.entity_type_manager().mesh_entity_manager.mesh_entity_types();
        for( const MeshEntityType& cur_entity_type : all_mesh_entity_types ) {
            index_t nb_cur_type_entities = geomodel_.nb_mesh_entities(
                cur_entity_type );
            vertex_maps_.at( cur_entity_type )->clear();
            vertex_maps_.at( cur_entity_type )->resize( nb_cur_type_entities );
            for( index_t e : range( nb_cur_type_entities ) ) {
                bind_vertex_map( { cur_entity_type, e } );
            }
        }
    }

    template< index_t DIMENSION >
    GEO::Attribute< index_t >&
    GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::bind_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        ringmesh_assert(
            mesh_entity_id.index() < vertex_maps_[mesh_entity_id.type()]->size() );
        if( geomodel_vertices_.is_initialized() ) {
            vertex_maps_.at( mesh_entity_id.type() )->bind_one_attribute(
                mesh_entity_id.index(),
                mesh_entity_vertex_attribute_manager( mesh_entity_id ),
                vertex_map_name() );
            vertex_map( mesh_entity_id ).fill( NO_ID );
        }
        return vertex_map( mesh_entity_id );
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::update_mesh_entity_maps_and_gmes(
        const std::vector< index_t >& old2new )
    {
        const std::vector< MeshEntityType >& all_mesh_entity_types =
            geomodel_.entity_type_manager().mesh_entity_manager.mesh_entity_types();
        for( const MeshEntityType& cur_entity_type : all_mesh_entity_types ) {
            for( index_t e : range( geomodel_.nb_mesh_entities( cur_entity_type ) ) ) {
                const GeoModelMeshEntity< DIMENSION >& E = geomodel_.mesh_entity(
                    cur_entity_type, e );
                gmme_id id = E.gmme();
                for( index_t v : range( E.nb_vertices() ) ) {
                    index_t old_m_id = geomodel_vertex_index( id, v );
                    index_t new_m_id = old2new[old_m_id];
                    set_vertex_map_value( id, v, new_m_id );

                    // Merge gme_vertices information
                    if( find( gme_vertices_[new_m_id], GMEVertex( id, v ) )
                        == NO_ID ) {
                        gme_vertices_[new_m_id].emplace_back( id, v );
                    }
                }
            }
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::unbind_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        resize_all_mesh_entity_vertex_maps( mesh_entity_id.type() );
        if( vertex_maps_.at( mesh_entity_id.type() )->is_attribute_bound(
            mesh_entity_id.index() ) ) {
            vertex_maps_.at( mesh_entity_id.type() )->unbind(
                mesh_entity_id.index() );
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::initialize_mesh_entity_vertex_map(
        const gmme_id& mesh_entity_id )
    {

        GEO::Attribute< index_t >& mesh_entity_vertex_map = bind_vertex_map(
            mesh_entity_id );

        const GeoModelMeshEntity< DIMENSION >& E = geomodel_.mesh_entity(
            mesh_entity_id );
        for( index_t v : range( E.nb_vertices() ) ) {
            mesh_entity_vertex_map[v] =
                geomodel_vertices_.nn_search().get_closest_neighbor( E.vertex( v ) );
        }
    }

    template< index_t DIMENSION >
    bool GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::test_and_initialize_mesh_entity_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        resize_all_mesh_entity_vertex_maps( mesh_entity_id.type() );
        if( !is_mesh_entity_vertex_map_initialized( mesh_entity_id ) ) {
            initialize_mesh_entity_vertex_map( mesh_entity_id );
            return false;
        }
        return true;
    }

    template< index_t DIMENSION >
    bool GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::is_mesh_entity_vertex_map_initialized(
        const gmme_id& mesh_entity_id ) const
    {
        return ( vertex_maps_.find( mesh_entity_id.type() )->second )->is_attribute_bound(
            mesh_entity_id.index() );
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::clear_all_mesh_entity_vertex_map()
    {
        for( auto& vertex_map : vertex_maps_ ) {
            for( index_t e : range( vertex_map.second->size() ) ) {
                vertex_map.second->unbind( e );
            }
            vertex_map.second->clear();
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::resize_all_mesh_entity_vertex_maps(
        const MeshEntityType& type )
    {
        vertex_maps_.at( type )->resize( geomodel_.nb_mesh_entities( type ),
            nullptr );
    }

    template< index_t DIMENSION >
    GEO::AttributesManager& GeoModelMeshVerticesBase< DIMENSION >::GeoModelVertexMapper::mesh_entity_vertex_attribute_manager(
        const gmme_id& mesh_entity_id ) const
    {
        const GeoModelMeshEntity< DIMENSION >& mesh_entity = geomodel_.mesh_entity(
            mesh_entity_id );
        return mesh_entity.vertex_attribute_manager();
    }

    template< index_t DIMENSION >
    GeoModelMeshVerticesBase< DIMENSION >::GeoModelMeshVerticesBase(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< PointSetMesh< DIMENSION > >& mesh )
        :
            GeoModelMeshCommon< DIMENSION >( gmm, gm ),
            mesh_( mesh ),
            vertex_mapper_( *this, gmm.geomodel() )
    {
        this->set_mesh( mesh_.get() );
    }

    template< index_t DIMENSION >
    bool GeoModelMeshVerticesBase< DIMENSION >::is_initialized() const
    {
        return mesh_->nb_vertices() > 0;
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshVerticesBase* >( this )->initialize();
        }
    }

    template< index_t DIMENSION >
    index_t nb_entity_vertices(
        const GeoModel< DIMENSION >& geomodel,
        const MeshEntityType& entity_type )
    {
        index_t count = 0;
        for( index_t i : range( geomodel.nb_mesh_entities( entity_type ) ) ) {
            count += geomodel.mesh_entity( entity_type, i ).nb_vertices();
        }
        return count;
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::fill_vertices_for_entity_type(
        const GeoModel< DIMENSION >& geomodel,
        const MeshEntityType& entity_type,
        index_t& count )
    {
        std::unique_ptr< PointSetMeshBuilder< DIMENSION > > mesh_builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        for( index_t i : range( geomodel.nb_mesh_entities( entity_type ) ) ) {
            GeoModelMeshEntity< DIMENSION >& E = const_cast< GeoModelMeshEntity<
                DIMENSION >& >( geomodel.mesh_entity( entity_type, i ) );
            if( E.nb_vertices() == 0 ) {
                continue;
            }

            // Map and vertex
            gmme_id id = E.gmme();
            for( index_t v : range( E.nb_vertices() ) ) {
                index_t local_count = count + v;
                mesh_builder->set_vertex( local_count, E.vertex( v ) );
                // Map from vertices of MeshEntities to GeoModelMeshVerticesBase
                vertex_mapper_.set_vertex_map_value( id, v, local_count );
                vertex_mapper_.add_to_gme_vertices( GMEVertex( id, v ),
                    local_count );
            }
            // Global vertex index increment
            count += E.nb_vertices();
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::nb_total_vertices() const
    {
        index_t nb = 0;
        nb += nb_entity_vertices( this->gm_,
            Corner< DIMENSION >::type_name_static() );
        nb += nb_entity_vertices( this->gm_, Line< DIMENSION >::type_name_static() );
        nb += nb_entity_vertices( this->gm_,
            Surface< DIMENSION >::type_name_static() );
        return nb;
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::initialize()
    {
        clear();

        // Total number of vertices in the
        // Corners, Lines, Surfaces and Regions of the GeoModel
        index_t nb = nb_total_vertices();

        // Get out if no vertices
        if( nb == 0 ) {
            return;
        }

        // Fill the vertices
        std::unique_ptr< PointSetMeshBuilder< DIMENSION > > builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        builder->create_vertices( nb );
        vertex_mapper_.clear_and_resize_geomodel_vertex_gmes( nb );
        vertex_mapper_.bind_all_mesh_entity_vertex_maps();

        fill_vertices();

        // Remove colocated vertices
        remove_colocated();
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::fill_vertices()
    {
        index_t count = 0;
        fill_vertices_for_entity_type( this->gm_,
            Corner< DIMENSION >::type_name_static(), count );
        fill_vertices_for_entity_type( this->gm_,
            Line< DIMENSION >::type_name_static(), count );
        fill_vertices_for_entity_type( this->gm_,
            Surface< DIMENSION >::type_name_static(), count );
        return count;
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::clear()
    {
        this->gmm_.polygons.clear();
        this->gmm_.edges.clear();
        vertex_mapper_.clear();

        std::unique_ptr< PointSetMeshBuilder< DIMENSION > > builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        builder->clear( true, false );
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::unbind_geomodel_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        vertex_mapper_.unbind_vertex_map( mesh_entity_id );
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::bind_geomodel_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        vertex_mapper_.bind_vertex_map( mesh_entity_id );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::nb() const
    {
        test_and_initialize();
        return mesh_->nb_vertices();
    }

    template< index_t DIMENSION >
    const vecn< DIMENSION >& GeoModelMeshVerticesBase< DIMENSION >::vertex(
        index_t v ) const
    {
        test_and_initialize();
        ringmesh_assert( v < nb() );
        return mesh_->vertex( v );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::index(
        const vecn< DIMENSION >& p ) const
    {
        test_and_initialize();
        const NNSearch< DIMENSION >& colocator = mesh_->vertex_nn_search();
        std::vector< index_t > vertices = colocator.get_neighbors( p,
            this->gm_.epsilon() );
        if( vertices.empty() ) {
            return NO_ID;
        } else {
            ringmesh_assert( vertices.size() == 1 );
            return vertices[0];
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::geomodel_vertex_id(
        const gmme_id& mesh_entity,
        index_t entity_vertex_index ) const
    {
        test_and_initialize();
        return vertex_mapper_.geomodel_vertex_index( mesh_entity,
            entity_vertex_index );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::geomodel_vertex_id(
        const gmme_id& mesh_entity,
        index_t entity_mesh_element_index,
        index_t vertex_local_index ) const
    {
        index_t entity_vertex_index =
            this->gm_.mesh_entity( mesh_entity ).mesh_element_vertex_index(
                entity_mesh_element_index, vertex_local_index );
        return geomodel_vertex_id( mesh_entity, entity_vertex_index );
    }

    template< index_t DIMENSION >
    std::vector< index_t > GeoModelMeshVerticesBase< DIMENSION >::mesh_entity_vertex_id(
        const gmme_id& mesh_entity,
        index_t geomodel_vertex_index ) const
    {
        test_and_initialize();
        ringmesh_assert( geomodel_vertex_index < nb() );
        return vertex_mapper_.mesh_entity_vertex_indices( geomodel_vertex_index,
            mesh_entity );
    }

    template< index_t DIMENSION >
    const std::vector< GMEVertex >& GeoModelMeshVerticesBase< DIMENSION >::gme_vertices(
        index_t v ) const
    {
        test_and_initialize();
        return vertex_mapper_.mesh_entity_vertex_indices( v );
    }

    template< index_t DIMENSION >
    std::vector< GMEVertex > GeoModelMeshVerticesBase< DIMENSION >::gme_type_vertices(
        const MeshEntityType& entity_type,
        index_t v ) const
    {
        test_and_initialize();
        return vertex_mapper_.mesh_entity_vertex_indices( v, entity_type );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::add_vertex(
        const vecn< DIMENSION >& point )
    {
        std::unique_ptr< PointSetMeshBuilder< DIMENSION > > builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        const index_t index = builder->create_vertex( point );
        vertex_mapper_.resize_geomodel_vertex_gmes( nb() );
        return index;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::add_vertices(
        const std::vector< vecn< DIMENSION > >& points )
    {
        ringmesh_assert( !points.empty() );
        std::unique_ptr< PointSetMeshBuilder< DIMENSION > > builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        const index_t start_index = builder->create_vertex( points[0] );
        for( size_t i = 1; i < points.size(); ++i ) {
            builder->create_vertex( points[i] );
        }
        vertex_mapper_.resize_geomodel_vertex_gmes( nb() );
        return start_index;
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::update_point(
        index_t v,
        const vecn< DIMENSION >& point )
    {
        test_and_initialize();
        ringmesh_assert( v < nb() );
        // Change the position of the unique_vertex
        std::unique_ptr< PointSetMeshBuilder< DIMENSION > > mesh_builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        mesh_builder->set_vertex( v, point );

        GeoModelBuilder< DIMENSION > builder( this->gm_ );

        const std::vector< GMEVertex >& gme_v = gme_vertices( v );
        for( const GMEVertex& info : gme_v ) {
            builder.geometry.set_mesh_entity_vertex( info.gmme, info.v_index, point,
                false );
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::update_vertex_mapping(
        const gmme_id& entity_id,
        index_t entity_vertex_index,
        index_t geomodel_vertex_index )
    {
        vertex_mapper_.set_vertex_map_value( entity_id, entity_vertex_index,
            geomodel_vertex_index );
        vertex_mapper_.add_to_gme_vertices(
            GMEVertex( entity_id, entity_vertex_index ), geomodel_vertex_index );
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::remove_colocated()
    {
        // Get out if nothing to do
        // and compute the points if they are not initialized yet
        if( nb() == 0 ) {
            return;
        }
        // Identify and invalidate colocated vertices
        std::vector< index_t > old2new;
        index_t nb_colocalised_vertices =
            mesh_->vertex_nn_search().get_colocated_index_mapping(
                this->gm_.epsilon(), old2new );
        if( nb_colocalised_vertices > 0 ) {
            erase_vertices( old2new );
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::erase_vertices(
        std::vector< index_t >& to_delete )
    {
        ringmesh_assert( to_delete.size() == nb() );

        // For mesh vertices deletion
        std::vector< bool > to_delete_bool( nb(), false );

        // Fill the delete information for geogram
        // Recycle the to_delete vertex to get the mapping between
        // new and old points. This is implemented to be the same
        // as what is done in the delete_elements function in geogram
        index_t nb_todelete = 0;
        index_t cur = 0;
        for( index_t v : range( nb() ) ) {
            if( to_delete[v] != v ) {
                to_delete_bool[v] = true;
                nb_todelete++;
                if( to_delete[v] != NO_ID ) {
                    ringmesh_assert( to_delete[v] < v );
                    to_delete[v] = to_delete[to_delete[v]];
                }
            } else {
                to_delete[v] = cur;
                ++cur;
            }
        }
        if( nb_todelete == 0 ) {
            return;
        }
        if( nb_todelete == nb() ) {
            // Clear everything
            clear();
            return;
        }

        // Empty the gme_vertices_ of the deleted vertices and erase them
        for( index_t v : range( nb() ) ) {
            vertex_mapper_.clear_geomodel_vertex_gmes( v );
        }

        // Delete the vertices - false is to not remove
        // isolated vertices (here all the vertices)
        PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ )->delete_vertices(
            to_delete_bool );

        vertex_mapper_.update_mesh_entity_maps_and_gmes( to_delete );
    }

    template< index_t DIMENSION >
    GeoModelMeshVertices< DIMENSION >::GeoModelMeshVertices(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< PointSetMesh< DIMENSION > >& mesh )
        : GeoModelMeshVerticesBase< DIMENSION >( gmm, gm, mesh )
    {
    }

    GeoModelMeshVertices< 3 >::GeoModelMeshVertices(
        GeoModelMesh< 3 >& gmm,
        GeoModel< 3 >& gm,
        std::unique_ptr< PointSetMesh< 3 > >& mesh )
        : GeoModelMeshVerticesBase< 3 >( gmm, gm, mesh )
    {
    }

    void GeoModelMeshVertices< 3 >::clear()
    {
        this->gmm_.cells.clear();
        GeoModelMeshVerticesBase< 3 >::clear();
    }

    index_t GeoModelMeshVertices< 3 >::nb_total_vertices() const
    {
        index_t nb = GeoModelMeshVerticesBase< 3 >::nb_total_vertices();
        nb += nb_entity_vertices( this->gm_, Region< 3 >::type_name_static() );
        return nb;
    }

    index_t GeoModelMeshVertices< 3 >::fill_vertices()
    {
        index_t count = GeoModelMeshVerticesBase< 3 >::fill_vertices();
        fill_vertices_for_entity_type( this->gm_, Region< 3 >::type_name_static(),
            count );
        return count;
    }

    template< >
    GeoModelMeshVerticesBase< 3 >::GeoModelVertexMapper::GeoModelVertexMapper(
        GeoModelMeshVerticesBase& geomodel_vertices,
        const GeoModel< 3 >& geomodel )
        : geomodel_vertices_( geomodel_vertices ), geomodel_( geomodel )
    {
        vertex_maps_[Corner< 3 >::type_name_static()] = &corner_vertex_maps_;
        vertex_maps_[Line< 3 >::type_name_static()] = &line_vertex_maps_;
        vertex_maps_[Surface< 3 >::type_name_static()] = &surface_vertex_maps_;
        vertex_maps_[Region< 3 >::type_name_static()] = &region_vertex_maps_;
    }

    /*******************************************************************************/

    template< index_t DIMENSION >
    GeoModelMeshCells< DIMENSION >::GeoModelMeshCells(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< VolumeMesh< DIMENSION > >& mesh )
        :
            GeoModelMeshCommon< DIMENSION >( gmm, gm ),
            mesh_( mesh ),
            nb_tet_( 0 ),
            nb_hex_( 0 ),
            nb_prism_( 0 ),
            nb_pyramid_( 0 ),
            nb_connector_( 0 ),
            mode_( NONE )
    {
        this->set_mesh( mesh_.get() );
    }

    template< index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::is_initialized() const
    {
        return mesh_->nb_cells() > 0;
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshCells* >( this )->initialize();
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::initialize()
    {
        this->gmm_.vertices.test_and_initialize();
        std::unique_ptr< VolumeMeshBuilder< DIMENSION > > mesh_builder =
            VolumeMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != this->gmm_.vertices.nb() ) {
            copy_vertices( mesh_builder.get(), *this->gmm_.vertices.mesh_ );
        }

        region_cell_ptr_.resize(
            this->gm_.nb_regions() * to_underlying_type( CellType::UNDEFINED ) + 1,
            0 );

        // Total number of  cells
        std::vector< index_t > nb_cells_per_type(
            to_underlying_type( CellType::UNDEFINED ), 0 );
        index_t nb = 0;
        for( const auto& region : this->gm_.regions() ) {
            nb += region.nb_mesh_elements();
        }

        // Get out if no cells
        if( nb == 0 ) {
            return;
        }

        // Compute the number of cell per type and per region
        for( const auto& region : this->gm_.regions() ) {
            index_t r = region.index();
            for( index_t c : range( region.nb_mesh_elements() ) ) {
                CellType cur_cell_type = region.cell_type( c );
                switch( cur_cell_type ) {
                    case CellType::TETRAHEDRON:
                        nb_cells_per_type[to_underlying_type( CellType::TETRAHEDRON )]++;
                        region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                            * r + to_underlying_type( CellType::TETRAHEDRON ) + 1]++;
                        break;
                    case CellType::HEXAHEDRON:
                        nb_cells_per_type[to_underlying_type( CellType::HEXAHEDRON )]++;
                        region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                            * r + to_underlying_type( CellType::HEXAHEDRON ) + 1]++;
                        break;
                    case CellType::PRISM:
                        nb_cells_per_type[to_underlying_type( CellType::PRISM )]++;
                        region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                            * r + to_underlying_type( CellType::PRISM ) + 1]++;
                        break;
                    case CellType::PYRAMID:
                        nb_cells_per_type[to_underlying_type( CellType::PYRAMID )]++;
                        region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                            * r + to_underlying_type( CellType::PYRAMID ) + 1]++;
                        break;
                    case CellType::UNCLASSIFIED:
                        nb_cells_per_type[to_underlying_type(
                            CellType::UNCLASSIFIED )]++;
                        region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                            * r + to_underlying_type( CellType::UNCLASSIFIED ) + 1]++;
                        break;
                    default:
                        ringmesh_assert_not_reached;
                        break;
                }
            }
        }

        // Compute the cell offsets
        std::vector< index_t > cells_offset_per_type(
            to_underlying_type( CellType::UNDEFINED ), 0 );
        for( index_t t : range( to_underlying_type( CellType::TETRAHEDRON ) + 1,
            to_underlying_type( CellType::UNDEFINED ) ) ) {
            cells_offset_per_type[t] += cells_offset_per_type[t - 1];
            cells_offset_per_type[t] += nb_cells_per_type[t - 1];
        }

        for( index_t i : range( 1, region_cell_ptr_.size() - 1 ) ) {
            region_cell_ptr_[i + 1] += region_cell_ptr_[i];
        }

        // Create "empty" tet, hex, pyr and prism
        for( index_t i : range( to_underlying_type( CellType::UNDEFINED ) ) ) {
            mesh_builder->create_cells( nb_cells_per_type[i], CellType( i ) );
        }

        // Fill the cells with vertices
        bind_attribute();
        std::vector< index_t > cur_cell_per_type(
            to_underlying_type( CellType::UNDEFINED ), 0 );
        const GeoModelMeshVerticesBase< DIMENSION >& geomodel_vertices =
            this->gmm_.vertices;
        for( const auto& region : this->gm_.regions() ) {
            for( index_t c : range( region.nb_mesh_elements() ) ) {
                CellType cur_cell_type = region.cell_type( c );
                index_t cur_cell = cells_offset_per_type[to_underlying_type(
                    cur_cell_type )]
                    + cur_cell_per_type[to_underlying_type( cur_cell_type )]++;
                for( index_t v : range( mesh_->nb_cell_vertices( cur_cell ) ) ) {
                    index_t region_vertex_index = region.mesh_element_vertex_index(
                        c, v );
                    index_t global_vertex_id = geomodel_vertices.geomodel_vertex_id(
                        region.gmme(), region_vertex_index );
                    mesh_builder->set_cell_vertex( cur_cell, v, global_vertex_id );
                }
                region_id_[cur_cell] = region.index();
                cell_id_[cur_cell] = c;
            }
        }

        // Retrieve the adjacencies
        mesh_builder->connect_cells();

        // Permute cells to sort them per region and per type
        std::vector< index_t > sorted_indices( mesh_->nb_cells() );
        std::iota( sorted_indices.begin(), sorted_indices.end(), 0 );
        GeoModelMeshCellsSort< DIMENSION > action( *mesh_, region_id_ );
        std::sort( sorted_indices.begin(), sorted_indices.end(), action );
        mesh_builder->permute_cells( sorted_indices );

        // Cache some values
        nb_tet_ = nb_cells_per_type[to_underlying_type( CellType::TETRAHEDRON )];
        nb_hex_ = nb_cells_per_type[to_underlying_type( CellType::HEXAHEDRON )];
        nb_prism_ = nb_cells_per_type[to_underlying_type( CellType::PRISM )];
        nb_pyramid_ = nb_cells_per_type[to_underlying_type( CellType::PYRAMID )];
        nb_connector_ = nb_cells_per_type[to_underlying_type(
            CellType::UNCLASSIFIED )];
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::bind_attribute()
    {
        if( !region_id_.is_bound() ) {
            region_id_.bind( attribute_manager(), region_att_name );
        }
        if( !cell_id_.is_bound() ) {
            cell_id_.bind( attribute_manager(), cell_region_att_name );
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::unbind_attribute()
    {
        if( region_id_.is_bound() ) {
            region_id_.unbind();
        }
        if( cell_id_.is_bound() ) {
            cell_id_.unbind();
        }
        if( polygon_id_.is_bound() ) {
            polygon_id_.unbind();
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb() const
    {
        test_and_initialize();
        return mesh_->nb_cells();
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_vertices( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return mesh_->nb_cell_vertices( c );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::vertex( index_t c, index_t v ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        ringmesh_assert( v < mesh_->nb_cell_vertices( c ) );
        return mesh_->cell_vertex( c, v );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_edges( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return mesh_->nb_cell_edges( c );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_facets( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return mesh_->nb_cell_facets( c );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_facet_vertices(
        index_t c,
        index_t lf ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        ringmesh_assert( lf < mesh_->nb_cell_facets( c ) );
        return mesh_->nb_cell_facet_vertices( c, lf );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::facet_vertex(
        index_t c,
        index_t lf,
        index_t lv ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        ringmesh_assert( lf < mesh_->nb_cell_facets( c ) );
        return mesh_->cell_facet_vertex( c, lf, lv );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::edge_vertex(
        index_t c,
        index_t le,
        index_t lv ) const
    {
        geo_debug_assert( le < nb_edges( c ) );
        geo_debug_assert( lv < 2 );
        return mesh_->cell_edge_vertex( c, le, lv );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::adjacent( index_t c, index_t f ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        ringmesh_assert( f < mesh_->nb_cell_facets( c ) );
        return mesh_->cell_adjacent( c, f );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::region( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return region_id_[c];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::index_in_region( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return cell_id_[c];
    }

    template< index_t DIMENSION >
    CellType GeoModelMeshCells< DIMENSION >::type( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return mesh_->cell_type( c );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_cells( CellType type ) const
    {
        test_and_initialize();
        switch( type ) {
            case CellType::TETRAHEDRON:
                return nb_tet();
            case CellType::HEXAHEDRON:
                return nb_hex();
            case CellType::PRISM:
                return nb_prism();
            case CellType::PYRAMID:
                return nb_pyramid();
            case CellType::UNCLASSIFIED:
                return nb_connector();
            case CellType::UNDEFINED:
                return nb();
            default:
                ringmesh_assert_not_reached;
                return 0;
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_cells(
        index_t r,
        CellType type ) const
    {
        test_and_initialize();
        switch( type ) {
            case CellType::TETRAHEDRON:
                return nb_tet( r );
            case CellType::HEXAHEDRON:
                return nb_hex( r );
            case CellType::PRISM:
                return nb_prism( r );
            case CellType::PYRAMID:
                return nb_pyramid( r );
            case CellType::UNCLASSIFIED:
                return nb_connector( r );
            case CellType::UNDEFINED:
                ringmesh_assert(
                    region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                        * ( r + 1 )]
                        - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                            * r] == this->gm_.region( r ).nb_mesh_elements() );
                return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                    * ( r + 1 )]
                    - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r];
            default:
                ringmesh_assert_not_reached;
                return 0;
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::cell(
        index_t r,
        index_t c,
        CellType type ) const
    {
        test_and_initialize();
        switch( type ) {
            case CellType::TETRAHEDRON:
                return tet( r, c );
            case CellType::HEXAHEDRON:
                return hex( r, c );
            case CellType::PRISM:
                return prism( r, c );
            case CellType::PYRAMID:
                return pyramid( r, c );
            case CellType::UNCLASSIFIED:
                return connector( r, c );
            case CellType::UNDEFINED:
                return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r]
                    + c;
            default:
                ringmesh_assert_not_reached;
                return 0;
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_tet() const
    {
        test_and_initialize();
        return nb_tet_;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_tet( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::TETRAHEDRON ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::TETRAHEDRON )];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::tet( index_t r, index_t t ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::TETRAHEDRON )] + t;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_hex() const
    {
        test_and_initialize();
        return nb_hex_;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_hex( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::HEXAHEDRON ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::HEXAHEDRON )];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::hex( index_t r, index_t h ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::HEXAHEDRON )] + h;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_prism() const
    {
        test_and_initialize();
        return nb_prism_;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_prism( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::PRISM ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::PRISM )];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::prism( index_t r, index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::PRISM )] + p;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_pyramid() const
    {
        test_and_initialize();
        return nb_pyramid_;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_pyramid( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::PYRAMID ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::PYRAMID )];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::pyramid( index_t r, index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::PYRAMID )] + p;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_connector() const
    {
        test_and_initialize();
        return nb_connector_;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_connector( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::UNCLASSIFIED ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::UNCLASSIFIED )];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::connector( index_t r, index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( r < this->gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::UNCLASSIFIED )] + c;
    }

    template< index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::is_duplication_initialized() const
    {
        return mode_ == this->gmm_.duplicate_mode();
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::test_and_initialize_duplication() const
    {
        if( !is_duplication_initialized() ) {
            const_cast< GeoModelMeshCells< DIMENSION >* >( this )->initialize_duplication();
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::initialize_duplication()
    {
        test_and_initialize();
        clear_duplication();

        /// 1. Get all the corner vertices (a lot of duplicated vertices)
        std::vector< vec3 > corner_vertices(
            mesh_->cell_end( mesh_->nb_cells() - 1 ) );
        for( index_t c : range( mesh_->nb_cells() ) ) {
            index_t begin = mesh_->cell_begin( c );
            for( index_t v : range( mesh_->nb_cell_vertices( c ) ) ) {
                corner_vertices[begin + v] = mesh_->vertex(
                    mesh_->cell_vertex( c, v ) );
            }
        }

        /// 2. Tag all corners to duplicate (vertices on a surface to duplicate)
        std::vector< ActionOnSurface > actions_on_surfaces( this->gm_.nb_surfaces(),
            SKIP );
        std::vector< bool > is_vertex_to_duplicate( corner_vertices.size(), false );
        {
            NNSearch< DIMENSION > nn_search( corner_vertices, false );
            for( const auto& surface : this->gm_.surfaces() ) {
                if( !is_surface_to_duplicate( surface.index() ) ) continue;
                actions_on_surfaces[surface.index()] = TO_PROCESS;
                for( index_t v : range( surface.nb_vertices() ) ) {
                    std::vector< index_t > colocated_corners =
                        nn_search.get_neighbors( surface.vertex( v ),
                            this->gm_.epsilon() );
                    for( index_t co : colocated_corners ) {
                        is_vertex_to_duplicate[co] = true;
                    }
                }
            }
        }
        // Free some memory
        corner_vertices.clear();

        /// 3. Duplicate the corners
        /* The goal is to visit the corners of the GeoModelMesh
         * that are on one side of a Surface. We propagate through the cells
         * that have one vertex on a Surface without crossing the Surface.
         * All the corners visited during this propagation around the vertex
         * are duplicated if needed.
         */
        this->gmm_.polygons.test_and_initialize();
        for( index_t c : range( mesh_->nb_cells() ) ) {
            for( index_t v : range( mesh_->nb_cell_vertices( c ) ) ) {
                // get the index of the corner inside cell_corners_
                index_t co = mesh_->cell_begin( c ) + v;

                if( !is_vertex_to_duplicate[co] ) continue;
                // The vertex is on a surface to duplicate

                // Propagate on the cells around the corresponding vertex.
                // The propagation process cannot cross any surface.
                index_t vertex_id = mesh_->cell_vertex( c, v );

                // all the cell corners resulting of the propagation
                std::vector< index_t > corner_used;

                // all cells used during the propagation, used to provide
                // adding the same cell several times into the stack
                std::vector< index_t > cell_added;

                // all the surfaces encountered during the propagation
                // and which side stopped the propagation
                std::vector< action_on_surface > surfaces;

                // stack of the front of cells
                std::stack< index_t > S;
                S.push( c );
                cell_added.push_back( c );
                do {
                    index_t cur_c = S.top();
                    S.pop();
                    // Find which corner of the current cell matches vertex_id
                    index_t cur_co = mesh_->find_cell_corner( cur_c, vertex_id );
                    ringmesh_assert( cur_co != NO_ID );
                    is_vertex_to_duplicate[cur_co] = false;
                    corner_used.push_back( cur_co );

                    // Find the cell facets including the vertex
                    std::vector< index_t > facets = cell_facets_around_vertex(
                        *mesh_, cur_c, vertex_id );
                    for( index_t cur_f : facets ) {
                        // Find if the facet is on a surface or inside the domain
                        index_t polygon = NO_ID;
                        bool side;
                        if( is_cell_facet_on_surface( cur_c, cur_f, polygon,
                            side ) ) {
                            index_t surface_id = this->gmm_.polygons.surface(
                                polygon );
                            surfaces.push_back(
                                action_on_surface( surface_id,
                                    ActionOnSurface( side ) ) );
                        } else {
                            // The cell facet is not on a surface.
                            // Add the adjacent cell to the stack if it exists
                            // and has not already been processed or added into the stack
                            index_t cur_adj = mesh_->cell_adjacent( cur_c, cur_f );
                            if( cur_adj != GEO::NO_CELL
                                && !contains( cell_added, cur_adj ) ) {
                                cell_added.push_back( cur_adj );
                                S.push( cur_adj );
                            }
                        }
                    }
                } while( !S.empty() );

                // Remove redundant occurrences and sort the remaining ones
                sort_unique( surfaces );

                // Determine if the corners should be duplicated or not because
                // we need to duplicate only one side of the surface
                if( are_corners_to_duplicate( surfaces, actions_on_surfaces ) ) {
                    // Add a new duplicated vertex and its associated vertex

                    /* @todo Review : Use the total_nb_vertices function [JP]
                     * why mm_.vertices.nb_vertices() and not nb_vertices() ?
                     * Please help the reader !! same thing 2 lines below [JP]
                     */
                    index_t duplicated_vertex_id =
                        this->gmm_.vertices.nb()
                            + static_cast< index_t >( duplicated_vertex_indices_.size() );
                    duplicated_vertex_indices_.push_back( vertex_id );

                    // Update all the cell corners on this side of the surface
                    // to the new duplicated vertex index
                    std::unique_ptr< VolumeMeshBuilder< DIMENSION > > mesh_builder =
                        VolumeMeshBuilder< DIMENSION >::create_builder( *mesh_ );
                    for( index_t cur_co : corner_used ) {
                        mesh_builder->set_cell_corner_vertex_index( cur_co,
                            duplicated_vertex_id );
                    }
                }
            }
        }
        mode_ = this->gmm_.duplicate_mode();
    }

    template< index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::is_cell_facet_on_surface(
        index_t c,
        index_t f,
        index_t& polygon,
        bool& side ) const
    {
        test_and_initialize_cell_facet();
        polygon = polygon_id_[mesh_->cell_facet( c, f )];
        if( polygon != NO_ID ) {
            vec3 facet_normal = this->gmm_.polygons.normal( polygon );
            vec3 cell_facet_normal = mesh_->cell_facet_normal( c, f );
            side = dot( facet_normal, cell_facet_normal ) > 0;
        }
        return polygon != NO_ID;
    }

    template< index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::are_corners_to_duplicate(
        const std::vector< action_on_surface >& surfaces,
        std::vector< ActionOnSurface >& info )
    {
        // Temporary vector, it is equal to surfaces
        std::vector< action_on_surface > temp_surfaces;
        temp_surfaces.reserve( temp_surfaces.size() );

        // Find if a free border was found, if so we encountered the
        // two sides of the surface during the propagation around a vertex
        for( index_t i : range( 1, surfaces.size() ) ) {
            if( surfaces[i - 1].first == surfaces[i].first ) {
                // Found free border -> skip
                ringmesh_assert( surfaces[i - 1].second != surfaces[i].second );
                i++; // skip the action_on_surface (i-1) and i
            } else {
                temp_surfaces.push_back( surfaces[i - 1] );
            }
        }
        temp_surfaces.push_back( surfaces.back() );

        for( const action_on_surface& action : temp_surfaces ) {
            index_t s = action.first;
            switch( info[s] ) {
                case SKIP:
                    break;
                case TO_PROCESS:
                    // First time we encounter this surface, do not duplicate
                    // this side but wait to see if we encounter the other.
                    // In the case of surfaces in the VOI, it is encountered only once
                    ringmesh_assert( action.second > TO_PROCESS );
                    info[s] = ActionOnSurface( !action.second );
                    break;
                default:
                    // If the side matches -> duplicate
                    if( info[s] == action.second ) {
                        return true;
                    }
            }
        }

        return false;
    }

    template< index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::is_surface_to_duplicate(
        index_t surface_id ) const
    {
        const Surface< DIMENSION >& cur_surface = this->gm_.surface( surface_id );
        if( cur_surface.is_on_voi() ) return false;
        switch( this->gmm_.duplicate_mode() ) {
            case ALL:
                return true;
            case FAULT: {
                gmge_id parent_interface = cur_surface.parent_gmge(
                    Interface< DIMENSION >::type_name_static() );
                if( parent_interface.is_defined() ) {
                    typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE feature =
                        this->gm_.geological_entity( parent_interface ).geological_feature();
                    return GeoModelGeologicalEntity< DIMENSION >::is_fault( feature );
                }
                return false;
            }
            default:
                return false;
        }
        return false;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_duplicated_vertices() const
    {
        test_and_initialize_duplication();
        return static_cast< index_t >( duplicated_vertex_indices_.size() );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_total_vertices() const
    {
        return nb_duplicated_vertices() + mesh_->nb_vertices();
    }

    template< index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::is_corner_duplicated(
        index_t c,
        index_t v,
        index_t& duplicate_vertex_index ) const
    {
        test_and_initialize_duplication();
        ringmesh_assert( c < mesh_->nb_cells() );
        ringmesh_assert( v < mesh_->nb_cell_vertices( c ) );
        index_t corner_value = mesh_->cell_vertex( c, v );
        if( corner_value < mesh_->nb_vertices() ) {
            return false;
        } else {
            duplicate_vertex_index = corner_value - mesh_->nb_vertices();
            return true;
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::duplicated_vertex(
        index_t duplicate_vertex_index ) const
    {
        test_and_initialize_duplication();
        ringmesh_assert(
            duplicate_vertex_index < duplicated_vertex_indices_.size() );
        return duplicated_vertex_indices_[duplicate_vertex_index];
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::clear()
    {
        std::unique_ptr< VolumeMeshBuilder< DIMENSION > > mesh_builder =
            VolumeMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        mesh_builder->clear( true, false );
        region_cell_ptr_.clear();
        nb_tet_ = 0;
        nb_hex_ = 0;
        nb_prism_ = 0;
        nb_pyramid_ = 0;
        nb_connector_ = 0;

        mode_ = NONE;
        duplicated_vertex_indices_.clear();
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::clear_duplication()
    {
        std::unique_ptr< VolumeMeshBuilder< DIMENSION > > mesh_builder =
            VolumeMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        for( index_t c : range( mesh_->nb_cells() ) ) {
            for( index_t v : range( mesh_->nb_cell_vertices( c ) ) ) {
                index_t index = NO_ID;
                if( is_corner_duplicated( c, v, index ) ) {
                    mesh_builder->set_cell_corner_vertex_index( c,
                        duplicated_vertex( index ) );
                }
            }
        }

        mode_ = NONE;
        duplicated_vertex_indices_.clear();
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::test_and_initialize_cell_facet() const
    {
        if( !polygon_id_.is_bound() ) {
            const_cast< GeoModelMeshCells* >( this )->initialize_cell_facet();
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::initialize_cell_facet()
    {
        this->gmm_.polygons.test_and_initialize();

        polygon_id_.bind( mesh_->cell_facet_attribute_manager(), "polygon_id" );
        polygon_id_.fill( NO_ID );
        const NNSearch< DIMENSION >& nn_search = this->gmm_.polygons.nn_search();
        for( index_t c : range( mesh_->nb_cells() ) ) {
            for( index_t f : range( mesh_->nb_cell_facets( c ) ) ) {
                std::vector< index_t > result = nn_search.get_neighbors(
                    mesh_->cell_facet_barycenter( c, f ), this->gm_.epsilon() );
                if( !result.empty() ) {
                    polygon_id_[mesh_->cell_facet( c, f )] = result[0];
                    // If there are more than 1 matching facet, this is WRONG
                    // and the vertex indices should be checked too [Jeanne]
                    ringmesh_assert( result.size() == 1 );
                }
            }
        }
    }

    template< index_t DIMENSION >
    vec3 GeoModelMeshCells< DIMENSION >::barycenter( index_t c ) const
    {
        test_and_initialize();
        return mesh_->cell_barycenter( c );
    }

    template< index_t DIMENSION >
    double GeoModelMeshCells< DIMENSION >::volume( index_t c ) const
    {
        test_and_initialize();
        return mesh_->cell_volume( c );
    }

    template< index_t DIMENSION >
    const VolumeAABBTree< DIMENSION >& GeoModelMeshCells< DIMENSION >::aabb() const
    {
        test_and_initialize();
        return mesh_->cell_aabb();
    }

    /*******************************************************************************/

    template< index_t DIMENSION >
    GeoModelMeshPolygonsBase< DIMENSION >::GeoModelMeshPolygonsBase(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< SurfaceMesh< DIMENSION > >& mesh )
        :
            GeoModelMeshCommon< DIMENSION >( gmm, gm ),
            mesh_( mesh ),
            nb_triangle_( 0 ),
            nb_quad_( 0 ),
            nb_unclassified_polygon_( 0 )
    {
        this->set_mesh( mesh_.get() );
    }

    template< index_t DIMENSION >
    GeoModelMeshPolygonsBase< DIMENSION >::~GeoModelMeshPolygonsBase()
    {
        unbind_attribute();
    }

    template< index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::bind_attribute()
    {
        if( !surface_id_.is_bound() ) {
            surface_id_.bind( attribute_manager(), surface_att_name );
        }
        if( !polygon_id_.is_bound() ) {
            polygon_id_.bind( attribute_manager(), polygon_surface_att_name );
        }

    }

    template< index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::unbind_attribute()
    {
        if( surface_id_.is_bound() ) {
            surface_id_.unbind();
        }
        if( polygon_id_.is_bound() ) {
            polygon_id_.unbind();
        }

    }

    template< index_t DIMENSION >
    bool GeoModelMeshPolygonsBase< DIMENSION >::is_initialized() const
    {
        return mesh_->nb_polygons() > 0;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb() const
    {
        test_and_initialize();
        return mesh_->nb_polygons();
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_vertices( index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        return mesh_->nb_polygon_vertices( p );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::vertex(
        index_t p,
        index_t v ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        ringmesh_assert( v < mesh_->nb_polygon_vertices( p ) );
        return mesh_->polygon_vertex( p, v );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::adjacent(
        index_t p,
        index_t e ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        ringmesh_assert( e < mesh_->nb_polygon_vertices( p ) );
        return mesh_->polygon_adjacent( p, e );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::surface( index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        return surface_id_[p];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::index_in_surface(
        index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        return polygon_id_[p];
    }

    template< index_t DIMENSION >
    PolygonType GeoModelMeshPolygonsBase< DIMENSION >::type(
        index_t p,
        index_t& index ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        index_t polygon = index_in_surface( p );
        index_t s = surface( p );
        for( index_t t : range( to_underlying_type( PolygonType::TRIANGLE ),
            to_underlying_type( PolygonType::UNDEFINED ) ) ) {
            PolygonType T = static_cast< PolygonType >( t );
            if( polygon < nb_polygons( s, T ) ) {
                index = polygon;
                return T;
            }
            polygon -= nb_polygons( s, T );
        }
        index = NO_ID;
        ringmesh_assert_not_reached;
        return PolygonType::UNDEFINED;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_polygons(
        PolygonType type ) const
    {
        test_and_initialize();
        switch( type ) {
            case PolygonType::TRIANGLE:
                return nb_triangle();
            case PolygonType::QUAD:
                return nb_quad();
            case PolygonType::UNCLASSIFIED:
                return nb_unclassified_polygon();
            case PolygonType::UNDEFINED:
                return nb();
            default:
                ringmesh_assert_not_reached;
                return 0;
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_polygons(
        index_t s,
        PolygonType type ) const
    {
        test_and_initialize();
        ringmesh_assert( s < this->gm_.nb_surfaces() );
        switch( type ) {
            case PolygonType::TRIANGLE:
                return nb_triangle( s );
            case PolygonType::QUAD:
                return nb_quad( s );
            case PolygonType::UNCLASSIFIED:
                return nb_unclassified_polygon( s );
            case PolygonType::UNDEFINED:
                return surface_polygon_ptr_[to_underlying_type(
                    PolygonType::UNDEFINED ) * ( s + 1 )]
                    - surface_polygon_ptr_[to_underlying_type(
                        PolygonType::UNDEFINED ) * s];
            default:
                ringmesh_assert_not_reached;
                return 0;
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::polygon(
        index_t s,
        index_t p,
        PolygonType type ) const
    {
        test_and_initialize();
        ringmesh_assert( s < this->gm_.nb_surfaces() );
        switch( type ) {
            case PolygonType::TRIANGLE:
                return triangle( s, p );
            case PolygonType::QUAD:
                return quad( s, p );
            case PolygonType::UNCLASSIFIED:
                return unclassified_polygon( s, p );
            case PolygonType::UNDEFINED:
                return surface_polygon_ptr_[to_underlying_type(
                    PolygonType::UNDEFINED ) * s] + p;
            default:
                ringmesh_assert_not_reached;
                return 0;
        }
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_triangle() const
    {
        test_and_initialize();
        return nb_triangle_;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_triangle( index_t s ) const
    {
        test_and_initialize();
        ringmesh_assert( s < this->gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + ( to_underlying_type( PolygonType::TRIANGLE ) + 1 )]
            - surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
                + to_underlying_type( PolygonType::TRIANGLE )];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::triangle(
        index_t s,
        index_t t ) const
    {
        test_and_initialize();
        ringmesh_assert( s < this->gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + to_underlying_type( PolygonType::TRIANGLE )] + t;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_quad() const
    {
        test_and_initialize();
        return nb_quad_;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_quad( index_t s ) const
    {
        test_and_initialize();
        ringmesh_assert( s < this->gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + ( to_underlying_type( PolygonType::QUAD ) + 1 )]
            - surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
                + to_underlying_type( PolygonType::QUAD )];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::quad( index_t s, index_t q ) const
    {
        test_and_initialize();
        ringmesh_assert( s < this->gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + to_underlying_type( PolygonType::QUAD )] + q;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_unclassified_polygon() const
    {
        test_and_initialize();
        return nb_unclassified_polygon_;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_unclassified_polygon(
        index_t s ) const
    {
        test_and_initialize();
        ringmesh_assert( s < this->gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + ( to_underlying_type( PolygonType::UNCLASSIFIED ) + 1 )]
            - surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
                + to_underlying_type( PolygonType::UNCLASSIFIED )];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::unclassified_polygon(
        index_t s,
        index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( s < this->gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + to_underlying_type( PolygonType::UNCLASSIFIED )] + p;
    }

    template< index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::clear()
    {
        surface_polygon_ptr_.clear();
        nb_triangle_ = 0;
        nb_quad_ = 0;
        std::unique_ptr< SurfaceMeshBuilder< DIMENSION > > mesh_builder =
            SurfaceMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        mesh_builder->clear( true, false );
    }

    template< index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshPolygonsBase* >( this )->initialize();
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::initialize()
    {
        this->gmm_.vertices.test_and_initialize();
        clear();
        surface_polygon_ptr_.resize(
            this->gm_.nb_surfaces() * to_underlying_type( PolygonType::UNDEFINED )
                + 1, 0 );
        std::unique_ptr< SurfaceMeshBuilder< DIMENSION > > mesh_builder =
            SurfaceMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != this->gmm_.vertices.nb() ) {
            copy_vertices( mesh_builder.get(), *this->gmm_.vertices.mesh_ );
        }

        // Compute the total number of polygons per type and per surface
        std::map< PolygonType, index_t > nb_polygon_per_type = {
            { PolygonType::TRIANGLE, 0 }, { PolygonType::QUAD, 0 }, {
                PolygonType::UNCLASSIFIED, 0 } };
        for( index_t s : range( this->gm_.nb_surfaces() ) ) {
            const Surface< DIMENSION >& surface = this->gm_.surface( s );
            if( surface.is_simplicial() ) {
                nb_polygon_per_type[PolygonType::TRIANGLE] +=
                    surface.nb_mesh_elements();
                surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
                    + to_underlying_type( PolygonType::TRIANGLE ) + 1] +=
                    surface.nb_mesh_elements();
            } else {
                for( index_t p : range( surface.nb_mesh_elements() ) ) {
                    switch( surface.nb_mesh_element_vertices( p ) ) {
                        case 3:
                            nb_polygon_per_type[PolygonType::TRIANGLE]++;
                            surface_polygon_ptr_[to_underlying_type(
                                PolygonType::UNDEFINED ) * s
                                + to_underlying_type( PolygonType::TRIANGLE ) + 1]++;
                            break;
                        case 4:
                            nb_polygon_per_type[PolygonType::QUAD]++;
                            surface_polygon_ptr_[to_underlying_type(
                                PolygonType::UNDEFINED ) * s
                                + to_underlying_type( PolygonType::QUAD ) + 1]++;
                            break;
                        default:
                            nb_polygon_per_type[PolygonType::UNCLASSIFIED]++;
                            surface_polygon_ptr_[to_underlying_type(
                                PolygonType::UNDEFINED ) * s
                                + to_underlying_type( PolygonType::UNCLASSIFIED ) + 1]++;
                            break;
                    }
                }
            }
        }

        // Get out if no polygons
        index_t nb_total_polygons = nb_polygon_per_type[PolygonType::TRIANGLE]
            + nb_polygon_per_type[PolygonType::QUAD]
            + nb_polygon_per_type[PolygonType::UNCLASSIFIED];

        if( nb_total_polygons == 0 ) {
            return;
        }

        // Create triangles and quads, the polygons will be handle later
        if( nb_polygon_per_type[PolygonType::TRIANGLE] ) {
            mesh_builder->create_triangles(
                nb_polygon_per_type[PolygonType::TRIANGLE] );
        }
        if( nb_polygon_per_type[PolygonType::QUAD] ) {
            mesh_builder->create_quads( nb_polygon_per_type[PolygonType::QUAD] );
        }

        // Compute the polygon offset
        std::map< PolygonType, index_t > polygon_offset_per_type = {
            { PolygonType::TRIANGLE, 0 }, { PolygonType::QUAD, 0 }, {
                PolygonType::UNCLASSIFIED, 0 } };
        polygon_offset_per_type[PolygonType::QUAD] =
            nb_polygon_per_type[PolygonType::TRIANGLE];
        polygon_offset_per_type[PolygonType::UNCLASSIFIED] =
            nb_polygon_per_type[PolygonType::TRIANGLE]
                + nb_polygon_per_type[PolygonType::QUAD];

        for( index_t i : range( 1, surface_polygon_ptr_.size() - 1 ) ) {
            surface_polygon_ptr_[i + 1] += surface_polygon_ptr_[i];
        }

        // Fill the triangles and quads created above
        // Create and fill polygons
        bind_attribute();
        const GeoModelMeshVerticesBase< DIMENSION >& geomodel_vertices =
            this->gmm_.vertices;
        std::vector< index_t > cur_polygon_per_type(
            to_underlying_type( PolygonType::UNDEFINED ), 0 );
        for( index_t s : range( this->gm_.nb_surfaces() ) ) {
            const Surface< DIMENSION >& surface = this->gm_.surface( s );
            gmme_id surface_id = surface.gmme();
            for( index_t p : range( surface.nb_mesh_elements() ) ) {
                index_t nb_vertices = surface.nb_mesh_element_vertices( p );
                index_t cur_polygon = NO_ID;
                if( nb_vertices < 5 ) {
                    PolygonType T = static_cast< PolygonType >( nb_vertices - 3 );
                    cur_polygon = polygon_offset_per_type[T]
                        + cur_polygon_per_type[to_underlying_type( T )]++;
                    for( index_t v : range( nb_vertices ) ) {
                        index_t v_id = geomodel_vertices.geomodel_vertex_id(
                            surface_id, p, v );
                        ringmesh_assert( v_id != NO_ID );
                        mesh_builder->set_polygon_vertex( cur_polygon, v, v_id );
                    }
                } else {
                    std::vector< index_t > vertices( nb_vertices );
                    for( index_t v : range( nb_vertices ) ) {
                        vertices[v] = geomodel_vertices.geomodel_vertex_id(
                            surface_id, p, v );
                    }
                    cur_polygon = mesh_builder->create_polygon( vertices );
                }
                surface_id_[cur_polygon] = s;
                polygon_id_[cur_polygon] = p;
            }
        }

        // Permute polygons to sort them per surface and per type
        // Example for a mesh with two surfaces and only triangles and quads
        // [TRGL,TRGL, .. , QUAD, QUAD .. , TRGL, TRGL, ... , QUAD, QUAD ..]
        // |          surface 0           |             surface 1           |
        std::vector< index_t > sorted_indices( mesh_->nb_polygons() );
        std::iota( sorted_indices.begin(), sorted_indices.end(), 0 );
        GeoModelMeshPolygonsBaseSort< DIMENSION > action( *mesh_, surface_id_ );
        std::sort( sorted_indices.begin(), sorted_indices.end(), action );
        mesh_builder->permute_polygons( sorted_indices );

        // Compute polygon adjacencies
        mesh_builder->connect_polygons();
        disconnect_along_lines();

        // Cache some values
        nb_triangle_ = nb_polygon_per_type[PolygonType::TRIANGLE];
        nb_quad_ = nb_polygon_per_type[PolygonType::QUAD];
        nb_unclassified_polygon_ = nb_polygon_per_type[PolygonType::UNCLASSIFIED];
    }

    template< index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::disconnect_along_lines()
    {
        std::unique_ptr< SurfaceMeshBuilder< DIMENSION > > mesh_builder =
            SurfaceMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        for( index_t s : range( this->gm_.nb_surfaces() ) ) {
            const Surface< DIMENSION >& surface = this->gm_.surface( s );
            for( index_t p : range( nb_polygons( s ) ) ) {
                index_t polygon_id = polygon( s, p );
                index_t surface_polygon_id = index_in_surface( polygon_id );
                for( index_t v : range( nb_vertices( polygon_id ) ) ) {
                    index_t adj = surface.polygon_adjacent_index( surface_polygon_id,
                        v );
                    if( adj == NO_ID ) {
                        mesh_builder->set_polygon_adjacent( polygon_id, v, NO_ID );
                    }
                }
            }
        }
    }

    template< index_t DIMENSION >
    vecn< DIMENSION > GeoModelMeshPolygonsBase< DIMENSION >::center(
        index_t p ) const
    {
        test_and_initialize();
        return mesh_->polygon_barycenter( p );
    }

    template< index_t DIMENSION >
    double GeoModelMeshPolygonsBase< DIMENSION >::area( index_t p ) const
    {
        test_and_initialize();
        return mesh_->polygon_area( p );
    }

    template< index_t DIMENSION >
    const SurfaceAABBTree< DIMENSION >& GeoModelMeshPolygonsBase< DIMENSION >::aabb() const
    {
        test_and_initialize();
        return mesh_->polygon_aabb();
    }

    template< index_t DIMENSION >
    GeoModelMeshPolygons< DIMENSION >::GeoModelMeshPolygons(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< SurfaceMesh< DIMENSION > >& mesh )
        : GeoModelMeshPolygonsBase< DIMENSION >( gmm, gm, mesh )
    {
    }

    GeoModelMeshPolygons< 3 >::GeoModelMeshPolygons(
        GeoModelMesh< 3 >& gmm,
        GeoModel< 3 >& gm,
        std::unique_ptr< SurfaceMesh< 3 > >& mesh )
        : GeoModelMeshPolygonsBase< 3 >( gmm, gm, mesh )
    {
    }

    vec3 GeoModelMeshPolygons< 3 >::normal( index_t p ) const
    {
        test_and_initialize();
        return mesh_->polygon_normal( p );
    }

    /*******************************************************************************/

    template< index_t DIMENSION >
    GeoModelMeshEdges< DIMENSION >::GeoModelMeshEdges(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< LineMesh< DIMENSION > >& mesh )
        : GeoModelMeshCommon< DIMENSION >( gmm, gm ), mesh_( mesh )
    {
        this->set_mesh( mesh_.get() );
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::nb_wells() const
    {
        test_and_initialize();
        return this->gm_.wells() ? this->gm_.wells()->nb_wells() : 0;
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::nb_edges() const
    {
        test_and_initialize();
        return mesh_->nb_edges();
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::nb_edges( index_t w ) const
    {
        test_and_initialize();
        return well_ptr_[w + 1] - well_ptr_[w];
    }

    template< index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::vertex(
        index_t w,
        index_t e,
        index_t v ) const
    {
        test_and_initialize();
        return mesh_->edge_vertex( well_ptr_[w] + e, v );
    }

    template< index_t DIMENSION >
    void GeoModelMeshEdges< DIMENSION >::clear()
    {
        std::unique_ptr< LineMeshBuilder< DIMENSION > > mesh_builder =
            LineMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        mesh_builder->clear( true, false );
        well_ptr_.clear();
    }

    template< index_t DIMENSION >
    bool GeoModelMeshEdges< DIMENSION >::is_initialized() const
    {
        return mesh_->nb_edges() > 0;
    }

    template< index_t DIMENSION >
    void GeoModelMeshEdges< DIMENSION >::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshEdges* >( this )->initialize();
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshEdges< DIMENSION >::initialize()
    {
        if( !this->gm_.wells() ) return;
        this->gmm_.vertices.test_and_initialize();
        clear();
        std::unique_ptr< LineMeshBuilder< DIMENSION > > mesh_builder =
            LineMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != this->gmm_.vertices.nb() ) {
            copy_vertices( mesh_builder.get(), *this->gmm_.vertices.mesh_ );
        }

        // Compute the total number of edge per well
        const WellGroup< DIMENSION >& wells = *this->gm_.wells();
        well_ptr_.resize( wells.nb_wells() + 1, 0 );
        index_t nb_edges = 0;
        for( index_t w : range( wells.nb_wells() ) ) {
            nb_edges += wells.well( w ).nb_edges();
            well_ptr_[w + 1] = nb_edges;
        }

        // Compute the edge offset
        for( index_t i : range( 1, well_ptr_.size() - 1 ) ) {
            well_ptr_[i + 1] += well_ptr_[i];
        }

        // Create edges
        mesh_builder->create_edges( well_ptr_.back() );

        // Fill edges
        index_t cur_edge = 0;
        for( index_t w : range( 0, wells.nb_wells() ) ) {
            const Well< DIMENSION >& well = wells.well( w );
            for( index_t p : range( well.nb_parts() ) ) {
                for( index_t e : range( well.part( p ).nb_edges() ) ) {
                    const vecn< DIMENSION >& e0 = well.part( p ).edge_vertex( e, 0 );
                    mesh_builder->set_edge_vertex( cur_edge, 0,
                        this->gmm_.vertices.index( e0 ) );
                    const vecn< DIMENSION >& e1 = well.part( p ).edge_vertex( e, 1 );
                    mesh_builder->set_edge_vertex( cur_edge, 1,
                        this->gmm_.vertices.index( e1 ) );
                    cur_edge++;
                }
            }
        }
    }

    template< index_t DIMENSION >
    const LineAABBTree< DIMENSION >& GeoModelMeshEdges< DIMENSION >::aabb() const
    {
        test_and_initialize();
        return mesh_->edge_aabb();
    }

    /*******************************************************************************/

    template< index_t DIMENSION >
    GeoModelMeshBase< DIMENSION >::GeoModelMeshBase(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& geomodel )
        :
            geomodel_( geomodel ),
            vertices( gmm, geomodel, mesh_set_.point_set_mesh ),
            edges( gmm, geomodel, mesh_set_.line_mesh ),
            polygons( gmm, geomodel, mesh_set_.surface_mesh )
    {
    }

    template< index_t DIMENSION >
    GeoModelMeshBase< DIMENSION >::~GeoModelMeshBase()
    {
        polygons.unbind_attribute();
    }

    template< index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::remove_colocated_vertices()
    {
        vertices.remove_colocated();
    }

    template< index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::erase_vertices(
        std::vector< index_t >& to_delete )
    {
        vertices.erase_vertices( to_delete );
    }

    template< index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::change_point_set_mesh_data_structure(
        const MeshType& type )
    {
        if( mesh_set_.point_set_mesh->type_name() != type ) {
            vertices.clear();
            mesh_set_.create_point_set_mesh( type );
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::change_line_mesh_data_structure(
        const MeshType& type )
    {
        if( mesh_set_.line_mesh->type_name() != type ) {
            edges.clear();
            mesh_set_.create_line_mesh( type );
        }
    }

    template< index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::change_surface_mesh_data_structure(
        const MeshType& type )
    {
        if( mesh_set_.surface_mesh->type_name() != type ) {
            polygons.clear();
            mesh_set_.create_surface_mesh( type );
        }
    }

    template< index_t DIMENSION >
    GeoModelMesh< DIMENSION >::GeoModelMesh( GeoModel< DIMENSION >& geomodel )
        : GeoModelMeshBase< DIMENSION >( *this, geomodel )
    {
    }

    GeoModelMesh< 3 >::GeoModelMesh( GeoModel< 3 >& geomodel )
        :
            GeoModelMeshBase< 3 >( *this, geomodel ),
            cells( *this, geomodel, mesh_set_.volume_mesh )
    {
    }

    GeoModelMesh< 3 >::~GeoModelMesh()
    {
        cells.unbind_attribute();
    }

    void GeoModelMesh< 3 >::change_volume_mesh_data_structure( const MeshType& type )
    {
        if( mesh_set_.volume_mesh->type_name() != type ) {
            cells.clear();
            mesh_set_.create_volume_mesh( type );
        }
    }

    void GeoModelMesh< 3 >::transfer_attributes() const
    {
        transfer_vertex_attributes();
        transfer_cell_attributes();
    }

    void GeoModelMesh< 3 >::transfer_cell_attributes() const
    {

        GEO::vector< std::string > att_c_names;
        cells.attribute_manager().list_attribute_names( att_c_names );

        const NNSearch< 3 >& nn_search = cells.cell_nn_search();

        for( const std::string& att_c : att_c_names ) {
            if( !GEO::Attribute< double >::is_defined( cells.attribute_manager(),
                att_c ) ) {
                continue;
            }
            GEO::Attribute< double > cur_att_on_geomodel_mesh(
                cells.attribute_manager(), att_c );
            index_t att_dim = cur_att_on_geomodel_mesh.dimension();

            for( const auto& region : geomodel_.regions() ) {
                if( region.cell_attribute_manager().is_defined( att_c ) ) {
                    Logger::warn( "Transfer attribute", "The attribute ", att_c,
                        " already exists on the ", region.gmme() );
                    continue;
                }
                GEO::Attribute< double > cur_att_on_geomodel_mesh_entity;
                cur_att_on_geomodel_mesh_entity.create_vector_attribute(
                    region.cell_attribute_manager(), att_c, att_dim );
                for( index_t c : range( region.nb_mesh_elements() ) ) {
                    vec3 center = region.mesh_element_barycenter( c );
                    std::vector< index_t > c_in_geom_model_mesh =
                        nn_search.get_neighbors( center, geomodel_.epsilon() );
                    ringmesh_assert( c_in_geom_model_mesh.size() == 1 );
                    for( index_t att_e : range( att_dim ) ) {
                        cur_att_on_geomodel_mesh_entity[c * att_dim + att_e] =
                            cur_att_on_geomodel_mesh[c_in_geom_model_mesh[0]
                                * att_dim + att_e];
                    }
                }
            }
        }
    }

    void GeoModelMesh< 3 >::transfer_vertex_attributes() const
    {
        GEO::vector< std::string > att_v_names;
        std::vector< std::string > att_v_double_names;
        vertices.attribute_manager().list_attribute_names( att_v_names );
        for( index_t att_v : range( vertices.attribute_manager().nb() ) ) {
            if( !GEO::Attribute< double >::is_defined( vertices.attribute_manager(),
                att_v_names[att_v] ) ) {
                continue;
            }
            att_v_double_names.push_back( att_v_names[att_v] );
            for( const auto& region : geomodel_.regions() ) {
                if( region.vertex_attribute_manager().is_defined(
                    att_v_names[att_v] ) ) {
                    Logger::warn( "Transfer attribute", "The attribute ",
                        att_v_names[att_v], " already exists on the ",
                        region.gmme() );
                    continue;
                }
                GEO::Attribute< double > cur_v_att;
                cur_v_att.create_vector_attribute( region.vertex_attribute_manager(),
                    att_v_names[att_v],
                    vertices.attribute_manager().find_attribute_store(
                        att_v_names[att_v] )->dimension() );
            }
        }
        for( const std::string& att_v : att_v_double_names ) {
            GEO::Attribute< double > cur_att_on_geomodelmesh(
                vertices.attribute_manager(), att_v );
            index_t att_dim = cur_att_on_geomodelmesh.dimension();

            AttributeVector< double > att_on_regions( geomodel_.nb_regions() );

            for( const auto& region : geomodel_.regions() ) {
                att_on_regions.bind_one_attribute( region.index(),
                    region.vertex_attribute_manager(), att_v );
            }

            for( index_t v : range( vertices.nb() ) ) {
                std::vector< GMEVertex > vertices_on_geomodel_region =
                    vertices.gme_type_vertices( Region< 3 >::type_name_static(), v );
                for( const GMEVertex& cur_vertex_on_geomodel : vertices_on_geomodel_region ) {
                    for( index_t att_e : range( att_dim ) ) {
                        att_on_regions[cur_vertex_on_geomodel.gmme.index()][cur_vertex_on_geomodel.v_index
                            * att_dim + att_e] = cur_att_on_geomodelmesh[v * att_dim
                            + att_e];
                    }
                }

            }
        }
    }

    template class RINGMESH_API GeoModelMeshBase< 2 > ;
    template class RINGMESH_API GeoModelMesh< 2 > ;
    template class RINGMESH_API GeoModelMeshVerticesBase< 2 > ;
    template class RINGMESH_API GeoModelMeshEdges< 2 > ;
    template class RINGMESH_API GeoModelMeshPolygonsBase< 2 > ;

    template class RINGMESH_API GeoModelMeshBase< 3 > ;
    template class RINGMESH_API GeoModelMesh< 3 > ;
    template class RINGMESH_API GeoModelMeshVerticesBase< 3 > ;
    template class RINGMESH_API GeoModelMeshEdges< 3 > ;
    template class RINGMESH_API GeoModelMeshPolygonsBase< 3 > ;
    template class RINGMESH_API GeoModelMeshCells< 3 > ;

}
