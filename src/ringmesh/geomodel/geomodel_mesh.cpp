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

    class GeoModelMeshPolygonsSort {
    public:
        GeoModelMeshPolygonsSort(
            const SurfaceMesh& mesh,
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
        const SurfaceMesh& mesh_;
        const GEO::Attribute< index_t >& surface_id_;
    };

    class GeoModelMeshCellsSort {
    public:
        GeoModelMeshCellsSort(
            const VolumeMesh& mesh,
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
        const VolumeMesh& mesh_;
        const GEO::Attribute< index_t >& region_id_;
    };

    std::string vertex_map_name()
    {
        return "model_vertex_map";
    }

    void cell_facets_around_vertex(
        const VolumeMesh& mesh,
        index_t cell,
        index_t vertex_id,
        std::vector< index_t >& facets )
    {
        facets.reserve( mesh.nb_cell_facets( cell ) );
        for( index_t f = 0; f < mesh.nb_cell_facets( cell ); f++ ) {
            for( index_t v = 0; v < mesh.nb_cell_facet_vertices( cell, f ); v++ ) {
                if( mesh.cell_facet_vertex( cell, f, v ) == vertex_id ) {
                    facets.push_back( f );
                    break;
                }
            }
        }
    }

    void copy_vertices( MeshBaseBuilder* builder, const MeshBase& mesh )
    {
        builder->clear_vertices( true, false );
        builder->create_vertices( mesh.nb_vertices() );
        for( index_t v = 0; v < mesh.nb_vertices(); v++ ) {
            builder->set_vertex( v, mesh.vertex( v ) );
        }
    }

    bool is_attribute_a_double(
        GEO::AttributesManager& att_manager,
        const std::string& att_name )
    {
        return GEO::Attribute< double >::is_defined( att_manager, att_name ) ;
    }
}

namespace RINGMesh {

    GeoModelMeshBase::GeoModelMeshBase( GeoModelMesh& gmm, GeoModel& gm )
        : gmm_( gmm ), gm_( gm ), mesh_base_( nullptr )
    {
    }

    GeoModelMeshVertices::GeoModelVertexMapper::GeoModelVertexMapper(
        GeoModelMeshVertices& geomodel_vertices,
        const GeoModel& geomodel )
        : geomodel_vertices_( geomodel_vertices ), geomodel_( geomodel )
    {
        vertex_maps_[Corner::type_name_static()] = &corner_vertex_maps_;
        vertex_maps_[Line::type_name_static()] = &line_vertex_maps_;
        vertex_maps_[Surface::type_name_static()] = &surface_vertex_maps_;
        vertex_maps_[Region::type_name_static()] = &region_vertex_maps_;
    }

    void GeoModelMeshVertices::GeoModelVertexMapper::clear()
    {
        gme_vertices_.clear();
        clear_all_mesh_entity_vertex_map();
    }

    index_t GeoModelMeshVertices::GeoModelVertexMapper::geomodel_vertex_index(
        const gmme_id& mesh_entity_id,
        index_t mesh_entity_vertex_index ) const
    {
        ringmesh_assert(
            mesh_entity_vertex_index
                < geomodel_.mesh_entity( mesh_entity_id ).nb_vertices() );

        return vertex_map( mesh_entity_id )[mesh_entity_vertex_index];
    }

    const std::vector< GMEVertex >& GeoModelMeshVertices::GeoModelVertexMapper::mesh_entity_vertex_indices(
        index_t v ) const
    {
        ringmesh_assert( v < gme_vertices_.size() );
        return gme_vertices_[v];
    }

    std::vector< GMEVertex > GeoModelMeshVertices::GeoModelVertexMapper::mesh_entity_vertex_indices(
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

    std::vector< index_t > GeoModelMeshVertices::GeoModelVertexMapper::mesh_entity_vertex_indices(
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

    const GEO::Attribute< index_t >&
    GeoModelMeshVertices::GeoModelVertexMapper::vertex_map(
        const gmme_id& mesh_entity_id ) const
    {
        return ( *vertex_maps_.at( mesh_entity_id.type() ) )[mesh_entity_id.index()];
    }

    GEO::Attribute< index_t >& GeoModelMeshVertices::GeoModelVertexMapper::vertex_map(
        const gmme_id& mesh_entity_id )
    {
        return ( *vertex_maps_[mesh_entity_id.type()] )[mesh_entity_id.index()];
    }

    void GeoModelMeshVertices::GeoModelVertexMapper::set_vertex_map_value(
        const gmme_id& mesh_entity_id,
        index_t mesh_entity_vertex_index,
        index_t geomodel_entity_vertex_index )
    {
        test_and_initialize_mesh_entity_vertex_map( mesh_entity_id );
        vertex_map( mesh_entity_id )[mesh_entity_vertex_index] =
            geomodel_entity_vertex_index;
    }

    void GeoModelMeshVertices::GeoModelVertexMapper::add_to_gme_vertices(
        const GMEVertex& gme_vertex,
        index_t geomodel_vertex_index )
    {
        gme_vertices_[geomodel_vertex_index].push_back( gme_vertex );
    }

    void GeoModelMeshVertices::GeoModelVertexMapper::bind_all_mesh_entity_vertex_maps()
    {
        const std::vector< MeshEntityType >& all_mesh_entity_types =
            MeshEntityTypeManager::mesh_entity_types();
        for( const MeshEntityType& cur_entity_type : all_mesh_entity_types ) {
            index_t nb_cur_type_entities = geomodel_.nb_mesh_entities(
                cur_entity_type );
            vertex_maps_.at( cur_entity_type )->clear();
            vertex_maps_.at( cur_entity_type )->resize( nb_cur_type_entities );
            for( index_t e = 0; e < nb_cur_type_entities; e++ ) {
                const gmme_id cur_entity( cur_entity_type, e );
                bind_vertex_map( cur_entity );
            }
        }
    }

    GEO::Attribute< index_t >&
    GeoModelMeshVertices::GeoModelVertexMapper::bind_vertex_map(
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

    void GeoModelMeshVertices::GeoModelVertexMapper::update_mesh_entity_maps_and_gmes(
        const std::vector< index_t >& old2new )
    {
        const std::vector< MeshEntityType >& all_mesh_entity_types =
            MeshEntityTypeManager::mesh_entity_types();
        for( const MeshEntityType& cur_entity_type : all_mesh_entity_types ) {
            for( index_t e = 0; e < geomodel_.nb_mesh_entities( cur_entity_type );
                e++ ) {
                const GeoModelMeshEntity& E = geomodel_.mesh_entity( cur_entity_type,
                    e );
                gmme_id id = E.gmme();
                for( index_t v = 0; v < E.nb_vertices(); v++ ) {
                    index_t old_m_id = geomodel_vertex_index( id, v );
                    index_t new_m_id = old2new[old_m_id];
                    set_vertex_map_value( id, v, new_m_id );

                    // Merge gme_vertices information
                    if( std::find( gme_vertices_[new_m_id].begin(),
                        gme_vertices_[new_m_id].end(), GMEVertex( id, v ) )
                        == gme_vertices_[new_m_id].end() ) {
                        gme_vertices_[new_m_id].push_back( GMEVertex( id, v ) );
                    }
                }
            }
        }
    }

    void GeoModelMeshVertices::GeoModelVertexMapper::unbind_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        resize_all_mesh_entity_vertex_maps( mesh_entity_id.type() );
        if( vertex_maps_.at( mesh_entity_id.type() )->is_attribute_bound(
            mesh_entity_id.index() ) ) {
            vertex_maps_.at( mesh_entity_id.type() )->unbind(
                mesh_entity_id.index() );
        }
    }

    void GeoModelMeshVertices::GeoModelVertexMapper::initialize_mesh_entity_vertex_map(
        const gmme_id& mesh_entity_id )
    {

        GEO::Attribute< index_t >& mesh_entity_vertex_map = bind_vertex_map(
            mesh_entity_id );

        const GeoModelMeshEntity& E = geomodel_.mesh_entity( mesh_entity_id );
        for( index_t v = 0; v < E.nb_vertices(); v++ ) {
            mesh_entity_vertex_map[v] =
                geomodel_vertices_.nn_search().get_closest_neighbor( E.vertex( v ) );
        }
    }

    bool GeoModelMeshVertices::GeoModelVertexMapper::test_and_initialize_mesh_entity_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        resize_all_mesh_entity_vertex_maps( mesh_entity_id.type() );
        if( !is_mesh_entity_vertex_map_initialized( mesh_entity_id ) ) {
            initialize_mesh_entity_vertex_map( mesh_entity_id );
            return false;
        }
        return true;
    }

    bool GeoModelMeshVertices::GeoModelVertexMapper::is_mesh_entity_vertex_map_initialized(
        const gmme_id& mesh_entity_id ) const
    {
        return ( vertex_maps_.find( mesh_entity_id.type() )->second )->is_attribute_bound(
            mesh_entity_id.index() );
    }

    void GeoModelMeshVertices::GeoModelVertexMapper::clear_all_mesh_entity_vertex_map()
    {
        for( index_t t = 0; t < MeshEntityTypeManager::nb_mesh_entity_types();
            t++ ) {
            const MeshEntityType& cur_type =
                MeshEntityTypeManager::mesh_entity_types()[t];
            for( index_t e = 0; e < vertex_maps_[cur_type]->size(); e++ ) {
                vertex_maps_[cur_type]->unbind( e );
            }
            vertex_maps_[cur_type]->clear();
        }
    }

    void GeoModelMeshVertices::GeoModelVertexMapper::resize_all_mesh_entity_vertex_maps(
        const MeshEntityType& type )
    {
        vertex_maps_.at( type )->resize( geomodel_.nb_mesh_entities( type ),
            nullptr );
    }

    GEO::AttributesManager& GeoModelMeshVertices::GeoModelVertexMapper::mesh_entity_vertex_attribute_manager(
        const gmme_id& mesh_entity_id ) const
    {
        const GeoModelMeshEntity& mesh_entity = geomodel_.mesh_entity(
            mesh_entity_id );
        return mesh_entity.vertex_attribute_manager();
    }

    GeoModelMeshVertices::GeoModelMeshVertices(
        GeoModelMesh& gmm,
        GeoModel& gm,
        std::unique_ptr< PointSetMesh >& mesh )
        :
            GeoModelMeshBase( gmm, gm ),
            mesh_( mesh ),
            vertex_mapper_( *this, gmm.geomodel() )
    {
        set_mesh( mesh_.get() );
    }

    GeoModelMeshVertices::~GeoModelMeshVertices()
    {
    }

    bool GeoModelMeshVertices::is_initialized() const
    {
        return mesh_->nb_vertices() > 0;
    }

    void GeoModelMeshVertices::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshVertices* >( this )->initialize();
        }
    }

    index_t nb_entity_vertices(
        const GeoModel& M,
        const MeshEntityType& entity_type )
    {
        index_t count = 0;
        for( index_t i = 0; i < M.nb_mesh_entities( entity_type ); ++i ) {
            count += M.mesh_entity( entity_type, i ).nb_vertices();
        }
        return count;
    }

    void GeoModelMeshVertices::fill_vertices(
        const GeoModel& M,
        const MeshEntityType& entity_type,
        index_t& count )
    {
        std::unique_ptr< PointSetMeshBuilder > mesh_builder =
            PointSetMeshBuilder::create_builder( *mesh_ );
        for( index_t i = 0; i < M.nb_mesh_entities( entity_type ); ++i ) {
            GeoModelMeshEntity& E = const_cast< GeoModelMeshEntity& >( M.mesh_entity(
                entity_type, i ) );
            if( E.nb_vertices() == 0 ) {
                continue;
            }

            // Map and vertex
            gmme_id id = E.gmme();
            for( index_t v = 0; v < E.nb_vertices(); v++ ) {
                index_t local_count = count + v;
                mesh_builder->set_vertex( local_count, E.vertex( v ) );
                // Map from vertices of MeshEntities to GeoModelMeshVertices
                vertex_mapper_.set_vertex_map_value( id, v, local_count );
                vertex_mapper_.add_to_gme_vertices( GMEVertex( id, v ),
                    local_count );
            }
            // Global vertex index increment
            count += E.nb_vertices();
        }
    }

    void GeoModelMeshVertices::initialize()
    {
        clear();

        // Total number of vertices in the
        // Corners, Lines, Surfaces and Regions of the GeoModel
        index_t nb = 0;
        nb += nb_entity_vertices( gm_, Corner::type_name_static() );
        nb += nb_entity_vertices( gm_, Line::type_name_static() );
        nb += nb_entity_vertices( gm_, Surface::type_name_static() );
        nb += nb_entity_vertices( gm_, Region::type_name_static() );

        // Get out if no vertices
        if( nb == 0 ) {
            return;
        }

        // Fill the vertices
        std::unique_ptr< PointSetMeshBuilder > builder =
            PointSetMeshBuilder::create_builder( *mesh_ );
        builder->create_vertices( nb );
        vertex_mapper_.clear_and_resize_geomodel_vertex_gmes( nb );
        vertex_mapper_.bind_all_mesh_entity_vertex_maps();

        index_t count = 0;
        fill_vertices( gm_, Corner::type_name_static(), count );
        fill_vertices( gm_, Line::type_name_static(), count );
        fill_vertices( gm_, Surface::type_name_static(), count );
        fill_vertices( gm_, Region::type_name_static(), count );

        // Remove colocated vertices
        remove_colocated();
    }

    void GeoModelMeshVertices::clear()
    {
        gmm_.cells.clear();
        gmm_.polygons.clear();
        gmm_.edges.clear();
        vertex_mapper_.clear();

        std::unique_ptr< PointSetMeshBuilder > builder =
            PointSetMeshBuilder::create_builder( *mesh_ );
        builder->clear( true, false );
    }

    void GeoModelMeshVertices::unbind_geomodel_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        vertex_mapper_.unbind_vertex_map( mesh_entity_id );
    }

    void GeoModelMeshVertices::bind_geomodel_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        vertex_mapper_.bind_vertex_map( mesh_entity_id );
    }

    index_t GeoModelMeshVertices::nb() const
    {
        test_and_initialize();
        return mesh_->nb_vertices();
    }

    const vec3& GeoModelMeshVertices::vertex( index_t v ) const
    {
        test_and_initialize();
        ringmesh_assert( v < nb() );
        return mesh_->vertex( v );
    }

    index_t GeoModelMeshVertices::index( const vec3& p ) const
    {
        test_and_initialize();
        const NNSearch& colocator = mesh_->vertex_nn_search();
        std::vector< index_t > vertices = colocator.get_neighbors( p,
            gm_.epsilon() );
        if( vertices.empty() ) {
            return NO_ID;
        } else {
            ringmesh_assert( vertices.size() == 1 );
            return vertices[0];
        }
    }

    index_t GeoModelMeshVertices::geomodel_vertex_id(
        const gmme_id& mesh_entity,
        index_t entity_vertex_index ) const
    {
        test_and_initialize();
        return vertex_mapper_.geomodel_vertex_index( mesh_entity,
            entity_vertex_index );
    }

    index_t GeoModelMeshVertices::geomodel_vertex_id(
        const gmme_id& mesh_entity,
        index_t entity_mesh_element_index,
        index_t vertex_local_index ) const
    {
        index_t entity_vertex_index =
            gm_.mesh_entity( mesh_entity ).mesh_element_vertex_index(
                entity_mesh_element_index, vertex_local_index );
        return geomodel_vertex_id( mesh_entity, entity_vertex_index );
    }

    std::vector< index_t > GeoModelMeshVertices::mesh_entity_vertex_id(
        const gmme_id& mesh_entity,
        index_t geomodel_vertex_index ) const
    {
        test_and_initialize();
        ringmesh_assert( geomodel_vertex_index < nb() );
        return vertex_mapper_.mesh_entity_vertex_indices( geomodel_vertex_index,
            mesh_entity );
    }

    const std::vector< GMEVertex >& GeoModelMeshVertices::gme_vertices(
        index_t v ) const
    {
        test_and_initialize();
        return vertex_mapper_.mesh_entity_vertex_indices( v );
    }

    std::vector< GMEVertex > GeoModelMeshVertices::gme_type_vertices(
        const MeshEntityType& entity_type,
        index_t v ) const
    {
        test_and_initialize();
        return vertex_mapper_.mesh_entity_vertex_indices( v, entity_type );
    }

    index_t GeoModelMeshVertices::add_vertex( const vec3& point )
    {
        std::unique_ptr< PointSetMeshBuilder > builder =
            PointSetMeshBuilder::create_builder( *mesh_ );
        const index_t index = builder->create_vertex( point );
        vertex_mapper_.resize_geomodel_vertex_gmes( nb() );
        return index;
    }

    index_t GeoModelMeshVertices::add_vertices( const std::vector< vec3 >& points )
    {
        ringmesh_assert( !points.empty() );
        std::unique_ptr< PointSetMeshBuilder > builder =
            PointSetMeshBuilder::create_builder( *mesh_ );
        const index_t start_index = builder->create_vertex( points[0] );
        for( size_t i = 1; i < points.size(); ++i ) {
            builder->create_vertex( points[i] );
        }
        vertex_mapper_.resize_geomodel_vertex_gmes( nb() );
        return start_index;
    }

    void GeoModelMeshVertices::update_point( index_t v, const vec3& point )
    {
        test_and_initialize();
        ringmesh_assert( v < nb() );
        // Change the position of the unique_vertex
        std::unique_ptr< PointSetMeshBuilder > mesh_builder =
            PointSetMeshBuilder::create_builder( *mesh_ );
        mesh_builder->set_vertex( v, point );

        GeoModelBuilder builder( gm_ );

        const std::vector< GMEVertex >& gme_v = gme_vertices( v );
        for( const GMEVertex& info : gme_v ) {
            builder.geometry.set_mesh_entity_vertex( info.gmme, info.v_index, point,
                false );
        }
    }

    void GeoModelMeshVertices::update_vertex_mapping(
        const gmme_id& entity_id,
        index_t entity_vertex_index,
        index_t geomodel_vertex_index )
    {
        vertex_mapper_.set_vertex_map_value( entity_id, entity_vertex_index,
            geomodel_vertex_index );
        vertex_mapper_.add_to_gme_vertices(
            GMEVertex( entity_id, entity_vertex_index ), geomodel_vertex_index );
    }

    void GeoModelMeshVertices::remove_colocated()
    {
        // Get out if nothing to do
        // and compute the points if they are not initialized yet
        if( nb() == 0 ) {
            return;
        }
        // Identify and invalidate colocated vertices
        std::vector< index_t > old2new;
        index_t nb_colocalised_vertices =
            mesh_->vertex_nn_search().get_colocated_index_mapping( gm_.epsilon(),
                old2new );
        if( nb_colocalised_vertices > 0 ) {
            erase_vertices( old2new );
        }
    }

    void GeoModelMeshVertices::erase_vertices( std::vector< index_t >& to_delete )
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
        for( index_t v = 0; v < nb(); ++v ) {
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
        for( index_t v = 0; v < nb(); ++v ) {
            vertex_mapper_.clear_geomodel_vertex_gmes( v );
        }

        // Delete the vertices - false is to not remove
        // isolated vertices (here all the vertices)
        PointSetMeshBuilder::create_builder( *mesh_ )->delete_vertices(
            to_delete_bool );

        vertex_mapper_.update_mesh_entity_maps_and_gmes( to_delete );
    }

    /*******************************************************************************/

    GeoModelMeshCells::GeoModelMeshCells(
        GeoModelMesh& gmm,
        GeoModel& gm,
        std::unique_ptr< VolumeMesh >& mesh )
        :
            GeoModelMeshBase( gmm, gm ),
            mesh_( mesh ),
            nb_tet_( 0 ),
            nb_hex_( 0 ),
            nb_prism_( 0 ),
            nb_pyramid_( 0 ),
            nb_connector_( 0 ),
            mode_( NONE )
    {
        set_mesh( mesh_.get() );
    }

    bool GeoModelMeshCells::is_initialized() const
    {
        return mesh_->nb_cells() > 0;
    }

    void GeoModelMeshCells::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshCells* >( this )->initialize();
        }
    }

    void GeoModelMeshCells::initialize()
    {
        gmm_.vertices.test_and_initialize();
        clear();
        std::unique_ptr< VolumeMeshBuilder > mesh_builder =
            VolumeMeshBuilder::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != gmm_.vertices.nb() ) {
            copy_vertices( mesh_builder.get(), *gmm_.vertices.mesh_ );
        }

        region_cell_ptr_.resize(
            gm_.nb_regions() * to_underlying_type( CellType::UNDEFINED ) + 1, 0 );

        // Total number of  cells
        std::vector< index_t > nb_cells_per_type(
            to_underlying_type( CellType::UNDEFINED ), 0 );
        index_t nb = 0;

        for( index_t r = 0; r < gm_.nb_regions(); ++r ) {
            nb += gm_.region( r ).nb_mesh_elements();
        }

        // Get out if no cells
        if( nb == 0 ) {
            return;
        }

        // Compute the number of cell per type and per region
        for( index_t r = 0; r < gm_.nb_regions(); ++r ) {
            const Region& cur_region = gm_.region( r );
            for( index_t c = 0; c < gm_.region( r ).nb_mesh_elements(); ++c ) {
                CellType cur_cell_type = cur_region.cell_type( c );
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
        for( index_t t = to_underlying_type( CellType::TETRAHEDRON ) + 1;
            t < to_underlying_type( CellType::UNDEFINED ); t++ ) {
            cells_offset_per_type[t] += cells_offset_per_type[t - 1];
            cells_offset_per_type[t] += nb_cells_per_type[t - 1];
        }

        for( index_t i = 1; i < region_cell_ptr_.size() - 1; i++ ) {
            region_cell_ptr_[i + 1] += region_cell_ptr_[i];
        }

        // Create "empty" tet, hex, pyr and prism
        for( index_t i = 0; i < to_underlying_type( CellType::UNDEFINED ); ++i ) {
            mesh_builder->create_cells( nb_cells_per_type[i], CellType( i ) );
        }

        // Fill the cells with vertices
        bind_attribute();
        std::vector< index_t > cur_cell_per_type(
            to_underlying_type( CellType::UNDEFINED ), 0 );
        const GeoModelMeshVertices& geomodel_vertices = gmm_.vertices;
        for( index_t r = 0; r < gm_.nb_regions(); ++r ) {
            const Region& cur_region = gm_.region( r );
            for( index_t c = 0; c < cur_region.nb_mesh_elements(); ++c ) {
                CellType cur_cell_type = cur_region.cell_type( c );
                index_t cur_cell = cells_offset_per_type[to_underlying_type(
                    cur_cell_type )]
                    + cur_cell_per_type[to_underlying_type( cur_cell_type )]++;
                for( index_t v = 0; v < mesh_->nb_cell_vertices( cur_cell ); v++ ) {
                    index_t region_vertex_index =
                        cur_region.mesh_element_vertex_index( c, v );
                    index_t global_vertex_id = geomodel_vertices.geomodel_vertex_id(
                        cur_region.gmme(), region_vertex_index );
                    mesh_builder->set_cell_vertex( cur_cell, v, global_vertex_id );
                }
                region_id_[cur_cell] = r;
                cell_id_[cur_cell] = c;
            }
        }

        // Retrieve the adjacencies
        mesh_builder->connect_cells();

        // Permute cells to sort them per region and per type
        std::vector< index_t > sorted_indices( mesh_->nb_cells() );
        for( index_t i = 0; i < mesh_->nb_cells(); i++ ) {
            sorted_indices[i] = i;
        }
        GeoModelMeshCellsSort action( *mesh_, region_id_ );
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

    void GeoModelMeshCells::bind_attribute()
    {
        if( !region_id_.is_bound() ) {
            region_id_.bind( attribute_manager(), region_att_name );
        }
        if( !cell_id_.is_bound() ) {
            cell_id_.bind( attribute_manager(), cell_region_att_name );
        }
    }

    void GeoModelMeshCells::unbind_attribute()
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

    index_t GeoModelMeshCells::nb() const
    {
        test_and_initialize();
        return mesh_->nb_cells();
    }

    index_t GeoModelMeshCells::nb_vertices( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return mesh_->nb_cell_vertices( c );
    }

    index_t GeoModelMeshCells::vertex( index_t c, index_t v ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        ringmesh_assert( v < mesh_->nb_cell_vertices( c ) );
        return mesh_->cell_vertex( c, v );
    }

    index_t GeoModelMeshCells::nb_edges( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return mesh_->nb_cell_edges( c );
    }

    index_t GeoModelMeshCells::nb_facets( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return mesh_->nb_cell_facets( c );
    }

    index_t GeoModelMeshCells::nb_facet_vertices( index_t c, index_t lf ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        ringmesh_assert( lf < mesh_->nb_cell_facets( c ) );
        return mesh_->nb_cell_facet_vertices( c, lf );
    }

    index_t GeoModelMeshCells::facet_vertex(
        index_t c,
        index_t lf,
        index_t lv ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        ringmesh_assert( lf < mesh_->nb_cell_facets( c ) );
        return mesh_->cell_facet_vertex( c, lf, lv );
    }

    index_t GeoModelMeshCells::edge_vertex( index_t c, index_t le, index_t lv ) const
    {
        geo_debug_assert( le < nb_edges( c ) );
        geo_debug_assert( lv < 2 );
        return mesh_->cell_edge_vertex( c, le, lv );
    }

    index_t GeoModelMeshCells::adjacent( index_t c, index_t f ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        ringmesh_assert( f < mesh_->nb_cell_facets( c ) );
        return mesh_->cell_adjacent( c, f );
    }

    index_t GeoModelMeshCells::region( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return region_id_[c];
    }

    index_t GeoModelMeshCells::index_in_region( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return cell_id_[c];
    }

    CellType GeoModelMeshCells::type( index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( c < mesh_->nb_cells() );
        return mesh_->cell_type( c );
    }

    index_t GeoModelMeshCells::nb_cells( CellType type ) const
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

    index_t GeoModelMeshCells::nb_cells( index_t r, CellType type ) const
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
                            * r] == gm_.region( r ).nb_mesh_elements() );
                return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                    * ( r + 1 )]
                    - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r];
            default:
                ringmesh_assert_not_reached;
                return 0;
        }
    }

    index_t GeoModelMeshCells::cell( index_t r, index_t c, CellType type ) const
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

    index_t GeoModelMeshCells::nb_tet() const
    {
        test_and_initialize();
        return nb_tet_;
    }

    index_t GeoModelMeshCells::nb_tet( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::TETRAHEDRON ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::TETRAHEDRON )];
    }

    index_t GeoModelMeshCells::tet( index_t r, index_t t ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::TETRAHEDRON )] + t;
    }

    index_t GeoModelMeshCells::nb_hex() const
    {
        test_and_initialize();
        return nb_hex_;
    }

    index_t GeoModelMeshCells::nb_hex( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::HEXAHEDRON ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::HEXAHEDRON )];
    }

    index_t GeoModelMeshCells::hex( index_t r, index_t h ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::HEXAHEDRON )] + h;
    }

    index_t GeoModelMeshCells::nb_prism() const
    {
        test_and_initialize();
        return nb_prism_;
    }

    index_t GeoModelMeshCells::nb_prism( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::PRISM ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::PRISM )];
    }

    index_t GeoModelMeshCells::prism( index_t r, index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::PRISM )] + p;
    }

    index_t GeoModelMeshCells::nb_pyramid() const
    {
        test_and_initialize();
        return nb_pyramid_;
    }

    index_t GeoModelMeshCells::nb_pyramid( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::PYRAMID ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::PYRAMID )];
    }

    index_t GeoModelMeshCells::pyramid( index_t r, index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::PYRAMID )] + p;
    }

    index_t GeoModelMeshCells::nb_connector() const
    {
        test_and_initialize();
        return nb_connector_;
    }

    index_t GeoModelMeshCells::nb_connector( index_t r ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + ( to_underlying_type( CellType::UNCLASSIFIED ) + 1 )]
            - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
                + to_underlying_type( CellType::UNCLASSIFIED )];
    }

    index_t GeoModelMeshCells::connector( index_t r, index_t c ) const
    {
        test_and_initialize();
        ringmesh_assert( r < gm_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED ) * r
            + to_underlying_type( CellType::UNCLASSIFIED )] + c;
    }

    bool GeoModelMeshCells::is_duplication_initialized() const
    {
        return mode_ == gmm_.duplicate_mode();
    }

    void GeoModelMeshCells::test_and_initialize_duplication() const
    {
        if( !is_duplication_initialized() ) {
            const_cast< GeoModelMeshCells* >( this )->initialize_duplication();
        }
    }

    void GeoModelMeshCells::initialize_duplication()
    {
        test_and_initialize();
        clear_duplication();

        /// 1. Get all the corner vertices (a lot of duplicated vertices)
        std::vector< vec3 > corner_vertices(
            mesh_->cell_end( mesh_->nb_cells() - 1 ) );
        for( index_t c = 0; c < mesh_->nb_cells(); c++ ) {
            index_t begin = mesh_->cell_begin( c );
            for( index_t v = 0; v < mesh_->nb_cell_vertices( c ); v++ ) {
                corner_vertices[begin + v] = mesh_->vertex(
                    mesh_->cell_vertex( c, v ) );
            }
        }

        /// 2. Tag all corners to duplicate (vertices on a surface to duplicate)
        std::vector< ActionOnSurface > actions_on_surfaces( gm_.nb_surfaces(),
            SKIP );
        std::vector< bool > is_vertex_to_duplicate( corner_vertices.size(), false );
        {
            NNSearch nn_search( corner_vertices, false );
            for( index_t s = 0; s < gm_.nb_surfaces(); s++ ) {
                if( !is_surface_to_duplicate( s ) ) continue;
                actions_on_surfaces[s] = TO_PROCESS;
                const Surface& surface = gm_.surface( s );
                for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
                    std::vector< index_t > colocated_corners =
                        nn_search.get_neighbors( surface.vertex( v ),
                            gm_.epsilon() );
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
        gmm_.polygons.test_and_initialize();
        for( index_t c = 0; c < mesh_->nb_cells(); c++ ) {
            for( index_t v = 0; v < mesh_->nb_cell_vertices( c ); v++ ) {
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
                    std::vector< index_t > facets;
                    cell_facets_around_vertex( *mesh_, cur_c, vertex_id, facets );
                    for( index_t cur_f = 0; cur_f < facets.size(); cur_f++ ) {
                        // Find if the facet is on a surface or inside the domain
                        index_t polygon = NO_ID;
                        bool side;
                        if( is_cell_facet_on_surface( cur_c, cur_f, polygon,
                            side ) ) {
                            index_t surface_id = gmm_.polygons.surface( polygon );
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
                        gmm_.vertices.nb()
                            + static_cast< index_t >( duplicated_vertex_indices_.size() );
                    duplicated_vertex_indices_.push_back( vertex_id );

                    // Update all the cell corners on this side of the surface
                    // to the new duplicated vertex index
                    std::unique_ptr< VolumeMeshBuilder > mesh_builder =
                        VolumeMeshBuilder::create_builder( *mesh_ );
                    for( index_t cur_co : corner_used ) {
                        mesh_builder->set_cell_corner_vertex_index( cur_co,
                            duplicated_vertex_id );
                    }
                }
            }
        }
        mode_ = gmm_.duplicate_mode();
    }

    bool GeoModelMeshCells::is_cell_facet_on_surface(
        index_t c,
        index_t f,
        index_t& polygon,
        bool& side ) const
    {
        test_and_initialize_cell_facet();
        polygon = polygon_id_[mesh_->cell_facet( c, f )];
        if( polygon != NO_ID ) {
            vec3 facet_normal = gmm_.polygons.normal( polygon );
            vec3 cell_facet_normal = mesh_->cell_facet_normal( c, f );
            side = dot( facet_normal, cell_facet_normal ) > 0;
        }
        return polygon != NO_ID;
    }

    bool GeoModelMeshCells::are_corners_to_duplicate(
        const std::vector< action_on_surface >& surfaces,
        std::vector< ActionOnSurface >& info )
    {
        // Temporary vector, it is equal to surfaces
        std::vector< action_on_surface > temp_surfaces;
        temp_surfaces.reserve( temp_surfaces.size() );

        // Find if a free border was found, if so we encountered the
        // two sides of the surface during the propagation around a vertex
        for( index_t i = 1; i < surfaces.size(); i++ ) {
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

    bool GeoModelMeshCells::is_surface_to_duplicate( index_t surface_id ) const
    {
        const Surface& cur_surface = gm_.surface( surface_id );
        if( cur_surface.is_on_voi() ) return false;
        switch( gmm_.duplicate_mode() ) {
            case ALL:
                return true;
            case FAULT: {
                gmge_id parent_interface = cur_surface.parent_gmge(
                    Interface::type_name_static() );
                if( parent_interface.is_defined() ) {
                    GeoModelGeologicalEntity::GEOL_FEATURE feature =
                        gm_.geological_entity( parent_interface ).geological_feature();
                    return GeoModelGeologicalEntity::is_fault( feature );
                }
                return false;
            }
            default:
                return false;
        }
        return false;
    }

    index_t GeoModelMeshCells::nb_duplicated_vertices() const
    {
        test_and_initialize_duplication();
        return static_cast< index_t >( duplicated_vertex_indices_.size() );
    }

    index_t GeoModelMeshCells::nb_total_vertices() const
    {
        return nb_duplicated_vertices() + mesh_->nb_vertices();
    }

    bool GeoModelMeshCells::is_corner_duplicated(
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

    index_t GeoModelMeshCells::duplicated_vertex(
        index_t duplicate_vertex_index ) const
    {
        test_and_initialize_duplication();
        ringmesh_assert(
            duplicate_vertex_index < duplicated_vertex_indices_.size() );
        return duplicated_vertex_indices_[duplicate_vertex_index];
    }

    void GeoModelMeshCells::clear()
    {
        std::unique_ptr< VolumeMeshBuilder > mesh_builder =
            VolumeMeshBuilder::create_builder( *mesh_ );
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

    void GeoModelMeshCells::clear_duplication()
    {
        std::unique_ptr< VolumeMeshBuilder > mesh_builder =
            VolumeMeshBuilder::create_builder( *mesh_ );
        for( index_t c = 0; c < mesh_->nb_cells(); c++ ) {
            for( index_t v = 0; v < mesh_->nb_cell_vertices( c ); v++ ) {
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

    void GeoModelMeshCells::test_and_initialize_cell_facet() const
    {
        if( !polygon_id_.is_bound() ) {
            const_cast< GeoModelMeshCells* >( this )->initialize_cell_facet();
        }
    }

    void GeoModelMeshCells::initialize_cell_facet()
    {
        gmm_.polygons.test_and_initialize();

        polygon_id_.bind( mesh_->cell_facet_attribute_manager(), "polygon_id" );
        polygon_id_.fill( NO_ID );
        const NNSearch& nn_search = gmm_.polygons.nn_search();
        for( index_t c = 0; c < mesh_->nb_cells(); c++ ) {
            for( index_t f = 0; f < mesh_->nb_cell_facets( c ); f++ ) {
                std::vector< index_t > result = nn_search.get_neighbors(
                    mesh_->cell_facet_barycenter( c, f ), gm_.epsilon() );
                if( !result.empty() ) {
                    polygon_id_[mesh_->cell_facet( c, f )] = result[0];
                    // If there are more than 1 matching facet, this is WRONG
                    // and the vertex indices should be checked too [Jeanne]
                    ringmesh_assert( result.size() == 1 );
                }
            }
        }
    }

    vec3 GeoModelMeshCells::barycenter( index_t c ) const
    {
        test_and_initialize();
        return mesh_->cell_barycenter( c );
    }

    double GeoModelMeshCells::volume( index_t c ) const
    {
        test_and_initialize();
        return mesh_->cell_volume( c );
    }

    const VolumeAABBTree& GeoModelMeshCells::aabb() const
    {
        test_and_initialize();
        return mesh_->cell_aabb();
    }

    /*******************************************************************************/

    GeoModelMeshPolygons::GeoModelMeshPolygons(
        GeoModelMesh& gmm,
        GeoModel& gm,
        std::unique_ptr< SurfaceMesh >& mesh )
        :
            GeoModelMeshBase( gmm, gm ),
            mesh_( mesh ),
            nb_triangle_( 0 ),
            nb_quad_( 0 ),
            nb_unclassified_polygon_( 0 )
    {
        set_mesh( mesh_.get() );
    }

    GeoModelMeshPolygons::~GeoModelMeshPolygons()
    {
        unbind_attribute();
    }

    void GeoModelMeshPolygons::bind_attribute()
    {
        if( !surface_id_.is_bound() ) {
            surface_id_.bind( attribute_manager(), surface_att_name );
        }
        if( !polygon_id_.is_bound() ) {
            polygon_id_.bind( attribute_manager(), polygon_surface_att_name );
        }

    }

    void GeoModelMeshPolygons::unbind_attribute()
    {
        if( surface_id_.is_bound() ) {
            surface_id_.unbind();
        }
        if( polygon_id_.is_bound() ) {
            polygon_id_.unbind();
        }

    }

    bool GeoModelMeshPolygons::is_initialized() const
    {
        return mesh_->nb_polygons() > 0;
    }

    index_t GeoModelMeshPolygons::nb() const
    {
        test_and_initialize();
        return mesh_->nb_polygons();
    }

    index_t GeoModelMeshPolygons::nb_vertices( index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        return mesh_->nb_polygon_vertices( p );
    }

    index_t GeoModelMeshPolygons::vertex( index_t p, index_t v ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        ringmesh_assert( v < mesh_->nb_polygon_vertices( p ) );
        return mesh_->polygon_vertex( p, v );
    }

    index_t GeoModelMeshPolygons::adjacent( index_t p, index_t e ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        ringmesh_assert( e < mesh_->nb_polygon_vertices( p ) );
        return mesh_->polygon_adjacent( p, e );
    }

    index_t GeoModelMeshPolygons::surface( index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        return surface_id_[p];
    }

    index_t GeoModelMeshPolygons::index_in_surface( index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        return polygon_id_[p];
    }

    PolygonType GeoModelMeshPolygons::type( index_t p, index_t& index ) const
    {
        test_and_initialize();
        ringmesh_assert( p < mesh_->nb_polygons() );
        index_t polygon = index_in_surface( p );
        index_t s = surface( p );
        for( index_t t = to_underlying_type( PolygonType::TRIANGLE );
            t < to_underlying_type( PolygonType::UNDEFINED ); t++ ) {
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

    index_t GeoModelMeshPolygons::nb_polygons( PolygonType type ) const
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

    index_t GeoModelMeshPolygons::nb_polygons( index_t s, PolygonType type ) const
    {
        test_and_initialize();
        ringmesh_assert( s < gm_.nb_surfaces() );
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

    index_t GeoModelMeshPolygons::polygon(
        index_t s,
        index_t p,
        PolygonType type ) const
    {
        test_and_initialize();
        ringmesh_assert( s < gm_.nb_surfaces() );
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

    index_t GeoModelMeshPolygons::nb_triangle() const
    {
        test_and_initialize();
        return nb_triangle_;
    }

    index_t GeoModelMeshPolygons::nb_triangle( index_t s ) const
    {
        test_and_initialize();
        ringmesh_assert( s < gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + ( to_underlying_type( PolygonType::TRIANGLE ) + 1 )]
            - surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
                + to_underlying_type( PolygonType::TRIANGLE )];
    }

    index_t GeoModelMeshPolygons::triangle( index_t s, index_t t ) const
    {
        test_and_initialize();
        ringmesh_assert( s < gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + to_underlying_type( PolygonType::TRIANGLE )] + t;
    }

    index_t GeoModelMeshPolygons::nb_quad() const
    {
        test_and_initialize();
        return nb_quad_;
    }

    index_t GeoModelMeshPolygons::nb_quad( index_t s ) const
    {
        test_and_initialize();
        ringmesh_assert( s < gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + ( to_underlying_type( PolygonType::QUAD ) + 1 )]
            - surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
                + to_underlying_type( PolygonType::QUAD )];
    }

    index_t GeoModelMeshPolygons::quad( index_t s, index_t q ) const
    {
        test_and_initialize();
        ringmesh_assert( s < gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + to_underlying_type( PolygonType::QUAD )] + q;
    }

    index_t GeoModelMeshPolygons::nb_unclassified_polygon() const
    {
        test_and_initialize();
        return nb_unclassified_polygon_;
    }

    index_t GeoModelMeshPolygons::nb_unclassified_polygon( index_t s ) const
    {
        test_and_initialize();
        ringmesh_assert( s < gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + ( to_underlying_type( PolygonType::UNCLASSIFIED ) + 1 )]
            - surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
                + to_underlying_type( PolygonType::UNCLASSIFIED )];
    }

    index_t GeoModelMeshPolygons::unclassified_polygon( index_t s, index_t p ) const
    {
        test_and_initialize();
        ringmesh_assert( s < gm_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
            + to_underlying_type( PolygonType::UNCLASSIFIED )] + p;
    }

    void GeoModelMeshPolygons::clear()
    {
        surface_polygon_ptr_.clear();
        nb_triangle_ = 0;
        nb_quad_ = 0;
        std::unique_ptr< SurfaceMeshBuilder > mesh_builder =
            SurfaceMeshBuilder::create_builder( *mesh_ );
        mesh_builder->clear( true, false );
    }

    void GeoModelMeshPolygons::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshPolygons* >( this )->initialize();
        }
    }

    void GeoModelMeshPolygons::initialize()
    {
        gmm_.vertices.test_and_initialize();
        clear();
        surface_polygon_ptr_.resize(
            gm_.nb_surfaces() * to_underlying_type( PolygonType::UNDEFINED ) + 1,
            0 );
        std::unique_ptr< SurfaceMeshBuilder > mesh_builder =
            SurfaceMeshBuilder::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != gmm_.vertices.nb() ) {
            copy_vertices( mesh_builder.get(), *gmm_.vertices.mesh_ );
        }

        // Compute the total number of polygons per type and per surface
        std::map< PolygonType, index_t > nb_polygon_per_type = {
            { PolygonType::TRIANGLE, 0 }, { PolygonType::QUAD, 0 }, {
                PolygonType::UNCLASSIFIED, 0 } };
        for( index_t s = 0; s < gm_.nb_surfaces(); s++ ) {
            const Surface& surface = gm_.surface( s );
            if( surface.is_simplicial() ) {
                nb_polygon_per_type[PolygonType::TRIANGLE] +=
                    surface.nb_mesh_elements();
                surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED ) * s
                    + to_underlying_type( PolygonType::TRIANGLE ) + 1] +=
                    surface.nb_mesh_elements();
            } else {
                for( index_t p = 0; p < surface.nb_mesh_elements(); p++ ) {
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

        for( index_t i = 1; i < surface_polygon_ptr_.size() - 1; i++ ) {
            surface_polygon_ptr_[i + 1] += surface_polygon_ptr_[i];
        }

        // Fill the triangles and quads created above
        // Create and fill polygons
        bind_attribute();
        const GeoModelMeshVertices& geomodel_vertices = gmm_.vertices;
        std::vector< index_t > cur_polygon_per_type(
            to_underlying_type( PolygonType::UNDEFINED ), 0 );
        for( index_t s = 0; s < gm_.nb_surfaces(); s++ ) {
            const Surface& surface = gm_.surface( s );
            gmme_id surface_id = surface.gmme();
            for( index_t p = 0; p < surface.nb_mesh_elements(); p++ ) {
                index_t nb_vertices = surface.nb_mesh_element_vertices( p );
                index_t cur_polygon = NO_ID;
                if( nb_vertices < 5 ) {
                    PolygonType T = static_cast< PolygonType >( nb_vertices - 3 );
                    cur_polygon = polygon_offset_per_type[T]
                        + cur_polygon_per_type[to_underlying_type( T )]++;
                    for( index_t v = 0; v < nb_vertices; v++ ) {
                        index_t v_id = geomodel_vertices.geomodel_vertex_id(
                            surface_id, p, v );
                        ringmesh_assert( v_id != NO_ID );
                        mesh_builder->set_polygon_vertex( cur_polygon, v, v_id );
                    }
                } else {
                    std::vector< index_t > vertices( nb_vertices );
                    for( index_t v = 0; v < nb_vertices; v++ ) {
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
        for( index_t i = 0; i < mesh_->nb_polygons(); i++ ) {
            sorted_indices[i] = i;
        }
        GeoModelMeshPolygonsSort action( *mesh_, surface_id_ );
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

    void GeoModelMeshPolygons::disconnect_along_lines()
    {
        std::unique_ptr< SurfaceMeshBuilder > mesh_builder =
            SurfaceMeshBuilder::create_builder( *mesh_ );
        for( index_t s = 0; s < gm_.nb_surfaces(); s++ ) {
            const Surface& surface = gm_.surface( s );
            for( index_t p = 0; p < nb_polygons( s ); p++ ) {
                index_t polygon_id = polygon( s, p );
                index_t surface_polygon_id = index_in_surface( polygon_id );
                for( index_t v = 0; v < nb_vertices( polygon_id ); v++ ) {
                    index_t adj = surface.polygon_adjacent_index( surface_polygon_id,
                        v );
                    if( adj == NO_ID ) {
                        mesh_builder->set_polygon_adjacent( polygon_id, v, NO_ID );
                    }
                }
            }
        }
    }

    vec3 GeoModelMeshPolygons::center( index_t p ) const
    {
        test_and_initialize();
        return mesh_->polygon_barycenter( p );
    }

    vec3 GeoModelMeshPolygons::normal( index_t p ) const
    {
        test_and_initialize();
        return mesh_->polygon_normal( p );
    }

    double GeoModelMeshPolygons::area( index_t p ) const
    {
        test_and_initialize();
        return mesh_->polygon_area( p );
    }

    const SurfaceAABBTree& GeoModelMeshPolygons::aabb() const
    {
        test_and_initialize();
        return mesh_->polygon_aabb();
    }
    /*******************************************************************************/

    GeoModelMeshEdges::GeoModelMeshEdges(
        GeoModelMesh& gmm,
        GeoModel& gm,
        std::unique_ptr< LineMesh >& mesh )
        : GeoModelMeshBase( gmm, gm ), mesh_( mesh )
    {
        set_mesh( mesh_.get() );
    }

    GeoModelMeshEdges::~GeoModelMeshEdges()
    {
    }

    index_t GeoModelMeshEdges::nb_wells() const
    {
        test_and_initialize();
        return gm_.wells() ? gm_.wells()->nb_wells() : 0;
    }

    index_t GeoModelMeshEdges::nb_edges() const
    {
        test_and_initialize();
        return mesh_->nb_edges();
    }

    index_t GeoModelMeshEdges::nb_edges( index_t w ) const
    {
        test_and_initialize();
        return well_ptr_[w + 1] - well_ptr_[w];
    }

    index_t GeoModelMeshEdges::vertex( index_t w, index_t e, index_t v ) const
    {
        test_and_initialize();
        return mesh_->edge_vertex( well_ptr_[w] + e, v );
    }

    void GeoModelMeshEdges::clear()
    {
        std::unique_ptr< LineMeshBuilder > mesh_builder =
            LineMeshBuilder::create_builder( *mesh_ );
        mesh_builder->clear( true, false );
        well_ptr_.clear();
    }

    bool GeoModelMeshEdges::is_initialized() const
    {
        return mesh_->nb_edges() > 0;
    }

    void GeoModelMeshEdges::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshEdges* >( this )->initialize();
        }
    }

    void GeoModelMeshEdges::initialize()
    {
        if( !gm_.wells() ) return;
        gmm_.vertices.test_and_initialize();
        clear();
        std::unique_ptr< LineMeshBuilder > mesh_builder =
            LineMeshBuilder::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != gmm_.vertices.nb() ) {
            copy_vertices( mesh_builder.get(), *gmm_.vertices.mesh_ );
        }

        // Compute the total number of edge per well
        const WellGroup& wells = *gm_.wells();
        well_ptr_.resize( wells.nb_wells() + 1, 0 );
        index_t nb_edges = 0;
        for( index_t w = 0; w < wells.nb_wells(); w++ ) {
            nb_edges += wells.well( w ).nb_edges();
            well_ptr_[w + 1] = nb_edges;
        }

        // Compute the edge offset
        for( index_t i = 1; i < well_ptr_.size() - 1; i++ ) {
            well_ptr_[i + 1] += well_ptr_[i];
        }

        // Create edges
        mesh_builder->create_edges( well_ptr_.back() );

        // Fill edges
        index_t cur_edge = 0;
        for( index_t w = 0; w < wells.nb_wells(); w++ ) {
            const Well& well = wells.well( w );
            for( index_t p = 0; p < well.nb_parts(); p++ ) {
                for( index_t e = 0; e < well.part( p ).nb_edges(); e++ ) {
                    const vec3& e0 = well.part( p ).edge_vertex( e, 0 );
                    mesh_builder->set_edge_vertex( cur_edge, 0,
                        gmm_.vertices.index( e0 ) );
                    const vec3& e1 = well.part( p ).edge_vertex( e, 1 );
                    mesh_builder->set_edge_vertex( cur_edge, 1,
                        gmm_.vertices.index( e1 ) );
                    cur_edge++;
                }
            }
        }
    }

    const LineAABBTree& GeoModelMeshEdges::aabb() const
    {
        test_and_initialize();
        return mesh_->edge_aabb();
    }

    /*******************************************************************************/

    GeoModelMesh::GeoModelMesh( GeoModel& geomodel )
        :
            geomodel_( geomodel ),
            mode_( GeoModelMeshCells::NONE ),
            vertices( *this, geomodel, mesh_set_.point_set_mesh ),
            edges( *this, geomodel, mesh_set_.line_mesh ),
            polygons( *this, geomodel, mesh_set_.surface_mesh ),
            cells( *this, geomodel, mesh_set_.volume_mesh )
    {
    }

    void GeoModelMesh::transfer_attributes_from_gmm_to_gm_regions() const
    {
        transfer_vertex_attributes_from_gmm_to_gm_regions();
        transfer_cell_attributes_from_gmm_to_gm_regions();
    }

    void GeoModelMesh::transfer_attributes_from_gm_regions_to_gmm() const
    {
        transfer_vertex_attributes_from_gm_regions_to_gmm() ;
        transfer_cell_attributes_from_gm_regions_to_gmm() ;
    }

    void GeoModelMesh::transfer_vertex_attributes_from_gmm_to_gm_regions() const
    {
        GEO::vector< std::string > att_v_names;
        std::vector< std::string > att_v_double_names;
        vertices.attribute_manager().list_attribute_names( att_v_names );
        for( index_t att_v = 0; att_v < vertices.attribute_manager().nb();
            att_v++ ) {

            if( !GEO::Attribute< double >::is_defined( vertices.attribute_manager(),
                att_v_names[att_v] ) ) {
                continue;
            }
            att_v_double_names.push_back( att_v_names[att_v] );
            for( index_t reg = 0; reg < geomodel_.nb_regions(); reg++ ) {

                if( geomodel_.region( reg ).vertex_attribute_manager().is_defined(
                    att_v_names[att_v] ) ) {
                    Logger::warn( "Transfer attribute", "The attribute ",
                        att_v_names[att_v], " already exists on the region ", reg );
                    continue;
                }
                GEO::Attribute< double > cur_v_att;
                cur_v_att.create_vector_attribute(
                    geomodel_.region( reg ).vertex_attribute_manager(),
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

            for( index_t reg = 0; reg < geomodel_.nb_regions(); reg++ ) {
                att_on_regions.bind_one_attribute( reg,
                    geomodel_.region( reg ).vertex_attribute_manager(), att_v );
            }

            for( index_t v = 0; v < vertices.nb(); v++ ) {
                std::vector< GMEVertex > vertices_on_geomodel_region =
                    vertices.gme_type_vertices( Region::type_name_static(), v );
                for( const GMEVertex& cur_vertex_on_geomodel : vertices_on_geomodel_region ) {
                    for( index_t att_e = 0; att_e < att_dim; att_e++ ) {
                        att_on_regions[cur_vertex_on_geomodel.gmme.index()][cur_vertex_on_geomodel.v_index
                            * att_dim + att_e] = cur_att_on_geomodelmesh[v * att_dim
                            + att_e];
                    }
                }

            }
        }
    }

    void GeoModelMesh::transfer_vertex_attributes_from_gm_regions_to_gmm() const
    {
        for( index_t reg_itr = 0; reg_itr < geomodel().nb_regions(); ++reg_itr ) {
            GEO::vector< std::string > att_v_names;
            const Region& cur_reg = geomodel().region( reg_itr );
            GEO::AttributesManager& reg_vertex_attr_mgr =
                cur_reg.vertex_attribute_manager();
            reg_vertex_attr_mgr.list_attribute_names( att_v_names );
            for( const std::string& cur_attr_name : att_v_names ) {

                // It is not necessary to copy the coordinates. There are already there.
                if( cur_attr_name == "point" ) {
                    continue;
                }

                const GEO::AttributeStore* cur_v_att_store_in_reg =
                    reg_vertex_attr_mgr.find_attribute_store( cur_attr_name );
                ringmesh_assert( cur_v_att_store_in_reg != nullptr );
                index_t dim = cur_v_att_store_in_reg->dimension();
                GEO::AttributeStore* cur_v_att_store = nullptr;
                if( !vertices.attribute_manager().is_defined( cur_attr_name ) ) {
                    const std::string cur_type_name =
                        GEO::AttributeStore::element_type_name_by_element_typeid_name(
                            cur_v_att_store_in_reg->element_typeid_name() );
                    ringmesh_assert(
                        GEO::AttributeStore::element_type_name_is_known(
                            cur_type_name ) );
                    cur_v_att_store =
                        GEO::AttributeStore::create_attribute_store_by_element_type_name(
                            cur_type_name, dim );
                    vertices.attribute_manager().bind_attribute_store( cur_attr_name,
                        cur_v_att_store );
                } else {
                    cur_v_att_store =
                        vertices.attribute_manager().find_attribute_store(
                            cur_attr_name );
                }
                ringmesh_assert( cur_v_att_store != nullptr );
                for( index_t v_in_reg_itr = 0; v_in_reg_itr < cur_reg.nb_vertices();
                    ++v_in_reg_itr ) {
                    index_t v_id_in_gmm = vertices.geomodel_vertex_id(
                        cur_reg.gmme(), v_in_reg_itr );
                    GEO::Memory::copy( (GEO::Memory::pointer) cur_v_att_store->data() + v_id_in_gmm * dim,
                        (GEO::Memory::pointer) cur_v_att_store_in_reg->data() + v_in_reg_itr * dim, dim);
                    /*for( index_t dim_itr = 0; dim_itr < dim; ++dim_itr ) {
                     cur_v_att[v_id_in_gmm * dim + dim_itr] =
                     cur_v_att_in_reg[v_in_reg_itr * dim + dim_itr] ;
                     }*/
                }
            }
        }
    }

    void GeoModelMesh::transfer_cell_attributes_from_gmm_to_gm_regions() const
    {

        GEO::vector< std::string > att_c_names;
        cells.attribute_manager().list_attribute_names( att_c_names );

        const NNSearch& nn_search = cells.cell_nn_search();

        for( const std::string& att_c : att_c_names ) {
            if( !GEO::Attribute< double >::is_defined( cells.attribute_manager(),
                att_c ) ) {
                continue;
            }
            GEO::Attribute< double > cur_att_on_geomodel_mesh(
                cells.attribute_manager(), att_c );
            index_t att_dim = cur_att_on_geomodel_mesh.dimension();

            for( index_t reg = 0; reg < geomodel_.nb_regions(); reg++ ) {
                if( geomodel_.region( reg ).cell_attribute_manager().is_defined(
                    att_c ) ) {
                    Logger::warn( "Transfer attribute", "The attribute ", att_c,
                        " already exists on the region ", reg );
                    continue;
                }
                GEO::Attribute< double > cur_att_on_geomodel_mesh_entity;
                cur_att_on_geomodel_mesh_entity.create_vector_attribute(
                    geomodel_.region( reg ).cell_attribute_manager(), att_c,
                    att_dim );
                for( index_t c = 0; c < geomodel_.region( reg ).nb_mesh_elements();
                    c++ ) {
                    vec3 center = geomodel_.region( reg ).mesh_element_barycenter(
                        c );
                    std::vector< index_t > c_in_geom_model_mesh =
                        nn_search.get_neighbors( center, geomodel_.epsilon() );
                    ringmesh_assert( c_in_geom_model_mesh.size() == 1 );
                    for( index_t att_e = 0; att_e < att_dim; att_e++ ) {
                        cur_att_on_geomodel_mesh_entity[c * att_dim + att_e] =
                            cur_att_on_geomodel_mesh[c_in_geom_model_mesh[0]
                                * att_dim + att_e];
                    }
                }
            }
        }
    }

    void GeoModelMesh::transfer_cell_attributes_from_gm_regions_to_gmm() const
    {
        cells.test_and_initialize(); /// TODO to keep?
        const NNSearch& nn_search = cells.cell_nn_search();
        for( index_t reg_itr = 0; reg_itr < geomodel().nb_regions(); ++reg_itr ) {
            GEO::vector< std::string > att_c_names;
            const Region& cur_reg = geomodel().region( reg_itr );
            GEO::AttributesManager& reg_cell_attr_mgr =
                cur_reg.cell_attribute_manager();
            reg_cell_attr_mgr.list_attribute_names( att_c_names );
            for( const std::string& cur_attr_name : att_c_names ) {

                /// TODO code tres similaire que pour les vertices... to merge ?
                const GEO::AttributeStore* cur_c_att_store_in_reg =
                    reg_cell_attr_mgr.find_attribute_store( cur_attr_name );
                ringmesh_assert( cur_c_att_store_in_reg != nullptr );
                index_t dim = cur_c_att_store_in_reg->dimension();
                GEO::AttributeStore* cur_c_att_store = nullptr;
                if( !cells.attribute_manager().is_defined( cur_attr_name ) ) {
                    const std::string cur_type_name =
                        GEO::AttributeStore::element_type_name_by_element_typeid_name(
                            cur_c_att_store_in_reg->element_typeid_name() );
                    ringmesh_assert(
                        GEO::AttributeStore::element_type_name_is_known(
                            cur_type_name ) );
                    cur_c_att_store =
                        GEO::AttributeStore::create_attribute_store_by_element_type_name(
                            cur_type_name, dim );
                    cells.attribute_manager().bind_attribute_store( cur_attr_name,
                        cur_c_att_store );
                } else {
                    cur_c_att_store = cells.attribute_manager().find_attribute_store(
                        cur_attr_name );
                }
                ringmesh_assert( cur_c_att_store != nullptr );
                ringmesh_assert(
                    cur_c_att_store->element_size()
                        == cur_c_att_store_in_reg->element_size() );
                /// TODO put other assert

                for( index_t c_in_reg_itr = 0;
                    c_in_reg_itr < cur_reg.nb_mesh_elements(); ++c_in_reg_itr ) {
                    vec3 center =
                        geomodel_.region( reg_itr ).mesh_element_barycenter(
                            c_in_reg_itr );
                    std::vector< index_t > c_in_geom_model_mesh =
                        nn_search.get_neighbors( center, geomodel_.epsilon() );
                    ringmesh_assert( c_in_geom_model_mesh.size() == 1 );
                    for( index_t dim_itr = 0; dim_itr < dim; ++dim_itr ) {
                        GEO::Memory::copy(
                            (GEO::Memory::pointer) cur_c_att_store->data()
                                + c_in_geom_model_mesh[0] * dim
                                    * cur_c_att_store->element_size(),
                            (GEO::Memory::pointer) cur_c_att_store_in_reg->data()
                                + c_in_reg_itr * dim
                                    * cur_c_att_store_in_reg->element_size(),
                            dim * cur_c_att_store->element_size() );
                        /*cur_c_att[c_in_geom_model_mesh[0] * dim + dim_itr] =
                         cur_c_att_in_reg[c_in_reg_itr * dim + dim_itr];*/
                    }
                }
            }
        }
    }

    GeoModelMesh::~GeoModelMesh()
    {
        polygons.unbind_attribute();
        cells.unbind_attribute();
    }

    void GeoModelMesh::remove_colocated_vertices()
    {
        vertices.remove_colocated();
    }

    void GeoModelMesh::erase_vertices( std::vector< index_t >& to_delete )
    {
        vertices.erase_vertices( to_delete );
    }

    void GeoModelMesh::change_point_set_mesh_data_structure( const MeshType& type )
    {
        if( mesh_set_.point_set_mesh->type_name() != type ) {
            vertices.clear();
            mesh_set_.create_point_set_mesh( type );
        }
    }

    void GeoModelMesh::change_line_mesh_data_structure( const MeshType& type )
    {
        if( mesh_set_.line_mesh->type_name() != type ) {
            edges.clear();
            mesh_set_.create_line_mesh( type );
        }
    }

    void GeoModelMesh::change_surface_mesh_data_structure( const MeshType& type )
    {
        if( mesh_set_.surface_mesh->type_name() != type ) {
            polygons.clear();
            mesh_set_.create_surface_mesh( type );
        }
    }

    void GeoModelMesh::change_volume_mesh_data_structure( const MeshType& type )
    {
        if( mesh_set_.volume_mesh->type_name() != type ) {
            cells.clear();
            mesh_set_.create_volume_mesh( type );
        }
    }

}
