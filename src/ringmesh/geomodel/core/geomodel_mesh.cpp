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

#include <ringmesh/geomodel/core/geomodel_mesh.h>

#include <numeric>
#include <stack>

#include <geogram/basic/algorithm.h>

#include <geogram/basic/permutation.h>
#include <geogram/mesh/mesh_geometry.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/pimpl_impl.h>
#include <ringmesh/geogram_extension/geogram_extension.h>
#include <ringmesh/geogram_extension/geogram_mesh.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/core/well.h>

#include <ringmesh/mesh/mesh_builder.h>

/*!
 * @author Arnaud Botella - Jeanne Pellerin - Antoine Mazuyer
 */

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    class GeoModelMeshPolygonsBaseSort
    {
    public:
        GeoModelMeshPolygonsBaseSort( const SurfaceMesh< DIMENSION >& mesh,
            const std::vector< index_t >& surface_id )
            : mesh_( mesh ), surface_id_( surface_id )
        {
        }

        bool operator()( index_t i, index_t j ) const
        {
            if( surface_id_[i] != surface_id_[j] )
            {
                return surface_id_[i] < surface_id_[j];
            }
            return mesh_.nb_polygon_vertices( i )
                   < mesh_.nb_polygon_vertices( j );
        }

    private:
        const SurfaceMesh< DIMENSION >& mesh_;
        const std::vector< index_t >& surface_id_;
    };

    template < index_t DIMENSION >
    class GeoModelMeshCellsSort
    {
    public:
        GeoModelMeshCellsSort( const VolumeMesh< DIMENSION >& mesh,
            const std::vector< index_t >& region_id )
            : mesh_( mesh ), region_id_( region_id )
        {
        }

        bool operator()( index_t i, index_t j ) const
        {
            if( region_id_[i] != region_id_[j] )
            {
                return region_id_[i] < region_id_[j];
            }
            return mesh_.cell_type( i ) < mesh_.cell_type( j );
        }

    private:
        const VolumeMesh< DIMENSION >& mesh_;
        const std::vector< index_t >& region_id_;
    };

    template < index_t DIMENSION >
    std::vector< index_t > cell_facets_around_vertex(
        const VolumeMesh< DIMENSION >& mesh, index_t cell, index_t vertex_id )
    {
        std::vector< index_t > facets;
        facets.reserve( mesh.nb_cell_facets( cell ) );
        for( auto f : range( mesh.nb_cell_facets( cell ) ) )
        {
            for( auto v : range( mesh.nb_cell_facet_vertices(
                     CellLocalFacet( cell, f ) ) ) )
            {
                if( mesh.cell_facet_vertex( CellLocalFacet( cell, f ), v )
                    == vertex_id )
                {
                    facets.push_back( f );
                    break;
                }
            }
        }
        return facets;
    }

    template < index_t DIMENSION >
    void copy_vertices( MeshBaseBuilder< DIMENSION >* builder,
        const MeshBase< DIMENSION >& mesh )
    {
        builder->clear_vertices( true, false );
        builder->create_vertices( mesh.nb_vertices() );
        for( auto v : range( mesh.nb_vertices() ) )
        {
            builder->set_vertex( v, mesh.vertex( v ) );
        }
    }
} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    GeoModelMeshCommon< DIMENSION >::GeoModelMeshCommon(
        GeoModelMesh< DIMENSION >& gmm, GeoModel< DIMENSION >& geomodel )
        : gmm_( gmm ), geomodel_( geomodel ), mesh_base_( nullptr )
    {
    }
    template < index_t DIMENSION >
    class GeoModelMeshVerticesBase< DIMENSION >::Impl
    {
    public:
        Impl( GeoModelMeshVerticesBase& geomodel_vertices,
            const GeoModel< DIMENSION >& geomodel )
            : geomodel_vertices_( geomodel_vertices ), geomodel_( geomodel )
        {
            vertex_maps_[Corner< DIMENSION >::type_name_static()] =
                &corner_vertex_maps_;
            vertex_maps_[Line< DIMENSION >::type_name_static()] =
                &line_vertex_maps_;
            vertex_maps_[Surface< DIMENSION >::type_name_static()] =
                &surface_vertex_maps_;
        }

        ~Impl() = default;

        /*!
         * \name Query
         * @{
         */

        /*!
         * @brief Returns the index of a GeoModelMeshEntity vertex in the
         * geomodel
         * global indexing
         * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
         * @param[in] mesh_entity_vertex_index Index of query vertex in the
         * GeoModelMeshEntity indexing
         * @returns Model index of the GeoModelMeshEntity vertex
         */
        index_t geomodel_vertex_index( const gmme_id& mesh_entity_id,
            index_t mesh_entity_vertex_index ) const
        {
            ringmesh_assert(
                mesh_entity_vertex_index
                < geomodel_.mesh_entity( mesh_entity_id ).nb_vertices() );

            return vertex_map( mesh_entity_id )[mesh_entity_vertex_index];
        }

        /*!
         * @brief Returns all the corresponding vertices in
         * GeoModelMeshEntities
         * to a given geomodel vertex
         * @param[in] vertex Model vertex index
         * @returns All the corresponding vertices in their local indexing
         */
        const std::vector< GMEVertex >& mesh_entity_vertex_indices(
            index_t v ) const
        {
            ringmesh_assert( v < gme_vertices_.size() );
            return gme_vertices_[v];
        }

        /*!
         * @brief Returns all the corresponding vertices in
         * GeoModelMeshEntities
         * of a specific type to a given geomodel vertex
         * @param[in] vertex Model vertex index
         * @param[in] mesh_entity_type Type of GeoModelMeshEntity
         * @return corresponding vertices in GeoModelMeshEntities
         * of a specific type
         */
        std::vector< GMEVertex > mesh_entity_vertex_indices(
            index_t v, const MeshEntityType& mesh_entity_type ) const
        {
            const auto& all_gmes = mesh_entity_vertex_indices( v );
            std::vector< GMEVertex > result;
            result.reserve( all_gmes.size() );
            for( const auto& vertex : all_gmes )
            {
                if( vertex.gmme.type() == mesh_entity_type )
                {
                    result.push_back( vertex );
                }
            }
            return result;
        }

        /*!
         * @brief Returns all the corresponding vertices to a geomodel
         * vertex
         * in a specific GeoModelMeshEntities
         * @param[in] vertex Model vertex index
         * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
         * @return corresponding vertices in the GeoModelMeshEntity
         * @returns All the corresponding vertices in their local indexing
         */
        std::vector< index_t > mesh_entity_vertex_indices(
            index_t v, const gmme_id& mesh_entity_id ) const
        {
            std::vector< index_t > result;
            auto all_gmes = mesh_entity_vertex_indices( v );
            for( const auto& vertex : all_gmes )
            {
                if( vertex.gmme == mesh_entity_id )
                {
                    result.push_back( vertex.v_index );
                }
            }
            return result;
        }

        std::vector< index_t >& vertex_map(
            const gmme_id& mesh_entity_id ) const
        {
            return ( *vertex_maps_.at(
                mesh_entity_id.type() ) )[mesh_entity_id.index()];
        }

        /*! @}
         * \name Updating
         * @{
         */

        /*!
         * @brief Sets the geomodel vertex mapping value of a given vertex
         * in a GeoModelMeshEntity
         * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
         * @param[in] mesh_entity_vertex_index Index of query vertex in the
         * GeoModelMeshEntity indexing
         * @param[in] geomodel_entity_vertex_index Model vertex index to map
         * with
         */
        void set_vertex_map_value( const gmme_id& mesh_entity_id,
            index_t mesh_entity_vertex_index,
            index_t geomodel_entity_vertex_index ) const
        {
            test_and_initialize_mesh_entity_vertex_map( mesh_entity_id );
            vertex_map( mesh_entity_id )[mesh_entity_vertex_index] =
                geomodel_entity_vertex_index;
        }

        void add_to_gme_vertices(
            const GMEVertex& gme_vertex, index_t geomodel_vertex_index ) const
        {
            gme_vertices_[geomodel_vertex_index].push_back( gme_vertex );
        }

        /*!
         * @brief Updates all the vertex maps with regards to the global
         * indexing
         * changes
         * @param[in] old2new Map between actual geomodel indexing and
         * wanted
         * geomodel indexing. Its size is equal to the number of geomodel
         * vertices.
         */
        void update_mesh_entity_maps_and_gmes(
            const std::vector< index_t >& old2new ) const
        {
            const auto& all_mesh_entity_types =
                geomodel_.entity_type_manager()
                    .mesh_entity_manager.mesh_entity_types();
            for( const auto& cur_entity_type : all_mesh_entity_types )
            {
                for( auto e :
                    range( geomodel_.nb_mesh_entities( cur_entity_type ) ) )
                {
                    const auto& E = geomodel_.mesh_entity( cur_entity_type, e );
                    auto id = E.gmme();
                    for( auto v : range( E.nb_vertices() ) )
                    {
                        auto old_m_id = geomodel_vertex_index( id, v );
                        auto new_m_id = old2new[old_m_id];
                        set_vertex_map_value( id, v, new_m_id );

                        // Merge gme_vertices information
                        if( find( gme_vertices_[new_m_id], GMEVertex( id, v ) )
                            == NO_ID )
                        {
                            gme_vertices_[new_m_id].emplace_back( id, v );
                        }
                    }
                }
            }
        }

        /*! @}
         * \name Initialization
         * @{
         */

        /*!
         * @brief Resizes the GME_Vertex vectors
         * @param[in] nb Size of the vector
         */
        void resize_geomodel_vertex_gmes( const index_t nb ) const
        {
            gme_vertices_.resize( nb );
        }

        /*!
         * @brief Clears and resizes the GME_Vertex vectors
         * @param[in] nb Size of the vector
         */
        void clear_and_resize_geomodel_vertex_gmes( const index_t nb ) const
        {
            gme_vertices_.clear();
            resize_geomodel_vertex_gmes( nb );
        }

        void bind_all_mesh_entity_vertex_maps() const
        {
            const auto& all_mesh_entity_types =
                geomodel_.entity_type_manager()
                    .mesh_entity_manager.mesh_entity_types();
            for( const auto& cur_entity_type : all_mesh_entity_types )
            {
                auto nb_cur_type_entities =
                    geomodel_.nb_mesh_entities( cur_entity_type );
                vertex_maps_.at( cur_entity_type )->clear();
                vertex_maps_.at( cur_entity_type )
                    ->resize( nb_cur_type_entities );
                for( auto e : range( nb_cur_type_entities ) )
                {
                    resize_vertex_map( { cur_entity_type, e } );
                }
            }
        }

        /*! @}
         * \name Clearing
         * @{
         */

        /*!
         * @brief Clears all the information about vertex mapping (vector
         * maps
         * and vectors of GME_Vertices
         */
        void clear() const
        {
            gme_vertices_.clear();
            clear_all_mesh_entity_vertex_map();
        }

        /*!
         * @brief Clears the GME_Vertices about one geomodel vertex
         */
        void clear_geomodel_vertex_gmes( index_t v ) const
        {
            ringmesh_assert( v < gme_vertices_.size() );
            gme_vertices_[v].clear();
        }

        void clear_vertex_map( const gmme_id& mesh_entity_id )
        {
            // This if statement is a quick dirty fix for MacOS X
            // since there is a different behavior of the destructor
            // of std::vector. The problem occurs during the deletion
            // of the mesh entities of a GeoModel (when the GeoModel
            // destructor is called). The size of the mesh entity vector
            // decreases during the mesh enity deletion in MacOS X instead
            // of when all the mesh entities have been deleted (as in Linux
            // or Windows). This fix is temporary and will be removed during
            // the attribute refactoring.
            if( vertex_maps_.at( mesh_entity_id.type() )->empty() )
            {
                resize_all_mesh_entity_vertex_maps( mesh_entity_id.type() );
            }
            ringmesh_assert(
                mesh_entity_id.index()
                < vertex_maps_.at( mesh_entity_id.type() )->size() );
            if( !vertex_maps_.at( mesh_entity_id.type() )
                     ->at( mesh_entity_id.index() )
                     .empty() )
            {
                vertex_maps_.at( mesh_entity_id.type() )
                    ->at( mesh_entity_id.index() )
                    .clear();
            }
        }

        std::vector< index_t >& resize_vertex_map(
            const gmme_id& mesh_entity_id ) const
        {
            ringmesh_assert( mesh_entity_id.index()
                             < vertex_maps_[mesh_entity_id.type()]->size() );
            if( geomodel_vertices_.is_initialized() )
            {
                const auto& mesh_entity =
                    geomodel_.mesh_entity( mesh_entity_id );
                vertex_maps_.at( mesh_entity_id.type() )
                    ->at( mesh_entity_id.index() )
                    .resize( mesh_entity.nb_vertices(), NO_ID );
            }
            return vertex_map( mesh_entity_id );
        }

        /*!
         * @}
         */

    private:
        /*!
         * @brief Initializes the given GeoModelMeshEntity vertex map
         * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
         */
        void initialize_mesh_entity_vertex_map(
            const gmme_id& mesh_entity_id ) const
        {
            auto& mesh_entity_vertex_map = resize_vertex_map( mesh_entity_id );
            const auto& E = geomodel_.mesh_entity( mesh_entity_id );
            for( auto v : range( E.nb_vertices() ) )
            {
                mesh_entity_vertex_map[v] =
                    geomodel_vertices_.nn_search().get_closest_neighbor(
                        E.vertex( v ) );
            }
        }

        /*!
         * @brief Tests if the given GeoModelMeshEntity vertex map is
         * initialized.
         * If not, initializes it.
         * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
         * @return True is the map was initialized, false if not.
         */
        bool test_and_initialize_mesh_entity_vertex_map(
            const gmme_id& mesh_entity_id ) const
        {
            resize_all_mesh_entity_vertex_maps( mesh_entity_id.type() );
            if( !is_mesh_entity_vertex_map_initialized( mesh_entity_id ) )
            {
                initialize_mesh_entity_vertex_map( mesh_entity_id );
                return false;
            }
            return true;
        }

        /*!
         * @brief Tests if the given GeoModelMeshEntity vertex map exists.
         * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
         * @return True is the map exists, false if not.
         */
        bool is_mesh_entity_vertex_map_initialized(
            const gmme_id& mesh_entity_id ) const
        {
            return !vertex_maps_.find( mesh_entity_id.type() )
                        ->second->at( mesh_entity_id.index() )
                        .empty();
        }

        /*!
         * @brief Unbinds all the GeoModelMeshEntity vertex maps
         */
        void clear_all_mesh_entity_vertex_map() const
        {
            for( auto& vertex_map : vertex_maps_ )
            {
                for( auto e : range( vertex_map.second->size() ) )
                {
                    vertex_map.second->at( e ).clear();
                }
                vertex_map.second->clear();
            }
        }

        void resize_all_mesh_entity_vertex_maps(
            const MeshEntityType& type ) const
        {
            vertex_maps_.at( type )->resize(
                geomodel_.nb_mesh_entities( type ) );
        }

        /*!
         * @brief Returns the vertex attribute of a GeoModelMeshEntity
         * @param[in] mesh_entity_id Unique id to a GeoModelMeshEntity
         */
        AttributesManager& mesh_entity_vertex_attribute_manager(
            const gmme_id& mesh_entity_id ) const
        {
            const auto& mesh_entity = geomodel_.mesh_entity( mesh_entity_id );
            return mesh_entity.vertex_attribute_manager();
        }

    private:
        GeoModelMeshVerticesBase< DIMENSION >& geomodel_vertices_;
        const GeoModel< DIMENSION >& geomodel_;

        /// Vertex maps
        std::vector< std::vector< index_t > > corner_vertex_maps_;
        std::vector< std::vector< index_t > > line_vertex_maps_;
        std::vector< std::vector< index_t > > surface_vertex_maps_;
        std::vector< std::vector< index_t > > region_vertex_maps_;
        mutable std::map< MeshEntityType,
            std::vector< std::vector< index_t > >* >
            vertex_maps_;

        /// GeoModelEntity Vertices for each geomodel vertex
        mutable std::vector< std::vector< GMEVertex > > gme_vertices_;
    };

    template < index_t DIMENSION >
    GeoModelMeshVerticesBase< DIMENSION >::~GeoModelMeshVerticesBase()
    {
    }

    template < index_t DIMENSION >
    GeoModelMeshVerticesBase< DIMENSION >::GeoModelMeshVerticesBase(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< PointSetMesh< DIMENSION > >& mesh )
        : GeoModelMeshCommon< DIMENSION >( gmm, gm ),
          mesh_( mesh ),
          impl_( *this, gmm.geomodel() )
    {
        this->set_mesh( mesh_.get() );
    }

    template < index_t DIMENSION >
    AttributesManager&
        GeoModelMeshVerticesBase< DIMENSION >::attribute_manager() const
    {
        return mesh_->vertex_attribute_manager();
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::test_and_initialize() const
    {
        if( !this->is_initialized() )
        {
            initialize();
        }
    }

    template < index_t DIMENSION >
    index_t nb_entity_vertices( const GeoModel< DIMENSION >& geomodel,
        const MeshEntityType& entity_type )
    {
        index_t count{ 0 };
        for( auto i : range( geomodel.nb_mesh_entities( entity_type ) ) )
        {
            count += geomodel.mesh_entity( entity_type, i ).nb_vertices();
        }
        return count;
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >& GeoModelMeshVerticesBase< DIMENSION >::nn_search() const
    {
        test_and_initialize();
        return mesh_->vertex_nn_search();
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::fill_vertices_for_entity_type(
        const GeoModel< DIMENSION >& geomodel,
        const MeshEntityType& entity_type,
        index_t& count ) const
    {
        auto mesh_builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        for( auto i : range( geomodel.nb_mesh_entities( entity_type ) ) )
        {
            const auto& E = geomodel.mesh_entity( entity_type, i );
            if( E.nb_vertices() == 0 )
            {
                continue;
            }

            // Map and vertex
            auto id = E.gmme();
            for( auto v : range( E.nb_vertices() ) )
            {
                auto local_count = count + v;
                mesh_builder->set_vertex( local_count, E.vertex( v ) );
                // Map from vertices of MeshEntities to GeoModelMeshVerticesBase
                impl_->set_vertex_map_value( id, v, local_count );
                impl_->add_to_gme_vertices( GMEVertex( id, v ), local_count );
            }
            // Global vertex index increment
            count += E.nb_vertices();
        }
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::nb_total_vertices() const
    {
        index_t nb{ 0 };
        nb += nb_entity_vertices(
            this->geomodel_, Corner< DIMENSION >::type_name_static() );
        nb += nb_entity_vertices(
            this->geomodel_, Line< DIMENSION >::type_name_static() );
        nb += nb_entity_vertices(
            this->geomodel_, Surface< DIMENSION >::type_name_static() );
        return nb;
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::initialize() const
    {
        this->set_is_initialized( true );

        // Total number of vertices in the
        // Corners, Lines, Surfaces and Regions of the GeoModel
        auto nb = nb_total_vertices();

        // Get out if no vertices
        if( nb == 0 )
        {
            return;
        }

        // Fill the vertices
        auto builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        builder->create_vertices( nb );
        impl_->clear_and_resize_geomodel_vertex_gmes( nb );
        impl_->bind_all_mesh_entity_vertex_maps();

        fill_vertices();

        // Remove colocated vertices
        remove_colocated();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::fill_vertices() const
    {
        index_t count{ 0 };
        fill_vertices_for_entity_type(
            this->geomodel_, Corner< DIMENSION >::type_name_static(), count );
        fill_vertices_for_entity_type(
            this->geomodel_, Line< DIMENSION >::type_name_static(), count );
        fill_vertices_for_entity_type(
            this->geomodel_, Surface< DIMENSION >::type_name_static(), count );
        return count;
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::clear() const
    {
        this->set_is_initialized( false );

        this->gmm_.polygons.clear();
        this->gmm_.wells.clear();
        impl_->clear();

        auto builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        builder->clear( true, false );
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::unbind_geomodel_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        impl_->clear_vertex_map( mesh_entity_id );
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::bind_geomodel_vertex_map(
        const gmme_id& mesh_entity_id )
    {
        impl_->resize_vertex_map( mesh_entity_id );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::nb() const
    {
        test_and_initialize();
        return mesh_->nb_vertices();
    }

    template < index_t DIMENSION >
    const vecn< DIMENSION >& GeoModelMeshVerticesBase< DIMENSION >::vertex(
        index_t v ) const
    {
        test_and_initialize();
        ringmesh_assert( v < nb() );
        return mesh_->vertex( v );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::index(
        const vecn< DIMENSION >& p ) const
    {
        test_and_initialize();
        const auto& colocator = mesh_->vertex_nn_search();
        auto vertices = colocator.get_neighbors( p, this->geomodel_.epsilon() );
        if( vertices.empty() )
        {
            return NO_ID;
        }
        ringmesh_assert( vertices.size() == 1 );
        return vertices[0];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::geomodel_vertex_id(
        const gmme_id& mesh_entity, index_t entity_vertex_index ) const
    {
        test_and_initialize();
        return impl_->geomodel_vertex_index( mesh_entity, entity_vertex_index );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshVerticesBase< DIMENSION >::geomodel_vertex_id(
        const gmme_id& mesh_entity,
        const ElementLocalVertex& element_local_vertex ) const
    {
        auto entity_vertex_index =
            this->geomodel_.mesh_entity( mesh_entity )
                .mesh_element_vertex_index( element_local_vertex );
        return geomodel_vertex_id( mesh_entity, entity_vertex_index );
    }

    template < index_t DIMENSION >
    std::vector< index_t >
        GeoModelMeshVerticesBase< DIMENSION >::mesh_entity_vertex_id(
            const gmme_id& mesh_entity, index_t geomodel_vertex_id ) const
    {
        test_and_initialize();
        ringmesh_assert( geomodel_vertex_id < nb() );
        return impl_->mesh_entity_vertex_indices(
            geomodel_vertex_id, mesh_entity );
    }

    template < index_t DIMENSION >
    const std::vector< GMEVertex >&
        GeoModelMeshVerticesBase< DIMENSION >::gme_vertices( index_t v ) const
    {
        test_and_initialize();
        return impl_->mesh_entity_vertex_indices( v );
    }

    template < index_t DIMENSION >
    std::vector< GMEVertex >
        GeoModelMeshVerticesBase< DIMENSION >::gme_type_vertices(
            const MeshEntityType& entity_type, index_t vertex ) const
    {
        test_and_initialize();
        return impl_->mesh_entity_vertex_indices( vertex, entity_type );
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::set_point(
        index_t v, const vecn< DIMENSION >& point )
    {
        test_and_initialize();
        ringmesh_assert( v < nb() );
        // Change the position of the unique_vertex
        auto mesh_builder =
            PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        mesh_builder->set_vertex( v, point );
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::update_vertex_mapping(
        const gmme_id& entity_id,
        index_t entity_vertex_index,
        index_t geomodel_vertex_index )
    {
        impl_->set_vertex_map_value(
            entity_id, entity_vertex_index, geomodel_vertex_index );
        impl_->add_to_gme_vertices( GMEVertex( entity_id, entity_vertex_index ),
            geomodel_vertex_index );
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::remove_colocated() const
    {
        // Get out if nothing to do
        // and compute the points if they are not initialized yet
        if( nb() == 0 )
        {
            return;
        }
        // Identify and invalidate colocated vertices
        index_t nb_colocalised_vertices{ NO_ID };
        std::vector< index_t > old2new;
        std::tie( nb_colocalised_vertices, old2new ) =
            mesh_->vertex_nn_search().get_colocated_index_mapping(
                this->geomodel_.epsilon() );
        if( nb_colocalised_vertices > 0 )
        {
            erase_vertices( old2new );
        }
    }

    template < index_t DIMENSION >
    void GeoModelMeshVerticesBase< DIMENSION >::erase_vertices(
        std::vector< index_t >& to_delete ) const
    {
        ringmesh_assert( to_delete.size() == nb() );

        // For mesh vertices deletion
        std::vector< bool > to_delete_bool( nb(), false );

        // Fill the delete information for geogram
        // Recycle the to_delete vertex to get the mapping between
        // new and old points. This is implemented to be the same
        // as what is done in the delete_elements function in geogram
        index_t nb_todelete{ 0 };
        index_t cur{ 0 };
        for( auto v : range( nb() ) )
        {
            if( to_delete[v] != v )
            {
                to_delete_bool[v] = true;
                nb_todelete++;
                if( to_delete[v] != NO_ID )
                {
                    ringmesh_assert( to_delete[v] < v );
                    to_delete[v] = to_delete[to_delete[v]];
                }
            }
            else
            {
                to_delete[v] = cur;
                ++cur;
            }
        }
        if( nb_todelete == 0 )
        {
            return;
        }
        if( nb_todelete == nb() )
        {
            // Clear everything
            clear();
            return;
        }

        // Empty the gme_vertices_ of the deleted vertices and erase them
        for( auto v : range( nb() ) )
        {
            impl_->clear_geomodel_vertex_gmes( v );
        }

        // Delete the vertices - false is to not remove
        // isolated vertices (here all the vertices)
        PointSetMeshBuilder< DIMENSION >::create_builder( *mesh_ )
            ->delete_vertices( to_delete_bool );

        impl_->update_mesh_entity_maps_and_gmes( to_delete );
    }

    template < index_t DIMENSION >
    GeoModelMeshVertices< DIMENSION >::GeoModelMeshVertices(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< PointSetMesh< DIMENSION > >& mesh )
        : GeoModelMeshVerticesBase< DIMENSION >( gmm, gm, mesh )
    {
    }

    GeoModelMeshVertices< 3 >::GeoModelMeshVertices( GeoModelMesh3D& gmm,
        GeoModel3D& gm,
        std::unique_ptr< PointSetMesh3D >& mesh )
        : GeoModelMeshVerticesBase< 3 >( gmm, gm, mesh )
    {
    }

    void GeoModelMeshVertices< 3 >::clear() const
    {
        this->set_is_initialized( false );
        this->gmm_.cells.clear();
        GeoModelMeshVerticesBase3D::clear();
    }

    index_t GeoModelMeshVertices< 3 >::nb_total_vertices() const
    {
        auto nb = GeoModelMeshVerticesBase3D::nb_total_vertices();
        nb +=
            nb_entity_vertices( this->geomodel_, Region3D::type_name_static() );
        return nb;
    }

    index_t GeoModelMeshVertices< 3 >::fill_vertices() const
    {
        auto count = GeoModelMeshVerticesBase3D::fill_vertices();
        fill_vertices_for_entity_type(
            this->geomodel_, Region3D::type_name_static(), count );
        return count;
    }

    template <>
    GeoModelMeshVerticesBase< 3 >::Impl::Impl(
        GeoModelMeshVerticesBase& geomodel_vertices,
        const GeoModel3D& geomodel )
        : geomodel_vertices_( geomodel_vertices ), geomodel_( geomodel )
    {
        vertex_maps_[Corner3D::type_name_static()] = &corner_vertex_maps_;
        vertex_maps_[Line3D::type_name_static()] = &line_vertex_maps_;
        vertex_maps_[Surface3D::type_name_static()] = &surface_vertex_maps_;
        vertex_maps_[Region3D::type_name_static()] = &region_vertex_maps_;
    }

    /*******************************************************************************/

    template < index_t DIMENSION >
    GeoModelMeshCells< DIMENSION >::GeoModelMeshCells(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< VolumeMesh< DIMENSION > >& mesh )
        : GeoModelMeshCommon< DIMENSION >( gmm, gm ), mesh_( mesh )
    {
        this->set_mesh( mesh_.get() );
    }

    template < index_t DIMENSION >
    AttributesManager& GeoModelMeshCells< DIMENSION >::attribute_manager() const
    {
        return mesh_->cell_attribute_manager();
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::test_and_initialize() const
    {
        if( !this->is_initialized() )
        {
            const_cast< GeoModelMeshCells* >( this )->initialize();
        }
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::initialize()
    {
        this->set_is_initialized( true );

        this->gmm_.vertices.test_and_initialize();
        auto mesh_builder =
            VolumeMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != this->gmm_.vertices.nb() )
        {
            copy_vertices( mesh_builder.get(), *this->gmm_.vertices.mesh_ );
        }

        region_cell_ptr_.resize(
            this->geomodel_.nb_regions()
                    * to_underlying_type( CellType::UNDEFINED )
                + 1,
            0 );

        // Total number of  cells
        std::vector< index_t > nb_cells_per_type(
            to_underlying_type( CellType::UNDEFINED ), 0 );
        index_t nb{ 0 };

        for( const auto& region : this->geomodel_.regions() )
        {
            nb += region.nb_mesh_elements();
        }

        // Get out if no cells
        if( nb == 0 )
        {
            return;
        }

        // Compute the number of cell per type and per region
        for( const auto& region : this->geomodel_.regions() )
        {
            auto r = region.index();
            for( auto c : range( region.nb_mesh_elements() ) )
            {
                CellType cur_cell_type = region.cell_type( c );
                switch( cur_cell_type )
                {
                case CellType::TETRAHEDRON:
                    nb_cells_per_type[to_underlying_type(
                        CellType::TETRAHEDRON )]++;
                    region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                         * r
                                     + to_underlying_type(
                                           CellType::TETRAHEDRON )
                                     + 1]++;
                    break;
                case CellType::HEXAHEDRON:
                    nb_cells_per_type[to_underlying_type(
                        CellType::HEXAHEDRON )]++;
                    region_cell_ptr_
                        [to_underlying_type( CellType::UNDEFINED ) * r
                            + to_underlying_type( CellType::HEXAHEDRON ) + 1]++;
                    break;
                case CellType::PRISM:
                    nb_cells_per_type[to_underlying_type( CellType::PRISM )]++;
                    region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                         * r
                                     + to_underlying_type( CellType::PRISM )
                                     + 1]++;
                    break;
                case CellType::PYRAMID:
                    nb_cells_per_type[to_underlying_type(
                        CellType::PYRAMID )]++;
                    region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                         * r
                                     + to_underlying_type( CellType::PYRAMID )
                                     + 1]++;
                    break;
                case CellType::UNCLASSIFIED:
                    nb_cells_per_type[to_underlying_type(
                        CellType::UNCLASSIFIED )]++;
                    region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                         * r
                                     + to_underlying_type(
                                           CellType::UNCLASSIFIED )
                                     + 1]++;
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
        for( auto t : range( to_underlying_type( CellType::TETRAHEDRON ) + 1,
                 to_underlying_type( CellType::UNDEFINED ) ) )
        {
            cells_offset_per_type[t] += cells_offset_per_type[t - 1];
            cells_offset_per_type[t] += nb_cells_per_type[t - 1];
        }

        for( auto i : range( 1, region_cell_ptr_.size() - 1 ) )
        {
            region_cell_ptr_[i + 1] += region_cell_ptr_[i];
        }

        // Create "empty" tet, hex, pyr and prism
        for( auto i : range( to_underlying_type( CellType::UNDEFINED ) ) )
        {
            mesh_builder->create_cells( nb_cells_per_type[i], CellType( i ) );
        }

        // Fill the cells with vertices
        resize_cell_data();
        std::vector< index_t > cur_cell_per_type(
            to_underlying_type( CellType::UNDEFINED ), 0 );
        const auto& geomodel_vertices = this->gmm_.vertices;
        for( const auto& region : this->geomodel_.regions() )
        {
            for( auto c : range( region.nb_mesh_elements() ) )
            {
                CellType cur_cell_type = region.cell_type( c );
                index_t cur_cell =
                    cells_offset_per_type[to_underlying_type( cur_cell_type )]
                    + cur_cell_per_type[to_underlying_type( cur_cell_type )]++;
                for( auto v : range( mesh_->nb_cell_vertices( cur_cell ) ) )
                {
                    auto region_vertex_index =
                        region.mesh_element_vertex_index( { c, v } );
                    auto global_vertex_id =
                        geomodel_vertices.geomodel_vertex_id(
                            region.gmme(), region_vertex_index );
                    mesh_builder->set_cell_vertex(
                        ElementLocalVertex( cur_cell, v ), global_vertex_id );
                }
                region_id_[cur_cell] = region.index();
                cell_id_[cur_cell] = c;
            }
        }

        sort_cells();

        // Retrieve the adjacencies
        mesh_builder->connect_cells();

        // Cache some values
        nb_tets_ =
            nb_cells_per_type[to_underlying_type( CellType::TETRAHEDRON )];
        nb_hexs_ =
            nb_cells_per_type[to_underlying_type( CellType::HEXAHEDRON )];
        nb_prisms_ = nb_cells_per_type[to_underlying_type( CellType::PRISM )];
        nb_pyramids_ =
            nb_cells_per_type[to_underlying_type( CellType::PYRAMID )];
        nb_connectors_ =
            nb_cells_per_type[to_underlying_type( CellType::UNCLASSIFIED )];
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::sort_cells()
    {
        auto mesh_builder =
            VolumeMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        std::vector< index_t > sorted_indices( mesh_->nb_cells() );
        std::iota( sorted_indices.begin(), sorted_indices.end(), 0 );
        GeoModelMeshCellsSort< DIMENSION > action( *mesh_, region_id_ );
        std::sort( sorted_indices.begin(), sorted_indices.end(), action );
        mesh_builder->permute_cells( sorted_indices );

        auto sorted_indices_geo =
            copy_std_vector_to_geo_vector( sorted_indices );
        GEO::Permutation::apply(
            region_id_.data(), sorted_indices_geo, sizeof( index_t ) );
        GEO::Permutation::apply(
            cell_id_.data(), sorted_indices_geo, sizeof( index_t ) );
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::resize_cell_data()
    {
        region_id_.resize( mesh_->nb_cells(), NO_ID );
        cell_id_.resize( mesh_->nb_cells(), NO_ID );
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::clear_cell_data()
    {
        region_id_.clear();
        cell_id_.clear();
        polygon_id_.clear();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb() const
    {
        test_and_initialize();
        return mesh_->nb_cells();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_vertices( index_t cell ) const
    {
        test_and_initialize();
        ringmesh_assert( cell < mesh_->nb_cells() );
        return mesh_->nb_cell_vertices( cell );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::vertex(
        const ElementLocalVertex& cell_local_vertex ) const
    {
        test_and_initialize();
        ringmesh_assert( cell_local_vertex.element_id < mesh_->nb_cells() );
        ringmesh_assert(
            cell_local_vertex.local_vertex_id
            < mesh_->nb_cell_vertices( cell_local_vertex.element_id ) );
        return mesh_->cell_vertex( cell_local_vertex );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_edges( index_t cell ) const
    {
        test_and_initialize();
        ringmesh_assert( cell < mesh_->nb_cells() );
        return mesh_->nb_cell_edges( cell );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_facets( index_t cell ) const
    {
        test_and_initialize();
        ringmesh_assert( cell < mesh_->nb_cells() );
        return mesh_->nb_cell_facets( cell );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_facet_vertices(
        const CellLocalFacet& cell_local_facet ) const
    {
        test_and_initialize();
        ringmesh_assert( cell_local_facet.cell_id < mesh_->nb_cells() );
        ringmesh_assert( cell_local_facet.local_facet_id
                         < mesh_->nb_cell_facets( cell_local_facet.cell_id ) );
        return mesh_->nb_cell_facet_vertices( cell_local_facet );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::facet_vertex(
        const CellLocalFacet& cell_local_facet, index_t local_vertex ) const
    {
        test_and_initialize();
        ringmesh_assert( cell_local_facet.cell_id < mesh_->nb_cells() );
        ringmesh_assert( cell_local_facet.local_facet_id
                         < mesh_->nb_cell_facets( cell_local_facet.cell_id ) );
        return mesh_->cell_facet_vertex( cell_local_facet, local_vertex );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::edge_vertex(
        index_t cell, index_t local_edge, index_t local_vertex ) const
    {
        geo_debug_assert( local_edge < nb_edges( cell ) );
        geo_debug_assert( local_vertex < 2 );
        return mesh_->cell_edge_vertex( cell, local_edge, local_vertex );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::adjacent(
        index_t cell, index_t facet ) const
    {
        test_and_initialize();
        ringmesh_assert( cell < mesh_->nb_cells() );
        ringmesh_assert( facet < mesh_->nb_cell_facets( cell ) );
        return mesh_->cell_adjacent( CellLocalFacet( cell, facet ) );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::region( index_t cell ) const
    {
        test_and_initialize();
        ringmesh_assert( cell < mesh_->nb_cells() );
        return region_id_[cell];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::index_in_region(
        index_t cell ) const
    {
        test_and_initialize();
        ringmesh_assert( cell < mesh_->nb_cells() );
        return cell_id_[cell];
    }

    template < index_t DIMENSION >
    CellType GeoModelMeshCells< DIMENSION >::type( index_t cell ) const
    {
        test_and_initialize();
        ringmesh_assert( cell < mesh_->nb_cells() );
        return mesh_->cell_type( cell );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_cells( CellType type ) const
    {
        test_and_initialize();
        switch( type )
        {
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

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_cells(
        index_t region, CellType type ) const
    {
        test_and_initialize();
        switch( type )
        {
        case CellType::TETRAHEDRON:
            return nb_tet( region );
        case CellType::HEXAHEDRON:
            return nb_hex( region );
        case CellType::PRISM:
            return nb_prism( region );
        case CellType::PYRAMID:
            return nb_pyramid( region );
        case CellType::UNCLASSIFIED:
            return nb_connector( region );
        case CellType::UNDEFINED:
            ringmesh_assert(
                region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                 * ( region + 1 )]
                    - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                       * region]
                == this->geomodel_.region( region ).nb_mesh_elements() );
            return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * ( region + 1 )]
                   - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                      * region];
        default:
            ringmesh_assert_not_reached;
            return 0;
        }
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::cell(
        index_t region, index_t cell, CellType type ) const
    {
        test_and_initialize();
        switch( type )
        {
        case CellType::TETRAHEDRON:
            return tet( region, cell );
        case CellType::HEXAHEDRON:
            return hex( region, cell );
        case CellType::PRISM:
            return prism( region, cell );
        case CellType::PYRAMID:
            return pyramid( region, cell );
        case CellType::UNCLASSIFIED:
            return connector( region, cell );
        case CellType::UNDEFINED:
            return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region]
                   + cell;
        default:
            ringmesh_assert_not_reached;
            return 0;
        }
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_tet() const
    {
        test_and_initialize();
        return nb_tets_;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_tet( index_t region ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + ( to_underlying_type( CellType::TETRAHEDRON )
                                      + 1 )]
               - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                      * region
                                  + to_underlying_type(
                                        CellType::TETRAHEDRON )];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::tet(
        index_t region, index_t tet ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + to_underlying_type( CellType::TETRAHEDRON )]
               + tet;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_hex() const
    {
        test_and_initialize();
        return nb_hexs_;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_hex( index_t region ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + ( to_underlying_type( CellType::HEXAHEDRON )
                                      + 1 )]
               - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                      * region
                                  + to_underlying_type( CellType::HEXAHEDRON )];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::hex(
        index_t region, index_t hex ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + to_underlying_type( CellType::HEXAHEDRON )]
               + hex;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_prism() const
    {
        test_and_initialize();
        return nb_prisms_;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_prism( index_t region ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + ( to_underlying_type( CellType::PRISM ) + 1 )]
               - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                      * region
                                  + to_underlying_type( CellType::PRISM )];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::prism(
        index_t region, index_t prism ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + to_underlying_type( CellType::PRISM )]
               + prism;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_pyramid() const
    {
        test_and_initialize();
        return nb_pyramids_;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_pyramid( index_t region ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + ( to_underlying_type( CellType::PYRAMID )
                                      + 1 )]
               - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                      * region
                                  + to_underlying_type( CellType::PYRAMID )];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::pyramid(
        index_t region, index_t pyramid ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + to_underlying_type( CellType::PYRAMID )]
               + pyramid;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_connector() const
    {
        test_and_initialize();
        return nb_connectors_;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_connector( index_t region ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + ( to_underlying_type( CellType::UNCLASSIFIED )
                                      + 1 )]
               - region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                      * region
                                  + to_underlying_type(
                                        CellType::UNCLASSIFIED )];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::connector(
        index_t region, index_t connector ) const
    {
        test_and_initialize();
        ringmesh_assert( region < this->geomodel_.nb_regions() );
        return region_cell_ptr_[to_underlying_type( CellType::UNDEFINED )
                                    * region
                                + to_underlying_type( CellType::UNCLASSIFIED )]
               + connector;
    }

    template < index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::is_duplication_initialized() const
    {
        return mode_ == this->gmm_.duplicate_mode();
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::test_and_initialize_duplication() const
    {
        if( !is_duplication_initialized() )
        {
            const_cast< GeoModelMeshCells< DIMENSION >* >( this )
                ->initialize_duplication();
        }
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::initialize_duplication()
    {
        test_and_initialize();

        /// 1. Get all the corner vertices (a lot of duplicated vertices)
        std::vector< vec3 > corner_vertices(
            mesh_->cell_end( mesh_->nb_cells() - 1 ) );
        for( auto c : range( mesh_->nb_cells() ) )
        {
            index_t begin = mesh_->cell_begin( c );
            for( auto v : range( mesh_->nb_cell_vertices( c ) ) )
            {
                corner_vertices[begin + v] = mesh_->vertex(
                    mesh_->cell_vertex( ElementLocalVertex( c, v ) ) );
            }
        }

        /// 2. Tag all corners to duplicate (vertices on a surface to duplicate)
        std::vector< ActionOnSurface > actions_on_surfaces(
            this->geomodel_.nb_surfaces(), SKIP );
        std::vector< bool > is_vertex_to_duplicate(
            corner_vertices.size(), false );
        {
            NNSearch< DIMENSION > nn_search( corner_vertices, false );
            for( const auto& surface : this->geomodel_.surfaces() )
            {
                if( !is_surface_to_duplicate( surface.index() ) )
                {
                    continue;
                }
                actions_on_surfaces[surface.index()] = TO_PROCESS;
                for( auto v : range( surface.nb_vertices() ) )
                {
                    auto colocated_corners = nn_search.get_neighbors(
                        surface.vertex( v ), this->geomodel_.epsilon() );
                    for( auto co : colocated_corners )
                    {
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
        for( auto c : range( mesh_->nb_cells() ) )
        {
            for( auto v : range( mesh_->nb_cell_vertices( c ) ) )
            {
                // get the index of the corner inside cell_corners_
                index_t co = mesh_->cell_begin( c ) + v;

                if( !is_vertex_to_duplicate[co] )
                {
                    continue;
                }
                // The vertex is on a surface to duplicate

                // Propagate on the cells around the corresponding vertex.
                // The propagation process cannot cross any surface.
                auto vertex_id =
                    mesh_->cell_vertex( ElementLocalVertex( c, v ) );

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
                do
                {
                    auto cur_c = S.top();
                    S.pop();
                    // Find which corner of the current cell matches vertex_id
                    auto cur_co = mesh_->find_cell_corner( cur_c, vertex_id );
                    ringmesh_assert( cur_co != NO_ID );
                    is_vertex_to_duplicate[cur_co] = false;
                    corner_used.push_back( cur_co );

                    // Find the cell facets including the vertex
                    auto facets =
                        cell_facets_around_vertex( *mesh_, cur_c, vertex_id );
                    for( auto cur_f : facets )
                    {
                        // Find if the facet is on a surface or inside the
                        // domain
                        index_t polygon{ NO_ID };
                        bool side;
                        if( is_cell_facet_on_surface(
                                cur_c, cur_f, polygon, side ) )
                        {
                            auto surface_id =
                                this->gmm_.polygons.surface( polygon );
                            surfaces.emplace_back(
                                surface_id, ActionOnSurface( side ) );
                        }
                        else
                        {
                            // The cell facet is not on a surface.
                            // Add the adjacent cell to the stack if it exists
                            // and has not already been processed or added into
                            // the stack
                            auto cur_adj =
                                mesh_->cell_adjacent( { cur_c, cur_f } );
                            if( cur_adj != NO_ID
                                && !contains( cell_added, cur_adj ) )
                            {
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
                if( are_corners_to_duplicate( surfaces, actions_on_surfaces ) )
                {
                    // Add a new duplicated vertex and its associated vertex

                    /* @todo Review : Use the total_nb_vertices function [JP]
                     * why mm_.vertices.nb_vertices() and not nb_vertices() ?
                     * Please help the reader !! same thing 2 lines below [JP]
                     */
                    auto duplicated_vertex_id =
                        this->gmm_.vertices.nb()
                        + static_cast< index_t >(
                              duplicated_vertex_indices_.size() );
                    duplicated_vertex_indices_.push_back( vertex_id );

                    // Update all the cell corners on this side of the surface
                    // to the new duplicated vertex index
                    auto mesh_builder =
                        VolumeMeshBuilder< DIMENSION >::create_builder(
                            *mesh_ );
                    for( auto cur_co : corner_used )
                    {
                        mesh_builder->set_cell_corner_vertex_index(
                            cur_co, duplicated_vertex_id );
                    }
                }
            }
        }
        mode_ = this->gmm_.duplicate_mode();
    }

    template < index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::is_cell_facet_on_surface( index_t cell,
        index_t facet_index,
        index_t& colocated_facet_index,
        bool& side ) const
    {
        test_and_initialize_cell_facet();
        colocated_facet_index =
            polygon_id_[mesh_->cell_facet( { cell, facet_index } )];
        if( colocated_facet_index != NO_ID )
        {
            auto facet_normal =
                this->gmm_.polygons.normal( colocated_facet_index );
            auto cell_facet_normal =
                mesh_->cell_facet_normal( { cell, facet_index } );
            side = dot( facet_normal, cell_facet_normal ) > 0;
        }
        return colocated_facet_index != NO_ID;
    }

    template < index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::are_corners_to_duplicate(
        const std::vector< action_on_surface >& surfaces,
        std::vector< ActionOnSurface >& info )
    {
        // Temporary vector, it is equal to surfaces
        std::vector< action_on_surface > temp_surfaces;
        temp_surfaces.reserve( temp_surfaces.size() );

        // Find if a free border was found, if so we encountered the
        // two sides of the surface during the propagation around a vertex
        for( auto i : range( 1, surfaces.size() ) )
        {
            if( surfaces[i - 1].first == surfaces[i].first )
            {
                // Found free border -> skip
                ringmesh_assert( surfaces[i - 1].second != surfaces[i].second );
                i++; // skip the action_on_surface (i-1) and i
            }
            else
            {
                temp_surfaces.push_back( surfaces[i - 1] );
            }
        }
        temp_surfaces.push_back( surfaces.back() );

        for( const auto& action : temp_surfaces )
        {
            auto s = action.first;
            switch( info[s] )
            {
            case SKIP:
                break;
            case TO_PROCESS:
                // First time we encounter this surface, do not duplicate
                // this side but wait to see if we encounter the other.
                // In the case of surfaces in the VOI, it is encountered only
                // once
                ringmesh_assert( action.second > TO_PROCESS );
                info[s] = ActionOnSurface( !action.second );
                break;
            default:
                // If the side matches -> duplicate
                if( info[s] == action.second )
                {
                    return true;
                }
            }
        }

        return false;
    }

    template < index_t DIMENSION >
    bool GeoModelMeshCells< DIMENSION >::is_surface_to_duplicate(
        index_t surface_id ) const
    {
        const auto& cur_surface = this->geomodel_.surface( surface_id );
        if( cur_surface.is_on_voi() )
        {
            return false;
        }
        switch( this->gmm_.duplicate_mode() )
        {
        case ALL:
            return true;
        case FAULT:
        {
            auto parent_interface = cur_surface.parent_gmge(
                Interface< DIMENSION >::type_name_static() );
            if( parent_interface.is_defined() )
            {
                auto feature =
                    this->geomodel_.geological_entity( parent_interface )
                        .geological_feature();
                return GeoModelGeologicalEntity< DIMENSION >::is_fault(
                    feature );
            }
            return false;
        }
        default:
            return false;
        }
        return false;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_duplicated_vertices() const
    {
        test_and_initialize_duplication();
        return static_cast< index_t >( duplicated_vertex_indices_.size() );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::nb_total_vertices() const
    {
        return nb_duplicated_vertices() + mesh_->nb_vertices();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::duplicated_corner_index(
        const ElementLocalVertex& cell_local_vertex ) const
    {
        test_and_initialize_duplication();
        ringmesh_assert( cell_local_vertex.element_id < mesh_->nb_cells() );
        ringmesh_assert(
            cell_local_vertex.local_vertex_id
            < mesh_->nb_cell_vertices( cell_local_vertex.element_id ) );
        auto corner_value = mesh_->cell_vertex( cell_local_vertex );
        if( corner_value < mesh_->nb_vertices() )
        {
            return NO_ID;
        }
        return corner_value - mesh_->nb_vertices();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshCells< DIMENSION >::duplicated_vertex(
        index_t duplicate_vertex_index ) const
    {
        test_and_initialize_duplication();
        ringmesh_assert(
            duplicate_vertex_index < duplicated_vertex_indices_.size() );
        return duplicated_vertex_indices_[duplicate_vertex_index];
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::clear()
    {
        this->set_is_initialized( false );

        auto mesh_builder =
            VolumeMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        mesh_builder->clear( true, false );
        region_cell_ptr_.clear();
        nb_tets_ = 0;
        nb_hexs_ = 0;
        nb_prisms_ = 0;
        nb_pyramids_ = 0;
        nb_connectors_ = 0;

        mode_ = NONE;
        duplicated_vertex_indices_.clear();
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::clear_duplication()
    {
        auto mesh_builder =
            VolumeMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        for( auto c : range( mesh_->nb_cells() ) )
        {
            for( auto v : range( mesh_->nb_cell_vertices( c ) ) )
            {
                auto index = duplicated_corner_index( { c, v } );
                if( index != NO_ID )
                {
                    mesh_builder->set_cell_corner_vertex_index(
                        c, duplicated_vertex( index ) );
                }
            }
        }

        mode_ = NONE;
        duplicated_vertex_indices_.clear();
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::test_and_initialize_cell_facet() const
    {
        if( polygon_id_.empty() )
        {
            const_cast< GeoModelMeshCells* >( this )->initialize_cell_facet();
        }
    }

    template < index_t DIMENSION >
    void GeoModelMeshCells< DIMENSION >::initialize_cell_facet()
    {
        this->gmm_.polygons.test_and_initialize();

        polygon_id_.resize( mesh_->nb_cell_facets(), NO_ID );
        const auto& nn_search = this->gmm_.polygons.nn_search();
        for( auto c : range( mesh_->nb_cells() ) )
        {
            for( auto f : range( mesh_->nb_cell_facets( c ) ) )
            {
                auto result = nn_search.get_neighbors(
                    mesh_->cell_facet_barycenter( { c, f } ),
                    this->geomodel_.epsilon() );
                if( !result.empty() )
                {
                    polygon_id_[mesh_->cell_facet( { c, f } )] = result[0];
                    // If there are more than 1 matching facet, this is WRONG
                    // and the vertex indices should be checked too [Jeanne]
                    ringmesh_assert( result.size() == 1 );
                }
            }
        }
    }

    template < index_t DIMENSION >
    vec3 GeoModelMeshCells< DIMENSION >::barycenter( index_t cell ) const
    {
        test_and_initialize();
        return mesh_->cell_barycenter( cell );
    }

    template < index_t DIMENSION >
    double GeoModelMeshCells< DIMENSION >::volume( index_t cell ) const
    {
        test_and_initialize();
        return mesh_->cell_volume( cell );
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >& GeoModelMeshCells< DIMENSION >::cell_nn_search() const
    {
        test_and_initialize();
        return mesh_->cell_nn_search();
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >& GeoModelMeshCells< DIMENSION >::cell_facet_nn_search() const
    {
        test_and_initialize();
        return mesh_->cell_facet_nn_search();
    }

    template < index_t DIMENSION >
    const VolumeAABBTree< DIMENSION >&
        GeoModelMeshCells< DIMENSION >::aabb() const
    {
        test_and_initialize();
        return mesh_->cell_aabb();
    }

    /*******************************************************************************/
    template < index_t DIMENSION >
    const std::string GeoModelMeshEdges< DIMENSION >::line_att_name = "line";
    template < index_t DIMENSION >
    const std::string GeoModelMeshEdges< DIMENSION >::edge_line_att_name =
        "edge_line";

    template < index_t DIMENSION >
    GeoModelMeshEdges< DIMENSION >::GeoModelMeshEdges(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< LineMesh< DIMENSION > >& mesh )
        : GeoModelMeshCommon< DIMENSION >( gmm, gm ), mesh_( mesh )
    {
        this->set_mesh( mesh_.get() );
    }

    template < index_t DIMENSION >
    GeoModelMeshEdges< DIMENSION >::~GeoModelMeshEdges()
    {
        clear_edge_data();
    }

    template < index_t DIMENSION >
    void GeoModelMeshEdges< DIMENSION >::resize_edge_data()
    {
        line_id_.resize( mesh_->nb_edges(), NO_ID );
        edge_id_.resize( mesh_->nb_edges(), NO_ID );
    }

    template < index_t DIMENSION >
    void GeoModelMeshEdges< DIMENSION >::clear_edge_data()
    {
        line_id_.clear();
        edge_id_.clear();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::nb() const
    {
        test_and_initialize();
        return mesh_->nb_edges();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::vertex(
        const ElementLocalVertex& edge_local_vertex ) const
    {
        test_and_initialize();
        ringmesh_assert( edge_local_vertex.element_id < mesh_->nb_edges() );
        ringmesh_assert( edge_local_vertex.local_vertex_id < 2 );
        return mesh_->edge_vertex( edge_local_vertex );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::line( index_t edge ) const
    {
        test_and_initialize();
        ringmesh_assert( edge < mesh_->nb_edges() );
        return line_id_[edge];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::index_in_line( index_t edge ) const
    {
        test_and_initialize();
        ringmesh_assert( edge < mesh_->nb_edges() );
        return edge_id_[edge];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::nb_edges( index_t line ) const
    {
        test_and_initialize();
        ringmesh_assert( line < this->geomodel_.nb_lines() );
        return line_edge_ptr_[line + 1];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshEdges< DIMENSION >::edge(
        index_t line, index_t edge ) const
    {
        test_and_initialize();
        ringmesh_assert( line < this->geomodel_.nb_lines() );
        return line_edge_ptr_[line] + edge;
    }

    template < index_t DIMENSION >
    void GeoModelMeshEdges< DIMENSION >::clear()
    {
        line_edge_ptr_.clear();
        nb_edges_ = 0;
        auto mesh_builder =
            LineMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        mesh_builder->clear( true, false );
    }

    template < index_t DIMENSION >
    void GeoModelMeshEdges< DIMENSION >::test_and_initialize() const
    {
        if( !this->is_initialized() )
        {
            const_cast< GeoModelMeshEdges* >( this )->initialize();
        }
    }

    template < index_t DIMENSION >
    void GeoModelMeshEdges< DIMENSION >::initialize()
    {
        this->set_is_initialized( true );

        this->gmm_.vertices.test_and_initialize();
        line_edge_ptr_.resize( this->geomodel_.nb_lines() + 1, 0 );
        auto mesh_builder =
            LineMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != this->gmm_.vertices.nb() )
        {
            copy_vertices( mesh_builder.get(), *this->gmm_.vertices.mesh_ );
        }

        // Compute the total number of edges per line
        for( auto l : range( this->geomodel_.nb_lines() ) )
        {
            const Line< DIMENSION >& line = this->geomodel_.line( l );
            line_edge_ptr_[l + 1] = line_edge_ptr_[l] + line.nb_mesh_elements();
            nb_edges_ += line.nb_mesh_elements();
        }

        // Create  edges
        mesh_builder->create_edges( nb_edges_ );
        resize_edge_data();
        const auto& geomodel_vertices = this->gmm_.vertices;
        index_t cur_edge{ 0 };
        for( auto l : range( this->geomodel_.nb_lines() ) )
        {
            const auto& line = this->geomodel_.line( l );
            auto line_id = line.gmme();
            for( auto e : range( line.nb_mesh_elements() ) )
            {
                for( auto v : range( 2 ) )
                {
                    auto v_id = geomodel_vertices.geomodel_vertex_id(
                        line_id, ElementLocalVertex( e, v ) );
                    ringmesh_assert( v_id != NO_ID );
                    mesh_builder->set_edge_vertex(
                        ElementLocalVertex( cur_edge, v ), v_id );
                }
                line_id_[cur_edge] = l;
                edge_id_[cur_edge] = e;
                cur_edge++;
            }
        }
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > GeoModelMeshEdges< DIMENSION >::center(
        index_t edge ) const
    {
        test_and_initialize();
        return mesh_->edge_barycenter( edge );
    }

    template < index_t DIMENSION >
    double GeoModelMeshEdges< DIMENSION >::length( index_t edge ) const
    {
        test_and_initialize();
        return mesh_->edge_length( edge );
    }

    template < index_t DIMENSION >
    AttributesManager& GeoModelMeshEdges< DIMENSION >::attribute_manager() const
    {
        return mesh_->edge_attribute_manager();
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >& GeoModelMeshEdges< DIMENSION >::nn_search() const
    {
        test_and_initialize();
        return mesh_->edge_nn_search();
    }

    template < index_t DIMENSION >
    const LineAABBTree< DIMENSION >&
        GeoModelMeshEdges< DIMENSION >::aabb() const
    {
        test_and_initialize();
        return mesh_->edge_aabb();
    }

    /*******************************************************************************/
    template < index_t DIMENSION >
    const std::string GeoModelMeshPolygonsBase< DIMENSION >::surface_att_name =
        "surface";
    template < index_t DIMENSION >
    const std::string
        GeoModelMeshPolygonsBase< DIMENSION >::polygon_surface_att_name =
            "polygon_surface";

    template < index_t DIMENSION >
    GeoModelMeshPolygonsBase< DIMENSION >::GeoModelMeshPolygonsBase(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< SurfaceMesh< DIMENSION > >& mesh )
        : GeoModelMeshCommon< DIMENSION >( gmm, gm ), mesh_( mesh )
    {
        this->set_mesh( mesh_.get() );
    }

    template < index_t DIMENSION >
    GeoModelMeshPolygonsBase< DIMENSION >::~GeoModelMeshPolygonsBase()
    {
        clear_polygon_data();
    }

    template < index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::resize_polygon_data(
        index_t nb_polygons )
    {
        surface_id_.resize( nb_polygons, NO_ID );
        polygon_id_.resize( nb_polygons, NO_ID );
    }

    template < index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::clear_polygon_data()
    {
        surface_id_.clear();
        polygon_id_.clear();
    }

    template < index_t DIMENSION >
    AttributesManager& GeoModelMeshPolygonsBase< DIMENSION >::attribute_manager() const
    {
        return mesh_->polygon_attribute_manager();
    }

    template < index_t DIMENSION >
    const NNSearch< DIMENSION >& GeoModelMeshPolygonsBase< DIMENSION >::nn_search() const
    {
        test_and_initialize();
        return mesh_->polygon_nn_search();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb() const
    {
        test_and_initialize();
        return mesh_->nb_polygons();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_vertices(
        index_t polygon ) const
    {
        test_and_initialize();
        ringmesh_assert( polygon < mesh_->nb_polygons() );
        return mesh_->nb_polygon_vertices( polygon );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::vertex(
        const ElementLocalVertex& polygon_local_vertex ) const
    {
        test_and_initialize();
        ringmesh_assert(
            polygon_local_vertex.element_id < mesh_->nb_polygons() );
        ringmesh_assert(
            polygon_local_vertex.local_vertex_id
            < mesh_->nb_polygon_vertices( polygon_local_vertex.element_id ) );
        return mesh_->polygon_vertex( polygon_local_vertex );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::adjacent(
        const PolygonLocalEdge& polygon_local_edge ) const
    {
        test_and_initialize();
        ringmesh_assert( polygon_local_edge.polygon_id < mesh_->nb_polygons() );
        ringmesh_assert(
            polygon_local_edge.local_edge_id
            < mesh_->nb_polygon_vertices( polygon_local_edge.polygon_id ) );
        return mesh_->polygon_adjacent( polygon_local_edge );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::surface(
        index_t polygon ) const
    {
        test_and_initialize();
        ringmesh_assert( polygon < mesh_->nb_polygons() );
        return surface_id_[polygon];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::index_in_surface(
        index_t polygon ) const
    {
        test_and_initialize();
        ringmesh_assert( polygon < mesh_->nb_polygons() );
        return polygon_id_[polygon];
    }

    template < index_t DIMENSION >
    std::tuple< PolygonType, index_t >
        GeoModelMeshPolygonsBase< DIMENSION >::type( index_t polygon ) const
    {
        test_and_initialize();
        ringmesh_assert( polygon < mesh_->nb_polygons() );
        auto local_polygon = index_in_surface( polygon );
        auto s = surface( polygon );
        for( auto t : range( to_underlying_type( PolygonType::TRIANGLE ),
                 to_underlying_type( PolygonType::UNDEFINED ) ) )
        {
            auto T = static_cast< PolygonType >( t );
            if( local_polygon < nb_polygons( s, T ) )
            {
                return std::make_tuple( T, local_polygon );
            }
            local_polygon -= nb_polygons( s, T );
        }
        ringmesh_assert_not_reached;
        return std::make_tuple( PolygonType::UNDEFINED, NO_ID );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_polygons(
        PolygonType type ) const
    {
        test_and_initialize();
        switch( type )
        {
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

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_polygons(
        index_t surface, PolygonType type ) const
    {
        test_and_initialize();
        ringmesh_assert( surface < this->geomodel_.nb_surfaces() );
        switch( type )
        {
        case PolygonType::TRIANGLE:
            return nb_triangle( surface );
        case PolygonType::QUAD:
            return nb_quad( surface );
        case PolygonType::UNCLASSIFIED:
            return nb_unclassified_polygon( surface );
        case PolygonType::UNDEFINED:
            return surface_polygon_ptr_[to_underlying_type(
                                            PolygonType::UNDEFINED )
                                        * ( surface + 1 )]
                   - surface_polygon_ptr_[to_underlying_type(
                                              PolygonType::UNDEFINED )
                                          * surface];
        default:
            ringmesh_assert_not_reached;
            return 0;
        }
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::polygon(
        index_t surface, index_t polygon, PolygonType type ) const
    {
        test_and_initialize();
        ringmesh_assert( surface < this->geomodel_.nb_surfaces() );
        switch( type )
        {
        case PolygonType::TRIANGLE:
            return triangle( surface, polygon );
        case PolygonType::QUAD:
            return quad( surface, polygon );
        case PolygonType::UNCLASSIFIED:
            return unclassified_polygon( surface, polygon );
        case PolygonType::UNDEFINED:
            return surface_polygon_ptr_[to_underlying_type(
                                            PolygonType::UNDEFINED )
                                        * surface]
                   + polygon;
        default:
            ringmesh_assert_not_reached;
            return 0;
        }
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_triangle() const
    {
        test_and_initialize();
        return nb_triangles_;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_triangle(
        index_t surface ) const
    {
        test_and_initialize();
        ringmesh_assert( surface < this->geomodel_.nb_surfaces() );
        return surface_polygon_ptr_
                   [to_underlying_type( PolygonType::UNDEFINED ) * surface
                       + ( to_underlying_type( PolygonType::TRIANGLE ) + 1 )]
               - surface_polygon_ptr_
                     [to_underlying_type( PolygonType::UNDEFINED ) * surface
                         + to_underlying_type( PolygonType::TRIANGLE )];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::triangle(
        index_t surface, index_t triangle ) const
    {
        test_and_initialize();
        ringmesh_assert( surface < this->geomodel_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED )
                                        * surface
                                    + to_underlying_type(
                                          PolygonType::TRIANGLE )]
               + triangle;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_quad() const
    {
        test_and_initialize();
        return nb_quads_;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_quad(
        index_t surface ) const
    {
        test_and_initialize();
        ringmesh_assert( surface < this->geomodel_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED )
                                        * surface
                                    + ( to_underlying_type( PolygonType::QUAD )
                                          + 1 )]
               - surface_polygon_ptr_
                     [to_underlying_type( PolygonType::UNDEFINED ) * surface
                         + to_underlying_type( PolygonType::QUAD )];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::quad(
        index_t surface, index_t quad ) const
    {
        test_and_initialize();
        ringmesh_assert( surface < this->geomodel_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED )
                                        * surface
                                    + to_underlying_type( PolygonType::QUAD )]
               + quad;
    }

    template < index_t DIMENSION >
    index_t
        GeoModelMeshPolygonsBase< DIMENSION >::nb_unclassified_polygon() const
    {
        test_and_initialize();
        return nb_unclassified_polygons_;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::nb_unclassified_polygon(
        index_t surface ) const
    {
        test_and_initialize();
        ringmesh_assert( surface < this->geomodel_.nb_surfaces() );
        return surface_polygon_ptr_
                   [to_underlying_type( PolygonType::UNDEFINED ) * surface
                       + ( to_underlying_type( PolygonType::UNCLASSIFIED )
                             + 1 )]
               - surface_polygon_ptr_
                     [to_underlying_type( PolygonType::UNDEFINED ) * surface
                         + to_underlying_type( PolygonType::UNCLASSIFIED )];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshPolygonsBase< DIMENSION >::unclassified_polygon(
        index_t surface, index_t polygon ) const
    {
        test_and_initialize();
        ringmesh_assert( surface < this->geomodel_.nb_surfaces() );
        return surface_polygon_ptr_[to_underlying_type( PolygonType::UNDEFINED )
                                        * surface
                                    + to_underlying_type(
                                          PolygonType::UNCLASSIFIED )]
               + polygon;
    }

    template < index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::clear()
    {
        this->set_is_initialized( false );
        surface_polygon_ptr_.clear();
        nb_triangles_ = 0;
        nb_quads_ = 0;
        auto mesh_builder =
            SurfaceMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        mesh_builder->clear( true, false );
    }

    template < index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::test_and_initialize() const
    {
        if( !this->is_initialized() )
        {
            const_cast< GeoModelMeshPolygonsBase* >( this )->initialize();
        }
    }

    template < index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::initialize()
    {
        this->set_is_initialized( true );

        this->gmm_.vertices.test_and_initialize();
        surface_polygon_ptr_.resize(
            this->geomodel_.nb_surfaces()
                    * to_underlying_type( PolygonType::UNDEFINED )
                + 1,
            0 );
        auto mesh_builder =
            SurfaceMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != this->gmm_.vertices.nb() )
        {
            copy_vertices( mesh_builder.get(), *this->gmm_.vertices.mesh_ );
        }

        // Compute the total number of polygons per type and per surface
        std::map< PolygonType, index_t > nb_polygon_per_type = {
            { PolygonType::TRIANGLE, 0 }, { PolygonType::QUAD, 0 },
            { PolygonType::UNCLASSIFIED, 0 }
        };
        for( auto s : range( this->geomodel_.nb_surfaces() ) )
        {
            const auto& surface = this->geomodel_.surface( s );
            if( surface.is_simplicial() )
            {
                nb_polygon_per_type[PolygonType::TRIANGLE] +=
                    surface.nb_mesh_elements();
                surface_polygon_ptr_
                    [to_underlying_type( PolygonType::UNDEFINED ) * s
                        + to_underlying_type( PolygonType::TRIANGLE ) + 1] +=
                    surface.nb_mesh_elements();
            }
            else
            {
                for( auto p : range( surface.nb_mesh_elements() ) )
                {
                    switch( surface.nb_mesh_element_vertices( p ) )
                    {
                    case 3:
                        nb_polygon_per_type[PolygonType::TRIANGLE]++;
                        surface_polygon_ptr_
                            [to_underlying_type( PolygonType::UNDEFINED ) * s
                                + to_underlying_type( PolygonType::TRIANGLE )
                                + 1]++;
                        break;
                    case 4:
                        nb_polygon_per_type[PolygonType::QUAD]++;
                        surface_polygon_ptr_
                            [to_underlying_type( PolygonType::UNDEFINED ) * s
                                + to_underlying_type( PolygonType::QUAD )
                                + 1]++;
                        break;
                    default:
                        nb_polygon_per_type[PolygonType::UNCLASSIFIED]++;
                        surface_polygon_ptr_[to_underlying_type(
                                                 PolygonType::UNDEFINED )
                                                 * s
                                             + to_underlying_type(
                                                   PolygonType::UNCLASSIFIED )
                                             + 1]++;
                        break;
                    }
                }
            }
        }

        // Get out if no polygons
        auto nb_total_polygons =
            nb_polygon_per_type[PolygonType::TRIANGLE]

            + nb_polygon_per_type[PolygonType::QUAD]
            + nb_polygon_per_type[PolygonType::UNCLASSIFIED];

        if( nb_total_polygons == 0 )
        {
            return;
        }

        // Create triangles and quads, the polygons will be handled later
        if( nb_polygon_per_type[PolygonType::TRIANGLE] != 0 )
        {
            mesh_builder->create_triangles(
                nb_polygon_per_type[PolygonType::TRIANGLE] );
        }
        if( nb_polygon_per_type[PolygonType::QUAD] != 0 )
        {
            mesh_builder->create_quads(
                nb_polygon_per_type[PolygonType::QUAD] );
        }

        // Compute the polygon offset
        std::map< PolygonType, index_t > polygon_offset_per_type = {
            { PolygonType::TRIANGLE, 0 }, { PolygonType::QUAD, 0 },
            { PolygonType::UNCLASSIFIED, 0 }
        };
        polygon_offset_per_type[PolygonType::QUAD] =
            nb_polygon_per_type[PolygonType::TRIANGLE];
        polygon_offset_per_type[PolygonType::UNCLASSIFIED] =
            nb_polygon_per_type[PolygonType::TRIANGLE]
            + nb_polygon_per_type[PolygonType::QUAD];

        for( auto i : range( 1, surface_polygon_ptr_.size() - 1 ) )
        {
            surface_polygon_ptr_[i + 1] += surface_polygon_ptr_[i];
        }

        // Fill the triangles and quads created above
        // Create and fill polygons
        resize_polygon_data( nb_total_polygons );
        const auto& geomodel_vertices = this->gmm_.vertices;
        std::vector< index_t > cur_polygon_per_type(
            to_underlying_type( PolygonType::UNDEFINED ), 0 );
        for( auto s : range( this->geomodel_.nb_surfaces() ) )
        {
            const auto& surface = this->geomodel_.surface( s );
            auto surface_id = surface.gmme();
            for( auto p : range( surface.nb_mesh_elements() ) )
            {
                auto nb_vertices = surface.nb_mesh_element_vertices( p );
                index_t cur_polygon{ NO_ID };
                if( nb_vertices < 5 )
                {
                    auto T = static_cast< PolygonType >( nb_vertices - 3 );
                    cur_polygon =
                        polygon_offset_per_type[T]
                        + cur_polygon_per_type[to_underlying_type( T )]++;
                    for( auto v : range( nb_vertices ) )
                    {
                        auto v_id = geomodel_vertices.geomodel_vertex_id(
                            surface_id, ElementLocalVertex( p, v ) );
                        ringmesh_assert( v_id != NO_ID );
                        mesh_builder->set_polygon_vertex(
                            ElementLocalVertex( cur_polygon, v ), v_id );
                    }
                }
                else
                {
                    std::vector< index_t > vertices( nb_vertices );
                    for( auto v : range( nb_vertices ) )
                    {
                        vertices[v] = geomodel_vertices.geomodel_vertex_id(
                            surface_id, ElementLocalVertex( p, v ) );
                    }
                    cur_polygon = mesh_builder->create_polygon( vertices );
                }
                surface_id_[cur_polygon] = s;
                polygon_id_[cur_polygon] = p;
            }
        }

        sort_polygons();

        // Compute polygon adjacencies
        mesh_builder->connect_polygons();
        disconnect_along_lines();

        // Cache some values
        nb_triangles_ = nb_polygon_per_type[PolygonType::TRIANGLE];
        nb_quads_ = nb_polygon_per_type[PolygonType::QUAD];
        nb_unclassified_polygons_ =
            nb_polygon_per_type[PolygonType::UNCLASSIFIED];
    }

    template < index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::sort_polygons()
    {
        auto mesh_builder =
            SurfaceMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        std::vector< index_t > sorted_indices( mesh_->nb_polygons() );
        std::iota( sorted_indices.begin(), sorted_indices.end(), 0 );
        GeoModelMeshPolygonsBaseSort< DIMENSION > action( *mesh_, surface_id_ );
        std::sort( sorted_indices.begin(), sorted_indices.end(), action );
        mesh_builder->permute_polygons( sorted_indices );

        auto sorted_indices_geo =
            copy_std_vector_to_geo_vector( sorted_indices );
        GEO::Permutation::apply(
            surface_id_.data(), sorted_indices_geo, sizeof( index_t ) );
        GEO::Permutation::apply(
            polygon_id_.data(), sorted_indices_geo, sizeof( index_t ) );
    }

    template < index_t DIMENSION >
    void GeoModelMeshPolygonsBase< DIMENSION >::disconnect_along_lines()
    {
        auto mesh_builder =
            SurfaceMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        for( auto s : range( this->geomodel_.nb_surfaces() ) )
        {
            const auto& surface = this->geomodel_.surface( s );
            for( auto p : range( nb_polygons( s ) ) )
            {
                auto polygon_id = polygon( s, p );
                auto surface_polygon_id = index_in_surface( polygon_id );
                for( auto v : range( nb_vertices( polygon_id ) ) )
                {
                    auto adj = surface.polygon_adjacent_index(
                        PolygonLocalEdge( surface_polygon_id, v ) );
                    if( adj == NO_ID )
                    {
                        mesh_builder->set_polygon_adjacent(
                            ElementLocalVertex( polygon_id, v ), NO_ID );
                    }
                }
            }
        }
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > GeoModelMeshPolygonsBase< DIMENSION >::center(
        index_t polygon ) const
    {
        test_and_initialize();
        return mesh_->polygon_barycenter( polygon );
    }

    template < index_t DIMENSION >
    double GeoModelMeshPolygonsBase< DIMENSION >::area( index_t polygon ) const
    {
        test_and_initialize();
        return mesh_->polygon_area( polygon );
    }

    template < index_t DIMENSION >
    const SurfaceAABBTree< DIMENSION >&
        GeoModelMeshPolygonsBase< DIMENSION >::aabb() const
    {
        test_and_initialize();
        return mesh_->polygon_aabb();
    }

    template < index_t DIMENSION >
    GeoModelMeshPolygons< DIMENSION >::GeoModelMeshPolygons(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< SurfaceMesh< DIMENSION > >& mesh )
        : GeoModelMeshPolygonsBase< DIMENSION >( gmm, gm, mesh )
    {
    }

    GeoModelMeshPolygons< 3 >::GeoModelMeshPolygons( GeoModelMesh< 3 >& gmm,
        GeoModel3D& gm,
        std::unique_ptr< SurfaceMesh3D >& mesh )
        : GeoModelMeshPolygonsBase< 3 >( gmm, gm, mesh )
    {
    }

    vec3 GeoModelMeshPolygons< 3 >::normal( index_t p ) const
    {
        test_and_initialize();
        return mesh_->polygon_normal( p );
    }

    /*******************************************************************************/

    template < index_t DIMENSION >
    GeoModelMeshWells< DIMENSION >::GeoModelMeshWells(
        GeoModelMesh< DIMENSION >& gmm,
        GeoModel< DIMENSION >& gm,
        std::unique_ptr< LineMesh< DIMENSION > >& mesh )
        : GeoModelMeshCommon< DIMENSION >( gmm, gm ), mesh_( mesh )
    {
        this->set_mesh( mesh_.get() );
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshWells< DIMENSION >::nb_wells() const
    {
        test_and_initialize();
        return this->geomodel_.wells() ? this->geomodel_.wells()->nb_wells()
                                       : 0;
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshWells< DIMENSION >::nb_edges() const
    {
        test_and_initialize();
        return mesh_->nb_edges();
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshWells< DIMENSION >::nb_edges( index_t well ) const
    {
        test_and_initialize();
        return well_ptr_[well + 1] - well_ptr_[well];
    }

    template < index_t DIMENSION >
    index_t GeoModelMeshWells< DIMENSION >::vertex(
        index_t well, index_t edge, index_t vertex ) const
    {
        test_and_initialize();
        return mesh_->edge_vertex(
            ElementLocalVertex( well_ptr_[well] + edge, vertex ) );
    }

    template < index_t DIMENSION >
    void GeoModelMeshWells< DIMENSION >::clear()
    {
        this->set_is_initialized( false );
        auto mesh_builder =
            LineMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        mesh_builder->clear( true, false );
        well_ptr_.clear();
    }

    template < index_t DIMENSION >
    void GeoModelMeshWells< DIMENSION >::test_and_initialize() const
    {
        if( !this->is_initialized() )
        {
            const_cast< GeoModelMeshWells* >( this )->initialize();
        }
    }

    template < index_t DIMENSION >
    void GeoModelMeshWells< DIMENSION >::initialize()
    {
        this->set_is_initialized( true );

        if( !this->geomodel_.wells() )
        {
            return;
        }
        this->gmm_.vertices.test_and_initialize();
        auto mesh_builder =
            LineMeshBuilder< DIMENSION >::create_builder( *mesh_ );
        if( mesh_->nb_vertices() != this->gmm_.vertices.nb() )
        {
            copy_vertices( mesh_builder.get(), *this->gmm_.vertices.mesh_ );
        }

        // Compute the total number of edge per well
        const auto& wells = *this->geomodel_.wells();
        well_ptr_.resize( wells.nb_wells() + 1, 0 );
        index_t nb_edges{ 0 };
        for( auto w : range( wells.nb_wells() ) )
        {
            nb_edges += wells.well( w ).nb_edges();
            well_ptr_[w + 1] = nb_edges;
        }

        // Compute the edge offset
        for( auto i : range( 1, well_ptr_.size() - 1 ) )
        {
            well_ptr_[i + 1] += well_ptr_[i];
        }

        // Create edges
        mesh_builder->create_edges( well_ptr_.back() );

        // Fill edges
        index_t cur_edge{ 0 };
        for( auto w : range( 0, wells.nb_wells() ) )
        {
            const Well< DIMENSION >& well = wells.well( w );
            for( auto p : range( well.nb_parts() ) )
            {
                for( auto e : range( well.part( p ).nb_edges() ) )
                {
                    const auto& e0 = well.part( p ).edge_vertex( { e, 0 } );
                    mesh_builder->set_edge_vertex(
                        ElementLocalVertex( cur_edge, 0 ),
                        this->gmm_.vertices.index( e0 ) );
                    const auto& e1 = well.part( p ).edge_vertex( { e, 1 } );
                    mesh_builder->set_edge_vertex(
                        ElementLocalVertex( cur_edge, 1 ),
                        this->gmm_.vertices.index( e1 ) );
                    cur_edge++;
                }
            }
        }
    }

    template < index_t DIMENSION >
    AttributesManager& GeoModelMeshWells< DIMENSION >::attribute_manager() const
    {
        return mesh_->edge_attribute_manager();
    }

    template < index_t DIMENSION >
    const LineAABBTree< DIMENSION >&
        GeoModelMeshWells< DIMENSION >::aabb() const
    {
        test_and_initialize();
        return mesh_->edge_aabb();
    }

    /*******************************************************************************/

    template < index_t DIMENSION >
    GeoModelMeshBase< DIMENSION >::GeoModelMeshBase(
        GeoModelMesh< DIMENSION >& gmm, GeoModel< DIMENSION >& geomodel )
        : geomodel_( geomodel ),
          mesh_set_( new MeshSet< DIMENSION > ),
          vertices( gmm, geomodel, mesh_set_->point_set_mesh ),
          edges( gmm, geomodel, mesh_set_->line_mesh ),
          wells( gmm, geomodel, mesh_set_->well_mesh ),
          polygons( gmm, geomodel, mesh_set_->surface_mesh )
    {
    }

    template < index_t DIMENSION >
    GeoModelMeshBase< DIMENSION >::~GeoModelMeshBase()
    {
        polygons.clear_polygon_data();
    }

    template < index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::remove_colocated_vertices()
    {
        vertices.remove_colocated();
    }

    template < index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::erase_vertices(
        std::vector< index_t >& to_delete )
    {
        vertices.erase_vertices( to_delete );
    }

    template < index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::change_point_set_mesh_data_structure(
        const MeshType& type )
    {
        if( mesh_set_->point_set_mesh->type_name() != type )
        {
            vertices.clear();
            mesh_set_->create_point_set_mesh( type );
        }
    }

    template < index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::change_line_mesh_data_structure(
        const MeshType& type )
    {
        if( mesh_set_->line_mesh->type_name() != type )
        {
            wells.clear();
            mesh_set_->create_line_mesh( type );
        }
    }

    template < index_t DIMENSION >
    void GeoModelMeshBase< DIMENSION >::change_surface_mesh_data_structure(
        const MeshType& type )
    {
        if( mesh_set_->surface_mesh->type_name() != type )
        {
            polygons.clear();
            mesh_set_->create_surface_mesh( type );
        }
    }

    template < index_t DIMENSION >
    GeoModelMesh< DIMENSION >::GeoModelMesh( GeoModel< DIMENSION >& geomodel )
        : GeoModelMeshBase< DIMENSION >( *this, geomodel )
    {
    }

    GeoModelMesh< 3 >::GeoModelMesh( GeoModel3D& geomodel )
        : GeoModelMeshBase< 3 >( *this, geomodel ),
          cells( *this, geomodel, mesh_set_->volume_mesh )
    {
    }

    GeoModelMesh< 3 >::~GeoModelMesh()
    {
    }

    void GeoModelMesh< 3 >::change_volume_mesh_data_structure(
        const MeshType& type )
    {
        if( mesh_set_->volume_mesh->type_name() != type )
        {
            cells.clear();
            mesh_set_->create_volume_mesh( type );
        }
    }

    template class geomodel_core_api GeoModelMeshBase< 2 >;
    template class geomodel_core_api GeoModelMesh< 2 >;
    template class geomodel_core_api GeoModelMeshVerticesBase< 2 >;
    template class geomodel_core_api GeoModelMeshWells< 2 >;
    template class geomodel_core_api GeoModelMeshEdges< 2 >;
    template class geomodel_core_api GeoModelMeshPolygonsBase< 2 >;

    template class geomodel_core_api GeoModelMeshBase< 3 >;
    template class geomodel_core_api GeoModelMeshVerticesBase< 3 >;
    template class geomodel_core_api GeoModelMeshWells< 3 >;
    template class geomodel_core_api GeoModelMeshEdges< 3 >;
    template class geomodel_core_api GeoModelMeshPolygonsBase< 3 >;
    template class geomodel_core_api GeoModelMeshCells< 3 >;

} // namespace RINGMesh
