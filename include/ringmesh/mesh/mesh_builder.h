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

#pragma once

#include <ringmesh/basic/common.h>

#include <geogram/mesh/mesh_repair.h>

#include <ringmesh/mesh/mesh.h>

namespace RINGMesh {
    class GeoModel ;
}

namespace RINGMesh {
    class RINGMESH_API MeshBaseBuilder: public GEO::Counted {
    ringmesh_disable_copy( MeshBaseBuilder ) ;

    public:
        virtual ~MeshBaseBuilder()
        {
        }
        /*!
         * \name general methods
         * @{
         */
        virtual void copy( const MeshBase& rhs, bool copy_attributes ) = 0 ;

        virtual void load_mesh( const std::string& filename ) = 0 ;
        /*!
         * @brief Removes all the entities and attributes of this mesh.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void clear( bool keep_attributes, bool keep_memory ) = 0 ;
        /**
         * \brief Fixes some defaults in a mesh.
         * \param[in] mode a combination of #MeshRepairMode flags.
         *  Combine them with the 'bitwise or' (|) operator.
         * \param[in] colocate_epsilon tolerance used to colocate vertices
         *  (if #MESH_REPAIR_COLOCATE is set in mode).
         */
        virtual void mesh_repair(
            GEO::MeshRepairMode mode,
            double colocate_epsilon ) = 0 ;
        /*!@}
         * \name Vertex related methods
         * @{
         */
        virtual void set_vertex( index_t v_id, const vec3& vertex ) = 0 ;
        /*!
         * @brief Creates a new vertex.
         * @return the index of the created vertex
         */
        virtual index_t create_vertex() = 0 ;
        /*!
         * @brief Creates a new vertex.
         * @param[in] coords a pointer to @function dimension() coordinate.
         * @return the index of the created vertex
         */
        virtual index_t create_vertex( const vec3& vertex )
        {
            index_t index = create_vertex() ;
            set_vertex( index, vertex ) ;
            return index ;
        }
        /*!
         * @brief Creates a contiguous chunk of vertices.
         * @param[in] nb number of sub-entities to create.
         * @return the index of the first created vertex
         */
        virtual index_t create_vertices( index_t nb ) = 0 ;
        /*!
         * @brief Deletes a set of vertices.
         * @param[in] to_delete     a vector of size @function nb(). If to_delete[e] is true,
         * then entity e will be destroyed, else it will be kept.
         */
        virtual void delete_vertices( const std::vector< bool >& to_delete ) = 0 ;
        /*!
         * @brief Removes all the vertices and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void clear_vertices( bool keep_attributes, bool keep_memory ) = 0 ;
        virtual void permute_vertices(
            const std::vector< index_t >& permutation ) = 0 ;
        /*!
         * @brief Deletes the NNSearch on vertices
         */
        virtual void clear_vertex_linked_objects() = 0 ;
        /*!@}
         */

        static MeshBaseBuilder* create_builder( MeshBase& mesh ) ;
    protected:
        MeshBaseBuilder()
        {
        }
        virtual void set_mesh( MeshBase& mesh ) = 0 ;
    } ;
    typedef GEO::SmartPointer< MeshBaseBuilder > MeshBaseBuilder_var ;
    typedef GEO::Factory0< MeshBaseBuilder > MeshBaseBuilderFactory ;
#define ringmesh_register_mesh_base_builder(type) \
    geo_register_creator(RINGMesh::MeshBaseBuilderFactory, type ## Builder, type::type_name_static())

    class RINGMESH_API Mesh0DBuilder: public virtual MeshBaseBuilder {
    ringmesh_disable_copy( Mesh0DBuilder ) ;
    public:
        virtual ~Mesh0DBuilder()
        {
        }

        virtual void set_mesh( Mesh0D& mesh ) = 0 ;

        static Mesh0DBuilder* create_builder( Mesh0D& mesh ) ;

        virtual void remove_isolated_vertices()
        {
            // All vertices are isolated in a Mesh0D
        }

    protected:
        Mesh0DBuilder()
            : MeshBaseBuilder()
        {
        }
    } ;
    typedef GEO::SmartPointer< Mesh0DBuilder > Mesh0DBuilder_var ;
    typedef GEO::Factory0< Mesh0DBuilder > Mesh0DBuilderFactory ;
#define ringmesh_register_mesh_0d_builder(type) \
    geo_register_creator(RINGMesh::Mesh0DBuilderFactory, type ## Builder, type::type_name_static())

    class RINGMESH_API Mesh1DBuilder: public virtual MeshBaseBuilder {
    ringmesh_disable_copy( Mesh1DBuilder ) ;
    public:
        virtual ~Mesh1DBuilder()
        {
        }

        virtual void set_mesh( Mesh1D& mesh ) = 0 ;

        static Mesh1DBuilder* create_builder( Mesh1D& mesh ) ;

        /*!
         * @brief Create a new edge.
         * @param[in] v1_id index of the starting vertex.
         * @param[in] v2_id index of the ending vertex.
         */
        virtual void create_edge( index_t v1_id, index_t v2_id ) = 0 ;
        /*!
         * \brief Creates a contiguous chunk of edges
         * \param[in] nb_edges number of edges to create
         * \return the index of the first edge
         */
        virtual index_t create_edges( index_t nb_edges ) = 0 ;
        /*!
         * @brief Sets a vertex of a facet by local vertex index.
         * @param[in] edge_id index of the edge, in 0..nb()-1.
         * @param[in] local_vertex_id index of the vertex in the facet.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of facet
         * \param of the facet facet_id. Index between 0 and @function nb() - 1.
         */
        virtual void set_edge_vertex(
            index_t edge_id,
            index_t local_vertex_id,
            index_t vertex_id ) = 0 ;
        /*!
         * @brief Deletes a set of edges.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it will be kept.
         * @param[in] remove_isolated_vertices if true, then the vertices
         * that are no longer incident to any entity are deleted.
         */
        virtual void delete_edges(
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) = 0 ;
        /*!
         * @brief Removes all the edges and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void clear_edges( bool keep_attributes, bool keep_memory ) = 0 ;
        virtual void permute_edges( const std::vector< index_t >& permutation ) = 0 ;

        virtual void clear_edge_linked_objects() = 0 ;
        /*!
         * @brief Remove vertices not connected to any mesh element
         */
        virtual void remove_isolated_vertices() = 0 ;
    protected:
        Mesh1DBuilder()
            : MeshBaseBuilder()
        {
        }
    } ;
    typedef GEO::SmartPointer< Mesh1DBuilder > Mesh1DBuilder_var ;
    typedef GEO::Factory0< Mesh1DBuilder > Mesh1DBuilderFactory ;
#define ringmesh_register_mesh_1d_builder(type) \
    geo_register_creator(RINGMesh::Mesh1DBuilderFactory, type ## Builder, type::type_name_static())

    class RINGMESH_API Mesh2DBuilder: public virtual MeshBaseBuilder {
    ringmesh_disable_copy( Mesh2DBuilder ) ;
    public:
        virtual ~Mesh2DBuilder()
        {
        }

        virtual void set_mesh( Mesh2D& mesh ) = 0 ;

        static Mesh2DBuilder* create_builder( Mesh2D& mesh ) ;

        /*!@}
         * \name Facet related methods
         * @{
         */
        /*!
         * brief create facet polygons
         * @param[in] facets is the vector of vertex index for each facet
         * @param[in] facet_ptr is the vector addressing the first facet vertex for each facet.
         */
        virtual void create_facet_polygons(
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) = 0 ;
        /*!
         * \brief Creates a polygonal facet
         * \param[in] vertices a const reference to a vector that
         *  contains the vertices
         * \return the index of the created facet
         */
        virtual index_t create_facet_polygon(
            const std::vector< index_t >& vertices ) = 0 ;
        /*!
         * \brief Creates a contiguous chunk of triangles
         * \param[in] nb_triangles number of triangles to create
         * \return the index of the first triangle
         */
        virtual index_t create_facet_triangles( index_t nb_triangles ) = 0 ;
        /*!
         * \brief Creates a contiguous chunk of quads
         * \param[in] nb_quads number of quads to create
         * \return the index of the first quad
         */
        virtual index_t create_facet_quads( index_t nb_quads ) = 0 ;
        /*!
         * @brief Sets a vertex of a facet by local vertex index.
         * @param[in] facet_id index of the facet, in 0.. @function nb() - 1.
         * @param[in] local_vertex_id index of the vertex in the facet.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of the
         * facet \param facet_id. Index between 0 and @function nb() - 1.
         */
        virtual void set_facet_vertex(
            index_t facet_id,
            index_t local_vertex_id,
            index_t vertex_id ) = 0 ;
        /*!
         * @brief Sets an adjacent facet by both its facet \param facet_id
         * and its local edge index \param edge_id.
         * @param[in] facet_id the facet index
         * @param[in] edge_id the local index of an edge in facet \p facet_id
         * @param[in] specifies the facet adjacent to \param facet_id along edge
         * \param edge_id or GEO::NO_FACET if the parameter \param edge_id is
         * on the border.
         */
        virtual void set_facet_adjacent(
            index_t facet_id,
            index_t edge_id,
            index_t specifies ) = 0 ;
        /*
         * \brief Copies a triangle mesh into this Mesh.
         * \details Facet adjacence are not computed.
         *   Facet and corner attributes are zeroed.
         * \param[in] triangles facet to vertex links
         * \param[in] steal_args if set, vertices and triangles
         * are 'stolen' from the arguments
         * (using vector::swap).
         */
        virtual void assign_facet_triangle_mesh(
            const std::vector< index_t >& triangles,
            bool steal_args ) = 0 ;
        /*!
         * @brief Removes all the facets and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void clear_facets( bool keep_attributes, bool keep_memory ) = 0 ;
        /*!
         * @brief Retrieve the adjacencies of facets
         */
        virtual void connect_facets() = 0 ;
        virtual void permute_facets(
            const std::vector< index_t >& permutation ) = 0 ;
        /*!
         * @brief Deletes a set of facets.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it will be kept.
         * @param[in] remove_isolated_vertices if true, then the vertices that are
         * no longer incident to any entity are deleted.
         */
        virtual void delete_facets(
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) = 0 ;

        virtual void clear_facet_linked_objects() = 0 ;
        /*!@}
         * \name Mesh2D algorithms
         * @{
         */
        /**
         * \brief Removes the connected components that have an area
         *  smaller than a given threshold.
         * \param[in] min_area the connected components with an
         *  area smaller than this threshold are removed
         * \param[in] min_facets the connected components with
         *  less than \param min_facets facets are removed
         */
        virtual void remove_small_connected_components(
            double min_area,
            index_t min_facets ) = 0 ;
        virtual void triangulate( const Mesh2D& surface_in ) = 0 ;
        /*!@}
         */
        /*!
         * @brief Remove vertices not connected to any mesh element
         */
        virtual void remove_isolated_vertices() = 0 ;
    protected:
        Mesh2DBuilder()
            : MeshBaseBuilder()
        {
        }
    } ;
    typedef GEO::SmartPointer< Mesh2DBuilder > Mesh2DBuilder_var ;
    typedef GEO::Factory0< Mesh2DBuilder > Mesh2DBuilderFactory ;
#define ringmesh_register_mesh_2d_builder(type) \
    geo_register_creator(RINGMesh::Mesh2DBuilderFactory, type ## Builder, type::type_name_static())

    class RINGMESH_API Mesh3DBuilder: public virtual MeshBaseBuilder {
    ringmesh_disable_copy( Mesh3DBuilder ) ;
    public:
        virtual ~Mesh3DBuilder()
        {
        }

        virtual void set_mesh( Mesh3D& mesh ) = 0 ;

        static Mesh3DBuilder* create_builder( Mesh3D& mesh ) ;

        /*!
         * @brief Creates a contiguous chunk of cells of the same type.
         * @param[in] nb_cells number of cells to create
         * @param[in] type type of the cells to create, one of GEO::MESH_TET, GEO::MESH_HEX,
         * GEO::MESH_PRISM, GEO::MESH_PYRAMID, GEO::MESH_CONNECTOR.
         * @return the first created cell.
         */
        virtual index_t create_cells(
            index_t nb_cells,
            GEO::MeshCellType type ) = 0 ;
        /*
         * \brief Copies a tets mesh into this Mesh.
         * \details Cells adjacence are not computed.
         *   cell and corner attributes are zeroed.
         * \param[in] tets cells to vertex links
         * \param[in] steal_args if set, vertices and tets
         * are 'stolen' from the arguments
         * (using vector::swap).
         */
        virtual void assign_cell_tet_mesh(
            const std::vector< index_t >& tets,
            bool steal_args ) = 0 ;
        /*!
         * @brief Sets a vertex of a cell by local vertex index.
         * @param[in] cell_id index of the cell, in 0.. @function nb() - 1.
         * @param[in] local_vertex_id index of the vertex in the cell.
         * Local index between 0 and @function nb_vertices(cell_id) - 1.
         * @param[in] vertex_id specifies the global index of the vertex \param
         * local_vertex_id in the cell \param cell_id. Index between 0 and
         * @function nb() - 1.
         */
        virtual void set_cell_vertex(
            index_t cell_id,
            index_t local_vertex_id,
            index_t vertex_id ) = 0 ;
        /*!
         * \brief Sets the vertex that a corner is incident to
         * \param[in] corner_index the corner, in 0.. @function nb() - 1
         * \param[in] vertex_index specifies the vertex that corner
         * \param corner_index is incident to
         */
        virtual void set_cell_corner_vertex_index(
            index_t corner_index,
            index_t vertex_index ) = 0 ;
        /*!
         * \brief Sets the cell adjacent
         * \param[in] cell_index index of the cell
         * \param[in] facet_index local index of the cell facet
         * \param[in] cell_adjacent adjacent value to set
         */
        virtual void set_cell_adjacent(
            index_t cell_index,
            index_t facet_index,
            index_t cell_adjacent ) = 0 ;

        /*!
         * @brief Retrieve the adjacencies
         */
        virtual void connect_cells() = 0 ;

        /*!
         * @brief Removes all the cells and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        virtual void clear_cells( bool keep_attributes, bool keep_memory ) = 0 ;
        /*!
         * @brief Applies a permutation to the entities and their attributes.
         * On exit, permutation is modified (used for internal bookkeeping).
         * Applying a permutation permutation is equivalent to:
         * <code>
         *  for( i = 0 ; i < permutation.size() ; i++) {
         *      data2[i] = data[permutation[i]]
         *       }
         *  data = data2 ;
         *  </code>
         */
        virtual void permute_cells( const std::vector< index_t >& permutation ) = 0 ;
        /*!
         * @brief Deletes a set of cells.
         * @param[in] to_delete     a vector of size @function nb().
         * If to_delete[e] is true, then entity e will be destroyed, else it will be kept.
         * @param[in] remove_isolated_vertices if true, then the vertices that are
         * no longer incident to any entity are deleted.
         */
        virtual void delete_cells(
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) = 0 ;

        virtual void clear_cell_linked_objects() = 0 ;
        /*!
         * @brief Remove vertices not connected to any mesh element
         */
        virtual void remove_isolated_vertices() = 0 ;
    protected:
        Mesh3DBuilder()
            : MeshBaseBuilder()
        {
        }
    } ;
    typedef GEO::SmartPointer< Mesh3DBuilder > Mesh3DBuilder_var ;
    typedef GEO::Factory0< Mesh3DBuilder > Mesh3DBuilderFactory ;
#define ringmesh_register_mesh_3d_builder(type) \
    geo_register_creator(RINGMesh::Mesh3DBuilderFactory, type ## Builder, type::type_name_static())

    class RINGMESH_API MeshAllDBuilder: public virtual Mesh0DBuilder,
        public virtual Mesh1DBuilder,
        public virtual Mesh2DBuilder,
        public virtual Mesh3DBuilder {
    ringmesh_disable_copy( MeshAllDBuilder ) ;
    public:
        virtual ~MeshAllDBuilder()
        {
        }

        virtual void set_mesh( MeshAllD& mesh ) = 0 ;

        virtual void remove_isolated_vertices() = 0 ;

        static MeshAllDBuilder* create_builder( MeshAllD& mesh ) ;

    protected:
        MeshAllDBuilder()
            : Mesh0DBuilder(), Mesh1DBuilder(), Mesh2DBuilder(), Mesh3DBuilder()
        {
        }
    } ;
    typedef GEO::SmartPointer< MeshAllDBuilder > MeshAllDBuilder_var ;
    typedef GEO::Factory0< MeshAllDBuilder > MeshAllDBuilderFactory ;
#define ringmesh_register_mesh_alld_builder(type) \
    geo_register_creator(RINGMesh::MeshAllDBuilderFactory, type ## Builder, type::type_name_static())
}
