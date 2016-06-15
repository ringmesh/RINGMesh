/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#ifndef __RINGMESH_MESH__
#define __RINGMESH_MESH__

#include <ringmesh/common.h>
#include <ringmesh/geometry.h>
#include <ringmesh/geogram_extension.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_AABB.h>

namespace RINGMesh {
    class GeoModel ;
}

namespace RINGMesh {

    /*!
     * @brief class to encapsulate mesh structure in order to provide an API
     * on which we base the RINGMesh algorithm 
     * @note For now, we encapsulate the GEO::Mesh class. We can develop the concept
     * using a factory to build several encapsulating classes. 
     */
    class RINGMESH_API Mesh {
    ringmesh_disable_copy( Mesh ) ;
        friend class MeshBuilder ;
    
    public:
        /*!
         * @brief Mesh constructor.
         * @param[in] dimension dimension of the vertices.
         * @param[in] single_precision if true, vertices are stored in single precision (float),
         * else they are stored as double precision (double)..
         */
        Mesh( const GeoModel& geo_model, index_t dimension, bool single_precision )
            : geo_model_( geo_model ), facets_aabb_( nil ), cells_aabb_( nil )
        {
            mesh_ = new GEO::Mesh( dimension, single_precision ) ;
            for( index_t i = 0; i < ColocaterANN::NB_LOCATION; i++ ) {
                ann_[i] = nil ;
            }
        }
        ~Mesh()
        {
            if( facets_aabb_ ) delete facets_aabb_ ;
            if( cells_aabb_ ) delete cells_aabb_ ;
            for( index_t i = 0; i < ColocaterANN::NB_LOCATION; i++ ) {
                if( ann_[i] ) delete ann_[i] ;
            }
            delete mesh_ ;
        }

        /*!
         * @brief Copy a mesh into this one.
         * @param[in] rhs a const reference to the mesh to be copied.
         * @param[in] copy_attributes if true, all attributes are copied.
         * @param[in] what a combination of MESH_VERTICES, MESH_EDGES, MESH_FACETS, MESH_CELLS flags. 
         * Set to MESH_ALL_ELEMENTS to copy everything (default). 
         * If MESH_VERTICES is not set, then the mesh is cleared.
         * @return a modifiable reference to the point that corresponds to the vertex.
         */
        void copy(
            const Mesh& rhs,
            bool copy_attributes,
            GEO::MeshElementsFlags what ) const
        {
            mesh_->copy( *rhs.mesh_, copy_attributes, what ) ;
        }

        void save_mesh(
            const std::string filename,
            const GEO::MeshIOFlags& ioflags ) const
        {
            GEO::mesh_save( *mesh_, filename, ioflags ) ;
        }
        /*!
         * @brief return the ColocaterANN at the given ColocaterANN::MeshLocation
         * @warning the ColocaterANN is destroy when calling the Mesh::facets_aabb() and Mesh::cells_aabb()
         */
        const ColocaterANN& colotater_ann(
            ColocaterANN::MeshLocation location ) const
        {
            if( ann_[location] == nil ) {
                ann_[location] = new ColocaterANN( *mesh_, location ) ;
            }
            return *ann_[location] ;
        }

        /*!
         * get access to GEO::MESH... only for GFX..
         * @todo Remove this function as soon as the GEO::MeshGFX is encapsulated
         */
        const GEO::Mesh& gfx_mesh() const
        {
            return *mesh_ ;
        }

        index_t nb_connected_components() const
        {
            return GEO::mesh_nb_connected_components( *mesh_ ) ;
        }

        //TODO maybe reimplement the function with a RINGMesh::Mesh??
        void print_mesh_bounded_attributes()
        {
            print_bounded_attributes( *mesh_ ) ;
        }

        /*!
         * \name Vertex methods
         * @{
         */
        /*
         * @brief Gets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @return const reference to the point that corresponds to the vertex.
         */
        const vec3& vertex( index_t v_id ) const
        {
            ringmesh_assert( v_id < nb_vertices() ) ;
            return mesh_->vertices.point( v_id ) ;
        }

        /*
         * @brief Gets the number of point in the Mesh.
         */
        index_t nb_vertices() const
        {
            return mesh_->vertices.nb() ;
        }
        GEO::AttributesManager& vertex_attribute_manager() const
        {
            return mesh_->vertices.attributes() ;
        }

        /*!
         *  @}
         * \name Edge methods
         * @{
         */
        /*
         * @brief Gets the index of an edge vertex.
         * @param[in] edge_id index of the edge.
         * @param[in] vertex_id local index of the vertex, in {0,1}
         * @return the global index of vertex \param vertex_id in edge \param edge_id.
         */
        index_t edge_vertex( index_t edge_id, index_t vertex_id ) const
        {
            return mesh_->edges.vertex( edge_id, vertex_id ) ;
        }
        /*!
         * @brief Gets the number of all the edges in the whole Mesh.
         */
        index_t nb_edges() const
        {
            return mesh_->edges.nb() ;
        }
        /*!
         * @brief Gets the length of the edge \param edge_id
         */
        double edge_length( index_t edge_id ) const
        {
            const vec3& e0 = vertex( edge_vertex( edge_id, 0 ) ) ;
            const vec3& e1 = vertex( edge_vertex( edge_id, 1 ) ) ;
            return GEO::Geom::distance( e0, e1 ) ;
        }
        GEO::AttributesManager& edge_attribute_manager() const
        {
            return mesh_->edges.attributes() ;
        }

        /*! @}
         * \name Facets methods
         * @{
         */
        /*!
         * @brief Gets the vertex index by facet index and local vertex index.
         * @param[in] facet_id the facet index.
         * @param[in] vertex_id the local edge index in \param facet_id.
         * @return the global facet index adjacent to the \param edge_id of the facet \param facet_id.
         * @precondition  \param edge_id < number of edge of the facet \param facet_id .
         */
        index_t facet_vertex( index_t facet_id, index_t vertex_id ) const
        {
            return mesh_->facets.vertex( facet_id, vertex_id ) ;
        }
        /*!
         * @brief return the vertex index of the corner \param corner_id
         */
        index_t facet_corner_vertex( index_t corner_id ) const
        {
            return mesh_->facet_corners.vertex( corner_id ) ;
        }
        /*!
         * @brief Gets the number of all facets in the whole Mesh.
         */
        index_t nb_facets() const
        {
            return mesh_->facets.nb() ;
        }
        /*!
         * @brief Gets the number of vertices in the facet \param facet_id.
         * @param[in] facet_id facet index
         */
        index_t nb_facet_vertices( index_t facet_id ) const
        {
            return mesh_->facets.nb_vertices( facet_id ) ;
        }
        /*!
         * Get the number of corners of all facets
         */
        index_t nb_facet_corners() const
        {
            return mesh_->facet_corners.nb() ;
        }
        /*!
         * @brief Get the first vertex index of a facet.
         * @param[in] facet_id facet index
         */
        index_t facet_begin( index_t facet_id ) const
        {
            return mesh_->facets.corners_begin( facet_id ) ;
        }
        /*!
         * @brief Get the last vertex index of a facet.
         * @param[in] facet_id facet index
         */
        index_t facet_end( index_t facet_id ) const
        {
            return mesh_->facets.corners_end( facet_id ) ;
        }
        /*!
         * @brief Gets the next vertex index in the facet \param facet_id.
         * @param[in] facet_id facet index
         * @param[in] vertex_id current index
         */
        index_t next_facet_vertex( index_t facet_id, index_t vertex_id ) const
        {
            ringmesh_assert( vertex_id < nb_facet_vertices( facet_id ) ) ;
            if( vertex_id != nb_facet_vertices( facet_id ) - 1 ) {
                return vertex_id + 1 ;
            } else {
                return 0 ;
            }
        }
        /*!
         * @brief Gets the previous vertex index in the facet \param facet_id.
         * @param[in] facet_id facet index
         * @param[in] vertex_id current index
         */
        index_t prev_facet_vertex( index_t facet_id, index_t vertex_id ) const
        {
            ringmesh_assert( vertex_id < nb_facet_vertices( facet_id ) ) ;
            if( vertex_id > 0 ) {
                return vertex_id - 1 ;
            } else {
                return nb_facet_vertices( facet_id ) - 1 ;
            }
        }
        /*!
         * @brief Gets an adjacent facet index by facet index and local edge index.
         * @param[in] facet_id the facet index.
         * @param[in] edge_id the local edge index in \param facet_id.
         * @return the global facet index adjacent to the \param edge_id of the facet \param facet_id.
         * @precondition  \param edge_id < number of edge of the facet \param facet_id .
         */
        index_t facet_adjacent( index_t facet_id, index_t edge_id ) const
        {
            return mesh_->facets.adjacent( facet_id, edge_id ) ;
        }
        GEO::AttributesManager& facet_attribute_manager() const
        {
            return mesh_->facets.attributes() ;
        }
        /*!
         * @brief Tests whether all the facets are triangles. when all the facets are triangles, storage and access is optimized.
         * @return True if all facets are triangles and False otherwise.
         */
        bool facets_are_simplicies() const
        {
            return mesh_->facets.are_simplices() ;
        }
        /*!
         * return true if the facet \param facet_id is a triangle
         */
        bool is_triangle( index_t facet_id ) const
        {
            return nb_facet_vertices( facet_id ) == 3 ;
        }
        /*!
         * @brief Create an AABB tree for a Mesh facets
         * @pre The GeoModelEntity must be simplicial
         * @warning SIDE EFFECTS: The mesh vertices are reordered.
         * @warning calling this function will destroy the ColocaterANN.
         */
        const GEO::MeshFacetsAABB& facets_aabb() const ;
        /*!
         * Computes the Mesh facet normal
         * @param[in] facet_id the facet index
         * @return the facet normal
         */
        vec3 facet_normal( index_t facet_id ) const
        {
            return normalize( GEO::Geom::mesh_facet_normal( *mesh_, facet_id ) ) ;
        }
        /*!
         * Computes the Mesh facet barycenter
         * @param[in] facet_id the facet index
         * @return the facet center
         */
        vec3 facet_barycenter( index_t facet_id ) const
        {
            return GEO::Geom::mesh_facet_center( *mesh_, facet_id ) ;
        }
        /*!
         * Computes the Mesh facet area
         * @param[in] facet_id the facet index
         * @return the facet area
         */
        double facet_area( index_t facet_id ) const
        {
            return GEO::Geom::mesh_facet_area( *mesh_, facet_id ) ;
        }
        /*! @}
         * \name Cells methods
         * @{
         */

        /*!
         * @brief Gets a vertex index by cell and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        index_t cell_vertex( index_t cell_id, index_t vertex_id ) const
        {
            return mesh_->cells.vertex( cell_id, vertex_id ) ;
        }
        /*!
         * @brief Gets a vertex index by cell and local edge and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] edge_id the local edge index in \param cell_id.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        index_t cell_edge_vertex(
            index_t cell_id,
            index_t edge_id,
            index_t vertex_id ) const
        {
            return mesh_->cells.edge_vertex( cell_id, edge_id, vertex_id ) ;
        }
        /*!
         * @brief Gets a vertex by cell facet and local vertex index.
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @param[in] vertex_id index of the vertex in the facet \param facet_id
         * @return the global vertex index.
         * @precondition vertex_id < number of vertices in the facet \param facet_id 
         * and facet_id number of facet in th cell \param cell_id
         */
        index_t cell_facet_vertex(
            index_t cell_id,
            index_t facet_id,
            index_t vertex_id ) const
        {
            return mesh_->cells.facet_vertex( cell_id, facet_id, vertex_id ) ;
        }
        /*!
         * @brief Gets a facet index by cell and local facet index.
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @return the global facet index.
         */
        index_t cell_facet( index_t cell_id, index_t facet_id ) const {
            return mesh_->cells.facet( cell_id,facet_id ) ;
        }
        /*!
         * Get the number of corners of all cells
         */
        index_t cell_corner_vertex( index_t corner_id ) const {
            return mesh_->cell_corners.vertex( corner_id ) ;
        }
        /*!
         * @brief Gets the number of facet in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        index_t nb_cell_facets( index_t cell_id ) const
        {
            return mesh_->cells.nb_facets( cell_id ) ;
        }
        /*!
         * @brief Gets the number of edges in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        index_t nb_cell_edges( index_t cell_id ) const
        {
            return mesh_->cells.nb_edges( cell_id ) ;
        }
        /*!
         * @brief Gets the number of vertices of a facet in a cell
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @return the number of vertices in the facet \param facet_id in the cell \param cell_id
         */
        index_t nb_cell_facet_vertices( index_t cell_id, index_t facet_id ) const
        {
            return mesh_->cells.facet_nb_vertices( cell_id, facet_id ) ;
        }
        /*!
         * @brief Gets the number of vertices of a cell
         * @param[in] cell_id index of the cell
         * @return the number of vertices in the cell \param cell_id
         */
        index_t nb_cell_vertices( index_t cell_id ) const
        {
            return mesh_->cells.nb_vertices( cell_id ) ;
        }
        /*!
         * @brief Gets the number of cells in the Mesh.
         */
        index_t nb_cells() const
        {
            return mesh_->cells.nb() ;
        }
        index_t nb_cell_corners()const {
            return mesh_->cell_corners.nb() ;
        }
        /*!
         * @return the index of the adjacent cell of \param cell_id along the facet \param facet_id
         */
        index_t cell_adjacent( index_t cell_id, index_t facet_id ) const
        {
            return mesh_->cells.adjacent( cell_id, facet_id ) ;
        }
        GEO::AttributesManager& cell_attribute_manager() const
        {
            return mesh_->cells.attributes() ;
        }
        /*!
         * @brief Gets the type of a cell.
         * @param[in] cell_id the cell index, in 0..nb()-1
         */
        GEO::MeshCellType cell_type( index_t cell_id ) const
        {
            return mesh_->cells.type( cell_id ) ;
        }
        /*!
         * @brief Tests whether all the cells are tetrahedra. when all the cells are tetrahedra, storage and access is optimized.
         * @return True if all cells are tetrahedra and False otherwise.
         */
        bool cells_are_simplicies() const
        {
            return mesh_->cells.are_simplices() ;
        }
        /*!
         * @brief Create an AABB tree for a Mesh cells
         * @pre The GeoModelEntity must be simplicial
         * @warning SIDE EFFECTS: The mesh vertices are reordered.
         * @warning calling this function will destroy the ColocaterANN.
         */
        const GEO::MeshCellsAABB& cells_aabb() const ;
        /*!
         * Computes the Mesh cell facet barycenter
         * @param[in] cell_id the cell index
         * @param[in] facet_id the facet index in the cell
         * @return the cell facet center
         */
        vec3 cell_facet_barycenter( index_t cell_id, index_t facet_id ) const
        {
            vec3 result( 0., 0., 0. ) ;
            index_t nb_vertices = nb_cell_facet_vertices( cell_id, facet_id ) ;
            for( index_t v = 0; v < nb_vertices; ++v ) {
                result += vertex( cell_facet_vertex( cell_id, facet_id, v ) ) ;
            }
            ringmesh_assert( nb_vertices > 0 ) ;

            return result / nb_vertices ;
        }
        /*!
         * Compute the non weighted barycenter of the \param cell_id
         */
        vec3 cell_barycenter( index_t cell_id ) const
        {
            return RINGMesh::mesh_cell_center( *mesh_, cell_id ) ;
        }
        /*!
         * Computes the Mesh cell facet normal
         * @param[in] cell_id the cell index
         * @param[in] facet_id the facet index in the cell
         * @return the cell facet normal
         */
        vec3 cell_facet_normal( index_t cell_id, index_t facet_id ) const
        {
            return GEO::mesh_cell_facet_normal( *mesh_, cell_id, facet_id ) ;
        }
        /*!
         * @brief compute the volume of the cell \param cell_id.
         */
        double cell_volume( index_t cell_id ) const
        {
            return RINGMesh::mesh_cell_volume( *mesh_, cell_id ) ;
        }

        index_t cell_begin(index_t cell_id )const{
            return mesh_->cells.corners_begin( cell_id ) ;
        }
        index_t cell_end( index_t cell_id )const{
            return mesh_->cells.corners_end( cell_id ) ;
        }
        index_t find_cell_corner( index_t cell_id, index_t vertex_id ) const {
            for( index_t v=0 ; v < nb_cell_vertices( cell_id ) ; ++v ){
                if( cell_vertex( cell_id, v ) == vertex_id ){
                    return cell_begin( cell_id ) + v ;
                }
            }
            return NO_ID ;
        }
        /*!
         * @}
         */

    private:
        mutable GEO::Mesh* mesh_ ;
        const GeoModel& geo_model_ ;

        mutable GEO::MeshFacetsAABB* facets_aabb_ ;
        mutable GEO::MeshCellsAABB* cells_aabb_ ;
        mutable ColocaterANN* ann_[ColocaterANN::NB_LOCATION] ;

    } ;

    class RINGMESH_API MeshBuilder {
    ringmesh_disable_copy( MeshBuilder ) ;

    public:
        MeshBuilder( Mesh& mesh )
            : mesh_( mesh )
        {
        }
        ~MeshBuilder() {} ;

        void load_mesh(
            const std::string& filename,
            const GEO::MeshIOFlags& ioflags )
        {
            GEO::mesh_load( filename, *mesh_.mesh_, ioflags ) ;
        }
        /*!
         * @brief Removes all the entities and attributes of this mesh.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear( bool keep_attributes, bool keep_memory )
        {
            mesh_.mesh_->clear( keep_attributes, keep_memory ) ;
        }
        /**
         * \brief Fixes some defaults in a mesh.
         * \param[in] mode a combination of #MeshRepairMode flags.
         *  Combine them with the 'bitwise or' (|) operator.
         * \param[in] colocate_epsilon tolerance used to colocate vertices
         *  (if #MESH_REPAIR_COLOCATE is set in mode).
         */
        void mesh_repair( GEO::MeshRepairMode mode, double colocate_epsilon )
        {
            GEO::mesh_repair( *mesh_.mesh_, mode, colocate_epsilon ) ;

        }
        /**
         * \brief Removes the connected components that have an area
         *  smaller than a given threshold.
         * \param[in] min_component_area the connected components with an
         *  area smaller than this threshold are removed
         * \param[in] min_component_facets the connected components with
         *  less than min_component_facets facets are removed
         */
        void remove_small_connected_components( double min_area, index_t min_facets )
        {
            GEO::remove_small_connected_components( *mesh_.mesh_, min_area,
                min_facets ) ;
        }

        /*!
         * \name Vertex methods
         * @{
         */
        /*!
         * @brief Gets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @return reference to the point that corresponds to the vertex.
         */
        vec3& vertex( index_t v_id )
        {
            return mesh_.mesh_->vertices.point( v_id ) ;
        }
        /*!
         * \brief Gets a point
         * \param[in] v_id the vertex, in 0..nb()-1
         * \return a pointer to the coordinates of the point
         *  that correspond to the vertex
         * \pre !single_precision()
         */
        double* point_ptr( index_t v_id )
        {
            return mesh_.mesh_->vertices.point_ptr( v_id ) ;
        }
        /*!
         * @brief Creates a new vertex.
         * @return the index of the created vertex
         */
        index_t create_vertex()
        {
            return mesh_.mesh_->vertices.create_vertex() ;
        }
        /*!
         * @brief Creates a new vertex.
         * @param[in] coords a pointer to dimension() coordinate.
         * @return the index of the created vertex
         */
        index_t create_vertex( const double* coords )
        {
            return mesh_.mesh_->vertices.create_vertex( coords ) ;
        }
        /*!
         * @brief Creates a contiguous chunk of vertices.
         * @param[in] nb number of sub-entities to create.
         * @return the index of the first created vertex
         */
        index_t create_vertices( index_t nb )
        {
            return mesh_.mesh_->vertices.create_vertices( nb ) ;
        }
        /*!
         * @brief Deletes a set of vertices.
         * @param[in] to_delete     a vector of size nb(). If to_delete[e] is different from 0,
         * then entity e will be destroyed, else it will be kept. On exit, to_delete is modified
         * (it is used for internal bookkeeping).
         * @param[in] remove_isolated_vertices if true, then the vertices that are no longer incident to any entity are deleted.
         */
        void delete_vertices(
            GEO::vector< index_t >& to_delete,
            bool remove_isolated_vertices )
        {
            mesh_.mesh_->vertices.delete_elements( to_delete,
                remove_isolated_vertices ) ;
            delete_vertex_colocater() ;
        }
        /*!
         * @brief Removes all the vertices and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_vertices( bool keep_attributes, bool keep_memory )
        {
            mesh_.mesh_->vertices.clear( keep_attributes, keep_memory ) ;
            delete_vertex_colocater() ;
        }
        /*!
         * @brief Deletes the ColocaterANN on vertices
         */
        void delete_vertex_colocater()
        {
            if( mesh_.ann_[ColocaterANN::VERTICES] ) {
                delete mesh_.ann_[ColocaterANN::VERTICES] ;
                mesh_.ann_[ColocaterANN::VERTICES] = nil ;
            }
        }

        /*!@}
         * \section Edge methods
         * @{
         */
        /*!
         * @brief Create a new edge.
         * @param[in] v1_id index of the starting vertex.
         * @param[in] v2_id index of the ending vertex.
         */
        void create_edge( index_t v1_id, index_t v2_id )
        {
            mesh_.mesh_->edges.create_edge( v1_id, v2_id ) ;
        }
        /*!
         * \brief Creates a contiguous chunk of edges
         * \param[in] nb_edges number of edges to create
         * \return the index of the first edge
         */
        index_t create_edges( index_t nb_edges )
        {
           return mesh_.mesh_->edges.create_edges( nb_edges ) ;
        }
        /*!
         * @brief Sets a vertex of a facet by local vertex index.
         * @param[in] edge_id index of the edge, in 0..nb()-1.
         * @param[in] local_vertex_id index of the vertex in the facet. Local index between 0 and nb_vertices(cell_id)-1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of facet \param of the facet facet_id. Index between 0 and nb()-1.
         */
        void set_edge_vertex(
            index_t edge_id,
            index_t local_vertex_id,
            index_t vertex_id )
        {
            mesh_.mesh_->edges.set_vertex( edge_id, local_vertex_id, vertex_id ) ;
        }
        /*!
         * @brief Deletes a set of edges.
         * @param[in] to_delete a vector of size nb(). If to_delete[e] is different from 0,
         * then entity e will be destroyed, else it will be kept. On exit, to_delete is modified
         * (it is used for internal bookkeeping).
         * @param[in] remove_isolated_vertices if true, then the vertices that are no longer incident to any entity are deleted.
         */
        void delete_edges(
            GEO::vector< index_t > to_delete,
            bool remove_isolated_vertices )
        {
            mesh_.mesh_->edges.delete_elements( to_delete,
                remove_isolated_vertices ) ;
            delete_edge_colocater() ;
        }
        /*!
         * @brief Removes all the edges and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_edges( bool keep_attributes, bool keep_memory ) {
            mesh_.mesh_->edges.clear( keep_attributes, keep_memory ) ;
            delete_edge_colocater() ;
        }
        /*!
         * @brief Deletes the ColocaterANN on edges
         */
        void delete_edge_colocater()
        {
            if( mesh_.ann_[ColocaterANN::EDGES] ) {
                delete mesh_.ann_[ColocaterANN::EDGES] ;
                mesh_.ann_[ColocaterANN::EDGES] = nil ;
            }
        }

        /*!@}
         * \section Facet methods
         * @{
         */
        /*!
         * brief create facet polygons
         * @param[in] facets is the vector of vertex index for each facet
         * @param[in] facet_ptr is the vector addressing the first facet vertex for each facet.
         */
        void create_facet_polygons(
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr )
        {
            for( index_t f = 0; f + 1 < facet_ptr.size(); f++ ) {
                index_t start = facet_ptr[f] ;
                index_t end = facet_ptr[f + 1] ;
                GEO::vector< index_t > facet_vertices ;
                copy_std_vector_to_geo_vector( facets, start, end, facet_vertices ) ;

                mesh_.mesh_->facets.create_polygon( facet_vertices ) ;
            }
        }
        /*!
         * \brief Creates a polygonal facet
         * \param[in] vertices a const reference to a vector that
         *  contains the vertices
         * \return the index of the created facet
         */
        index_t create_facet_polygon( const GEO::vector< index_t >& vertices ) {
            return mesh_.mesh_->facets.create_polygon( vertices ) ;
        }

        /*!
         * \brief Creates a contiguous chunk of triangles
         * \param[in] nb_triangles number of triangles to create
         * \return the index of the first triangle
         */
       index_t create_facet_triangles( index_t nb_triangles )
       {
           return mesh_.mesh_->facets.create_triangles( nb_triangles ) ;

       }
       /*!
        * \brief Creates a contiguous chunk of quads
        * \param[in] nb_quads number of quads to create
        * \return the index of the first quad
        */
       index_t create_facet_quads( index_t nb_quads )
       {
           return mesh_.mesh_->facets.create_quads( nb_quads ) ;
       }
        /*!
         * @brief Sets a vertex of a facet by local vertex index.
         * @param[in] facet_id index of the facet, in 0..nb()-1.
         * @param[in] local_vertex_id index of the vertex in the facet. Local index between 0 and nb_vertices(cell_id)-1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of facet \param of the facet facet_id. Index between 0 and nb()-1.
         */
        void set_facet_vertex(
            index_t facet_id,
            index_t local_vertex_id,
            index_t vertex_id )
        {
            mesh_.mesh_->facets.set_vertex( facet_id, local_vertex_id, vertex_id ) ;
        }
        /*!
         * @brief Sets a vertex of a facet by local vertex index.
         * @param[in] corner_id index of the corner.
         * @param[in] global_vertex_id.
         */
        void set_facet_corner( index_t corner_id, index_t vertex_id )
        {
            mesh_.mesh_->facet_corners.set_vertex( corner_id, vertex_id ) ;
        }
        /*!
         * @brief Sets an adjacent facet by corner index.
         * @param[in] corner_id the corner index starting edge.
         * @param[in] specifies the facet incident to f along edge le or GEO::NO_FACET if \p edge_id is on the border.
         */
        void set_facet_corners_adjacent( index_t corner_id, index_t specifies )
        {
            mesh_.mesh_->facet_corners.set_adjacent_facet( corner_id, specifies ) ;
        }
        /*!
         * @brief Sets an adjacent facet by facet and local edge index.
         * @param[in] facet_id the facet index
         * @param[in] edge_id the local index of an edge in facet \p facet_id
         * @param[in] specifies the facet incident to f along edge le or GEO::NO_FACET if \p edge_id is on the border.
         */
        void set_facet_adjacent(
            index_t facet_id,
            index_t edge_id,
            index_t specifies )
        {
            mesh_.mesh_->facets.set_adjacent( facet_id, edge_id, specifies ) ;
        }
        /*
         * \brief Copies a triangle mesh into this Mesh.
         * \details Facet adjacence are not computed.
         *   Facet and corner attributes are zeroed.
         * \param[in] triangles facet to vertex links
         * \param[in] steal_args if set, vertices and triangles
         * are 'stolen' from the arguments
         * (using vector::swap).
         */
        void assign_facet_triangle_mesh(
            const std::vector< index_t >& triangles,
            bool steal_args )
        {
            GEO::vector< index_t > copy ;
            copy_std_vector_to_geo_vector( triangles, copy ) ;
            mesh_.mesh_->facets.assign_triangle_mesh( copy, steal_args ) ;
        }
        /*!
         * @brief Removes all the facets and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_facets( bool keep_attributes, bool keep_memory ) {
            mesh_.mesh_->facets.clear( keep_attributes, keep_memory ) ;
        }
        void connect_facets() {
            mesh_.mesh_->facets.connect() ;
        }
        void permute_facets( GEO::vector<index_t>& permutation ) {
            mesh_.mesh_->facets.permute_elements( permutation ) ;
        }

        /*!@}
         * \section Cells methods
         * @{
         */
        /*!
         * @brief Creates a contiguous chunk of cells of the same type.
         * @param[in] nb_cells  number of cells to create
         * @param[in] type   type of the cells to create, one of GEO::MESH_TET, GEO::MESH_HEX,
         * GEO::MESH_PRISM, GEO::MESH_PYRAMID, GEO::MESH_CONNECTOR.
         * @return the first created cell.
         */
        index_t create_cells( index_t nb_cells, GEO::MeshCellType type )
        {
            return mesh_.mesh_->cells.create_cells( nb_cells, type ) ;
        }
        /*
         * \brief Copies a tets mesh into this Mesh.
         * \details Cells adjacence are not computed.
         *   cell and corner attributes are zeroed.
         * \param[in] tets cells to vertex links
         * \param[in] steal_args if set, vertices and tets
         * are 'stolen' from the arguments
         * (using vector::swap).
         */
        void assign_cell_tet_mesh(
            const std::vector< index_t >& tets,
            bool steal_args )
        {
            GEO::vector< index_t > copy ;
            copy_std_vector_to_geo_vector( tets, copy ) ;
            mesh_.mesh_->cells.assign_tet_mesh( copy, steal_args ) ;
        }

        /*!
         * @brief Sets a vertex of a cell by local vertex index.
         * @param[in] cell_id index of the cell, in 0..nb()-1.
         * @param[in] local_vertex_id index of the vertex in the cell. Local index between 0 and nb_vertices(cell_id)-1.
         * @param[in] vertex_id specifies the vertex \param local_vertex_id of cell \param of the cell cell_id. Index between 0 and nb()-1.
         */
        void set_cell_vertex(
            index_t cell_id,
            index_t local_vertex_id,
            index_t vertex_id )
        {
            mesh_.mesh_->cells.set_vertex( cell_id, local_vertex_id, vertex_id ) ;
        }
        /*!
         * \brief Sets the vertex that a corner is incident to
         * \param[in] corner_index the corner, in 0..nb()-1
         * \param[in] vertex_index specifies the vertex that corner \p c is incident to
         */
        void set_cell_corner_vertex_index( index_t corner_index, index_t vertex_index ) {
            mesh_.mesh_->cell_corners.set_vertex( corner_index, vertex_index ) ;
        }
        /*!
         * @brief Retrieve the adjacencies
         */
        void cells_connect()
        {
            mesh_.mesh_->cells.connect() ;
        }
        /*!
         * @brief Applies a permutation to the entities and their attributes.
         * On exit, permutation is modified (used for internal bookkeeping).
         * Applying a permutation permutation is equivalent to:
         * <code>
         *  for(i=0; i<permutation.size(); i++) {
         *      data2[i] = data[permutation[i]]
         *       }
         *  data = data2 ;
         *  </code>
         */
        void cells_permute_elements( GEO::vector< index_t >& permutation )
        {
            mesh_.mesh_->cells.permute_elements( permutation ) ;
        }
        /*!
         * @brief Removes all the cells and attributes.
         * @param[in] keep_attributes if true, then all the existing attribute
         * names / bindings are kept (but they are cleared). If false, they are destroyed.
         * @param[in] keep_memory if true, then memory is kept and can be reused
         * by subsequent mesh entity creations.
         */
        void clear_cells( bool keep_attributes, bool keep_memory ) {
            mesh_.mesh_->cells.clear( keep_attributes, keep_memory ) ;
        }
        void permute_cells( GEO::vector<index_t>& permutation ) {
            mesh_.mesh_->cells.permute_elements( permutation ) ;
        }
        /*!
         * @}
         */

    private:
        Mesh& mesh_ ;
    } ;
}

#endif
