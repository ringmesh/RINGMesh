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

#include <ringmesh/basic/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_compare.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_AABB.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geogram_extension/geogram_extension.h>

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

		/*
		 * \brief Compares the current mesh with an other.
		 *  Compares the current mesh with mesh \p M1 according to the comparison
		 *  flags \p flags (see #MeshCompareFlags). The function returns a
		 *  comparison status similar to the comparison \p flags:
		 *  - if a comparison failed for a test f in in \p flags, the same
		 *  flag f is set in the status
		 *  - otherwise the flag f is cleared in the status.
		 *  A status of zero indicates that all comparisons succeeded.
		 * \param[in] M1 the first input mesh
		 * \param[in] M2 the second input mesh
		 * \param[in] flags specifies which properties of the meshes should be
		 *  compared. By default it is set to #MESH_COMPARE_SURFACE_PROPS
		 *  information for the two meshes
		 * \param[in] tolerance relative tolerance used to compare floating point
		 * values (such as the mash areas)
		 * \param[in] verbose enables/disables the display of mesh information
		 *  for the two meshes, as well mesh comparison error messages.
		 * \retval #MESH_COMPARE_OK if meshes \p M1 and \p M2 are identical
		 *  according to the comparison criteria
		 * \retval the comparison status otherwise.
		 * \see MeshCompareFlags
		 * \see Geom::mesh_area()
		 * \see meshes_have_same_topology()
		 */
		 bool mesh_compare(
			const Mesh& M1, 
			GEO::MeshCompareFlags flags = GEO::MESH_COMPARE_SURFACE_PROPS,
			double tolerance = 0.0,
			bool verbose = true)
		 {
			 return GEO::mesh_compare(*mesh_, *M1.mesh_, flags, tolerance, verbose);
		 }

        void save_mesh(
            const std::string& filename,
            const GEO::MeshIOFlags& ioflags ) const
        {
            GEO::mesh_save( *mesh_, filename, ioflags ) ;
        }
        /*!
         * @brief return the ColocaterANN at the given ColocaterANN::MeshLocation
         * @warning the ColocaterANN is destroy when calling the Mesh::facets_aabb() and Mesh::cells_aabb()
         */
        const ColocaterANN& colocater_ann(
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
        void print_mesh_bounded_attributes() const
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
        GEO::AttributesManager& cell_facet_attribute_manager() const
        {
            return mesh_->cell_facets.attributes() ;
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
            return RINGMesh::mesh_cell_barycenter( *mesh_, cell_id ) ;
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
}

#endif
