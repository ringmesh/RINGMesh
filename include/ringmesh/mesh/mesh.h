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
     * class base class for encapsulating Mesh structure
     * @brief encapsulate adimensional mesh functionalities in order to provide an API
     * on which we base the RINGMesh algorithm
     * @note For now, we encapsulate the GEO::Mesh class.
     */
    class RINGMESH_API MeshBase {
    ringmesh_disable_copy( MeshBase ) ;
    public:

        virtual void save_mesh(
            const std::string& filename,
            const GEO::MeshIOFlags& ioflags ) const = 0 ;

        /*!
         * @brief return the ColocaterANN at the given ColocaterANN::MeshLocation
         * @warning the ColocaterANN is destroy when calling the Mesh::facets_aabb() and Mesh::cells_aabb()
         */
        virtual const ColocaterANN& colocater_ann(
            ColocaterANN::MeshLocation location ) const = 0 ;

        /*!
         * get access to GEO::MESH... only for GFX..
         * @todo Remove this function as soon as the GEO::MeshGFX is encapsulated
         */
        virtual const GEO::Mesh& gfx_mesh() const = 0 ;

        virtual index_t nb_connected_components() const = 0 ;

        //TODO maybe reimplement the function with a RINGMesh::Mesh??
        virtual void print_mesh_bounded_attributes() const = 0 ;
        /*!
         * \name Vertex methods
         * @{
         */
        /*!
         * @brief Gets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @return const reference to the point that corresponds to the vertex.
         */
        virtual const vec3& vertex( index_t v_id ) const = 0 ;

        /*
         * @brief Gets the number of point in the Mesh.
         */
        virtual index_t nb_vertices() const = 0 ;

        virtual GEO::AttributesManager& vertex_attribute_manager() const = 0 ;
        /*!
         * @}
         */
    protected:
        /*!
         * @brief MeshBase constructor.
         * @param[in] geo_model GeoModel from which the MeshBase depend on.
         * @param[in] dimension dimension of the vertices.
         * @param[in] single_precision if true, vertices are stored in single precision (float),
         * else they are stored as double precision (double)..
         */
        MeshBase( const GeoModel& geo_model )
            : geo_model_( geo_model )
        {
        }
        virtual ~MeshBase()
        {
        }

    protected:
        const GeoModel& geo_model_ ;

    } ;

    /*!
     * class for encapsulating isoated vertices structure
     */
    class RINGMESH_API Mesh0D: public virtual MeshBase {
    ringmesh_disable_copy( Mesh0D ) ;

    protected:
        /*!
         * @brief Mesh0D constructor.
         * @param[in] geo_model GeoModel from which the MeshBase depend on.
         * @param[in] dimension dimension of the vertices.
         * @param[in] single_precision if true, vertices are stored in single precision (float),
         * else they are stored as double precision (double)..
         */
        Mesh0D( const GeoModel& geo_model )
            : MeshBase( geo_model )
        {
        }
        virtual ~Mesh0D()
        {
        }

    } ;

    /*!
     * class for encapsulating 1D mesh component
     */
    class RINGMESH_API Mesh1D: public virtual MeshBase {
    ringmesh_disable_copy( Mesh1D ) ;
    public:
        /*
         * @brief Gets the index of an edge vertex.
         * @param[in] edge_id index of the edge.
         * @param[in] vertex_id local index of the vertex, in {0,1}
         * @return the global index of vertex \param vertex_id in edge \param edge_id.
         */
        virtual index_t edge_vertex( index_t edge_id, index_t vertex_id ) const = 0 ;

        /*!
         * @brief Gets the number of all the edges in the whole Mesh.
         */
        virtual index_t nb_edges() const =0 ;

        /*!
         * @brief Gets the length of the edge \param edge_id
         */
        virtual double edge_length( index_t edge_id ) const = 0 ;

        virtual GEO::AttributesManager& edge_attribute_manager() const = 0 ;

    protected:
        Mesh1D( const GeoModel& geo_model )
            : MeshBase( geo_model )
        {
        }
        virtual ~Mesh1D()
        {
        }
    } ;

    /*!
     * class for encapsulating 2D mesh component
     */
    class RINGMESH_API Mesh2D: public virtual MeshBase {
    ringmesh_disable_copy( Mesh2D ) ;
    public:
        /*!
         * @brief Gets the vertex index by facet index and local vertex index.
         * @param[in] facet_id the facet index.
         * @param[in] vertex_id the local edge index in \param facet_id.
         * @return the global facet index adjacent to the \param edge_id of the facet \param facet_id.
         * @precondition  \param edge_id < number of edge of the facet \param facet_id .
         */
        virtual index_t facet_vertex( index_t facet_id, index_t vertex_id ) const=0 ;
        /*!
         * @brief Gets the number of all facets in the whole Mesh.
         */
        virtual index_t nb_facets() const=0 ;
        /*!
         * @brief Gets the number of vertices in the facet \param facet_id.
         * @param[in] facet_id facet index
         */
        virtual index_t nb_facet_vertices( index_t facet_id ) const=0 ;
        /*!
         * @brief Get the first vertex index of a facet.
         * @param[in] facet_id facet index
         */
        virtual index_t facet_begin( index_t facet_id ) const=0 ;
        /*!
         * @brief Get the last vertex index of a facet.
         * @param[in] facet_id facet index
         */
        virtual index_t facet_end( index_t facet_id ) const=0 ;
        /*!
         * @brief Gets the next vertex index in the facet \param facet_id.
         * @param[in] facet_id facet index
         * @param[in] vertex_id current index
         */
        virtual index_t next_facet_vertex(
            index_t facet_id,
            index_t vertex_id ) const=0 ;
        /*!
         * @brief Gets the previous vertex index in the facet \param facet_id.
         * @param[in] facet_id facet index
         * @param[in] vertex_id current index
         */
        virtual index_t prev_facet_vertex(
            index_t facet_id,
            index_t vertex_id ) const=0 ;
        /*!
         * @brief Gets an adjacent facet index by facet index and local edge index.
         * @param[in] facet_id the facet index.
         * @param[in] edge_id the local edge index in \param facet_id.
         * @return the global facet index adjacent to the \param edge_id of the facet \param facet_id.
         * @precondition  \param edge_id < number of edge of the facet \param facet_id .
         */
        virtual index_t facet_adjacent( index_t facet_id, index_t edge_id ) const=0 ;
        virtual GEO::AttributesManager& facet_attribute_manager() const=0 ;
        /*!
         * @brief Tests whether all the facets are triangles. when all the facets are triangles, storage and access is optimized.
         * @return True if all facets are triangles and False otherwise.
         */
        virtual bool facets_are_simplicies() const=0 ;
        /*!
         * return true if the facet \param facet_id is a triangle
         */
        virtual bool is_triangle( index_t facet_id ) const=0 ;
        /*!
         * @brief Create an AABB tree for a Mesh facets
         * @pre The GeoModelEntity must be simplicial
         * @warning SIDE EFFECTS: The mesh vertices are reordered.
         * @warning calling this function will destroy the ColocaterANN.
         */
        virtual const GEO::MeshFacetsAABB& facets_aabb() const = 0 ;
        /*!
         * Computes the Mesh facet normal
         * @param[in] facet_id the facet index
         * @return the facet normal
         */
        virtual vec3 facet_normal( index_t facet_id ) const=0 ;
        /*!
         * Computes the Mesh facet barycenter
         * @param[in] facet_id the facet index
         * @return the facet center
         */
        virtual vec3 facet_barycenter( index_t facet_id ) const=0 ;
        /*!
         * Computes the Mesh facet area
         * @param[in] facet_id the facet index
         * @return the facet area
         */
        virtual double facet_area( index_t facet_id ) const=0 ;
    protected:
        Mesh2D( const GeoModel& geo_model )
            : MeshBase( geo_model )
        {
        }
        virtual ~Mesh2D()
        {
        }

    } ;

    /*!
     * class for encapsulating 3D mesh component
     */
    class RINGMESH_API Mesh3D: public virtual MeshBase {
    ringmesh_disable_copy( Mesh3D ) ;
    public:
        /*!
         * @brief Gets a vertex index by cell and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        virtual index_t cell_vertex( index_t cell_id, index_t vertex_id ) const = 0 ;

        /*!
         * @brief Gets a vertex index by cell and local edge and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] edge_id the local edge index in \param cell_id.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        virtual index_t cell_edge_vertex(
            index_t cell_id,
            index_t edge_id,
            index_t vertex_id ) const = 0 ;

        /*!
         * @brief Gets a vertex by cell facet and local vertex index.
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @param[in] vertex_id index of the vertex in the facet \param facet_id
         * @return the global vertex index.
         * @precondition vertex_id < number of vertices in the facet \param facet_id
         * and facet_id number of facet in th cell \param cell_id
         */
        virtual index_t cell_facet_vertex(
            index_t cell_id,
            index_t facet_id,
            index_t vertex_id ) const = 0 ;

        /*!
         * @brief Gets a facet index by cell and local facet index.
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @return the global facet index.
         */
        virtual index_t cell_facet( index_t cell_id, index_t facet_id ) const = 0 ;

        /*!
         * @brief Gets the number of facet in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        virtual index_t nb_cell_facets( index_t cell_id ) const = 0 ;

        /*!
         * @brief Gets the number of edges in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        virtual index_t nb_cell_edges( index_t cell_id ) const = 0 ;

        /*!
         * @brief Gets the number of vertices of a facet in a cell
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @return the number of vertices in the facet \param facet_id in the cell \param cell_id
         */
        virtual index_t nb_cell_facet_vertices(
            index_t cell_id,
            index_t facet_id ) const = 0 ;

        /*!
         * @brief Gets the number of vertices of a cell
         * @param[in] cell_id index of the cell
         * @return the number of vertices in the cell \param cell_id
         */
        virtual index_t nb_cell_vertices( index_t cell_id ) const = 0 ;

        /*!
         * @brief Gets the number of cells in the Mesh.
         */
        virtual index_t nb_cells() const = 0 ;

        /*!
         * @return the index of the adjacent cell of \param cell_id along the facet \param facet_id
         */
        virtual index_t cell_adjacent(
            index_t cell_id,
            index_t facet_id ) const = 0 ;

        virtual GEO::AttributesManager& cell_attribute_manager() const = 0 ;

        virtual GEO::AttributesManager& cell_facet_attribute_manager() const = 0 ;

        /*!
         * @brief Gets the type of a cell.
         * @param[in] cell_id the cell index, in 0..nb()-1
         */
        virtual GEO::MeshCellType cell_type( index_t cell_id ) const = 0 ;

        /*!
         * @brief Tests whether all the cells are tetrahedra. when all the cells are tetrahedra, storage and access is optimized.
         * @return True if all cells are tetrahedra and False otherwise.
         */
        virtual bool cells_are_simplicies() const = 0 ;

        /*!
         * @brief Create an AABB tree for a Mesh cells
         * @pre The GeoModelEntity must be simplicial
         * @warning SIDE EFFECTS: The mesh vertices are reordered.
         * @warning calling this function will destroy the ColocaterANN.
         */
        virtual const GEO::MeshCellsAABB& cells_aabb() const = 0 ;
        /*!
         * Computes the Mesh cell facet barycenter
         * @param[in] cell_id the cell index
         * @param[in] facet_id the facet index in the cell
         * @return the cell facet center
         */
        virtual vec3 cell_facet_barycenter(
            index_t cell_id,
            index_t facet_id ) const = 0 ;

        /*!
         * Compute the non weighted barycenter of the \param cell_id
         */
        virtual vec3 cell_barycenter( index_t cell_id ) const = 0 ;

        /*!
         * Computes the Mesh cell facet normal
         * @param[in] cell_id the cell index
         * @param[in] facet_id the facet index in the cell
         * @return the cell facet normal
         */
        virtual vec3 cell_facet_normal(
            index_t cell_id,
            index_t facet_id ) const = 0 ;

        /*!
         * @brief compute the volume of the cell \param cell_id.
         */
        virtual double cell_volume( index_t cell_id ) const = 0 ;

        virtual index_t cell_begin( index_t cell_id ) const = 0 ;

        virtual index_t cell_end( index_t cell_id ) const = 0 ;

        virtual index_t find_cell_corner(
            index_t cell_id,
            index_t vertex_id ) const = 0 ;

    protected:
        Mesh3D( const GeoModel& geo_model )
            : MeshBase( geo_model )
        {
        }
        virtual ~Mesh3D()
        {
        }

    } ;

    class RINGMESH_API MeshAllD: public Mesh0D,
        public Mesh1D,
        public Mesh2D,
        public Mesh3D {
    ringmesh_disable_copy( MeshAllD ) ;
    protected:
        MeshAllD( const GeoModel& geo_model )
            :
                MeshBase( geo_model ),
                Mesh0D( geo_model ),
                Mesh1D( geo_model ),
                Mesh2D( geo_model ),
                Mesh3D( geo_model )
        {
        }
        virtual ~MeshAllD()
        {
        }
    } ;

    /*!
     * @brief class to encapsulate mesh structure in order to provide an API
     * on which we base the RINGMesh algorithm 
     * @note For now, we encapsulate the GEO::Mesh class. We can develop the concept
     * using a factory to build several encapsulating classes. 
     */
    class RINGMESH_API GeogramMesh: public MeshAllD {
    ringmesh_disable_copy( GeogramMesh ) ;
        friend class GeogramMeshBuilder ;

    public:
        /*!
         * @brief Mesh constructor.
         * @param[in] dimension dimension of the vertices.
         * @param[in] single_precision if true, vertices are stored in single precision (float),
         * else they are stored as double precision (double)..
         */
        GeogramMesh(
            const GeoModel& geo_model,
            index_t dimension,
            bool single_precision )
            :
                MeshBase( geo_model ),
                MeshAllD( geo_model ),
                facets_aabb_( nil ),
                cells_aabb_( nil )
        {
            mesh_ = new GEO::Mesh( dimension, single_precision ) ;
            for( index_t i = 0; i < ColocaterANN::NB_LOCATION; i++ ) {
                ann_[i] = nil ;
            }
        }
        ~GeogramMesh()
        {
            if( facets_aabb_ ) delete facets_aabb_ ;
            if( cells_aabb_ ) delete cells_aabb_ ;
            for( index_t i = 0; i < ColocaterANN::NB_LOCATION; i++ ) {
                if( ann_[i] ) delete ann_[i] ;
            }
            delete mesh_ ;
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
        index_t cell_facet( index_t cell_id, index_t facet_id ) const
        {
            return mesh_->cells.facet( cell_id, facet_id ) ;
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

        index_t cell_begin( index_t cell_id ) const
        {
            return mesh_->cells.corners_begin( cell_id ) ;
        }
        index_t cell_end( index_t cell_id ) const
        {
            return mesh_->cells.corners_end( cell_id ) ;
        }
        index_t find_cell_corner( index_t cell_id, index_t vertex_id ) const
        {
            for( index_t v = 0; v < nb_cell_vertices( cell_id ); ++v ) {
                if( cell_vertex( cell_id, v ) == vertex_id ) {
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

        mutable GEO::MeshFacetsAABB* facets_aabb_ ;
        mutable GEO::MeshCellsAABB* cells_aabb_ ;
        mutable ColocaterANN* ann_[ColocaterANN::NB_LOCATION] ;

    } ;
}

#endif
