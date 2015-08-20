/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_MACRO_MESH__
#define __RINGMESH_MACRO_MESH__

#include <ringmesh/common.h>

#include <geogram/mesh/mesh.h>

namespace GEO {
    class MeshCellsAABB ;
    class MeshFacetsAABB ;
}

namespace RINGMesh {
    class BoundaryModel ;
    class MacroMesh ;
    class WellGroup ;
}

namespace RINGMesh {


    /*!
     * Several modes for vertex duplication algorithm:
     *  - NONE = no duplication
     *  - FAULT = duplication along faults
     *  - HORIZON = duplication along horizons
     *  - ALL = duplication along faults and horizons
     */
    enum DuplicateMode {
        NONE, FAULT, HORIZON, ALL
    } ;
    const static index_t NB_FACET_TYPES = 2 ;
    const static index_t NB_CELL_TYPES = 4 ;

    /*!
     * Optional storage of the MacroMesh vertices
     */
    class RINGMESH_API MacroMeshVertices {
        friend class MacroMesh ;
    public:
        MacroMeshVertices( const MacroMesh& mm )
            : mm_( mm )
        {
        }

        index_t nb_vertices() const ;
        index_t vertex_id( index_t mesh, index_t v ) const ;
        const vec3& vertex( index_t global_v ) const ;
        const vec3& vertex( index_t mesh, index_t v ) const ;
        const vec3& duplicated_vertex( index_t v ) const ;

        bool vertex_id(
            index_t mesh,
            index_t cell_corner,
            index_t& vertex_id,
            index_t& duplicated_vertex_id = dummy_index_t ) const ;

        index_t nb_duplicated_vertices() const ;
        index_t nb_total_vertices() const ;

        bool is_surface_to_duplicate( index_t s ) const ;

        void clear() ;

    private:
        /// enum to characterize the action to do concerning a surface
        enum SurfaceAction {
            SKIP = -2,          /// do nothing
            TO_PROCESS = -1,    /// need to be duplicated (don't know which side yet)
            NEG_SIDE = 0,       /// need to duplicate the side opposite to the facet normal
            POS_SIDE = 1        /// need to duplicate the side following the facet normal
        } ;
        /// Action to do according a surface index
        typedef std::pair< index_t, SurfaceAction > surface_side ;

        void initialize() ;
        void initialize_duplication() ;

        bool duplicate_corner(
            const std::set< surface_side >& surfaces,
            std::vector< SurfaceAction >& info ) ;

    private:
        /// Attached MaroMesh
        const MacroMesh& mm_ ;

        /// Vector of the vertices with different coordinates without duplication
        std::vector< vec3 > vertices_ ;
        /*!
         * @brief Mapping of each mesh vertex id into the global macromesh
         * vertex storage id
         * @details Let Vmi denote the vertex index value of ith vertex
         * of the mth mesh in the macromesh. The vector is stored like this :
         * [V11, V12, V13 ... , Vi1, Vi2, Vi3, ...]
         */
        std::vector< index_t > global_vertex_indices_ ;
        /*!
         * Vector of size mm_.nb_meshes(). Each value is the
         * number of all the number of vertices of the previous meshes.
         * This is equivalent to this code:
         * vertex2mesh_[m] = 0 ;
         * for( index_t i = 0; i < m; i++ ) {
         *     vertex2mesh_[m] += mm_.mesh(i).vertices.nb() ;
         * }
         */
        std::vector< index_t > vertex2mesh_ ;

        /*!
         * @brief Vector of the cell corners
         * @details A corner denote the point index of a cell vertex.
         * There is one corner per vertex per cell, so for a mesh of 3 tetrahedra
         * there are 3 * 4 = 12 corners. This vector stores continuously all the
         * cell corners of all the meshes.
         */
        std::vector< index_t > cell_corners_ ;
        /*!
         * Vector storing the index of where to start reading the cell_corners_
         * vector for a given mesh. Vector size is  mm_.nb_meshes() + 1.
         * Ex. The first cell corner of the mesh m is cell_corners_[mesh_cell_corner_ptr_[m]],
         * the second cell_corners_[mesh_cell_corner_ptr_[m]+1], ...
         */
        std::vector< index_t > mesh_cell_corner_ptr_ ;
        /*!
         * @brief Vector of duplicated vertices
         * @details Each value is a duplicated vertex, the index corresponds to
         * vertex index in vertices_.
         */
        std::vector< index_t > duplicated_vertex_indices_ ;
    } ;

    /*!
     * Optional storage of the MacroMesh facets
     */
    class RINGMESH_API MacroMeshFacets {
    public:
        MacroMeshFacets( const MacroMesh& mm )
            : mm_( mm ), nb_facets_( 0 ), nb_triangle_( 0 ), nb_quad_( 0 )
        {
        }

        index_t facet( index_t s, index_t f ) const ;
        index_t mesh( index_t s ) const ;
        index_t nb_facets() const ;
        index_t nb_facets( index_t s ) const ;

        index_t nb_triangle() const ;
        index_t nb_triangle( index_t s ) const ;
        index_t triangle_id( index_t s, index_t t ) const ;

        index_t nb_quad() const ;
        index_t nb_quad( index_t s ) const ;
        index_t quad_id( index_t s, index_t q ) const ;

        void clear() ;

    private:
        /*!
         * Tests if the MacroMeshFacets needs to be initialized and initialize it
         */
        void test_initialize() const {
            if( surface_facets_.empty() ) {
                const_cast< MacroMeshFacets* >( this )->initialize() ;
            }
        }
        /*!
         * Id where to start reading the vector surface_facets_ for a given surface
         * @param[in] s id of the surface
         * @return the corresponding id
         */
        index_t surface_begin( index_t s ) const
        {
            return surface_facet_ptr_[NB_FACET_TYPES * s] ;
        }
        /*!
         * Id where to stop reading the vector surface_facets_ for a given surface
         * @param[in] s id of the surface
         * @return the corresponding id
         */
        index_t surface_end( index_t s ) const
        {
            return surface_facet_ptr_[NB_FACET_TYPES * ( s + 1 )] ;
        }
        /*!
         * Accessor for the surface_facets_ vector
         * @param[in] global_f the id to read
         * @return the value in the vector
         */
        index_t facet( index_t global_f ) const
        {
            return surface_facets_[global_f] ;
        }

        void initialize() ;

    private:
        /// Attached MaroMesh
        const MacroMesh& mm_ ;

        /*!
         * @brief  Vector of the facet indices
         * @details This vector stores the facet indices sorted by mesh and by type.
         * Let Tsi denote the ith triangle index of the sth surface and
         * Qsi  denote the ith quad index of the sth surface. The vector storage is:
         * [T11, T12, .... , Q11, Q12 ... , T21, T22, ... , Q21, Q22 ... ]
         */
        std::vector< index_t > surface_facets_ ;
        /*!
         * Vector storing the index of where to start reading the surface_facets_
         * vector for a given surface and a given facet type.
         * For example:
         *    the 2nd quad index of the surface index S will be found here:
         *    surface_facets_[surface_facet_ptr_[NB_FACET_TYPES*S + 1] + 2]
         */
        std::vector< index_t > surface_facet_ptr_ ;
        /*!
         * @brief Mapping between a surface id and a mesh id of the MacroMesh
         * @details Since the same surfaces can be found twice, one in each adjacent
         * region. This vector stores the mesh indices to use corresponding to each
         * surface index.
         */
        std::vector< index_t > surface2mesh_ ;

        /// Number of facets in the MacroMesh
        index_t nb_facets_ ;
        /// Number of triangles in the MacroMesh
        index_t nb_triangle_ ;
        /// Number of quad in the MacroMesh
        index_t nb_quad_ ;
    } ;

    /*!
     * Optional storage of the MacroMesh cells
     */
    class RINGMESH_API MacroMeshCells {
    public:
        MacroMeshCells( const MacroMesh& mm )
            :
                mm_( mm ),
                nb_cells_( 0 ),
                nb_tet_( 0 ),
                nb_pyramid_( 0 ),
                nb_prism_( 0 ),
                nb_hex_( 0 )
        {
        }

        index_t cell_adjacent( index_t mesh, index_t c, index_t f ) const ;

        index_t cell_index_in_mesh( index_t global_index, index_t& mesh_id ) const ;

        index_t nb_cells() const ;
        index_t nb_cells( index_t r ) const ;

        index_t nb_tet() const ;
        index_t nb_tet( index_t r ) const ;
        index_t tet_id( index_t r, index_t t ) const ;

        index_t nb_pyramid() const ;
        index_t nb_pyramid( index_t r ) const ;
        index_t pyramid_id( index_t r, index_t p ) const ;

        index_t nb_prism() const ;
        index_t nb_prism( index_t r ) const ;
        index_t prism_id( index_t r, index_t p ) const ;

        index_t nb_hex() const ;
        index_t nb_hex( index_t r ) const ;
        index_t hex_id( index_t r, index_t h ) const ;

        void clear() ;

    private:
        /*!
         * Tests if the MacroMeshCells needs to be initialized and initialize it
         */
        void test_initialize() const
        {
            if( nb_cells_ == 0 ) {
                const_cast< MacroMeshCells* >( this )->initialize() ;
            }
        }
        void initialize() ;
        /*!
         * Id where to start reading the vector mesh_cell_ptr_ for a given mesh
         * @param[in] mesh id of the mesh
         * @return the corresponding id
         */
        index_t mesh_begin( index_t mesh ) const
        {
            return mesh_cell_ptr_[NB_CELL_TYPES * mesh] ;
        }
        /*!
         * Id where to stop reading the vector mesh_cell_ptr_ for a given mesh
         * @param[in] mesh id of the mesh
         * @return the corresponding id
         */
        index_t mesh_end( index_t mesh ) const
        {
            return mesh_cell_ptr_[NB_CELL_TYPES * mesh] ;
        }

    private:
        /// Attached MaroMesh
        const MacroMesh& mm_ ;

        /// Vector of the cell ids in the corresponding GEO::Mesh
        std::vector< index_t > cells_ ;
        /// Mapping between mesh id and cell elements in cells_
        std::vector< index_t > mesh_cell_ptr_ ;
        /// Vector of the adjacent cell ids in the MacroMesh
        std::vector< index_t > cell_adjacents_ ;
        /// Mapping between mesh id and cell elements in cell_adjacents_
        std::vector< index_t > mesh_cell_adjacent_ptr_ ;

        /// Number of cells in the MacroMesh
        index_t nb_cells_ ;
        /// Number of tetrahedra in the MacroMesh
        index_t nb_tet_ ;
        /// Number of pyramids in the MacroMesh
        index_t nb_pyramid_ ;
        /// Number of prisms in the MacroMesh
        index_t nb_prism_ ;
        /// Number of hexahedra in the MacroMesh
        index_t nb_hex_ ;
    } ;

    /*!
     * Optional storage of the MacroMesh tools
     */
    class RINGMESH_API MacroMeshTools {
    public:
        MacroMeshTools( MacroMesh& mm ) ;
        ~MacroMeshTools() ;

        const GEO::MeshFacetsAABB& facet_aabb( index_t region ) const ;
        const GEO::MeshCellsAABB& cell_aabb( index_t region ) const ;

        void clear() ;

    private:
        void init_facet_aabb( index_t region ) const ;
        void init_cell_aabb( index_t region ) const ;

    private:
        /// Attached MaroMesh
        MacroMesh& mm_ ;

        /// Storage of the AABB trees on the facets
        std::vector< GEO::MeshFacetsAABB* > facet_aabb_ ;
        /// Storage of the AABB trees on the cells
        std::vector< GEO::MeshCellsAABB* > cell_aabb_ ;
    } ;

    class RINGMESH_API MacroMeshOrder {
    public:
        MacroMeshOrder( MacroMesh& mm ) ;
        ~MacroMeshOrder() ;
        const index_t nb_total_vertices() const ;
        const index_t nb_vertices() const ;
        void clear() ;
        const vec3 point( const index_t id ) const ;
        void move_point( const index_t id, const vec3& u ) ;
    private:
        void initialize() ;
        /*!
         * Tests if the MacroMeshOrder needs to be initialized and initialize it
         */
        void test_initialize() const
        {
            if( nb_vertices_ == 0 ) {
                const_cast< MacroMeshOrder* >( this )->initialize() ;
            }
        }
        void test_point_list_initialize() ;
    private:
        /// Attached MacroMesh
        const MacroMesh& mm_ ;
        /// Total number of vertices + new nodes on cell edges
        index_t nb_vertices_ ;
        /// New points
        std::vector< vec3 > points_ ;

    } ;

    class RINGMESH_API MacroMeshEdges {
    public:
        MacroMeshEdges( MacroMesh& mm )
            : mm_( mm )
        {
        }

        index_t nb_wells() const ;
        index_t nb_edges() const ;
        index_t nb_edges( index_t w ) const ;
        index_t vertex_id( index_t w, index_t e, index_t v ) const ;

    private:
        /*!
         * Tests if the MacroMeshCells needs to be initialized and initialize it
         */
        void test_initialize() const
        {
            if( edges_.empty() ) {
                const_cast< MacroMeshEdges* >( this )->initialize() ;
            }
        }
        void initialize() ;
    private:
        /// Attached MaroMesh
        const MacroMesh& mm_ ;

        /// Vector of edge vertex id in the MacroMesh
        std::vector< index_t > edges_ ;
        /// Mapping between well id and edges_
        std::vector< index_t > well_ptr_ ;

    } ;

    static std::vector< std::vector< vec3 > > empty_vertices ;
    class RINGMESH_API MacroMesh {
    ringmesh_disable_copy( MacroMesh ) ;
    public:
        const static index_t ALL_REGIONS = index_t( -1 ) ;

        MacroMesh( const BoundaryModel& model ) ;
        MacroMesh() ;
        virtual ~MacroMesh() ;
        void copy( const MacroMesh& mm, bool copy_attributes = true ) ;

        void compute_tetmesh(
            const std::string& method,
            index_t region_id = ALL_REGIONS,
            bool add_steiner_points = true,
            std::vector< std::vector< vec3 > >& internal_vertices = empty_vertices ) ;

        /*!
         * Access to a GEO::Mesh of the MacroMesh
         * @param[in] region id of mesh/region
         * @return a reference to the GEO::Mesh
         */
        GEO::Mesh& mesh( index_t region ) const
        {
            return *meshes_[region] ;
        }

        /*!
         * Get the number of meshes/regions
         * @return the corresponding number
         */
        index_t nb_meshes() const
        {
            return meshes_.size() ;
        }

        void set_wells( const WellGroup* wells ) ;
        const WellGroup* wells() const
        {
            return wells_ ;
        }

        /*!
         * Access to the BoundaryModel attached to the MacroMesh
         * @return a const reference to the corresponding BoundaryModel
         */
        const BoundaryModel& model() const
        {
            ringmesh_debug_assert( model_ ) ;
            return *model_ ;
        }
        void set_model( const BoundaryModel& model ) ;

        /*!
         * Access the DuplicateMode
         * @return the current DuplicateMode
         */
        DuplicateMode duplicate_mode() const
        {
            return mode_ ;
        }
        /*!
         * Set a new DuplicateMode
         * @param[in] mode the new DuplicateMode for the MacroMesh
         */
        void set_duplicate_mode( const DuplicateMode& mode ) const
        {
            MacroMesh* not_const = const_cast< MacroMesh* >( this ) ;
            not_const->mode_ = mode ;
            not_const->vertices.cell_corners_.clear() ;
            not_const->vertices.duplicated_vertex_indices_.clear() ;
            not_const->vertices.mesh_cell_corner_ptr_.clear() ;
        }

        /*
         * Change the order of the mesh
         */
        void set_order( const index_t o )
        {
            if( o != order_ ) {
                order.clear() ;
            }
            order_ = o ;
        }

        /*
         * Gets the mesh elements order
         * @return the const order
         */
        const index_t get_order() const
        {
            return order_ ;
        }

        void translate( const vec3& translation_vector ) ;
        void rotate( const vec3& origin, const vec3& axis, float64 angle, bool degrees = false ) ;

    protected:
        /// BoundaryModel representing the structural information of the mesh
        const BoundaryModel* model_ ;
        /// Vector of meshes, one by region
        std::vector< GEO::Mesh* > meshes_ ;
        /// Optional duplication mode to compute the duplication of vertices on surfaces
        DuplicateMode mode_ ;

        /// Optional WellGroup associated with the model
        const WellGroup* wells_ ;

        /// Order of the mesh
        index_t order_ ;

    public:
        /// Optional storage of the MacroMesh vertices
        MacroMeshVertices vertices ;
        /// Optional storage of the MacroMesh edges
        MacroMeshEdges edges ;
        /// Optional storage of the MacroMesh facets
        MacroMeshFacets facets ;
        /// Optional storage of the MacroMesh cells
        MacroMeshCells cells ;
        /// Optional storage of tools
        MacroMeshTools tools ;
        /// Optional storage for high orders mesh
        MacroMeshOrder order ;

    } ;
}

#endif
