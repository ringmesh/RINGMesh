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
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */


#ifndef __RINGMESH_GEO_MODEL_MESH__
#define __RINGMESH_GEO_MODEL_MESH__

#include <ringmesh/common.h>
#include <ringmesh/geo_model_element.h>

#include <geogram/mesh/mesh.h>

namespace RINGMesh {
    class GeoModel ;
    class GeoModelMesh ;
}

namespace RINGMesh {


    class RINGMESH_API GeoModelMeshVertices {
        ringmesh_disable_copy( GeoModelMeshVertices ) ;
        friend class GeoModelMesh ;
    public:
        /*!
         * @brief Identification of a vertex in a GeoModelElement
         * @todo Je changerai bien ce nom moche au refactoring [JP]
         */
        struct VertexInGME {
            VertexInGME( GME::gme_t t, index_t vertex_id_in )
                : gme_id( t ), v_id( vertex_id_in )
            {
            }
            VertexInGME()
                : gme_id(), v_id( NO_ID )
            {
            }
            bool operator<( const VertexInGME& rhs ) const
            {
                if( gme_id != rhs.gme_id ) {
                    return gme_id < rhs.gme_id ;
                } else {
                    return v_id < rhs.v_id ;
                }
            }
            bool operator==( const VertexInGME& rhs ) const
            {
                return gme_id == rhs.gme_id && v_id == rhs.v_id ;
            }
            bool is_defined() const
            {
                return gme_id.is_defined() && v_id != NO_ID ;
            }
            /// Unique identifier of the associated GeoModelElement
            GME::gme_t gme_id ;
            /// Index of the vertex in the BME
            index_t v_id ;
        } ;

    public:
        GeoModelMeshVertices( GeoModelMesh& gmm, GeoModel& gm, GEO::Mesh& mesh ) ;
        ~GeoModelMeshVertices() ;

        /*!
         * Test if the mesh vertices are initialized
         */
        bool is_initialized() const ;
        /*!
         * Test if the mesh vertices need to be initialized,
         * if so initialize them.
         */
        void test_and_initialize() const ;

        /*!
         * @brief Number of vertices stored.
         */
        index_t nb() const ;

        /*!
         * @brief Coordinates of a vertex of the GeoModel
         * @pre v < nb()
         */
        const vec3& vertex( index_t v ) const ;

        /*!
         * @brief Returns the index of the given vertex in the model
         * @param[in] p input point coordinates
         * @return index of the vertex in the model if found
         * (distance < epsilon), otherwise NO_ID
         */
        index_t index( const vec3& p ) const ;

        /*!
         * @brief Get the vertices in GME corresponding to the given unique vertex
         */
        const std::vector< VertexInGME >& gme_vertices( index_t v ) const ;

        /*!
         * @brief To use when building the model by first adding its vertices
         * @warning The client is responsible for setting the mapping between the points
         * of the BME and the unique vertex
         */
        index_t add_vertex( const vec3& point ) ;

        /*!
         * @brief Add a vertex in a GeoModelElement
         *        corresponding to an existing vertex of the model
         */
        void add_to_bme( index_t v, const VertexInGME& v_gme ) ;

        /*!
         * @brief Change one of the GME vertex associated to a vertex
         * @param v Index of the vertex
         * @param i Index of the GME vertex
         * @param v_gme Id of GME and of the vertex in that GME
         */
        void set_gme( index_t v, index_t i, const VertexInGME& v_gme ) ;

        /*!
         * @brief Set the point coordinates of all the vertices that
         *        share this unique vertex, including the unique vertex itself.
         * @param[in] v Index of the vertex
         * @param[in] point New coordinates
         */
        void update_point( index_t v, const vec3& point ) ;

        /*!
         * @brief Clear the vertices - clear the gme_vertices_ -
         *        clear global vertex information in the all BMME
         * @warning Not stable - crashes if atributes are still bound
         */
        void clear() ;

    private:
        /*!
         * @brief Initialize the vertices from the vertices
         *        of the GeoModel Corners, Lines, and Surfaces
         * @details Fills the mesh_.vertices, gme_vertices_ and
         *         delete colocated vertices
         */
        void initialize() ;

        /*!
         * @brief Delete the KdTree and set the pointer to nil.
         */
        void clear_kdtree() ;

        /*!
         * Test if the kdtree need to be initialized,
         * if so initialize it.
         */
        void test_kdtree_and_initialize() const ;
        /*!
         * Initialize the kdtree with the mesh vertices
         */
        void initialize_kdtree() ;
        /*!
         * @brief Remove colocated vertices
         */
        void remove_colocated() ;

        /*!
         * @brief Delete vertices for which to_delete[i] != i
         * @detail The global vertices are deleted, gme_vertices_
         * is updated and the model_vertx_id in the GeoModelMeshElement
         * of the BoudnaryModel are updated too.
         *
         * @param[in,out] to_delete can be NO_ID or give the index of a
         *  kept vertex with wich information should be merged.
         *  It is recyled to give the mapping between old and new vertex indices
         * @pre to_delete[ v ] is either NO_ID, or is equal or inferior to v
         */
        void erase_vertices( std::vector< index_t >& to_delete ) ;

        /*!
         * @brief Remove all invalid VertexInGME and delete the vertices
         * that are not anymore in any GeoModelElement
         */
        void erase_invalid_vertices() ;

    private:
        /// Attached GeoModelMesh owning the vertices
        GeoModelMesh& gmm_ ;
        /// Attached GeoModel
        GeoModel& gm_ ;
        /// Attached Mesh
        GEO::Mesh& mesh_ ;


        /*!
         * Vertices in GeoModelElements corresponding to each vertex
         * @todo Change this extremely expensive storage !!!
         */
        std::vector< std::vector< VertexInGME > > gme_vertices_ ;
        /// Kd-tree of the model vertices
        ColocaterANN* kdtree_ ;

    } ;

    class RINGMESH_API GeoModelMeshFacets {
        ringmesh_disable_copy( GeoModelMeshFacets ) ;
        friend class GeoModelMesh ;
    public:
        enum FacetType {
            TRIANGLE, QUAD, POLYGON, ALL, NO_FACET
        };

    public:
        GeoModelMeshFacets( GeoModelMesh& gmm, GEO::Mesh& mesh ) ;
        ~GeoModelMeshFacets() ;

        /*!
         * Test if the mesh facets are initialized
         */
        bool is_initialized() const ;
        void test_and_initialize() const ;

        /*!
         * @brief Number of facets stored.
         */
        index_t nb() const ;

        /*!
         * Get the number of vertices in the facet
         * @param[in] f the facet index
         * @return the number of vertices
         */
        index_t nb_vertices( index_t f ) const ;
        /*!
         * Get the vertex index of a vertex in a facet
         * in the GeoModelMesh
         * @param[in] f the facet index
         * @param[in] v the local vertex index [0, nb_vertices_in_facet[
         * @return the vertex index
         */
        index_t vertex( index_t f, index_t v ) const ;
        /*!
         * Get the adjacent facet index in the GeoModelMesh
         * @param[in] f the facet index
         * @param[in] e the edge index
         * @return the adjacent facet index
         */
        index_t adjacent( index_t f, index_t e ) const ;
        /*!
         * Get the surface index in the GeoModel according the facet
         * index in the GeoModelMesh
         * @param[in] f the facet index
         * @return the surface index
         */
        index_t surface( index_t f ) const ;
        /*!
         * Get the facet index in the GeoModelMesh restricted to
         * the surface owing the facet
         * @param[in] f the facet index
         * @return the facet index varying from 0 to nb_facets
         * in the surface owing \p f
         */
        index_t index_in_surface( index_t f ) const ;
        /*!
         * Get the facet index in the GeoModelMesh restricted to
         * the surface owing the facet and its type
         * @param[in] f the facet index
         * @param[out] index the facet index varying from 0 to nb_facets
         * of the corresponding type of \p f in the owing surface
         * @return the type of the facet \p f
         */
        FacetType type( index_t f, index_t& index = dummy_index_t ) const ;

        /*!
         * Get the number of facets of the corresponding type
         * @param[in] type the corresponding type
         * @return the number of facets
         */
        index_t nb_facets( FacetType type = ALL ) const ;
        /*!
         * Get the number of facets of the corresponding type
         * in the given surface of the GeoModel
         * @param[in] s the surface index
         * @param[in] type the corresponding type
         * @return the number of facets
         */
        index_t nb_facets( index_t s, FacetType type = ALL ) const ;
        /*!
         * Get the facet index in the GeoModelMesh
         * @param[in] s the surface index owing the facet
         * @param[in] f the facet index varying from 0 to nb_facets in the surface
         * @param[in] type it can specify the facet type used. For example, if type = QUAD
         * then \p f represents the fth quad in the surface \p s and \p f can vary from 0
         * to nb_quads( s ).
         * @return the facet index
         */
        index_t facet( index_t s, index_t f, FacetType type = ALL ) const ;

        /*!
         * Get the number of triangles in the GeoModelMesh
         * @return the number of triangles
         */
        index_t nb_triangle() const ;
        /*!
         * Get the number of triangles in the given surface
         * @param[in] s the surface index
         * @return the number of triangles
         */
        index_t nb_triangle( index_t s ) const ;
        /*!
         * Get the facet index in the GeoModelMesh corresponding
         * to the asked triangle in the surface
         * @param[in] s the surface index
         * @param[in] t the tth triangle index varying from 0 to nb_triangles( s )
         * @return the facet index
         */
        index_t triangle( index_t s, index_t t ) const ;

        /*!
         * Get the number of quads in the GeoModelMesh
         * @return the number of quads
         */
        index_t nb_quad() const ;
        /*!
         * Get the number of quads in the given surface
         * @param[in] s the surface index
         * @return the number of quads
         */
        index_t nb_quad( index_t s ) const ;
        /*!
         * Get the facet index in the GeoModelMesh corresponding
         * to the asked quad in the surface
         * @param[in] s the surface index
         * @param[in] q the qth quad index varying from 0 to nb_quads( s )
         * @return the facet index
         */
        index_t quad( index_t s, index_t q ) const ;

        /*!
         * Get the number of polygons in the GeoModelMesh
         * @return the number of polygons
         */
        index_t nb_polygon() const ;
        /*!
         * Get the number of polygons in the given surface
         * @param[in] s the surface index
         * @return the number of polygons
         */
        index_t nb_polygon( index_t s ) const ;
        /*!
         * Get the facet index in the GeoModelMesh corresponding
         * to the asked polygon in the surface
         * @param[in] s the surface index
         * @param[in] p the pth polygon index varying from 0 to nb_polygons( s )
         * @return the facet index
         */
        index_t polygon( index_t s, index_t p ) const ;

        /*!
         * Clear the facets of the GeoModelMesh
         */
        void clear() ;

        /*!
         * Get the center of the given facet
         * @param[in] f the facet index
         */
        vec3 center( index_t f ) const ;
        /*!
         * Get the area of the facet
         * @param[in] f the facet index
         */
        double area( index_t f ) const ;

    private:
        /*!
         * Initialize the facets of the GeoModelMesh
         * and sort them per surface and facet type
         * Example for a mesh with two surfaces and only triangles and quads
         * [TRGL,TRGL, .. , QUAD, QUAD .. , TRGL, TRGL, ... , QUAD, QUAD ..]
         * |          surface 0           |             surface 1           |
         */
        void initialize() ;

        /*!
         * Bind attribute to the facets attribute manager
         */
        void bind_attribute() ;
        /*!
         * Unbind attribute to the facets attribute manager
         */
        void unbind_attribute() ;

    private:
        /// Attached GeoModelMesh owning the vertices
        GeoModelMesh& gmm_ ;
        /// Attached GeoModel
        const GeoModel& gm_ ;
        /// Attached Mesh
        GEO::Mesh& mesh_ ;

        /// Attribute storing the surface index per facet
        GEO::Attribute< index_t > surface_id_ ;
        /*!
         * Vector storing the index of the starting facet index
         * for a given surface and a given facet type.
         * For example:
         *    the 2nd quad index of the surface index S will be found here:
         *    surface_facet_ptr_[ALL*S + QUAD] + 2
         */
        std::vector< index_t > surface_facet_ptr_ ;

        /// Number of triangles in the GeoModelMesh
        index_t nb_triangle_ ;
        /// Number of quads in the GeoModelMesh
        index_t nb_quad_ ;
        /// Number of polygons in the GeoModelMesh
        index_t nb_polygon_ ;

    } ;

    class RINGMESH_API GeoModelMeshEdges {
        ringmesh_disable_copy( GeoModelMeshEdges ) ;
    public:
        GeoModelMeshEdges( GeoModelMesh& gmm, GEO::Mesh& mesh ) ;
        ~GeoModelMeshEdges() ;

        /*!
         * Test if the mesh edges are initialized
         */
        bool is_initialized() const ;
        /*!
         * Tests if the mesh edges needs to be initialized and initialize it
         */
        void test_and_initialize() const ;

        /*!
         * Gets the number of wells
         * @return the corresponding number
         */
        index_t nb_wells() const ;
        /*!
         * Gets the number of edges in the MacroMesh
         * @return the corresponding number
         */
        index_t nb_edges() const ;
        /*!
         * Gets the number of edges of a Well
         * @param[in] w the well index
         * @return the corresponding number
         */
        index_t nb_edges( index_t w ) const ;
        /*!
         * Gets the vertex index of the GeoModelMesh
         * @param[in] w the well index
         * @param[in] e the edge index in the well (from 0 to nb_edges in the well)
         * @param[in] v the vertex index of the edge (0 or 1 )
         * @return the global vertex index
         */
        index_t vertex( index_t w, index_t e, index_t v ) const ;
        /*!
         * Clear the mesh edges
         */
        void clear() ;

    private:
        /*!
         * Initialize the mesh edges
         */
        void initialize() ;

    private:
        /// Attached GeoModelMesh owning the vertices
        GeoModelMesh& gmm_ ;
        /// Attached GeoModel
        const GeoModel& gm_ ;
        /// Attached Mesh
        GEO::Mesh& mesh_ ;

        /*!
         * Vector storing the index of the starting edge index
         * for a given well
         */
        std::vector< index_t > well_ptr_ ;

    } ;

    class RINGMESH_API GeoModelMeshCells {
    ringmesh_disable_copy( GeoModelMeshCells ) ;
        friend class GeoModelMesh ;
    public:
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

    public:
        GeoModelMeshCells( GeoModelMesh& gmm, GEO::Mesh& mesh ) ;
        /*!
         * Test if the mesh cells are initialized
         */
        bool is_initialized() const ;
        /*!
         * Test if the mesh cells are duplicated
         */
        bool is_duplication_initialized() const ;

        /*!
         * Test if the mesh cells need to be initialized,
         * if so initialize them.
         */
        void test_and_initialize() const ;

        /*!
         * @brief Number of cells stored.
         */
        index_t nb() const ;
        /*!
         * Gets the number of duplicated points by the DuplicateMode algorithm
         * @return the corresponding number of duplications
         */
        index_t nb_duplicated_vertices() const ;
        /*!
         * Gets the total number of vertices (mesh.vertices.nb() + nb_duplicated_vertices())
         * @return the corresponding number of vertices
         */
        index_t nb_total_vertices() const ;
        /*!
         * Check if the corner in a cell is duplicated,
         * if so give the duplicated vertex index
         * @param[in] c the cell index in the GeoModelMesh
         * @param[in] v the local vertex index in the cell \p c (0 to nb_vertices( c ))
         * @param[out] duplicate_vertex_index the duplicated vertex index (0 to nb_duplicated_vertices())
         * @return true if the corner is duplicated
         */
        bool is_corner_duplicated(
            index_t c,
            index_t v,
            index_t& duplicate_vertex_index ) const ;
        /*!
         * Get the vertex index in the GeoModelMesh corresponding
         * to the given duplicated vertex index
         * @param[in] duplicate_vertex_index the duplicated vertex index
         * @return the vertex index
         */
        index_t duplicated_vertex( index_t duplicate_vertex_index ) const ;

        /*!
         * Get the number of vertices in the cell
         * @param[in] c the cell index
         * @return the number of vertices
         */
        index_t nb_vertices( index_t c ) const ;
        /*!
         * Get the vertex index of a vertex in a cell
         * in the GeoModelMesh
         * @param[in] c the cell index
         * @param[in] v the local vertex index [0, nb_vertices_in_cell[
         * @return the vertex index
         */
        index_t vertex( index_t c, index_t v ) const ;
        /*!
         * Get the number of facets in the cell
         * @param[in] c the cell index
         */
        index_t nb_facets( index_t c ) const ;
        /*!
         * Get the adjacent cell index in the GeoModelMesh
         * @param[in] c the cell index
         * @param[in] f the edge index
         * @return the adjacent cell index
         */
        index_t adjacent( index_t c, index_t f ) const ;
        /*!
         * Get the region index in the GeoModel according the cell
         * index in the GeoModelMesh
         * @param[in] c the cell index
         * @return the region index
         */
        index_t region( index_t c ) const ;
        /*!
         * Get the cell index in the GeoModelMesh restricted to
         * the region owing the cell
         * @param[in] c the cell index
         * @return the cell index varying from 0 to nb_cells
         * in the region owing \p f
         */
        index_t index_in_region( index_t c ) const ;
        /*!
         * Get the cell index in the GeoModelMesh restricted to
         * the region owing the cell and its type
         * @param[in] c the cell index
         * @param[out] index the cell index varying from 0 to nb_cells
         * of the corresponding type of \p c in the owing region
         * @return the type of the cell \p f
         */
        GEO::MeshCellType type( index_t c, index_t& index = dummy_index_t ) const ;

        /*!
         * Get the number of cells of the corresponding type
         * @param[in] type the corresponding type
         * @return the number of cells
         */
        index_t nb_cells( GEO::MeshCellType type = GEO::MESH_NB_CELL_TYPES ) const ;
        /*!
         * Get the number of cells of the corresponding type
         * in the given region of the GeoModel
         * @param[in] r the region index
         * @param[in] type the corresponding type
         * @return the number of cells
         */
        index_t nb_cells(
            index_t r,
            GEO::MeshCellType type = GEO::MESH_NB_CELL_TYPES ) const ;
        /*!
         * Get the cell index in the GeoModelMesh
         * @param[in] r the region index owing the cell
         * @param[in] c the cell index varying from 0 to nb_cells in the region
         * @param[in] type it can specify the cell type used. For example, if type = QUAD
         * then \p c represents the fth quad in the region \p s and \p c can vary from 0
         * to nb_quads( s ).
         * @return the cell index
         */
        index_t cell( index_t r, index_t c, GEO::MeshCellType type =
            GEO::MESH_NB_CELL_TYPES ) const ;

        /*!
         * Get the number of tets in the GeoModelMesh
         * @return the number of tets
         */
        index_t nb_tet() const ;
        /*!
         * Get the number of tets in the given region
         * @param[in] r the region index
         * @return the number of tets
         */
        index_t nb_tet( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked tet in the region
         * @param[in] r the region index
         * @param[in] t the tth tet index varying from 0 to nb_tet( r )
         * @return the cell index
         */
        index_t tet( index_t r, index_t t ) const ;

        /*!
         * Get the number of hexs in the GeoModelMesh
         * @return the number of hexs
         */
        index_t nb_hex() const ;
        /*!
         * Get the number of hexs in the given region
         * @param[in] r the region index
         * @return the number of hexs
         */
        index_t nb_hex( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked hex in the region
         * @param[in] r the region index
         * @param[in] h the hth hex index varying from 0 to nb_hex( r )
         * @return the cell index
         */
        index_t hex( index_t r, index_t h ) const ;

        /*!
         * Get the number of prisms in the GeoModelMesh
         * @return the number of prisms
         */
        index_t nb_prism() const ;
        /*!
         * Get the number of prisms in the given region
         * @param[in] r the region index
         * @return the number of prisms
         */
        index_t nb_prism( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked prism in the region
         * @param[in] r the region index
         * @param[in] p the pth prism index varying from 0 to nb_prism( r )
         * @return the cell index
         */
        index_t prism( index_t r, index_t p ) const ;


        /*!
         * Get the number of pyramids in the GeoModelMesh
         * @return the number of pyramids
         */
        index_t nb_pyramid() const ;
        /*!
         * Get the number of pyramids in the given region
         * @param[in] r the region index
         * @return the number of pyramids
         */
        index_t nb_pyramid( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked pyramid in the region
         * @param[in] r the region index
         * @param[in] p the pth pyramid index varying from 0 to nb_pyramid( r )
         * @return the cell index
         */
        index_t pyramid( index_t r, index_t p ) const ;

        /*!
         * Get the number of connectors in the GeoModelMesh
         * @return the number of connectors
         */
        index_t nb_connector() const ;
        /*!
         * Get the number of connectors in the given region
         * @param[in] r the region index
         * @return the number of connectors
         */
        index_t nb_connector( index_t r ) const ;
        /*!
         * Get the cell index in the GeoModelMesh corresponding
         * to the asked connector in the region
         * @param[in] r the region index
         * @param[in] c the cth connector index varying from 0 to nb_connector( r )
         * @return the cell index
         */
        index_t connector( index_t r, index_t c ) const ;

        /*!
         * Clear the mesh cells
         */
        void clear() ;
        /*!
         * Remove the duplication of the mesh cell facets
         */
        void clear_duplication() ;

        /*!
         * Determine if a cell facet is on a surface, if so fill the \p action
         * with the surface id and the surface side encountered
         * @param[in] c the cell index
         * @param[in] f the facet index
         * @param[out] facet the facet index colocalised with the cell facet
         * @param[out] side the side of the facet \p facet.
         * true = side of the facet normal, false = the other side
         * @return true is the cell facet is on a surface
         */
        bool is_cell_facet_on_surface(
            index_t c,
            index_t f,
            index_t& facet = dummy_index_t,
            bool& side = dummy_bool ) const ;

        /*!
         * Get the center of the given cell
         * @param[in] c the cell index
         */
        vec3 center( index_t c ) const ;
        /*!
         * Get the volume of the cell
         * @param[in] c the cell index
         */
        double volume( index_t c ) const ;

    private:
        /*!
         * @brief Initialize the  cells from the cells
         *        of the GeoModel Region cells
         * @details Fills the mesh_.cells
         */
        void initialize() ;

        /*!
         * Bind attribute to the cells attribute manager
         */
        void bind_attribute() ;
        /*!
         * Unbind attribute to the cells attribute manager
         */
        void unbind_attribute() ;

        /// enum to characterize the action to do concerning a surface
        /// Action concerns the vertices of a Surface and not the Surface
        enum ActionOnSurface {
            /// do nothing
            SKIP = -2,
            /// need to be duplicated (don't know which side yet)
            TO_PROCESS = -1,
            /// need to duplicate the side opposite to the facet normal
            NEG_SIDE = 0,
            /// need to duplicate the side following the facet normal
            POS_SIDE = 1
        } ;

        /// Action to do according a surface index
        typedef std::pair< index_t, ActionOnSurface > action_on_surface ;

        /*!
         * Test if the mesh cell are duplicated according
         * the duplication mode, if not duplicate them.
         */
        void test_and_initialize_duplication() const ;
        /*!
         * Duplicate the mesh cell along some surfaces defined
         * by the duplication mode
         */
        void initialize_duplication() ;
        /*!
         * Test if we need to duplicate mesh cell along the given
         * surface according the duplicate mode
         * @param[in] s the surface index in the GeoModel
         */
        bool is_surface_to_duplicate( index_t s ) const ;

        /*!
         * Determine the actions to do according the action_on_surfaces
         * encountered during the propagation around a vertex (initialize())
         * @param[in] surfaces the action_on_surfaces encountered
         * @param[in,out] info the global information on what to do for each surface.
         * This information is updated in this function according the encountered action_on_surfaces
         * @return true if the corners should be duplicated
         */
        bool are_corners_to_duplicate(
            const std::vector< action_on_surface >& surfaces,
            std::vector< ActionOnSurface >& info ) ;
        /*!
         * Test if the mesh cell facet attribute is filled with
         * the colocalised facet. If not fill it.
         */
        void test_and_initialize_cell_facet() const ;
        /*!
         * Initialize the mesh cell facet attribute of colocalised facet.
         */
        void initialize_cell_facet() ;

    private:
        /// Attached GeoModelMesh owning the vertices
        GeoModelMesh& gmm_ ;
        /// Attached GeoModel
        const GeoModel& gm_ ;
        /// Attached Mesh
        GEO::Mesh& mesh_ ;

        /// Attribute storing the region index per cell
        GEO::Attribute< index_t > region_id_ ;
        /*!
         * Vector storing the index of the starting cell index
         * for a given region and a given cell type.
         * For example:
         *    the 2nd hex index of the region index R will be found here:
         *    surface_facet_ptr_[GEO::MESH_NB_CELL_TYPES*R + HEX] + 2
         */
        std::vector< index_t > region_cell_ptr_ ;

        /// Number of tet in the GeoModelMesh
        index_t nb_tet_ ;
        /// Number of hex in the GeoModelMesh
        index_t nb_hex_ ;
        /// Number of prism in the GeoModelMesh
        index_t nb_prism_ ;
        /// Number of pyramid in the GeoModelMesh
        index_t nb_pyramid_ ;
        /// Number of connector in the GeoModelMesh
        index_t nb_connector_ ;

        /// Current duplicate mode applied on the mesh
        DuplicateMode mode_ ;
        /*!
         * @brief Vector of duplicated vertices
         * @details Each value is a duplicated vertex, the index corresponds to
         * vertex index in mesh.vertices.
         */
        std::vector< index_t > duplicated_vertex_indices_ ;

        /*!
         * @brief Attribute storing the colocalised facet index per cell facet
         * @detail If a cell facet is on a surface, the attribute is equal to
         * the index of the corresponding facet.
         */
        GEO::Attribute< index_t > facet_id_ ;
    } ;

    class RINGMESH_API GeoModelMeshOrder {
        ringmesh_disable_copy( GeoModelMeshOrder ) ;


    } ;

    class RINGMESH_API GeoModelMesh {
        ringmesh_disable_copy( GeoModelMesh ) ;
    public:
        GeoModelMesh( GeoModel& gm ) ;
        ~GeoModelMesh() ;

        const GeoModel& model() const
        {
            return gm_ ;
        }

        /*!
         * Copy the current GeoModelMesh into a Mesh
         * @param[out] mesh The mesh to fill
         */
        void copy_mesh( GEO::Mesh& mesh ) const
        {
            mesh.copy( *mesh_ ) ;
        }

        GEO::AttributesManager& vertex_attribute_manager() const
        {
            return mesh_->vertices.attributes() ;
        }
        GEO::AttributesManager& facet_attribute_manager() const
        {
            return mesh_->facets.attributes() ;
        }
        GEO::AttributesManager& cell_attribute_manager() const
        {
            return mesh_->cells.attributes() ;
        }

        /*!
         * Access the DuplicateMode
         * @return the current DuplicateMode
         */
        GeoModelMeshCells::DuplicateMode duplicate_mode() const
        {
            return mode_ ;
        }
        /*!
         * Set a new DuplicateMode
         * @param[in] mode the new DuplicateMode for the GeoModelMesh
         */
        void set_duplicate_mode( const GeoModelMeshCells::DuplicateMode& mode ) const
        {
            if( mode_ == mode ) return ;
            mode_ = mode ;
            const_cast< GeoModelMesh* >( this )->cells.clear_duplication() ;
        }

        /*!
         * @brief Remove colocated vertices
         */
        void remove_colocated_vertices() ;

        /*!
         * @brief Delete vertices for which to_delete[i] != i
         * @detail The global vertices are deleted, gme_vertices_
         * is updated and the model_vertx_id in the GeoModelMeshElement
         * of the BoudnaryModel are updated too.
         *
         * @param[in,out] to_delete can be NO_ID or give the index of a
         *  kept vertex with wich information should be merged.
         *  It is recyled to give the mapping between old and new vertex indices
         * @pre to_delete[ v ] is either NO_ID, or is equal or inferior to v
         */
        void erase_vertices( std::vector< index_t >& to_delete ) ;

        /*!
         * @brief Remove all invalid VertexInGME and delete the vertices
         * that are not anymore in any GeoModelElement
         */
        void erase_invalid_vertices() ;

    private:
        /// Attached GeoMode
        GeoModel& gm_ ;
        /*!
         * @brief Mesh storing the vertices, edges, facets and cells
         * that are not colocated/duplicated
         * @details Each mesh element is unique.
         * On these elements, attributes can be defined
         */
        GEO::Mesh* mesh_ ;
        /// Optional duplication mode to compute the duplication of cells on surfaces
        mutable GeoModelMeshCells::DuplicateMode mode_ ;

    public:
        GeoModelMeshVertices vertices ;
        GeoModelMeshEdges edges ;
        GeoModelMeshFacets facets ;
        GeoModelMeshCells cells ;
//        GeoModelMeshOrder order ;

    } ;

}

#endif
