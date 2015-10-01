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
        GeoModelMeshVertices( GeoModelMesh& gmm, GEO::Mesh& mesh ) ;
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
        const GeoModel& gm_ ;
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
        index_t facet_in_surface( index_t f ) const ;
        /*!
         * Get the facet index in the GeoModelMesh restricted to
         * the surface owing the facet and its type
         * @param[in] f the facet index
         * @param[out] index the facet index varying from 0 to nb_facets
         * of the corresponding type of \p f in the owing surface
         * @return the type of the facet \p f
         */
        FacetType facet_type( index_t f, index_t& index ) const ;

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

    class RINGMESH_API GeoModelMeshCells {
        ringmesh_disable_copy( GeoModelMeshCells ) ;


    } ;

    class RINGMESH_API GeoModelMeshOrder {
        ringmesh_disable_copy( GeoModelMeshOrder ) ;


    } ;

    class RINGMESH_API GeoModelMesh {
        ringmesh_disable_copy( GeoModelMesh ) ;
    public:
        GeoModelMesh( const GeoModel& gm ) ;
        ~GeoModelMesh() ;

        const GeoModel& model() const {
            return gm_ ;
        }

        GEO::AttributesManager& vertex_attribute_manager() const {
            return mesh_->vertices.attributes() ;
        }
        GEO::AttributesManager& facet_attribute_manager() const {
            return mesh_->facets.attributes() ;
        }
        GEO::AttributesManager& cell_attribute_manager() const {
            return mesh_->cells.attributes() ;
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
        const GeoModel& gm_ ;
        /*!
         * @brief Mesh storing the vertices, facets and cells that
         * are not colocated/duplicated
         * @details Each mesh element is unique.
         * On these elements, attributes can be defined
         */
        GEO::Mesh* mesh_ ;

    public:
        GeoModelMeshVertices vertices ;
        GeoModelMeshFacets facets ;
//        GeoModelMeshCells cells ;
//        GeoModelMeshOrder order ;

    } ;



}

#endif
