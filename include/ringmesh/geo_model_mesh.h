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
         * @brief Number of vertices stored.
         * @details Calls initialize() if vertices are not filled yet
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
         * Test if the mesh vertices need to be initialized,
         * if so initialize them.
         */
        void test_and_initialize() const ;
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


    } ;

    class RINGMESH_API GeoModelMeshCells {
        friend class GeoModelMesh ;
    public:
        GeoModelMeshCells( GeoModelMesh& gmm, GEO::Mesh& mesh ) ;
        /*!
         * Test if the mesh cells are initialized
         */
        bool is_initialized() const ;

    private:
        /*!
         * Test if the mesh cells need to be initialized,
         * if so initialize them.
         */
        void test_and_initialize() const ;
        /*!
         * @brief Initialize the  cells from the cells
         *        of the GeoModel Region cells
         * @details Fills the mesh_.cells
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
         * Vector storing the index of the starting cell index
         * for a given region and a given cell type.
         * For example:
         *    the 2nd hex index of the surface index S will be found here:
         *    surface_facet_ptr_[ALL*S + HEX] + 2
         */
        std::vector< index_t > region_cell_ptr_ ;

    } ;

    class RINGMESH_API GeoModelMeshOrder {


    } ;

    class RINGMESH_API GeoModelMesh {
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
        GeoModelMeshCells cells ;
        GeoModelMeshOrder order ;

    } ;



}

#endif
