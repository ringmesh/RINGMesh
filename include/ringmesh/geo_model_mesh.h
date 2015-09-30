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

        bool is_initialized() const ;
        void test_and_initialize() const ;

    private:
        /*!
         * @brief Initialize the vertices from the vertices
         *        of the GeoModel Corners, Lines, and Surfaces
         * @details Fills the mesh_.vertices, bme_vertices_ and
         *         delete colocated vertices
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
         * Vertices in GeoModelElements corresponding to each vertex
         * @todo Change this extremely expensive storage !!!
         */
        std::vector< std::vector< VertexInGME > > gme_vertices_ ;

    } ;

    class RINGMESH_API GeoModelMeshFacets {


    } ;

    class RINGMESH_API GeoModelMeshCells {


    } ;

    class RINGMESH_API GeoModelMeshOrder {


    } ;

    class RINGMESH_API GeoModelMesh {
    public:
        GeoModelMesh( const GeoModel& gm ) ;
        ~GeoModelMesh() ;

        const GEO::Mesh& mesh() const {
            return *mesh_ ;
        }
        const GeoModel& model() const {
            return gm_ ;
        }

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
