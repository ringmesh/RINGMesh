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

#ifndef __RINGMESH_GEO_MODEL_REPAIR__
#define __RINGMESH_GEO_MODEL_REPAIR__

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/geo_model_builder.h>

/*!
 * @file ringmesh/geo_model_repair.h
 * @brief Functions to repair GeoModel.
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    /*!
     * @brief Try repairing a supposedly invalid GeoModel
     * @details Remove colocated vertices in all GeoModelMeshEntity.
     *          Remove degenerated edges and facets in Surfaces and Lines.
     * @warning The Mesh of the model is deleted.
     *          This function will by no mean fix all errors in a GeoModel
     *          It has been tested on a very small number of models.
     *
     * @todo Convenience design to change. This allows an easy access for repair
     * to the internal meshes of the GeoModel. This class is otherwise artificial.
     */
    class RINGMESH_API GeoModelRepair: public GeoModelBuilder {
    public:
        GeoModelRepair( GeoModel& model )
            : GeoModelBuilder( model )
        {
        }
        virtual ~GeoModelRepair()
        {
        }
        /*!
         * @brief Detect and remove degenerated edges in a \param line.
         * @return the number of degenerated edges that have been removed from the line.
         */
        void geo_model_mesh_repair() ;
        /*!
         * @brief For all the lines in the geomodel, switch line boundaries
         * if the way of their indices do not follow the way of the vertex indices.
         */
        void repair_line_boundary_vertex_order() ;

    private:
        index_t repair_line_mesh( Line& line ) ;
        void mesh_detect_degenerate_edges(
            const Mesh1D& M,
            GEO::vector< index_t >& e_is_degenerate,
            GEO::vector< index_t >& colocated_vertices ) ;
        void mesh_detect_degenerate_facets(
            const Mesh2D& M,
            GEO::vector< index_t >& f_is_degenerate,
            GEO::vector< index_t >& colocated_vertices ) ;
        bool facet_is_degenerate(
            const Mesh2D& M,
            index_t f,
            GEO::vector< index_t >& colocated_vertices ) ;

        index_t detect_degenerate_facets( Mesh2D& M ) ;

        void remove_degenerate_facet_and_edges( std::set< gme_t >& to_remove ) ;

        void remove_colocated_element_vertices( std::set< gme_t >& to_remove ) ;

        void remove_colocated_entity_vertices( std::set< gme_t >& to_remove ) ;
        void vertices_on_inside_boundary(
            const gme_t& E_id,
            std::set< index_t >& vertices ) ;

        bool edge_is_degenerate(
            const Mesh1D& M,
            index_t e,
            GEO::vector< index_t >& colocated_vertices )
        {
            index_t v1 = colocated_vertices[M.edge_vertex( e, 0 )] ;
            index_t v2 = colocated_vertices[M.edge_vertex( e, 1 )] ;
            return v1 == v2 ;
        }
    } ;

} //namespace RINGMesh

#endif 
