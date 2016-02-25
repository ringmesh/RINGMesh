/*
*  Copyright (c) 2012-2014, Bruno Levy
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*
*  * Redistributions of source code must retain the above copyright notice,
*  this list of conditions and the following disclaimer.
*  * Redistributions in binary form must reproduce the above copyright notice,
*  this list of conditions and the following disclaimer in the documentation
*  and/or other materials provided with the distribution.
*  * Neither the name of the ALICE Project-Team nor the names of its
*  contributors may be used to endorse or promote products derived from this
*  software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
*  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
*  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
*  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
*  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
*  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
*  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
*  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
*  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
*  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: Bruno Levy
*
*     Bruno.Levy@inria.fr
*     http://www.loria.fr/~levy
*
*     ALICE Project
*     LORIA, INRIA Lorraine,
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX
*     FRANCE
*
*
*  MODIFIED BY 
* 
* Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
* All rights reserved.
*
*     http://www.ring-team.org
*
*     RING Project
*     Ecole Nationale Superieure de Geologie - GeoRessources
*     2 Rue du Doyen Marcel Roubault - TSA 70605
*     54518 VANDOEUVRE-LES-NANCY
*     FRANCE
*/
#ifndef __RINGMESH_GEOGRAM_MESH_REPAIR__
#define __RINGMESH_GEOGRAM_MESH_REPAIR__

#include <ringmesh/common.h>

#include <vector>

/* 
 * @file High level repair operations on GEO::Mesh modified from geogram/mesh/mesh_repair.cpp
 * @author Jeanne Pellerin
 */

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {
    
    /*!
     * @brief Connects the facets of a Mesh except along edges flagged so.
     * @details Resistent to non-manifold edges AND resistent to badly oriented facets
     */
    void RINGMESH_API connect_mesh_facets( GEO::Mesh& mesh,
        const std::vector<index_t>& border_edges ) ;
    
    /*!
    * @brief Connects the facets of a Mesh except along the mesh edges.
    * @details Resistent to non-manifold edges AND resistent to badly oriented facets
    */
    void RINGMESH_API connect_mesh_facets_except_on_mesh_edges( GEO::Mesh& mesh ) ;

    void RINGMESH_API connect_mesh_facets_except_on_mesh_edges( GEO::Mesh& mesh,
        std::vector<index_t>& non_manifold_edges ) ;

    /*!     
     * @brief Connects the facets of a Mesh with non-manifold edges 
     * @details Resistent to badly oriented facets     
     * @note For a Mesh with non-manifold edges (not the normal state of a Mesh) 
     *       mesh.facets.connect() fails big time without any warning.     
     */
    void RINGMESH_API connect_facets( GEO::Mesh& mesh ) ;

    /*!
    * @brief Returns true if there are colocated vertices in the Mesh
    * @details This is a wrapper around Geogram colocate functions.
    */
    bool RINGMESH_API has_mesh_colocate_vertices( const GEO::Mesh& M, double tolerance ) ;

    /*!
    * @brief Merges the vertices of a mesh that are at the same geometric location
    * @note Copied and modified from geogram/mes/mesh_repair.cpp. 
    *       No choice since BL will not give access to it.
    */
    void RINGMESH_API repair_colocate_vertices( GEO::Mesh& M, double tolerance ) ;
}

#endif