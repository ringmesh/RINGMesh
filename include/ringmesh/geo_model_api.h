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

#ifndef __RINGMESH_GEO_MODEL_API__
#define __RINGMESH_GEO_MODEL_API__

#include <ringmesh/common.h>

namespace RINGMesh {
    class GeoModel ;
    class GeoModelElement ;
}

namespace RINGMesh {

    /*!
     * @brief Compute the size (volume, area, length) of an Element
     *
     * @param[in] E Element to evaluate
     */
    double model_element_size( const GeoModelElement& E ) ;
    /*!
     * Compute the size (volume, area, length) of an Element cell (cell, facet, edge)
     * @param[in] E Element to evaluate
     * @param[in] c the cell index
     */
    double model_element_cell_size( const GeoModelElement& E, index_t c ) ;

    /*!
     * @brief Compute the center of a GeoModelElement
     *
     * @param[in] E Element to evaluate
     * @return The coordinates of the center
     */
    vec3 model_element_center( const GeoModelElement& E ) ;
    /*!
     * @brief Compute the center of a GeoModelElement cell (cell, facet, edge)
     *
     * @param[in] E Element to evaluate
     * @param[in] c the cell index
     * @return The coordinates of the center
     */
    vec3 model_element_cell_center( const GeoModelElement& E, index_t c  ) ;

    /*!
    * @brief Translates the boundary model by a vector.
    *
    * Every single mesh of the boundary model is translated:
    * corners, lines and surfaces.
    *
    * @param[in] M GeoModel on which compute the translation
    * @param[in] translation_vector vector of translation.
    */
    void translate( GeoModel& M, const vec3& ) ;

    /*!
    * \brief Rotate the boundary model.
    *
    * Applies a rotation about the line defined by the point
    * \p origin and the vector \p axis. The rotation angle is
    * \p angle. If \p degrees is true the angle is in degrees,
    * else in radians. All the vertices of the boundary model
    * undergo the rotation (each mesh inside the boundary model:
    * corners, lines and surfaces).
    *
    * @param[in] M GeoModel on which compute the rotation
    * @param[in] origin point in which passes the rotation axis.
    * @param[in] axis vector which defines the rotation axis.
    * @param[in] angle rotation angle (in radians or degrees).
    * @param[in] degrees true is \p theta is in degrees, false
    * if in radians.
    */
    void rotate(
        GeoModel& M, 
        const vec3& origin,
        const vec3& axis,
        float64 angle,
        bool degrees = false ) ;

    const static index_t ALL_REGIONS = index_t( -1 ) ;
    static std::vector< std::vector< vec3 > > empty_vertices ;
    /*!
     * Compute the tetrahedral mesh of the input structural model
     * @param[in] M GeoModel to tetrahedralize
     * @param[in] method Mesher used
     * @param[in] region_id Region to mesh, ALL_REGIONS for all
     * @param[in] add_steiner_points if true, the mesher will add some points inside the region
     * @param[in] internal_vertices points inside the domain to constrain during the
     * mesh generation. There is one vector per mesh.
     */
    void tetrahedralize(
        GeoModel& M,
        const std::string& method = "TetGen",
        index_t region_id = ALL_REGIONS,
        bool add_steiner_points = true,
        std::vector< std::vector< vec3 > >& internal_vertices = empty_vertices ) ;

}



#endif 
