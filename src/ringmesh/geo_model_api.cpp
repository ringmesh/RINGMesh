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

#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geometry.h>

namespace RINGMesh {

    /*!
    * @brief Translates the boundary model by a vector.
    *
    * Every single mesh of the boundary model is translated:
    * corners, lines and surfaces.
    *
    * @param[in] translation_vector vector of translation.
    *
    * @todo Review: Add documentation - Replace the return value by a gme_t [AB]
    */
    void translate(
        GeoModel& M,
        const vec3& translation_vector )
    {
        // Note: if the translation is null, do nothing.
        if( translation_vector == vec3( 0, 0, 0 ) ) {
            return ;
        }

        for( index_t v = 0; v < M.mesh.vertices.nb(); ++v ) {
            vec3 p = M.mesh.vertices.vertex( v ) ;
            for( index_t i = 0; i < 3; i++ ) {
                p[ i ] += translation_vector[ i ] ;
            }
            M.mesh.vertices.update_point( v, p ) ;
        }
    }

    /*!
    * \brief Rotate the boundary model.
    *
    * Applies a rotation about the line defined by the point
    * \p origin and the vector \p axis. The rotation angle is
    * \p theta. If \p degrees is true the angle is in degrees,
    * else in radians. All the vertices of the boundary model
    * undergo the rotation (each mesh inside the boundary model:
    * corners, lines and surfaces).
    *
    * @param origin point in which passes the rotation axis.
    *
    * @param axis vector which defines the rotation axis.
    *
    * @param theta rotation angle (in radians or degrees).
    *
    * @param degrees true is \p theta is in degrees, false
    * if in radians.
    */
    void rotate(
        GeoModel& M,
        const vec3& origin,
        const vec3& axis,
        float64 theta,
        bool degrees )
    {
        // Note: Rotation is impossible about an axis with null length.
        ringmesh_debug_assert( axis != vec3() ) ;
        if( theta == 0. ) {
            return ;
        }

        GEO::Matrix< float64, 4 > rot_mat ;
        rotation_matrix_about_arbitrary_axis(
            origin, axis, theta, degrees, rot_mat ) ;

        for( index_t v = 0; v < M.mesh.vertices.nb(); ++v ) {
            const vec3& p = M.mesh.vertices.vertex( v ) ;

            float64 old[ 4 ] = { p[ 0 ], p[ 1 ], p[ 2 ], 1. } ;
            float64 new_p[ 4 ] = { 0, 0, 0, 1. } ;
            GEO::mult( rot_mat, old, new_p ) ;
            ringmesh_debug_assert( new_p[ 3 ] == 1. ) ;

            M.mesh.vertices.update_point( v, vec3( new_p[ 0 ], new_p[ 1 ], new_p[ 2 ] ) ) ;
        }
    }


}