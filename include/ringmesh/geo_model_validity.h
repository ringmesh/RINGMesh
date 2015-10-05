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

/*! \author Jeanne Pellerin */

#ifndef __RINGMESH_GEO_MODEL_VALIDITY__
#define __RINGMESH_GEO_MODEL_VALIDITY__

#include <ringmesh/common.h>

#include <geogram/basic/file_system.h>

namespace RINGMesh {
    class GeoModel ;
    class GeoModelElement ;
}


namespace RINGMesh {
    
    
    /*!
    * @brief Check model validity
    * @details In debug mode problematic vertices, edges, elements are
    *          saved in the debug_directory_
    *
    * @param check_surface_intersections Optional expensive check of the
    *        intersections between the model surfaces
    *
    * @todo Check the consistency of index info for vertices -
    * gme_vertices model_vertex_id
    */
    bool is_geomodel_valid( 
        const GeoModel& GM, 
        bool check_surface_intersections = true 
    ) ;
    
    /*!
    * @brief Check the validity of all individual elements
    * @details Check that the elements belong to this model,
    *          call the check validity for each element
    *          For regions, check that their boundary is a one connected component
    *          manifold closed surface.
    *
    */
    bool are_geomodel_elements_valid( const GeoModel& GM ) ;


    /*!
    * @brief Check geological validity
    * @details Only a fault can have a free border and
    *          an stratigraphical interface can be on the boundary of maximum two layers
    *          See Building and Editing a Sealed Geological Model,
    *          Caumon et al. 2004
    */
    bool is_geomodel_geology_valid( const GeoModel& GM ) ;



    /// Do we need it to be in a class ? [JP]
    static std::string validity_errors_directory =
        GEO::FileSystem::get_current_working_directory() ;

    /*!
    * @brief Set the directory where debugging information shall be stored
    * @details Test that this directory exists, if not
    *          keep the previous value.
    *          The default directory is the executable directory.
    */
    static void set_debug_directory( const std::string& directory )
    {
        if( GEO::FileSystem::is_directory( directory ) ) {
            validity_errors_directory = directory ;
        }         
    }




}




#endif