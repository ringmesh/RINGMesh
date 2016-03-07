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

#ifndef __RINGMESH_GEO_MODEL_VALIDITY__
#define __RINGMESH_GEO_MODEL_VALIDITY__

#include <ringmesh/common.h>

#include <geogram/basic/file_system.h>

/*!
* @file ringmesh/geo_model_validity.h
* @brief Functions to check the validity of GeoModels
* @author Jeanne Pellerin
*/

namespace RINGMesh {
    class GeoModel ;
    class GeoModelElement ;
}

namespace RINGMesh {       
    /*! Set the default directory to store invalid elements of 
     *   models to be the current working directory
     */
    static std::string validity_errors_directory =
        GEO::FileSystem::get_current_working_directory() ;

    /*!
     * @brief Set the directory where debugging information on 
     * invalid elements shall be stored
     * @details If directory does not exist keep the previous value.
     */
    void set_validity_errors_directory( const std::string& directory ) ;

    /*!
    * @brief Check global model validity
    * @details In debug mode problematic vertices, edges, elements are
    *          saved in the validity_errors_directory
    * @param GM GeoModel to check
    * @param check_surface_intersections Optional expensive check of the
    *        intersections between the model surfaces
    * @todo Check the consistency of gme_vertices vs. model_vertex_id
    */
    bool RINGMESH_API is_geomodel_valid(
        const GeoModel& GM, 
        bool check_surface_intersections = true 
    ) ;
    
    /*!
    * @brief Check the validity of all individual elements
    * @details Check that the elements belong to this model,
    *          call the check validity for each element
    *          For regions, check that their boundary is a one connected component
    *          manifold closed surface.
    */
    bool RINGMESH_API are_geomodel_elements_valid( const GeoModel& GM ) ;

    /*!
    * @brief Check geological validity of a GeoModel
    * @details Only a fault can have a free border and
    *          an stratigraphical interface can be on the boundary of maximum two layers
    *          See Building and Editing a Sealed Geological Model, Caumon et al. 2004
    */
    bool RINGMESH_API is_geomodel_geology_valid( const GeoModel& GM ) ;
  
}




#endif
