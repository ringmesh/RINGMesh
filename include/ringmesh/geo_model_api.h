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
*     Ecole Nationale Superieure de Geologie - GeoRessources
*     2 Rue du Doyen Marcel Roubault - TSA 70605
*     54518 VANDOEUVRE-LES-NANCY
*     FRANCE
*/

#ifndef __RINGMESH_GEO_MODEL_API__
#define __RINGMESH_GEO_MODEL_API__

#include <ringmesh/common.h>

#include <geogram/basic/attributes.h>
#include <geogram/mesh/mesh.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_element.h>
#include <ringmesh/geogram_extension.h>


/*!
* @file ringmesh/geo_model_api.h
* @brief High level functions on GeoModel
* @author Jeanne Pellerin and Arnaud Botella
* @todo Encapsulate these functions in a namespace and TEST them.
*/

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {
    class GeoModel ;
    class GeoModelElement ;
    class GeoModelMeshElement ;
}

namespace RINGMesh {

    /*!
    * @brief Print in the console the model statistics
    * @details Output number of facets, vertices, and of the different element types.
    * @todo Implement a comparison of GeoModels for testing.
    */
    void RINGMESH_API print_model( const GeoModel& model ) ;

    
    bool RINGMESH_API are_geomodel_surface_meshes_simplicial( const GeoModel& geomodel );

    bool RINGMESH_API are_geomodel_region_meshes_simplicial( const GeoModel& geomodel );


    /*!
    * @brief Build a Mesh from the model non-duplicated vertices and its Surface facets.
    * @details Adjacencies are not set. Client should call mesh repair functions afterwards.
    * @todo Implement Implement it for meshed Regions.
    * Implement it for a set of GME.
    */
    void RINGMESH_API build_mesh_from_geomodel( const GeoModel& model, GEO::Mesh& M ) ;


    /*! 
     * @brief Bind named GEO::Attribute on the GeoModel element vertices
     * @warning It is up to the client to unbind the attribute    
     * @pre Elements of geomodel_element_type are GeoModelMeshElement
     */
    template< class T >
    void create_attributes_on_geomodel_element_facets(
        const GeoModel& geomodel,
        GeoModelElement::TYPE geomodel_element_type,
        const std::string& attribute_name,
        AttributeVector<T>& attributes )
    {
        index_t nb_elements = geomodel.nb_elements( geomodel_element_type ) ;
        attributes.resize( nb_elements ) ;
        for( index_t i = 0; i < nb_elements; ++i ) {
            const GeoModelMeshElement& E = geomodel.mesh_element( GME::gme_t( geomodel_element_type, i ) );
            GEO::AttributesManager& manager = E.facet_attribute_manager() ; 
            attributes.bind_one_attribute( i, manager, attribute_name ) ;
        }
    }

    /*!
    * @brief Bind named GEO::Attribute on theGeoModel elements vertices
    * @warning It is up to the client to unbind the attribute
    * @pre Elements of mesh_element_type are GeoModelMeshElement
    */
    template< class T >
    void create_attributes_on_geomodel_element_cells(
        const GeoModel& geomodel,
        GeoModelElement::TYPE geomodel_element_type,
        const std::string& attribute_name,
        AttributeVector<T>& attributes )
    {
        index_t nb_elements = geomodel.nb_elements( geomodel_element_type ) ;
        attributes.resize( nb_elements ) ;
        for( index_t i = 0; i < nb_elements; ++i ) {
            const GeoModelMeshElement& E = geomodel.mesh_element( GME::gme_t( geomodel_element_type, i ) );
            GEO::AttributesManager& manager = E.cell_attribute_manager() ;
            attributes.bind_one_attribute( i, manager, attribute_name ) ;
        }
    }


    /*!
    * Compute the tetrahedral mesh of the input structural model
    * @param[in] M GeoModel to tetrahedralize
    * @param[in] method External mesher used, Tetgen by default
    * @param[in] region_id Region to mesh. By default it set to NO_ID and all regions are meshed.
    * @param[in] add_steiner_points if true (default value), the mesher will add some points inside the region.
    */
    void RINGMESH_API tetrahedralize(
        GeoModel& M,
        const std::string& method = "TetGen",
        index_t region_id = NO_ID,
        bool add_steiner_points = true ) ;

    /*!
    * Compute the tetrahedral mesh of the input structural model
    * @param[in] M GeoModel to tetrahedralize
    * @param[in] method External mesher used
    * @param[in] region_id Region to mesh. If set to NO_ID and all regions are meshed.
    * @param[in] add_steiner_points if true, the mesher will add some points inside the region.
    * @param[in] internal_vertices points inside the domain to constrain mesh generation.
    * There is one vector per region.
    */
    void RINGMESH_API tetrahedralize(
        GeoModel& M,
        const std::string& method,
        index_t region_id,
        bool add_steiner_points,
        const std::vector< std::vector< vec3 > >& internal_vertices ) ;


    /*!
    * @brief Translates the boundary model by a vector.
    *
    * Every single mesh of the boundary model is translated:
    * corners, lines and surfaces.
    *
    * @param[in] M GeoModel on which compute the translation
    * @param[in] translation_vector vector of translation.
    */
    void RINGMESH_API translate( GeoModel& M, const vec3& ) ;

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
    * @param[in] degrees true is \p angle is in degrees, false
    * if in radians.
    */
    void RINGMESH_API rotate(
        GeoModel& M, 
        const vec3& origin,
        const vec3& axis,
        float64 angle,
        bool degrees = false ) ;


    /*-----------------------------------------------------------------------*/

    /*!
    * @brief Compute the size (volume, area, length) of an Element
    * @param[in] E Element to evaluate
    */
    double RINGMESH_API model_element_size( const GeoModelElement& E ) ;

    /*!
    * Compute the size (volume, area, length) of an Element cell (cell, facet, edge)
    * @param[in] E Element to evaluate
    * @param[in] c the cell index
    */
    double RINGMESH_API model_element_cell_size( const GeoModelElement& E, index_t c ) ;

    /*!
    * @brief Compute the center of a GeoModelElement
    *
    * @param[in] E Element to evaluate
    * @return The coordinates of the center
    */
    vec3 RINGMESH_API model_element_center( const GeoModelElement& E ) ;

    /*!
    * @brief Compute the centroid of a GeoModelMeshElement cell (cell, facet, edge)
    *
    * @param[in] E Element to evaluate
    * @param[in] c the cell index
    * @return The coordinates of the center
    *
    * @pre E has a valid mesh.
    */
    vec3 RINGMESH_API model_element_cell_center(
        const GeoModelMeshElement& E, index_t c ) ;

}



#endif 
