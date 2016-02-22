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
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_GEO_MODEL__
#define __RINGMESH_GEO_MODEL__

#include <ringmesh/common.h>

#include <ringmesh/geo_model_element.h>
#include <ringmesh/geo_model_mesh.h>

#include <vector>

/*!
 * @file ringmesh/geo_model.h
 * @brief Class representing a geological structural model: GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 */

namespace RINGMesh {
    class WellGroup ;
}

namespace RINGMesh {

    /*!
     * @brief The class to describe a geological model represented 
     * by its boundary surfaces and whose regions can be optionally meshed
     */
    class RINGMESH_API GeoModel {
    ringmesh_disable_copy( GeoModel ) ;
        friend class GeoModelBuilder ;
        friend class GeoModelEditor ;

    public:
        const static index_t NO_ID = index_t( -1 ) ;

        /*!
         * @brief Constructs an empty GeoModel
         */
        GeoModel() ;

        /*!
         * @brief Deletes all the GeoModelElements of the GeoModel
         */
        virtual ~GeoModel() ;

        /*
         * @todo Implement a real copy_constructor and operator= [JP]
         */
        void copy( const GeoModel& from ) ;

        /*!
         * @brief Name of the model
         */
        const std::string& name() const
        {
            return name_ ;
        }

        /*!
         * \name Generic GeoModelElement accessors
         * @{
         */

        /*!
         * @brief Returns the number of elements of the given type
         * @details Default value is 0.
         * @param[in] type the element type
         */
        index_t nb_elements( GME::TYPE type ) const
        {
            if( type < GME::NO_TYPE ) {
                return elements( type ).size() ;
            } else if( type == GME::ALL_TYPES ) {
                ringmesh_assert( !nb_elements_per_type_.empty() ) ;
                return nb_elements_per_type_.back() ;
            } else {
                ringmesh_debug_assert_not_reached;
            return 0 ;
        }
    }

    /*!
     * @brief Returns a const reference the identified GeoModelElement
     * @param[in] id Type and index of the element. For the
     * pair (Region, NO_ID) universe region is returned.
     * @pre Element identification is valid.
     */
    const GeoModelElement& element( GME::gme_t id ) const
    {
        return *element_ptr( id ) ;
    }

    /*!
     * Convenient overload of element( GME::gme_t id )
     */
    const GeoModelElement& element( GME::TYPE element_type, index_t element_index ) const
    {
        return element( GME::gme_t( element_type, element_index ) ) ;
    }

    /*!
     * @brief Generic access to a meshed element
     * @pre Type of the element is CORNER, LINE, SURFACE, or REGION
     */
    const GeoModelMeshElement& mesh_element( GME::gme_t id ) const
    {
        ringmesh_assert( GME::has_mesh( id.type ) ) ;
        return dynamic_cast< const GeoModelMeshElement& >( element( id ) ) ;
    }

    /*!
     * Convenient overload of element( GME::gme_t id )
     */
    const GeoModelMeshElement& mesh_element( GME::TYPE element_type, index_t element_index ) const
    {
        return mesh_element( GME::gme_t( element_type, element_index ) ) ;
    }

    /*! @}
     * \name Specialized accessors.
     * @{
     */
    index_t nb_corners() const
    {
        return nb_elements( GME::CORNER ) ;
    }
    index_t nb_lines() const
    {
        return nb_elements( GME::LINE ) ;
    }
    index_t nb_surfaces() const
    {
        return nb_elements( GME::SURFACE ) ;
    }
    index_t nb_regions() const
    {
        return nb_elements( GME::REGION ) ;
    }
    index_t nb_contacts() const
    {
        return nb_elements( GME::CONTACT ) ;
    }
    index_t nb_interfaces() const
    {
        return nb_elements( GME::INTERFACE ) ;
    }
    index_t nb_layers() const
    {
        return nb_elements( GME::LAYER ) ;
    }

    const Corner& corner( index_t index ) const
    {
        // Yes, we could use static_cast, but I trust nobody and check [JP]
        return dynamic_cast< const Corner& >( *corners_.at( index ) ) ;
    }

    const Line& line( index_t index ) const
    {
        return dynamic_cast< const Line& >( element( GME::LINE, index ) ) ;
    }

    const Surface& surface( index_t index ) const
    {
        return dynamic_cast< const Surface& >( element( GME::SURFACE, index ) ) ;
    }

    const Region& region( index_t index ) const
    {
        return dynamic_cast<const Region&>( element( GME::REGION, index ) ) ;
    }

    const GeoModelElement& contact( index_t index ) const
    {
        return element( GME::CONTACT, index ) ;
    }

    const GeoModelElement& one_interface( index_t index ) const
    {
        return element( GME::INTERFACE, index ) ;
    }

    const GeoModelElement& layer( index_t index ) const
    {
        return element( GME::LAYER, index ) ;
    }

    const Region& universe() const
    {
        return dynamic_cast<const Region&>( element( GME::REGION, NO_ID ) ) ;
    }

    /*!
     * @}
     */
    /* @todo Wells have nothing at all to do here [JP] */
    void set_wells( const WellGroup* wells ) ;
    const WellGroup* wells() const
    {
        return wells_ ;
    }

private:
    /*!
     * @brief Convert a global BME index into a typed index
     * @details Relies on the nb_elements_per_type_ vector that
     *          must be up to date at all times
     *          See the GeoModelBuilder::end_model() function
     * @param[in] global A BME id of TYPE - ALL_TYPES
     * @return A BME id of an element of the model, or a invalid one if nothing found
     */
    inline GME::gme_t global_to_typed_id( const GME::gme_t& global ) const
    {
        ringmesh_debug_assert( global.type == GME::ALL_TYPES ) ;

        index_t t = NO_ID ;
        for( index_t i = 1 ; i < nb_elements_per_type_.size() ; i++ ) {
            if( global.index >= nb_elements_per_type_[i - 1]
                && global.index < nb_elements_per_type_[i] ) {
                t = i - 1 ;
                break ;
            }
        }
        GME::TYPE T = static_cast< GME::TYPE >( t ) ;
        if( T < GME::NO_TYPE ) {
            index_t i = global.index - nb_elements_per_type_[t] ;
            return GME::gme_t( T, i ) ;
        } else {
            ringmesh_debug_assert_not_reached ;
            return GME::gme_t() ;
        }
    }

    /*!
     * @brief Generic accessor to the storage of elements of the given type
     * @pre The type must be valid NO_TYPE or ALL_TYPES will throw an assertion
     */
    inline std::vector< GME* >& modifiable_elements( GME::TYPE type )
    {
        return const_cast< std::vector< GME* >& >( elements( type ) ) ;
    }

    /*!
     * @brief Generic accessor to the storage of elements of the given type
     * @pre The type must be valid. NO_TYPE or ALL_TYPES will throw an assertion
     */
    const std::vector< GME* >& elements( GME::TYPE type ) const
    {
        switch( type ) {
            case GME::CORNER:
            return corners_ ;
            case GME::LINE:
            return lines_ ;
            case GME::SURFACE:
            return surfaces_ ;
            case GME::REGION:
            return regions_ ;
            case GME::CONTACT:
            return contacts_ ;
            case GME::INTERFACE:
            return interfaces_ ;
            case GME::LAYER:
            return layers_ ;
            default:
            ringmesh_assert_not_reached ;
            return surfaces_ ;
        }
    }

    /*!
     * @brief Modifiable pointer to an element of the model
     * @param[in] id Type and index of the element. For the
     * pair (Region, NO_ID) universe region is returned.
     *
     * @pre Element identification is valid.
     */
    GeoModelElement* element_ptr( const GME::gme_t& id ) const
    {
        if( id.type == GME::REGION && id.index == NO_ID ) {
            return const_cast< Region*> ( &universe_ ) ;
        } else {
            if( id.type < GME::NO_TYPE ) {
                ringmesh_assert( id.index < nb_elements( id.type ) ) ;
                return elements( id.type )[ id.index ] ;
            } else if( id.type == GME::ALL_TYPES ) {
                return element_ptr( global_to_typed_id( id ) ) ;
            } else {
                ringmesh_assert_not_reached ;
                return const_cast< Region*> ( &universe_ ) ;
            }
        }
    }

    /*!
     * @brief Reference to a modifiable element of the model
     * @pre The id must refer to a valid element of the model
     */
    GeoModelElement& modifiable_element(
        const GME::gme_t& id ) const
    {
        return *element_ptr( id ) ;
    }

    /*!
     * @brief Reference to a modifiable meshed element of the model
     * @pre Assert in debug model that the given id refers to a meshed element.
     *      The id must refer to a valid element.
     */
    inline GeoModelMeshElement& modifiable_mesh_element(
        const GME::gme_t& id ) const
    {
        ringmesh_debug_assert( GME::has_mesh( id.type ) ) ;
        return dynamic_cast<GeoModelMeshElement&>(
            modifiable_element( id ) ) ;
    }

    /*!
     * @brief Clears and fills the model nb_elements_per_type_ vector
     * @details See global element access with GeoModel::element( BME::TYPE, index_t )
     */
    void init_global_model_element_access()
    {
        nb_elements_per_type_.clear() ;

        index_t count = 0 ;
        nb_elements_per_type_.push_back( count ) ;
        for( index_t type = GME::CORNER ; type < GME::NO_TYPE ; type++ ) {
            count += nb_elements( ( GME::TYPE ) type ) ;
            nb_elements_per_type_.push_back( count ) ;
        }
    }

public:
    GeoModelMesh mesh ;

private:
    // Name of the model
    std::string name_ ;

    /*!
     * \name Mandatory elements of the model
     * @{
     */

    /*!
     * @brief Elements that are Corner
     */
    std::vector< GeoModelElement* > corners_ ;
    /*!
     * @brief Elements that are Line
     */
    std::vector< GeoModelElement* > lines_ ;
    /*!
     * @brief Elements that are Surface
     */
    std::vector< GeoModelElement* > surfaces_ ;
    /*!
     * @brief Elements that are Region
     */
    std::vector< GeoModelElement* > regions_ ;

    /*!
     * The Region defining the model extension
     * @todo Put it as the last item in regions_ and do not forget to create it in the Builder
     */
    Region universe_ ;

    /*! @}
     * \name Optional geological elements
     * @{
     */

    /*!
     * @brief Elements of type CONTACT
     */
    std::vector< GeoModelElement* > contacts_ ;
    /*!
     * @brief Elements of type INTERFACE
     */
    std::vector< GeoModelElement* > interfaces_ ;

    /*!
     * @brief Elements of type LAYER
     */
    std::vector< GeoModelElement* > layers_ ;

    /*!
     * @}
     */

    /*
     * @brief Global access to GeoModelElements.
     * It MUST be updated if one element is added.
     * @warning It must be up to date at all times
     */
    std::vector< index_t > nb_elements_per_type_ ;

    /*! Optional WellGroup associated with the model
     * @todo Give a more general name - this could be anything [JP]
     */
    const WellGroup* wells_ ;
} ;

}

#endif
