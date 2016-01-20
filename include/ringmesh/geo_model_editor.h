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

#ifndef __RINGMESH_GEO_MODEL_EDITOR__
#define __RINGMESH_GEO_MODEL_EDITOR__

#include <ringmesh/common.h>

#include <ringmesh/geo_model_element.h> 
#include <ringmesh/geo_model.h> 

#include <set>

namespace RINGMesh {

    /*!
     * @brief Basic edition of a GeoModel
     *        Topological level. NO GEOMETRY HERE.
     *        GeoModelBuilder deals with the dirty geometry
     */
    class RINGMESH_API GeoModelEditor {
    public:
        GeoModelEditor( GeoModel& M )
            :model_( M )
        {}

        ~GeoModelEditor()
        {}

        void copy( const GeoModel& from ) ;
        void copy_macro_topology( const GeoModel& from ) ;
        void copy_meshes( const GeoModel& from ) ;

        /*!
        *@brief The model under construction
        */
        const GeoModel& model() const
        {
            return model_ ;
        }

        /*!
        *@brief Set the name of the model
        */
        void set_model_name( const std::string& name )
        {
            model_.name_ = name ;
        }

        GME::gme_t create_element( GME::TYPE e_type ) ;

        void create_geomodel_elements( GME::TYPE element_type, index_t nb_elements )
        {
            for( index_t i = 0; i < nb_elements; ++i ) {
                create_element( element_type ) ;
            }
        }


        /*! @}
        * \name Creation - Deletion - Access to GeoModelElements.
        * @{
        */

        /*!
        * @brief Reference to a modifiable element of the model
        * @pre The id must refer to a valid element of the model
        * @note Stupid override of a function of GeoModel
        */
        GeoModelElement& element( const GME::gme_t& id ) const
        {
            return *model_.element_ptr( id ) ;
        }

        /*!
        * @brief Reference to a modifiable meshed element of the model
        * @pre The id must refer to a valid element.
        * @note Stupid override of a function of GeoModel
        */
        GeoModelMeshElement& mesh_element( const GME::gme_t& id ) const
        {
            ringmesh_debug_assert( GME::has_mesh( id.type ) ) ;
            return dynamic_cast<GeoModelMeshElement&>( element( id ) ) ;
        }

        /*!
        * @brief Reference to a modifiable meshed element of the model
        */
        GeoModelMeshElement& mesh_element( GeoModelElement::TYPE element_type,
                                           index_t element_index ) const
        {
            return mesh_element( GME::gme_t( element_type, element_index ) ) ;
        }


        /*! @}
        * \name Filling GeoModelElement attributes.
        * @{
        */

        /*!
        * @brief Complete missing information in GeoModelElements
        * boundaries - in_boundary - parent - children
        */
        bool complete_element_connectivity() ;

        /*!
        * @brief Fill the boundaries of all elements of the given type
        * @details If the boundary elements do not have any in_boundary
        * information, nothing is done.
        */
        void fill_elements_boundaries( GME::TYPE type ) ;

        /*!
        * @brief Fill the in_boundary vector of all elements of the given type
        * @details If the in_boundary elements do not have any boundary
        * information, nothing is done, and model construction will eventually fail.
        */
        void fill_elements_in_boundaries( GME::TYPE type ) ;

        /*!
        * @brief Fill the parent of all elements of the given type
        * @details If the parents do not have any child nothing is done.
        */
        void fill_elements_parent( GME::TYPE type ) ;

        /*!
        * @brief Fill the children of all elements of the given type
        * @details If the children elements do not have any parent information
        * nothing is done.
        */
        void fill_elements_children( GME::TYPE type ) ;


        void set_element_name( const GME::gme_t& t, const std::string& name ) const
        {
            element( t ).name_ = name ;
        }

        void set_element_geol_feature( const GME::gme_t& t, GME::GEOL_FEATURE geol ) const
        {
            element( t ).geol_feature_ = geol ;
        }

        void add_element_boundary( const GME::gme_t& t,
                                   const GME::gme_t& boundary,
                                   bool side = false ) const
        {
            ringmesh_debug_assert( boundary.is_defined() ) ;
            ringmesh_debug_assert( GME::boundary_type( t.type ) == boundary.type ) ;
            element( t ).boundaries_.push_back( boundary ) ;

            if( t.type == GME::REGION ) {
                dynamic_cast< Region& >( element( t ) ).sides_.push_back( side ) ;
            }
        }

        void set_element_boundary( const GME::gme_t& t,
                                   index_t id,
                                   const GME::gme_t& boundary,
                                   bool side = false ) const
        {
            /// No check on the validity of the index of the element boundary
            /// NO_ID is used to flag elements to delete
            ringmesh_debug_assert( GME::boundary_type( t.type ) == boundary.type ) ;
            ringmesh_debug_assert( id < element( t ).nb_boundaries() ) ;
            element( t ).boundaries_[id] = boundary ;

            if( t.type == GME::REGION ) {
                dynamic_cast< Region& >( element( t ) ).sides_[id] = side ;
            }
        }

        void add_element_in_boundary( const GME::gme_t& t,
                                      const GME::gme_t& in_boundary ) const
        {
            ringmesh_debug_assert( in_boundary.is_defined() ) ;
            ringmesh_debug_assert( GME::in_boundary_type( t.type ) == in_boundary.type ) ;
            element( t ).in_boundary_.push_back( in_boundary ) ;
        }

        void set_element_in_boundary( const GME::gme_t& t,
                                      index_t id,
                                      const GME::gme_t& in_boundary ) const
        {
            /// No check on the validity of the index of the element in_boundary
            /// NO_ID is used to flag elements to delete
            ringmesh_debug_assert( GME::in_boundary_type( t.type ) == in_boundary.type ) ;
            ringmesh_debug_assert( id < element( t ).nb_in_boundary() ) ;
            element( t ).in_boundary_[id] = in_boundary ;
        }

        void set_element_parent( const GME::gme_t& t, const GME::gme_t& parent_index ) const
        {
            ringmesh_debug_assert( GME::parent_type( t.type ) == parent_index.type ) ;
            element( t ).parent_ = parent_index ;
        }

        void add_element_child( const GME::gme_t& t, const GME::gme_t& child_index ) const
        {
            ringmesh_debug_assert( child_index.is_defined() ) ;
            ringmesh_debug_assert( GME::child_type( t.type ) == child_index.type ) ;
            element( t ).children_.push_back( child_index ) ;
        }

        void set_element_child( const GME::gme_t& t,
                                index_t id,
                                const GME::gme_t& child_index ) const
        {
            /// No check on the validity of the index of the element child_index
            /// NO_ID is used to flag elements to delete
            ringmesh_debug_assert( GME::child_type( t.type ) == child_index.type ) ;
            element( t ).children_[id] = child_index ;
        }

  
        /*!
        * @brief Set an element of the model.
        * @details It is on purpose that element validity is not checked.
        *          This way nil pointers can be set for a further element removal.
        * @param id Id card of the element to modify. The ownership of the previous element
        *           is given up by the GeoModel.
        * @param E Element to set. The ownership is transferred to the GeoModel.
        */
        void set_element( const GME::gme_t& id, GeoModelElement* E ) const
        {
            if( id.type < GME::NO_TYPE ) {
                model_.modifiable_elements( id.type )[ id.index ] = E ;
            } else {
                ringmesh_assert_not_reached;
            }
        }

        void remove_elements( const std::set< GME::gme_t >& elements ) ;

        /*!
         * @todo Could be moved in the API [JP]
         */
        bool get_dependent_elements( std::set< GME::gme_t >& elements ) const ;
        
        void remove_elements_and_dependencies(
            const std::set< GME::gme_t >& elements_to_remove ) ;

    protected:
        void delete_elements( std::vector< std::vector< index_t > >& to_erase ) ;

        void resize_elements( GME::TYPE type, index_t nb ) ;

        void erase_invalid_element_references( GeoModelElement& E ) ;

    private:
        void copy_element_topology( GeoModelElement& lhs,
            const GeoModelElement& rhs ) ;

    protected:
        GeoModel& model_ ;
    };
}

#endif
