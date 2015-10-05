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

        /*! @}
        * \name Creation - Deletion - Access to GeoModelElements.
        * @{
        */

        /*!
        * @brief Reference to a modifiable element of the model
        * @pre The id must refer to a valid element of the model
        * @note Stupid override of a function of GeoModel
        */
        GeoModelElement& element(
            const GME::gme_t& id ) const
        {
            return *model_.element_ptr( id ) ;
        }

        /*!
        * @brief Reference to a modifiable meshed element of the model
        * @pre Assert in debug model that the given id refers to a meshed element.
        *      The id must refer to a valid element.
        * @note Stupid override of a function of GeoModel
        */
        GeoModelMeshElement& mesh_element(
            const GME::gme_t& id ) const
        {
            ringmesh_debug_assert( GME::has_mesh( id.type ) ) ;
            return dynamic_cast<GeoModelMeshElement&>( element( id ) ) ;
        }

        /*! @}
        * \name Filling GeoModelElement attributes.
        * @{
        */
        void set_model(
            const GME::gme_t& t,
            GeoModel* m )
        {
            element( t ).set_model( m ) ;
        }

        void set_element_index(
            const GME::gme_t& t )
        {
            element( t ).set_id( t.index ) ;
        }

        void set_element_name(
            const GME::gme_t& t,
            const std::string& name )
        {
            element( t ).set_name( name ) ;
        }

        void set_element_geol_feature(
            const GME::gme_t& t,
            GME::GEOL_FEATURE geol )
        {
            element( t ).set_geological_feature( geol ) ;
        }

        void add_element_boundary(
            const GME::gme_t& t,
            const GME::gme_t& boundary,
            bool side = false )
        {
            if( t.type == GME::REGION ) {
                dynamic_cast< Region& >(
                    element( t ) ).add_boundary( boundary, side ) ;
            } else {
                element( t ).add_boundary( boundary ) ;
            }
        }

        void add_element_in_boundary(
            const GME::gme_t& t,
            const GME::gme_t& in_boundary )
        {
            element( t ).add_in_boundary( in_boundary ) ;
        }

        void set_parent(
            const GME::gme_t& t,
            const GME::gme_t& parent_index )
        {
            element( t ).set_parent( parent_index ) ;
        }

        void add_child(
            const GME::gme_t& t,
            const GME::gme_t& child_index )
        {
            element( t ).add_child( child_index ) ;
        }

        // Universe
        void set_universe( const std::vector<
                           std::pair< index_t, bool > >& boundaries ) ;


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


    protected:
        GeoModel& model_ ;
    };
}

#endif