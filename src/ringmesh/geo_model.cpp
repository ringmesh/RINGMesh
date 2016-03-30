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
 
/*!
 * @file Implementation of THE GeoModel
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geo_model.h>

#include <ringmesh/geo_model_builder.h>

namespace RINGMesh {

    typedef GME::gme_t gme_t ;

    GeoModel::GeoModel()
        :
            mesh( *this ),
            universe_( *this, NO_ID, "Universe", GME::NO_GEOL ),
            wells_( nil )
    {
        // Initialization for the next debug test
#ifdef RINGMESH_DEBUG
        for( index_t gme_type_itr = 0; gme_type_itr < GME::NO_TYPE;
            ++gme_type_itr ) {
            gme_modifiers_[gme_type_itr] = nil ;
        }
#endif
        /// @todo Review: This usage of this pointer in initialization list is a time bomb [JP]
        gme_modifiers_[GME::CORNER] = new CornerModifier( *this ) ;
        gme_modifiers_[GME::LINE] = new LineModifier( *this ) ;
        gme_modifiers_[GME::SURFACE] = new SurfaceModifier( *this ) ;
        gme_modifiers_[GME::REGION] = new RegionModifier( *this ) ;
        gme_modifiers_[GME::CONTACT] = new ContactModifier( *this ) ;
        gme_modifiers_[GME::INTERFACE] = new InterfaceModifier( *this ) ;
        gme_modifiers_[GME::LAYER] = new LayerModifier( *this ) ;

        // Test to check that no new GME::TYPE has been forgotten
#ifdef RINGMESH_DEBUG
        for( index_t gme_type_itr = 0; gme_type_itr < GME::NO_TYPE;
            ++gme_type_itr ) {
            ringmesh_assert( gme_modifiers_[gme_type_itr] != nil ) ;
        }
#endif
    }

    GeoModel::~GeoModel()
    {
        for( index_t t = GME::CORNER; t < GME::NO_TYPE; ++t ) {
            GME::TYPE T = (GME::TYPE) t ;
            const GeoModelElementModifier& gmef = gme_modifier( T ) ;
            for( index_t i = 0; i < nb_elements( T ); ++i ) {
                delete gmef.elements()[i] ;
            }
            delete &gmef ;
        }
    }

    /*!
     * Copies a GeoModel in another one
     * @param[in] from GeoModel to copy
     * 
     * @todo This shouln't be a member function because it does not do nothing
     * with what the class has.
     */
    void GeoModel::copy( const GeoModel& from )
    {
        GeoModelBuilder builder( *this ) ;
        builder.copy_macro_topology( from ) ;
        builder.copy_meshes( from ) ;
    }

    /*!
     * Associates a WellGroup to the GeoModel
     * @param[in] wells the WellGroup
     * @todo Review : What is this for ?
     * @todo Extend to other object types.
     */
    void GeoModel::set_wells( const WellGroup* wells )
    {
        wells_ = wells ;
    }

    GME* GeoModel::CornerModifier::new_element( index_t id ) const
    {
        return new Corner( model_, id ) ;
    }

    std::vector< GME* >& GeoModel::CornerModifier::elements()
    {
        return model_.corners_ ;
    }

    GME* GeoModel::LineModifier::new_element( index_t id ) const
    {
        return new Line( model_, id ) ;
    }

    std::vector< GME* >& GeoModel::LineModifier::elements()
    {
        return model_.lines_ ;
    }

    GME* GeoModel::SurfaceModifier::new_element( index_t id ) const
    {
        return new Surface( model_, id ) ;
    }

    std::vector< GME* >& GeoModel::SurfaceModifier::elements()
    {
        return model_.surfaces_ ;
    }

    GME* GeoModel::RegionModifier::new_element( index_t id ) const
    {
        return new Region( model_, id ) ;
    }

    std::vector< GME* >& GeoModel::RegionModifier::elements()
    {
        return model_.regions_ ;
    }

    std::vector< GME* >& GeoModel::ContactModifier::elements()
    {
        return model_.contacts_ ;
    }

    GME* GeoModel::ContactModifier::new_element( index_t id ) const
    {
        return new GeoModelElement( model_, GME::CONTACT, id ) ;
    }

    GME* GeoModel::InterfaceModifier::new_element( index_t id ) const
    {
        return new GeoModelElement( model_, GME::INTERFACE, id ) ;
    }

    std::vector< GME* >& GeoModel::InterfaceModifier::elements()
    {
        return model_.interfaces_ ;
    }

    GME* GeoModel::LayerModifier::new_element( index_t id ) const
    {
        return new GeoModelElement( model_, GME::LAYER, id ) ;
    }

    std::vector< GME* >& GeoModel::LayerModifier::elements()
    {
        return model_.layers_ ;
    }

}
