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
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*! \author Jeanne Pellerin and Arnaud Botella */

#include <ringmesh/geo_model.h>

namespace RINGMesh {

    typedef GME::gme_t gme_t ;

    GeoModel::GeoModel()
        :
            mesh( *this ),            
            wells_( nil )
    {
    }

    GeoModel::~GeoModel()
    {
        for( index_t t = GME::CORNER; t < GME::NO_TYPE; ++t ) {
            GME::TYPE T = (GME::TYPE) t ;
            for( index_t i = 0; i < nb_elements( T ); ++i ) {
                delete elements( T )[i] ;
            }
        }
    }
   

    /*!
     * Copies a GeoModel in another one
     * @param[in] from GeoModel to copy
     */
    void GeoModel::copy( const GeoModel& from )
    {
        copy_macro_topology( from ) ;
        copy_meshes( from ) ;
    }

    /*!
     * @brief Mini-factory. Create an empty element of the right type 
     */
    GME* create_element( GME::TYPE T )
    {

        if( T == GME::CORNER ) {
            return new Corner ;
        } else if( T == GME::LINE ) {
            return new Line ;
        } else if( T == GME::SURFACE ) {
            return new Surface ;
        } else if( T > GME::SURFACE && T < GME::NO_TYPE ) {
            return new GeoModelElement ;
        } else {
            return nil ;
        }
    }

    /*!
     * @brief Copy macro information from a model
     * @details Copy the all the model elements and their relationship ignoring their geometry
     *
     * @param[in] from Model to copy the information from
     */
    void GeoModel::copy_macro_topology( const GeoModel& from )
    {
        name_ = from.name_ ;
        for( index_t t = GME::CORNER; t < GME::NO_TYPE; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            std::vector< GME* >& store = modifiable_elements( T ) ;
            store.resize( from.nb_elements( T ), nil ) ;

            for( index_t e = 0; e < nb_elements( T ); ++e ) {
                store[e] = create_element( T ) ;
                ringmesh_debug_assert( store[ e ] != nil ) ;
            }
            RINGMESH_PARALLEL_LOOP
            for( index_t e = 0; e < nb_elements( T ); ++e ) {
                store[e]->copy_macro_topology( from.element( gme_t( T, e ) ),
                    *this ) ;
            }
        }
        universe_.copy_macro_topology( from.universe_, *this ) ;

        nb_elements_per_type_ = from.nb_elements_per_type_ ;
    }

    /*!
     * @brief Copy meshes from a model
     * @details Copy the all the element meshes
     *
     * @param[in] from Model to copy the meshes from
     *
     * @pre The two models must have the same number of elements
     */
    void GeoModel::copy_meshes( const GeoModel& from )
    {
        for( index_t t = GME::CORNER; t < GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            RINGMESH_PARALLEL_LOOP
            for( index_t e = 0; e < elements( T ).size(); ++e ) {
                GeoModelMeshElement* E =
                    dynamic_cast< GeoModelMeshElement* >( elements( T )[e] ) ;
                ringmesh_debug_assert( E != nil ) ;
                const GeoModelMeshElement& E_from =
                    dynamic_cast< const GeoModelMeshElement& >( from.element(
                        GME::gme_t( T, e ) ) ) ;

                E->unbind_attributes() ;
                E->mesh().copy( E_from.mesh() ) ;
                E->bind_attributes() ;
            }
        }
    }  

  
    /*!
     * Associates a WellGroup to the GeoModel
     * @param[in] wells the WellGroup
     */
    void GeoModel::set_wells( const WellGroup* wells )
    {
        wells_ = wells ;
    }

} // namespace
