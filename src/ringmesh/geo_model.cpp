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
#include <ringmesh/geo_model_mesh_entity.h>
#include <ringmesh/geo_model_entity.h>
#include <ringmesh/geo_model_geological_entity.h>

#include <ringmesh/algorithm.h>

namespace RINGMesh {

    typedef GME::gme_t gme_t ;

    // Not the smartest but hopefully compiles in C++98
    const std::vector< std::string > mesh_entity_types = {
        Corner::type_name_static(),
        Line::type_name_static(),
        Surface::type_name_static(),
        Region::type_name_static()} ;

    GeoModel::GeoModel()
        :
            mesh( *this ),
            universe_( *this ),
            wells_( nil ),
            mesh_entity_types_( mesh_entity_types )             
    {
    }

    GeoModel::~GeoModel()
    {        
        for( index_t i = 0; i < corners_.size(); ++i ) {
            delete corners_[i] ;
        }
        for( index_t i = 0; i < lines_.size(); ++i ) {
            delete lines_[i] ;
        }
        for( index_t i = 0; i < surfaces_.size(); ++i ) {
            delete surfaces_[i] ;
        }
        for( index_t i = 0; i < regions_.size(); ++i ) {
            delete regions_[i] ;
        }

        for( index_t i = 0 ; i < geological_entities_.size(); ++i ){
            for( index_t j = 0 ; j < geological_entities_[i].size(); ++j ) {
                delete geological_entities_[i][j] ;
            }
        }
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

    index_t GeoModel::geological_entity_type( const std::string& type ) const
    {
        return find( geological_entity_types_, type ) ;
    }

} // namespace
