/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
 * @file Implementation of all GeoModelEntities classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geomodel/geomodel_entity.h>

#include <ringmesh/geomodel/geomodel.h>

namespace RINGMesh {

    GeoModelEntity::~GeoModelEntity()
    {
    }

    Universe::Universe( const GeoModel& geomodel )
        : GeoModelEntity( geomodel, NO_ID )
    {
        name_ = universe_type_name();
    }

    GeoModelEntity::GEOL_FEATURE GeoModelEntity::determine_geological_type(
        const std::string& in )
    {
        if( in == "reverse_fault" ) {
            return REVERSE_FAULT;
        } else if( in == "normal_fault" ) {
            return NORMAL_FAULT;
        } else if( in == "fault" ) {
            return FAULT;
        } else if( in == "top" ) {
            return STRATI;
        } else if( in == "none" ) {
            // This might seem strange - but it seems that what's
            // Gocad is doing
            return STRATI;
        } else if( in == "topographic" ) {
            return STRATI;
        } else if( in == "unconformity" ) {
            return UNCONFORMITY;
        } else if( in == "boundary" ) {
            return VOI;
        } else {
            // Default case - no information
            return NO_GEOL;
        }
    }

    std::string GeoModelEntity::geol_name( GeoModelEntity::GEOL_FEATURE t )
    {
        switch( t ) {
            case STRATI:
                return "top";
            case FAULT:
                return "fault";
            case REVERSE_FAULT:
                return "reverse_fault";
            case NORMAL_FAULT:
                return "normal_fault";
            case UNCONFORMITY:
                return "unconformity";
            case VOI:
                return "boundary";
            case NO_GEOL:
                return "no_geological_feature";
            default:
                return "no_geological_feature";
                break;
        }
    }

    bool Universe::is_valid() const
    {
        return true;
    }

}
