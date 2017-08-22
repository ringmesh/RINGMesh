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

#pragma once

#include <ringmesh/basic/common.h>
#include <ringmesh/geomodel/entity_type_manager.h>

/*!
 * @brief Structures and classes used to index elements in a GeoModel,
 * in the meshes of its entities etc.
 * @author Jeanne Pellerin & Arnaud Botella
 */

namespace RINGMesh {
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelGeologicalEntity );
}

namespace RINGMesh {

    /*!
     * @brief Vertex in a GeoModelEntity
     */
    struct GMEVertex {
        GMEVertex( gmme_id t, index_t vertex_id_in )
            : gmme( std::move( t ) ), v_index( vertex_id_in )
        {
        }
        bool operator==( const GMEVertex& rhs ) const
        {
            return gmme == rhs.gmme && v_index == rhs.v_index;
        }
        bool is_defined() const
        {
            return gmme.is_defined() && v_index != NO_ID;
        }
        /// GeoModelEntity index in the GeoModel that owns it
        gmme_id gmme;
        /// Index of the vertex in the GeoModelEntity
        index_t v_index { NO_ID };
    };

    template< index_t DIMENSION >
    class RINGMESH_API entity_range: public range {
    protected:
        entity_range( const GeoModel< DIMENSION >& geomodel, index_t last )
            : range( last ), geomodel_( geomodel )
        {
        }

    protected:
        const GeoModel< DIMENSION >& geomodel_;
    };

    template< index_t DIMENSION >
    class RINGMESH_API corner_range: public entity_range< DIMENSION > {
    public:
        corner_range( const GeoModel< DIMENSION >& geomodel )
            : entity_range< DIMENSION >( geomodel, geomodel.nb_corners() )
        {
        }
        const corner_range< DIMENSION >& begin() const
        {
            return *this;
        }
        const corner_range< DIMENSION >& end() const
        {
            return *this;
        }
        const Corner< DIMENSION >& operator*() const
        {
            return this->geomodel_.corner( this->iter_ );
        }
    };

    template< index_t DIMENSION >
    class RINGMESH_API line_range: public entity_range< DIMENSION > {
    public:
        line_range( const GeoModel< DIMENSION >& geomodel )
            : entity_range< DIMENSION >( geomodel, geomodel.nb_lines() )
        {
        }
        const line_range< DIMENSION >& begin() const
        {
            return *this;
        }
        const line_range< DIMENSION >& end() const
        {
            return *this;
        }
        const Line< DIMENSION >& operator*() const
        {
            return this->geomodel_.line( this->iter_ );
        }
    };

    template< index_t DIMENSION >
    class RINGMESH_API surface_range: public entity_range< DIMENSION > {
    public:
        surface_range( const GeoModel< DIMENSION >& geomodel )
            : entity_range< DIMENSION >( geomodel, geomodel.nb_surfaces() )
        {
        }
        const surface_range< DIMENSION >& begin() const
        {
            return *this;
        }
        const surface_range< DIMENSION >& end() const
        {
            return *this;
        }
        const Surface< DIMENSION >& operator*() const
        {
            return this->geomodel_.surface( this->iter_ );
        }
    };

    template< index_t DIMENSION >
    class RINGMESH_API region_range: public entity_range< DIMENSION > {
    public:
        region_range( const GeoModel< DIMENSION >& geomodel )
            : entity_range< DIMENSION >( geomodel, geomodel.nb_regions() )
        {
        }
        const region_range< DIMENSION >& begin() const
        {
            return *this;
        }
        const region_range< DIMENSION >& end() const
        {
            return *this;
        }
        const Region< DIMENSION >& operator*() const
        {
            return this->geomodel_.region( this->iter_ );
        }
    };

    template< index_t DIMENSION >
    class RINGMESH_API geol_entity_range: public entity_range< DIMENSION > {
    public:
        geol_entity_range(
            const GeoModel< DIMENSION >& geomodel,
            const GeologicalEntityType& geological_entity_type )
            :
                entity_range< DIMENSION >( geomodel,
                    geomodel.nb_geological_entities( geological_entity_type ) ),
                type_( geological_entity_type )
        {
        }
        const geol_entity_range< DIMENSION >& begin() const
        {
            return *this;
        }
        const geol_entity_range< DIMENSION >& end() const
        {
            return *this;
        }
        const GeoModelGeologicalEntity< DIMENSION >& operator*() const
        {
            return this->geomodel_.geological_entity( this->type_, this->iter_ );
        }

    protected:
        const GeologicalEntityType& type_;
    };
}
