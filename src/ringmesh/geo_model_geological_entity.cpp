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
 * @file Implementation of all GeoModelEntities classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geo_model_geological_entity.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/io.h>

namespace RINGMesh {

    const GeoModelMeshEntity& GeoModelGeologicalEntity::child( index_t x ) const
    {
        return model().mesh_entity( child_id( x ) ) ;
    }

    bool GeoModelGeologicalEntity::is_on_voi() const
    {
        for( index_t i = 0; i < nb_children(); i++ ) {
            if( !child( i ).is_on_voi() ) return false ;
        }
        return true ;
    }
    class MSHIOHandler2: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from GMSH not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            /// @todo after implementing GMMOrder
            throw RINGMeshException( "I/O",
                "Saving of a GeoModel from GMSH not implemented yet" ) ;
//                gm.set_duplicate_mode( FAULT ) ;

            std::ofstream out( filename.c_str() ) ;
            out.precision( 16 ) ;

            out << "$MeshFormat" << std::endl ;
            out << "2.2 0 8" << std::endl ;
            out << "$EndMeshFormat" << std::endl ;

            out << "$Nodes" << std::endl ;
        }
    } ;

    void GeoModelGeologicalEntity::initialize()
    {
        ringmesh_register_GeoModelGeologicalEntity_creator( Contact ) ;
        ringmesh_register_GeoModelGeologicalEntity_creator( Interface ) ;
        ringmesh_register_GeoModelGeologicalEntity_creator( Layer ) ;
    }
}
