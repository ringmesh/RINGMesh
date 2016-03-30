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

#ifndef __RINGMESH_DUPLICATE_INTERFACE_BUILDER__
#define __RINGMESH_DUPLICATE_INTERFACE_BUILDER__

#include <ringmesh/common.h>
#include <ringmesh/geo_model_builder.h>

/*!
 * @file ringmesh/duplicate_interface_builder.h
 * @brief Class to duplicate GeoModel Interface to
 * enable sliding along them (faults) and unconformal
 * mesh generation.
 * @author Benjamin Chauvin
 */

namespace RINGMesh {
    class GeoModel ;
}

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {

    class RINGMESH_API DuplicateInterfaceBuilder: public GeoModelBuilder {
    ringmesh_disable_copy(DuplicateInterfaceBuilder) ;
    public:
        DuplicateInterfaceBuilder( GeoModel& model ) ;
        virtual ~DuplicateInterfaceBuilder() ;
        void duplicate_interface( index_t interface_id_to_duplicate ) ;

    public:
        void get_new_surfaces( index_t interface_id_to_duplicate ) ;
    private:
        void build_merged_and_bad_lines(
            const std::map< index_t, std::vector< index_t > >& surfaces_boundary_regions,
            const std::string& side_name,
            std::vector< std::vector< index_t > >& to_erase_by_type,
            const GME::gme_t& sided_interface_gme_t) ;
    } ;

}

#endif
