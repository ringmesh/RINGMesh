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

#include <memory>

#include <geogram/basic/line_stream.h>

#include <zlib/unzip.h>

#include <ringmesh/geomodel/geomodel_builder_file.h>

/*!
 * @file ringmesh/geomodel_builder_ringmesh.h
 * @brief Classes to build GeoModel from various inputs
 * @author 
 */

namespace RINGMesh {
    template< index_t DIMENSION > class GeoModelBuilderGMImpl;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelBuilderGM final : public GeoModelBuilderFile< 3 > {
    public:
        static const index_t NB_VERSION = 3;
        GeoModelBuilderGM( GeoModel< 3 >& geomodel, const std::string& filename );
        virtual ~GeoModelBuilderGM();

    private:
        void load_geological_entities( const std::string& geological_entity_file );

        /*!
         * @brief Load meshes of all the mesh entities from a zip file
         * @param[in] uz the zip file
         */
        void load_meshes( unzFile& uz );

        virtual void load_file() final;

        void load_mesh_entities( const std::string& mesh_entity_file );

    private:
        index_t file_version_;
        std::unique_ptr< GeoModelBuilderGMImpl< 3 > > version_impl_[NB_VERSION];
    };
}
