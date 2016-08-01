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

#ifndef __RINGMESH_GEO_MODEL_BUILDER_RINGMESH__
#define __RINGMESH_GEO_MODEL_BUILDER_RINGMESH__

#include <ringmesh/basic/common.h>

#include <vector>
#include <string>
#include <stack>

#include <geogram/basic/line_stream.h>

#include <third_party/zlib/unzip.h>

#include <ringmesh/geomodel/geo_model_builder.h>

#define MAX_FILENAME 512
#define READ_SIZE 8192

/*!
 * @file ringmesh/geo_model_builder_ringmesh.h
 * @brief Classes to build GeoModel from various inputs
 * @author 
 */


namespace RINGMesh {

    class RINGMESH_API GeoModelBuilderGM: public GeoModelBuilderFile {
    public:
        GeoModelBuilderGM( GeoModel& model, const std::string& filename )
            : GeoModelBuilderFile( model, filename )
        {
        }
        virtual ~GeoModelBuilderGM()
        {
        }

    private:
        void load_geological_entities( GEO::LineInput& file_line ) ;

        /*!
         * @brief Load meshes of mesh entities of one type from a zip file
         * @param[in] gme_t the GeoModelMeshEntity type
         * @param[in] uz the zip file
         */
        void load_meshes( const std::string& type, unzFile& uz ) ;

        void load_file() ;

        /*!
         * @brief Unzip a file in a zip file and set it to the current unZIP file
         */
        void unzip_one_file( unzFile& uz, const char filename[MAX_FILENAME] ) ;

        void load_mesh_entities( GEO::LineInput& file_line ) ;
    } ;
}

#endif
