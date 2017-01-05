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

#ifndef __RINGMESH_GEOMODEL_BUILDER_RINGMESH__
#define __RINGMESH_GEOMODEL_BUILDER_RINGMESH__

#include <ringmesh/basic/common.h>

#include <geogram/basic/line_stream.h>

#include <third_party/zlib/unzip.h>

#include <ringmesh/geomodel/geomodel_builder.h>

/*!
 * @file ringmesh/geomodel_builder_ringmesh.h
 * @brief Classes to build GeoModel from various inputs
 * @author 
 */

namespace RINGMesh {
    class GeoModelBuilderGMImpl ;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelBuilderGM: public GeoModelBuilderFile {
    public:
        static const index_t NB_VERSION = 2 ;
        GeoModelBuilderGM( GeoModel& geomodel, const std::string& filename ) ;
        virtual ~GeoModelBuilderGM() ;

    private:
        void load_geological_entities( GEO::LineInput& file_line ) ;

        /*!
         * @brief Load meshes of all the mesh entities from a zip file
         * @param[in] uz the zip file
         */
        void load_meshes( unzFile& uz ) ;

        void load_file() ;

        void load_mesh_entities( GEO::LineInput& file_line ) ;

    private:
        index_t file_version_ ;
        GeoModelBuilderGMImpl* version_impl_[NB_VERSION] ;
    } ;

    class RINGMESH_API OldGeoModelBuilderGM: public GeoModelBuilderFile {
    public:
        OldGeoModelBuilderGM( GeoModel& geomodel, const std::string& filename )
            : GeoModelBuilderFile( geomodel, filename )
        {
        }
        virtual ~OldGeoModelBuilderGM()
        {
        }

    private:

        /*!
         * @brief Load elements of one type from a zip file
         * @param[in] gme_t the GeoModelElement type
         * @param[in] uz the zip file
         */
        void load_entities( const std::string& old_type_name, unzFile& uz ) ;

        void load_file() ;
        /*!
         * @brief Load the topology. Topology is how corners, lines, surfaces and
         * regions are organized into contacts, interfaces and layers. It also contains
         * basics information on the GeoModel.
         */
        void load_topology( GEO::LineInput& file_line ) ;
        std::string match_nb_entities( const char* s ) const ;
        /*!
         * @brief Convert an old entity name to its new name.
         */
        EntityType type_name_old_to_new( const std::string& old_type_name ) const ;
        bool child_allowed( const char* s ) const ;
        /*!
         * @brief Load the connectivities. These are how corners are
         * connected to lines, lines connected to surfaces and surfaces
         * connected to regions
         */
        void load_connectivities( GEO::LineInput& file_line ) ;
    } ;
}

#endif
