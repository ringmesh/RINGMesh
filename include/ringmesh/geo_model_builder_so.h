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

#ifndef __RINGMESH_GEO_MODEL_BUILDER_SO__
#define __RINGMESH_GEO_MODEL_BUILDER_SO__

#include <ringmesh/common.h>
#include <ringmesh/geo_model_builder.h>

#include <geogram/basic/factory.h>

namespace RINGMesh {
    class GeoModelBuilderTSolid ;
    struct VertexMap ;
    struct TSolidLoadingStorage ;
}

namespace RINGMesh {
    void tsolid_import_factory_initialize() ;

    /*!
     * @brief Build a GeoModel (with meshed regions) from a Gocad TSolid (file.so)
     */
    class TSolidLineParser: public GEO::Counted {
    ringmesh_disable_copy(TSolidLineParser) ;
    public:
        static TSolidLineParser* create(
            const std::string& keyword,
            GeoModelBuilderTSolid& gm_builder,
            GeoModel& geomodel ) ;
        virtual void execute(
            const GEO::LineInput& line,
            TSolidLoadingStorage& load_storage ) = 0 ;

    protected:
        TSolidLineParser()
            : Counted(), builder_( nil ), geomodel_( nil )
        {
        }
        virtual ~TSolidLineParser()
        {
        }

        GeoModelBuilderTSolid& builder()
        {
            ringmesh_assert( builder_ != nil ) ;
            return *builder_ ;
        }

        GeoModel& geomodel()
        {
            ringmesh_assert( geomodel_ != nil ) ;
            return *geomodel_ ;
        }

        virtual void set_builder( GeoModelBuilderTSolid& builder )
        {
            builder_ = &builder ;
        }

        virtual void set_geomodel( GeoModel& geomodel )
        {
            geomodel_ = &geomodel ;
        }

    private:
        GeoModelBuilderTSolid* builder_ ;
        GeoModel* geomodel_ ;
    } ;

    typedef GEO::SmartPointer< TSolidLineParser > TSolidLineParser_var ;
    typedef GEO::Factory0< TSolidLineParser > TSolidLineParserFactory ;
#define ringmesh_register_TSolidLineParser_creator(type, name) \
                        geo_register_creator(TSolidLineParserFactory, type, name)

}

namespace RINGMesh {
    /*!
     * @brief Builds a meshed GeoModel from a Gocad TSolid (file.so)
     */
    class RINGMESH_API GeoModelBuilderTSolid: public GeoModelBuilderFile {
    public:
        GeoModelBuilderTSolid( GeoModel& model, const std::string& filename )
            :
                GeoModelBuilderFile( model, filename ), file_line_( filename )
        {
            if( !file_line_.OK() ) {
                throw RINGMeshException( "I/O", "Failed to open file " + filename ) ;
            }
        }
        virtual ~GeoModelBuilderTSolid()
        {
        }

    private:
        virtual void load_file() ;

        /*!
         * @brief Parses the file and loads the GeoModel
         * @details The GeoModel loaded by this function is not valid because
         * some computation are still not done (i.e., surface internal borders,
         * lines and corners computation, boundary links between region and
         * surface, contacts)
         */
        void read_file() ;

        /*!
         * @brief Reads the first word of the current line (keyword)
         * and executes the good action with the information of the line
         * @details Uses the TsolidLineParser factory
         */
        void read_line( TSolidLoadingStorage& load_utils ) ;

    private:
        GEO::LineInput file_line_ ;
        friend class RINGMesh::TSolidLineParser ;
    } ;
}

#endif
