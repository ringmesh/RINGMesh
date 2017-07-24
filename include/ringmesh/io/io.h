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

#include <zlib/zip.h>
#include <zlib/unzip.h>

#include <geogram/basic/factory.h>
#include <geogram/basic/string.h>

#define MAX_FILENAME 512
#define READ_SIZE 8192

const char COMMA = ',';
const char EOL = '\n';
const char SPACE = ' ';
const char TAB = '\t';

/*!
 * @file Global input - output functions of RINGMesh
 * @author Various
 */

namespace RINGMesh {
    class StratigraphicColumn;
    template< index_t DIMENSION > class GeoModel;
    template< index_t DIMENSION > class WellGroup;

    CLASS_DIMENSION_ALIASES( GeoModel );
    CLASS_DIMENSION_ALIASES( WellGroup );
}

namespace GEO {
    class MeshSubElementsStore;
}

namespace RINGMesh {

    /*!
     * Compares the contains of two files
     * @param[in] f1 the first filename
     * @param[in] f2 the second filename
     * @return return True if the files are identical
     */
    bool RINGMESH_API compare_files( const std::string& f1, const std::string& f2 );
    /*!
     * Loads a GeoModel from a file
     * @param[out] geomodel the geomodel to fill
     * @param[in] filename the file to load
     */
    template< index_t DIMENSION >
    bool geomodel_load(
        GeoModel< DIMENSION >& geomodel,
        const std::string& filename );
    /*!
     * Saves a GeoModel to a file
     * @param[in] geomodel the geomodel to save
     * @param[in] filename the file to save
     */
    template< index_t DIMENSION >
    void geomodel_save(
        const GeoModel< DIMENSION >& geomodel,
        const std::string& filename );
    /*!
     * Loads a WellGroup from a file
     * @param[in] filename the file to load
     * @param][in,out] wells the wells to fill
     */
    void RINGMESH_API well_load(
        const std::string& filename,
        WellGroup< 3 >& wells );

    /*!
     * Returns the dimension of the GeoModel in the \p filename
     */
    index_t RINGMESH_API find_geomodel_dimension( const std::string& filename );

    template< index_t DIMENSION >
    class GeoModelIOHandler: public GEO::Counted {
    public:
        virtual ~GeoModelIOHandler() = default;

        static void initialize_geomodel_output();

        static std::unique_ptr< GeoModelIOHandler< DIMENSION > > get_handler(
            const std::string& filename );

        bool load_geomodel(
            const std::string& filename,
            GeoModel< DIMENSION >& geomodel );

        void save_geomodel(
            const GeoModel< DIMENSION >& geomodel,
            const std::string& filename );

        virtual index_t dimension( const std::string& filename ) const
        {
            return DIMENSION;
        }

    protected:
        GeoModelIOHandler() = default;
        virtual void load(
            const std::string& filename,
            GeoModel< DIMENSION >& geomodel ) = 0;

        virtual void save(
            const GeoModel< DIMENSION >& geomodel,
            const std::string& filename ) = 0;

    private:
        static GeoModelIOHandler* create( const std::string& format );
    };

    CLASS_DIMENSION_ALIASES( GeoModelIOHandler );

    template< index_t DIMENSION >
    using GeoModelIOHandlerFactory = GEO::Factory0< GeoModelIOHandler< DIMENSION > >;

    using GeoModelIOHandlerFactory2D = GeoModelIOHandlerFactory< 2 >;
    using GeoModelIOHandlerFactory3D = GeoModelIOHandlerFactory< 3 >;

#define ringmesh_register_GeoModelIOHandler2D_creator( type, name ) \
    geo_register_creator( GeoModelIOHandlerFactory2D, type, name )

#define ringmesh_register_GeoModelIOHandler3D_creator( type, name ) \
    geo_register_creator( GeoModelIOHandlerFactory3D, type, name )

    /***************************************************************************/
    class RINGMESH_API WellGroupIOHandler: public GEO::Counted {
    public:
        virtual ~WellGroupIOHandler() = default;

        static void initialize();

        static std::unique_ptr< WellGroupIOHandler > get_handler(
            const std::string& filename );

        virtual void load( const std::string& filename, WellGroup3D& mesh ) = 0;

        virtual void save(
            const WellGroup3D& mesh,
            const std::string& filename ) = 0;

    protected:
        WellGroupIOHandler() = default;

    private:
        static WellGroupIOHandler* create( const std::string& format );
    };
    using WellGroupIOHandlerFactory = GEO::Factory0< WellGroupIOHandler >;

#define ringmesh_register_WellGroupIOHandler_creator( type, name ) \
    geo_register_creator( WellGroupIOHandlerFactory, type, name )

    /***************************************************************************/

    void RINGMESH_API mesh_initialize();

    void RINGMESH_API zip_file( zipFile zf, const std::string& name );

    void RINGMESH_API unzip_file( unzFile uz, const char filename[MAX_FILENAME] );

    void RINGMESH_API unzip_current_file(
        unzFile uz,
        const char filename[MAX_FILENAME] );

    /*********************************************************************************************/
    class RINGMESH_API StratigraphicColumnIOHandler: public GEO::Counted {
    public:
        virtual ~StratigraphicColumnIOHandler() = default;

        static void initialize();

        static StratigraphicColumnIOHandler* create( const std::string& format );

        static StratigraphicColumnIOHandler* get_handler(
            const std::string& filename );

        virtual void load(
            const std::string& filename,
            StratigraphicColumn& column,
            GeoModel3D& geomodel ) = 0;

        virtual void save(
            const StratigraphicColumn& column,
            const std::string& filename ) = 0;

    protected:
        StratigraphicColumnIOHandler() = default;

    };
    typedef GEO::SmartPointer< StratigraphicColumnIOHandler > StratigraphicColumnIOHandler_var;
    typedef GEO::Factory0< StratigraphicColumnIOHandler > StratigraphicColumnIOHandlerFactory;

#define ringmesh_register_StratigraphicColumnIOHandler_creator( type, name ) \
		geo_register_creator( StratigraphicColumnIOHandlerFactory, type, name )
}
