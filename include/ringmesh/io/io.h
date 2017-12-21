/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include <ringmesh/io/common.h>

#include <memory>

#include <ringmesh/basic/factory.h>

const char COMMA = ',';
const char EOL = '\n';
const char SPACE = ' ';
const char TAB = '\t';

/*!
 * @file Global input - output functions of RINGMesh
 * @author Various
 */

namespace RINGMesh
{
    class StratigraphicColumn;
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( WellGroup );

    ALIAS_3D( GeoModel );
    ALIAS_3D( WellGroup );
} // namespace RINGMesh

namespace GEO
{
    class MeshSubElementsStore;
} // namespace GEO

namespace RINGMesh
{
    /*!
     * Compares the contains of two files
     * @param[in] f1 the first filename
     * @param[in] f2 the second filename
     * @return return True if the files are identical
     */
    bool io_api compare_files( const std::string& f1, const std::string& f2 );
    /*!
     * Loads a GeoModel from a file
     * @param[out] geomodel the geomodel to fill
     * @param[in] filename the file to load
     */
    template < index_t DIMENSION >
    bool geomodel_load(
        GeoModel< DIMENSION >& geomodel, const std::string& filename );
    /*!
     * Saves a GeoModel to a file
     * @param[in] geomodel the geomodel to save
     * @param[in] filename the file to save
     */
    template < index_t DIMENSION >
    void geomodel_save(
        const GeoModel< DIMENSION >& geomodel, const std::string& filename );
    /*!
     * Loads a WellGroup from a file
     * @param[in] filename the file to load
     * @param][in,out] wells the wells to fill
     */
    void io_api well_load( const std::string& filename, WellGroup3D& wells );

    /*!
     * Returns the dimension of the GeoModel in the \p filename
     */
    index_t io_api find_geomodel_dimension( const std::string& filename );

    template < index_t DIMENSION >
    class io_api GeoModelInputHandler
    {
        ringmesh_disable_copy_and_move( GeoModelInputHandler );

    public:
        virtual ~GeoModelInputHandler() = default;

        static void initialize();

        static std::unique_ptr< GeoModelInputHandler< DIMENSION > > get_handler(
            const std::string& filename );

        bool load_geomodel(
            const std::string& filename, GeoModel< DIMENSION >& geomodel );

        virtual index_t dimension( const std::string& filename ) const
        {
            ringmesh_unused( filename );
            return DIMENSION;
        }

    protected:
        GeoModelInputHandler() = default;
        virtual void load(
            const std::string& filename, GeoModel< DIMENSION >& geomodel ) = 0;
    };

    ALIAS_2D_AND_3D( GeoModelInputHandler );

    template < index_t DIMENSION >
    using GeoModelInputHandlerFactory =
        Factory< std::string, GeoModelInputHandler< DIMENSION > >;

    ALIAS_2D_AND_3D( GeoModelInputHandlerFactory );

    template < index_t DIMENSION >
    class io_api GeoModelOutputHandler
    {
        ringmesh_disable_copy_and_move( GeoModelOutputHandler );

    public:
        virtual ~GeoModelOutputHandler() = default;

        static void initialize();

        static std::unique_ptr< GeoModelOutputHandler< DIMENSION > >
            get_handler( const std::string& filename );

        void save_geomodel( const GeoModel< DIMENSION >& geomodel,
            const std::string& filename );

    protected:
        GeoModelOutputHandler() = default;
        virtual void save( const GeoModel< DIMENSION >& geomodel,
            const std::string& filename ) = 0;
    };

    ALIAS_2D_AND_3D( GeoModelOutputHandler );

    template < index_t DIMENSION >
    using GeoModelOutputHandlerFactory =
        Factory< std::string, GeoModelOutputHandler< DIMENSION > >;

    ALIAS_2D_AND_3D( GeoModelOutputHandlerFactory );

    /***************************************************************************/
    class io_api WellGroupIOHandler
    {
        ringmesh_disable_copy_and_move( WellGroupIOHandler );

    public:
        virtual ~WellGroupIOHandler() = default;

        static void initialize();

        static std::unique_ptr< WellGroupIOHandler > get_handler(
            const std::string& filename );

        virtual void load( const std::string& filename, WellGroup3D& mesh ) = 0;

        virtual void save(
            const WellGroup3D& mesh, const std::string& filename ) = 0;

    protected:
        WellGroupIOHandler() = default;

    private:
        static std::unique_ptr< WellGroupIOHandler > create(
            const std::string& format );
    };
    using WellGroupIOHandlerFactory =
        Factory< std::string, WellGroupIOHandler >;

    /***************************************************************************/

    void io_api mesh_initialize();

    /*********************************************************************************************/
    class io_api StratigraphicColumnIOHandler
    {
        ringmesh_disable_copy_and_move( StratigraphicColumnIOHandler );

    public:
        virtual ~StratigraphicColumnIOHandler() = default;

        static void initialize();

        static std::unique_ptr< StratigraphicColumnIOHandler > get_handler(
            const std::string& filename );

        virtual void load( const std::string& filename,
            StratigraphicColumn& column,
            GeoModel3D& geomodel ) = 0;

        virtual void save( const StratigraphicColumn& column,
            const std::string& filename ) = 0;

    protected:
        StratigraphicColumnIOHandler() = default;

    private:
        static std::unique_ptr< StratigraphicColumnIOHandler > create(
            const std::string& format );
    };
    using StratigraphicColumnIOHandlerFactory =
        Factory< std::string, StratigraphicColumnIOHandler >;
} // namespace RINGMesh
