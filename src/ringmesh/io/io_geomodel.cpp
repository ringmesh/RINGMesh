/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/io/io.h>

#include <cctype>
#include <deque>
#include <iomanip>

#include <tinyxml2.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/task_handler.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/core/well.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/geomodel_adapter_resqml.h>
#include <ringmesh/io/geomodel_builder_gocad.h>
#include <ringmesh/io/geomodel_builder_resqml.h>

#include <ringmesh/io/zip_file.h>

#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>
#include <ringmesh/mesh/point_set_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>
/*!
 * @file Implementation of classes loading GeoModels
 * @author Arnaud Botella and Antoine Mazuyer
 */

namespace
{
    using namespace RINGMesh;

#include "geomodel/io_abaqus.hpp"
#include "geomodel/io_adeli.hpp"
#include "geomodel/io_aster.hpp"
#include "geomodel/io_csmp.hpp"
#include "geomodel/io_feflow.hpp"
#include "geomodel/io_gm.hpp"
#include "geomodel/io_gprs.hpp"
#include "geomodel/io_mfem.hpp"
#include "geomodel/io_model3d.hpp"
#include "geomodel/io_msh.hpp"
#include "geomodel/io_resqml.hpp"
#include "geomodel/io_smesh.hpp"
#include "geomodel/io_stl.hpp"
#include "geomodel/io_stradivarius.hpp"
#include "geomodel/io_svg.hpp"
#include "geomodel/io_tetgen.hpp"
#include "geomodel/io_tsolid.hpp"
#include "geomodel/io_vtk.hpp"

    template < typename Class, typename Factory >
    std::unique_ptr< Class > create_handler( const std::string& format )
    {
        auto handler = Factory::create( format );
        if( !handler )
        {
            Logger::err( "I/O", "Currently supported file formats are: " );
            for( const std::string& name : Factory::list_creators() )
            {
                Logger::err( "I/O", " ", name );
            }

            throw RINGMeshException(
                "I/O", "Unsupported file format: ", format );
        }
        return handler;
    }
}

namespace RINGMesh
{
    template <>
    void GeoModelOutputHandler< 2 >::initialize()
    {
        GeoModelOutputHandlerFactory2D::register_creator< GeoModelHandlerGM2D >(
            "gm" );
        GeoModelOutputHandlerFactory2D::register_creator< MFEMIOHandler2D >(
            "mfem" );
    }

    template <>
    void GeoModelInputHandler< 2 >::initialize()
    {
        GeoModelInputHandlerFactory2D::register_creator< GeoModelHandlerGM2D >(
            "gm" );
        GeoModelInputHandlerFactory2D::register_creator<
            StradivariusIOHandler >( "model" );
        GeoModelInputHandlerFactory2D::register_creator< SVGIOHandler >(
            "svg" );
    }

    /*
     * Initializes the possible handler for IO files
     */
    template <>
    void GeoModelOutputHandler< 3 >::initialize()
    {
        GeoModelOutputHandlerFactory3D::register_creator< TetGenIOHandler >(
            "tetgen" );
        GeoModelOutputHandlerFactory3D::register_creator< TSolidIOHandler >(
            "so" );
        GeoModelOutputHandlerFactory3D::register_creator< CSMPIOHandler >(
            "csmp" );
        GeoModelOutputHandlerFactory3D::register_creator< AsterIOHandler >(
            "mail" );
        GeoModelOutputHandlerFactory3D::register_creator< VTKIOHandler >(
            "vtk" );
        // todo GPRS export is not working for the moment [AB]
        //        GeoModelOutputHandlerFactory3D::register_creator<
        //        GPRSIOHandler >( "gprs" );
        GeoModelOutputHandlerFactory3D::register_creator< MSHIOHandler >(
            "msh" );
        GeoModelOutputHandlerFactory3D::register_creator< MFEMIOHandler3D >(
            "mfem" );
        GeoModelOutputHandlerFactory3D::register_creator< GeoModelHandlerGM3D >(
            "gm" );
        GeoModelOutputHandlerFactory3D::register_creator< AbaqusIOHandler >(
            "inp" );
        GeoModelOutputHandlerFactory3D::register_creator< AdeliIOHandler >(
            "adeli" );
        GeoModelOutputHandlerFactory3D::register_creator< FeflowIOHandler >(
            "fem" );
        GeoModelOutputHandlerFactory3D::register_creator< MLIOHandler >( "ml" );
        GeoModelOutputHandlerFactory3D::register_creator< SMESHIOHandler >(
            "smesh" );
        GeoModelOutputHandlerFactory3D::register_creator< STLIOHandler >(
            "stl" );
        GeoModelOutputHandlerFactory3D::register_creator< RESQMLIOHandler >(
            "epc" );
    }

    template <>
    void GeoModelInputHandler< 3 >::initialize()
    {
        GeoModelInputHandlerFactory3D::register_creator< GeoModelHandlerGM3D >(
            "gm" );
        GeoModelInputHandlerFactory3D::register_creator< MLIOHandler >( "ml" );
        GeoModelInputHandlerFactory3D::register_creator< TSolidIOHandler >(
            "so" );
        GeoModelInputHandlerFactory3D::register_creator< RESQMLIOHandler >(
            "epc" );
    }

    /***************************************************************************/

    template < index_t DIMENSION >
    std::unique_ptr< GeoModelInputHandler< DIMENSION > >
        GeoModelInputHandler< DIMENSION >::get_handler(
            const std::string& filename )
    {
        return create_handler< GeoModelInputHandler< DIMENSION >,
            GeoModelInputHandlerFactory< DIMENSION > >(
            GEO::FileSystem::extension( filename ) );
    }

    template < index_t DIMENSION >
    bool GeoModelInputHandler< DIMENSION >::load_geomodel(
        const std::string& filename, GeoModel< DIMENSION >& geomodel )
    {
        load( filename, geomodel );
        Logger::out(
            "I/O", " Loaded geomodel ", geomodel.name(), " from ", filename );

        return is_geomodel_valid( geomodel );
    }

    /***************************************************************************/

    template < index_t DIMENSION >
    std::unique_ptr< GeoModelOutputHandler< DIMENSION > >
        GeoModelOutputHandler< DIMENSION >::get_handler(
            const std::string& filename )
    {
        return create_handler< GeoModelOutputHandler< DIMENSION >,
            GeoModelOutputHandlerFactory< DIMENSION > >(
            GEO::FileSystem::extension( filename ) );
    }

    template < index_t DIMENSION >
    void GeoModelOutputHandler< DIMENSION >::save_geomodel(
        const GeoModel< DIMENSION >& geomodel, const std::string& filename )
    {
        save( geomodel, filename );
    }

    template class io_api GeoModelInputHandler< 2 >;
    template class io_api GeoModelOutputHandler< 2 >;
    template class io_api GeoModelOutputHandler< 3 >;
    template class io_api GeoModelInputHandler< 3 >;

} // namespace RINGMesh
