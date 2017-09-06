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

#include <ringmesh/io/io.h>

#include <cctype>
#include <iomanip>
#include <deque>

#include <tinyxml2/tinyxml2.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_api.h>
#include <ringmesh/geomodel/geomodel_builder.h>
#include <ringmesh/geomodel/geomodel_builder_gocad.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_validity.h>

#include <ringmesh/mesh/well.h>
#include <ringmesh/mesh/geogram_mesh.h>

#include <ringmesh/io/zip_file.h>

#include <minizip/unzip.h>
#include <minizip/zip.h>

/*!
 * @file Implementation of classes loading GeoModels
 * @author Arnaud Botella and Antoine Mazuyer
 */

namespace {
    using namespace RINGMesh;

#include "geomodel/io_abaqus.cpp"
#include "geomodel/io_adeli.cpp"
#include "geomodel/io_aster.cpp"
#include "geomodel/io_csmp.cpp"
#include "geomodel/io_feflow.cpp"
#include "geomodel/io_gm.cpp"
#include "geomodel/io_gprs.cpp"
#include "geomodel/io_mfem.cpp"
#include "geomodel/io_model3d.cpp"
#include "geomodel/io_msh.cpp"
#include "geomodel/io_smesh.cpp"
#include "geomodel/io_stl.cpp"
#include "geomodel/io_svg.cpp"
#include "geomodel/io_tetgen.cpp"
#include "geomodel/io_tsolid.cpp"
#include "geomodel/io_vtk.cpp"

#ifdef RINGMESH_WITH_GEOLOGYJS
#    include "boundary_geomodel/io_html.cpp"
#endif

}

namespace RINGMesh {

    template< >
    void GeoModelIOHandler< 2 >::initialize()
    {
        GeoModelIOHandlerFactory2D::register_creator< GeoModelHandlerGM2D >( "gm" );
        GeoModelIOHandlerFactory2D::register_creator< SVGIOHandler >( "svg" );
        GeoModelIOHandlerFactory2D::register_creator< MFEMIOHandler2D >( "mfem" );
    }

    /*
     * Initializes the possible handler for IO files
     */
    template< >
    void GeoModelIOHandler< 3 >::initialize()
    {
        GeoModelIOHandlerFactory3D::register_creator< TetGenIOHandler >( "tetgen" );
        GeoModelIOHandlerFactory3D::register_creator< TSolidIOHandler >( "so" );
        GeoModelIOHandlerFactory3D::register_creator< CSMPIOHandler >( "csmp" );
        GeoModelIOHandlerFactory3D::register_creator< AsterIOHandler >( "mail" );
        GeoModelIOHandlerFactory3D::register_creator< VTKIOHandler >( "vtk" );
        GeoModelIOHandlerFactory3D::register_creator< GPRSIOHandler >( "gprs" );
        GeoModelIOHandlerFactory3D::register_creator< MSHIOHandler >( "msh" );
        GeoModelIOHandlerFactory3D::register_creator< MFEMIOHandler3D >( "mfem" );
        GeoModelIOHandlerFactory3D::register_creator< GeoModelHandlerGM3D >( "gm" );
        GeoModelIOHandlerFactory3D::register_creator< AbaqusIOHandler >( "inp" );
        GeoModelIOHandlerFactory3D::register_creator< AdeliIOHandler >( "adeli" );
        GeoModelIOHandlerFactory3D::register_creator< FeflowIOHandler >( "fem" );
        GeoModelIOHandlerFactory3D::register_creator< MLIOHandler >( "ml" );
        GeoModelIOHandlerFactory3D::register_creator< SMESHIOHandler >( "smesh" );
        GeoModelIOHandlerFactory3D::register_creator< STLIOHandler >( "stl" );

#ifdef RINGMESH_WITH_GEOLOGYJS
        GeoModelIOHandlerFactory3D::register_creator< HTMLIOHandler >( "html" );
#endif
    }

} // namespace RINGMesh
