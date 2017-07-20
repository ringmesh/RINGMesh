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
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
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

#include <ringmesh/io/io.h>

#include <cctype>
#include <iomanip>
#include <deque>

#include <tinyxml2/tinyxml2.h>

#include <geogram/basic/command_line.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_api.h>
#include <ringmesh/geomodel/geomodel_builder_gocad.h>
#include <ringmesh/geomodel/geomodel_builder_ringmesh.h>
#include <ringmesh/geomodel/geomodel_validity.h>

#include <ringmesh/mesh/well.h>
#include <ringmesh/mesh/geogram_mesh.h>

/*!
 * @file Implementation of classes loading volumetric GeoModels
 * @author Arnaud Botella and Antoine Mazuyer
 */

namespace {
    using namespace RINGMesh;

#include "full_geomodel/io_abaqus.cpp"
#include "full_geomodel/io_adeli.cpp"
#include "full_geomodel/io_aster.cpp"
#include "full_geomodel/io_csmp.cpp"
#include "full_geomodel/io_feflow.cpp"
#include "full_geomodel/io_geomodel.cpp"
#include "full_geomodel/io_gprs.cpp"
#include "full_geomodel/io_mfem.cpp"
#include "full_geomodel/io_msh.cpp"
#include "full_geomodel/io_svg.cpp"
#include "full_geomodel/io_tetgen.cpp"
#include "full_geomodel/io_tsolid.cpp"
#include "full_geomodel/io_vtk.cpp"

}
namespace RINGMesh {

    template< >
    void GeoModelIOHandler< 2 >::initialize_full_geomodel_output()
    {
        ringmesh_register_GeoModelIOHandler2D_creator( GeoModelHandlerGM< 2 >, "gm" );
        ringmesh_register_GeoModelIOHandler2D_creator( SVGIOHandler, "svg" );
    }

    /*
     * Initializes the possible handler for IO files
     */
    template< >
    void GeoModelIOHandler< 3 >::initialize_full_geomodel_output()
    {
        ringmesh_register_GeoModelIOHandler3D_creator( TetGenIOHandler, "tetgen" );
		ringmesh_register_GeoModelIOHandler3D_creator( TSolidIOHandler, "so" );
		ringmesh_register_GeoModelIOHandler3D_creator( TSolidIOHandler, "lso" );
        ringmesh_register_GeoModelIOHandler3D_creator( CSMPIOHandler, "csmp" );
        ringmesh_register_GeoModelIOHandler3D_creator( AsterIOHandler, "mail" );
        ringmesh_register_GeoModelIOHandler3D_creator( VTKIOHandler, "vtk" );
        ringmesh_register_GeoModelIOHandler3D_creator( GPRSIOHandler, "gprs" );
        ringmesh_register_GeoModelIOHandler3D_creator( MSHIOHandler, "msh" );
        ringmesh_register_GeoModelIOHandler3D_creator( MFEMIOHandler, "mfem" );
        ringmesh_register_GeoModelIOHandler3D_creator( GeoModelHandlerGM< 3 >, "gm" );
        ringmesh_register_GeoModelIOHandler3D_creator( AbaqusIOHandler, "inp" );
        ringmesh_register_GeoModelIOHandler3D_creator( AdeliIOHandler, "adeli" );
        ringmesh_register_GeoModelIOHandler3D_creator( FeflowIOHandler, "fem" );
    }

}
