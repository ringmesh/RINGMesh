/*
* Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
*  Contacts:
*     Arnaud.Botella@univ-lorraine.fr
*     Antoine.Mazuyer@univ-lorraine.fr
*     Jeanne.Pellerin@wias-berlin.de
*
*     http://www.ring-team.org
*
*     RING Project
*     Ecole Nationale Superieure de Geologie - GeoRessources
*     2 Rue du Doyen Marcel Roubault - TSA 70605
*     54518 VANDOEUVRE-LES-NANCY
*     FRANCE
*/


#include <geogram/basic/logger.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/io.h>
#include <ringmesh/ringmesh_tests_config.h>

/*!
 * @file Test global tetrahedralization of a GeoModel
 * @author Jeanne Pellerin
 */

using GEO::Logger ;;
using RINGMesh::index_t ;
using RINGMesh::GeoModel ;
using RINGMesh::Surface ;


int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    // Set an output log file
    std::string log_file( ringmesh_test_output_path ) ;
    log_file += "log.txt" ;
    GEO::FileLogger* file_logger = new GEO::FileLogger( log_file ) ;
    GEO::Logger::instance()->register_client( file_logger ) ;

    std::string file_name( ringmesh_test_data_path ) ;
    file_name += "modelA6.ml" ;
    std::string result_file_name( ringmesh_test_output_path ) ;
    result_file_name += "geomodel_tet_mesh.mesh" ;

    GeoModel geomodel ;
    geomodel_surface_load( file_name, geomodel ) ;
    tetgen_tetrahedralize_geomodel_regions( geomodel ) ;
  
    // Save volumetric mesh with cell region attribute
    /// @todo Implement a function to do it
    GEO::Mesh geomodel_regions_mesh ;
    build_mesh_from_geomodel( geomodel, geomodel_regions_mesh ) ;
    GEO::MeshIOFlags mesh_io_flags ;
    mesh_io_flags.set_elements( GEO::MeshElementsFlags( GEO::MESH_VERTICES | GEO::MESH_CELLS ) ) ;
    mesh_io_flags.set_attribute( GEO::MESH_CELL_REGION ) ;
    GEO::mesh_save( geomodel_regions_mesh, result_file_name, mesh_io_flags ) ;

    return 0 ;
}

