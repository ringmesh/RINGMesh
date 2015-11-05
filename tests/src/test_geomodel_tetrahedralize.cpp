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
#include <ringmesh/io.h>
#include <ringmesh/ringmesh_tests_config.h>

/*!
 * @file Test global tetrahedralization of a GeoModel
 * @author Jeanne Pellerin
 */

using GEO::Logger ;
using RINGMesh::index_t ;

bool mesh_constrained_tetgen( GEO::Mesh& M, bool refine, double quality )
{
    if( !M.facets.are_simplices() ) {
        Logger::err( "TetMeshing" )
            << "Mesh is not triangulated"
            << std::endl;
        return false;
    }

    GEO::Delaunay_var delaunay = GEO::Delaunay::create( 3, "tetgen" );
    delaunay->set_refine( refine );
    delaunay->set_quality( quality );
    delaunay->set_constraints( &M );
    // Compute the tetrahedrons 
    delaunay->set_vertices( 0, nil ); 

    GEO::vector<double> pts( delaunay->nb_vertices() * 3 );
    GEO::vector<index_t> tet2v( delaunay->nb_cells() * 4 );
    for( index_t v = 0; v < delaunay->nb_vertices(); ++v ) {
        pts[ 3 * v ] = delaunay->vertex_ptr( v )[ 0 ];
        pts[ 3 * v + 1 ] = delaunay->vertex_ptr( v )[ 1 ];
        pts[ 3 * v + 2 ] = delaunay->vertex_ptr( v )[ 2 ];
    }
    for( index_t t = 0; t < delaunay->nb_cells(); ++t ) {
        tet2v[ 4 * t ] = index_t( delaunay->cell_vertex( t, 0 ) );
        tet2v[ 4 * t + 1 ] = index_t( delaunay->cell_vertex( t, 1 ) );
        tet2v[ 4 * t + 2 ] = index_t( delaunay->cell_vertex( t, 2 ) );
        tet2v[ 4 * t + 3 ] = index_t( delaunay->cell_vertex( t, 3 ) );
    }

    M.cells.assign_tet_mesh( 3, pts, tet2v, true );
    M.cells.connect();
    M.show_stats( "TetMeshing" );
    return true;
}

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    // Set an output log file
    std::string log_file( ringmesh_test_output_path ) ;
    log_file += "log.txt" ;
    GEO::FileLogger* file_logger = new GEO::FileLogger( log_file ) ;
    GEO::Logger::instance()->register_client( file_logger ) ;


    std::string file_name( ringmesh_test_data_path ) ;
    file_name += "split_cube.ml" ;
    std::string result_file_name( ringmesh_test_output_path ) ;
    result_file_name += "geomodel_tet_mesh.mesh" ;


    // 0. Load a Geomodel
    GeoModel geomodel ;
    geomodel_surface_load( file_name, geomodel ) ;


    // 1. Get the global Mesh of Surfaces of the GeoModel
    GEO::Mesh geomodel_surfaces_mesh ;
    geomodel.mesh.copy_mesh( geomodel_surfaces_mesh ) ;

    // 2. Tetrahedralize this Mesh
    mesh_constrained_tetgen( geomodel_surfaces_mesh, false, 2.0 ) ;

    GEO::mesh_save( geomodel_surfaces_mesh, result_file_name ) ;


    // 3. Assign pieces of the generated Mesh to the Regions 
    // of the model

    

}