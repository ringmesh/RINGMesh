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

#include <ringmesh/ringmesh_tests_config.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/io.h>
#include <ringmesh/utils.h>

#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_io.h>

/*!
* @file Test GeoModel building from a mesh loaded from a .so file
* @author Pierre Anquez
*/

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    GEO::Logger::out( "TEST" ) << "Test IO for a mesh GeoModel in .so" << std::endl ;

    std::string file_name( ringmesh_test_data_path ) ;
    file_name += "modelA4.so" ;
//    file_name = "/home/anquez/anquez.so" ;
//    file_name = "/home/anquez/anquez_fault_f4a.so" ;

    GeoModel model ;
    if( !geomodel_volume_load( file_name, model ) ) {
        return 1 ;
    }

    GeoModel out_model ;
    if( !geomodel_surface_load( "imported_tsolid_surf.bm", out_model ) ) {
        return 1 ;
    }

    bool res = true ;
    if ( model.nb_corners() != 52 ||
         model.nb_lines() != 98 ||
         model.nb_surfaces() != 55 ||
         model.nb_regions() != 8 ||
         model.nb_interfaces() != 11 ||
         model.mesh.vertices.nb() != 6691 ||
         model.mesh.facets.nb() != 10049 ||
         model.mesh.cells.nb() != 34540 ||
         out_model.mesh.vertices.nb() != 4465 ||
         out_model.mesh.facets.nb() != 10049 ||
         out_model.mesh.cells.nb() != 0) {
        res = false ;
    }
    if( res ) {
        GEO::Logger::out( "TEST" ) << "SUCCESS" << std::endl ;
    } else {
        GEO::Logger::out( "TEST" ) << "FAILED" << std::endl ;
    }

    return !res ;
}
