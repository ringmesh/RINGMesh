/* * Copyright (c)2018, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/ringmesh_tests_config.h>

#include <geogram/basic/line_stream.h>

#include <geogram/basic/command_line.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

/*!
 * @file Test GeoModel RESQML2 mixed elements in an UnstructuredGrid
 * @author Wan-Chiu Li
 */

using namespace RINGMesh;
using GMGE = GeoModelGeologicalEntity< 3 >;

void load_geomodel( GeoModel3D& geomodel, const std::string& file )
{
    // TODO
    GEO::CmdLine::set_arg( "validity:do_not_check", "A" );

    bool loaded_model_is_valid = geomodel_load( geomodel, file );
    if( !loaded_model_is_valid )
    {
        throw RINGMeshException(
            "RINGMesh Test", "Failed when loading model from file " + file
                                 + ": the loaded model is not valid." );
    }
}

/*
The unstructured grid in the following .epc is built with the code below
UnstructuredGridMixedElementsRepresentationTest.epc

const ULONG64 nodesCount = 10;
double pts[] = { 0, 0, 300, 700, 0, 350, 0, 150, 300, 700, 150, 350, 0, 0, 500,
    700, 0, 550, 0, 150, 500, 700, 150, 550, 50, 50, 800, 900, 80, 900 };

// getting the local depth 3d crs
LocalDepth3dCrsTest* crsTest = new LocalDepth3dCrsTest( this->epcDoc, true );
RESQML2_0_1_NS::LocalDepth3dCrs* crs =
    epcDoc->getResqmlAbstractObjectByUuid< RESQML2_0_1_NS::LocalDepth3dCrs >(
        LocalDepth3dCrsTest::defaultUuid );

// getting the hdf proxy
AbstractHdfProxy* hdfProxy = this->epcDoc->getHdfProxySet()[0];

// creating the ijk grid
RESQML2_0_1_NS::UnstructuredGridRepresentation* rep =
    this->epcDoc->createUnstructuredGridRepresentation( crs, uuid, title, 3 );

REQUIRE( rep != nullptr );

ULONG64 face_indices_per_cell[15] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
    13, 14 };
ULONG64 faceIndicesCumulativeCountPerCell[3] = { 6, 11, 15 };
ULONG64 nodeIndicesCumulativeCountPerFace[15] = { 4, 8, 12, 16, 20, 24, 28, 31,
    34, 37, 40, 43, 46, 49, 52 };
ULONG64 node_indices_per_face[52] = { 0, 1, 3, 2, 4, 5, 7, 6, 0, 1, 5, 4, 1, 3,
    7, 5, 3, 2, 6, 7, 0, 4, 6, 2, 4, 5, 7, 6, 8, 4, 5, 4, 8, 6, 6, 8, 7, 5, 7,
    8, 5, 8, 7, 9, 5, 7, 9, 7, 8, 9, 8, 5 };

unsigned char face_righthandness[15] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1 };

rep->setGeometry( face_righthandness,
    pts,
    nodesCount,
    this->epcDoc->getHdfProxySet()[0],
    face_indices_per_cell,
    faceIndicesCumulativeCountPerCell,
    15,
    node_indices_per_face,
    nodeIndicesCumulativeCountPerFace,
    gsoap_resqml2_0_1::resqml2__CellShape::resqml2__CellShape__polyhedral );
// cleaning
delete crsTest;

*/

void do_test()
{
    GeoModel3D origin_geomodel;
    load_geomodel( origin_geomodel,
        ringmesh_test_data_path
            + "UnstructuredGridMixedElementsRepresentationTest.epc" );

    if( origin_geomodel.nb_mesh_entities( Region3D::type_name_static() ) != 1 )
    {
        throw RINGMeshException( "TEST", "Only one region is expected" );
    }

    const Region3D& region = origin_geomodel.region( 0 );
    if( region.nb_mesh_elements() != 3 )
    {
        throw RINGMeshException( "TEST", "3 cells in the region are expected" );
    }
    if( region.cell_type( 0 ) != CellType::HEXAHEDRON )

    {
        throw RINGMeshException(
            "TEST", "First cell of type HEXAHEDRON is expected" );
    }
    if( region.cell_type( 1 ) != CellType::PYRAMID )

    {
        throw RINGMeshException(
            "TEST", "Second cell of type PYRAMID is expected" );
    }
    if( region.cell_type( 2 ) != CellType::TETRAHEDRON )

    {
        throw RINGMeshException(
            "TEST", "Third cell of type TETRAHEDRON is expected" );
    }
}

int main()
{
    try
    {
        Logger::out( "TEST", "GeoModel3D RESQML2 Round-trip" );
        do_test();
    }
    catch( const RINGMeshException& e )
    {
        Logger::err( e.category(), e.what() );
        return 1;
    }
    catch( const std::exception& e )
    {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
