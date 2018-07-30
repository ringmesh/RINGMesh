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
#include <ringmesh/basic/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>

#include <ringmesh/ringmesh_tests_config.h>

#include <vector>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <ringmesh/io/io.h>

using namespace RINGMesh;

// For the next section, read the documentation to understand
// the concept of Geological Entity and Mesh Entities
// Let's to a sum up of the GeoModel we try to build:
// For the Geological Entities (handle by the class
// GeoModelGeologicalEntity):
// 9 Contacts
const index_t nb_contacts = 9;
// 1 horizons + 6 boundaries = 7 Interfaces
const index_t nb_interfaces = 7;
// 2 Layers
const index_t nb_layers = 2;

// For the Meshed Entities, (handle by the class GeoModelMeshEntity)
// 5 Corners
const index_t nb_corners = 5;
// 9 Lines
const index_t nb_lines = 9;
// 7 Surfaces
const index_t nb_surfaces = 7;
// stergd Regions
const index_t nb_regions = 2;

bool test_build_geomodel( const GeoModel3D& geomodel )
{
    const index_t nb_surface_entities =
        geomodel.nb_mesh_entities( surface_type_name_static() );
    const index_t nb_region_entities =
        geomodel.nb_mesh_entities( region_type_name_static() );
    const index_t nb_corner_entities =
        geomodel.nb_mesh_entities( corner_type_name_static() );
    const index_t nb_line_entities =
        geomodel.nb_mesh_entities( line_type_name_static() );

    const index_t nb_contact_entities =
        geomodel.nb_geological_entities( Contact3D::type_name_static() );
    const index_t nb_interface_entities =
        geomodel.nb_geological_entities( Interface3D::type_name_static() );
    const index_t nb_layer_entities =
        geomodel.nb_geological_entities( Layer3D::type_name_static() );

    if( nb_contacts != nb_contact_entities )
    {
        throw RINGMeshException( "RINGMesh Test",
            "The numbers of contacts are not equal " + geomodel.name()
                + ": the loaded model is not valid." );
    }

    if( nb_interfaces != nb_interface_entities )
    {
        throw RINGMeshException( "RINGMesh Test",
            "The numbers of interfaces are not equal " + geomodel.name()
                + ": the loaded model is not valid." );
    }

    if( nb_layers != nb_layer_entities )
    {
        throw RINGMeshException( "RINGMesh Test",
            "The numbers of layers are not equal " + geomodel.name()
                + ": the loaded model is not valid." );
    }

    if( nb_surfaces != nb_surface_entities )
    {
        throw RINGMeshException( "RINGMesh Test",
            "The numbers of surfaces are not equal " + geomodel.name()
                + ": the loaded model is not valid." );
    }

    if( nb_regions != nb_region_entities )
    {
        throw RINGMeshException( "RINGMesh Test",
            "The numbers of regions are not equal " + geomodel.name()
                + ": the loaded model is not valid." );
    }

    if( nb_lines != nb_line_entities )
    {
        throw RINGMeshException( "RINGMesh Test",
            "The numbers of lines are not equal " + geomodel.name()
                + ": the loaded model is not valid." );
    }

    if( nb_corners != nb_corner_entities )
    {
        throw RINGMeshException( "RINGMesh Test",
            "The numbers of corners are not equal " + geomodel.name()
                + ": the loaded model is not valid." );
    }
    return true;
}
bool build_geomodel( GeoModel3D& geomodel )
{
    GeoModelBuilder3D builder( geomodel );

    // We first create the GeoModelGeoglogicalEntity
    // Create the contacts
    for( index_t contact = 0; contact < nb_contacts; contact++ )
    {
        builder.geology.create_geological_entity(
            Contact3D::type_name_static() );
        // the static method type_name_static() is available for each
        // GeoModelEntity. It returns an EntityType which is a string
        // corresponding to the Type of the entity.
    }

    // Create the Interfaces
    for( index_t interface_itr = 0; interface_itr < nb_interfaces;
         interface_itr++ )
    {
        builder.geology.create_geological_entity(
            Interface3D::type_name_static() );
    }

    // Create the Layers
    for( index_t layer = 0; layer < nb_layers; layer++ )
    {
        builder.geology.create_geological_entity( Layer3D::type_name_static() );
    }

    // Then we create the GeoModelMeshEntity
    // Create the Corners
    for( index_t corner = 0; corner < nb_corners; corner++ )
    {
        builder.topology.create_mesh_entity( Corner3D::type_name_static() );
    }

    // Create the Lines
    for( index_t lines = 0; lines < nb_lines; lines++ )
    {
        builder.topology.create_mesh_entity( Line3D::type_name_static() );
    }

    // Create the Surfaces
    for( index_t surface = 0; surface < nb_surfaces; surface++ )
    {
        builder.topology.create_mesh_entity( Surface3D::type_name_static() );
    }

    // Create the Regions
    for( index_t region = 0; region < nb_regions; region++ )
    {
        builder.topology.create_mesh_entity( Region3D::type_name_static() );
    }

    //#############################
    // Setting the Geometry       #
    //#############################

    // We declare the coordinates of the corners. We arrange the corner in a
    // table
    vec3 corners_table[nb_corners];
    corners_table[0] = vec3( 0, 0, 0 );
    corners_table[1] = vec3( 4, 0, 0 );
    corners_table[2] = vec3( 2, 2, 0 );
    corners_table[3] = vec3( 2, 1, 4 );
    corners_table[4] = vec3( 2, 1, -4 );

    // We associate the coordinates with the corners
    for( index_t corner = 0; corner < nb_corners; ++corner )
    {
        builder.geometry.set_corner( corner, corners_table[corner] );
    }

    // We associate the coordinates with the lines
    // We create a vector cur_coor_line containing the 2 vertices
    // for each line. Of course, you can have more vertices in a Line
    std::vector< vec3 > cur_coor_line( 2 );
    cur_coor_line[0] = corners_table[0];
    cur_coor_line[1] = corners_table[3];
    builder.geometry.set_line( 0, cur_coor_line );

    const GeoModelMeshEntity3D& line =
        geomodel.mesh_entity( line_type_name_static(), 0 );

    if( line.vertex( 0 ) != cur_coor_line[0]
        || line.vertex( 1 ) != cur_coor_line[1] )
    {
        throw RINGMeshException(
            "RINGMesh Test", "The line's corners are not equal" );
    }

    cur_coor_line[0] = corners_table[2];
    cur_coor_line[1] = corners_table[3];
    builder.geometry.set_line( 1, cur_coor_line );

    cur_coor_line[0] = corners_table[1];
    cur_coor_line[1] = corners_table[2];
    builder.geometry.set_line( 2, cur_coor_line );

    cur_coor_line[0] = corners_table[0];
    cur_coor_line[1] = corners_table[1];
    builder.geometry.set_line( 3, cur_coor_line );

    cur_coor_line[0] = corners_table[1];
    cur_coor_line[1] = corners_table[4];
    builder.geometry.set_line( 4, cur_coor_line );

    cur_coor_line[0] = corners_table[0];
    cur_coor_line[1] = corners_table[4];
    builder.geometry.set_line( 5, cur_coor_line );

    cur_coor_line[0] = corners_table[1];
    cur_coor_line[1] = corners_table[3];
    builder.geometry.set_line( 6, cur_coor_line );

    cur_coor_line[0] = corners_table[2];
    cur_coor_line[1] = corners_table[4];
    builder.geometry.set_line( 7, cur_coor_line );

    cur_coor_line[0] = corners_table[0];
    cur_coor_line[1] = corners_table[2];
    builder.geometry.set_line( 8, cur_coor_line );

    // We associate the coordinates with the Surfaces
    // We create a vector cur_coor_surface containing 4 vertices.
    // These 4 vertices delimits each surface so each surface
    // will contain one unique quad as a facet.
    // You can defined a more complicated mesh (for example a
    // triangular mesh) with these methods.
    std::vector< index_t > facet( 3, 0 );
    facet[0] = 0;
    facet[1] = 1;
    facet[2] = 2;
    std::vector< index_t > facet_ptr( 2 );
    facet_ptr[0] = 0;
    facet_ptr[1] = 3;
    std::vector< vec3 > cur_coor_surface( 3 );
    cur_coor_surface[0] = corners_table[0];
    cur_coor_surface[1] = corners_table[1];
    cur_coor_surface[2] = corners_table[2];
    builder.geometry.set_surface_geometry(
        0, cur_coor_surface, facet, facet_ptr );

    const GeoModelMeshEntity3D& surface =
        geomodel.mesh_entity( surface_type_name_static(), 0 );

    if( surface.vertex( 0 ) != cur_coor_surface[0]
        || surface.vertex( 1 ) != cur_coor_surface[1]
        || surface.vertex( 2 ) != cur_coor_surface[2] )
    {
        throw RINGMeshException(
            "RINGMesh Test", "The surface's corners are not equal" );
    }

    cur_coor_surface[0] = corners_table[3];
    cur_coor_surface[1] = corners_table[2];
    cur_coor_surface[2] = corners_table[1];
    builder.geometry.set_surface_geometry(
        1, cur_coor_surface, facet, facet_ptr );

    cur_coor_surface[0] = corners_table[3];
    cur_coor_surface[1] = corners_table[0];
    cur_coor_surface[2] = corners_table[1];
    builder.geometry.set_surface_geometry(
        2, cur_coor_surface, facet, facet_ptr );

    cur_coor_surface[0] = corners_table[3];
    cur_coor_surface[1] = corners_table[2];
    cur_coor_surface[2] = corners_table[0];
    builder.geometry.set_surface_geometry(
        3, cur_coor_surface, facet, facet_ptr );

    cur_coor_surface[0] = corners_table[0];
    cur_coor_surface[1] = corners_table[1];
    cur_coor_surface[2] = corners_table[4];
    builder.geometry.set_surface_geometry(
        4, cur_coor_surface, facet, facet_ptr );

    cur_coor_surface[0] = corners_table[0];
    cur_coor_surface[1] = corners_table[4];
    cur_coor_surface[2] = corners_table[2];
    builder.geometry.set_surface_geometry(
        5, cur_coor_surface, facet, facet_ptr );

    cur_coor_surface[0] = corners_table[1];
    cur_coor_surface[1] = corners_table[4];
    cur_coor_surface[2] = corners_table[2];
    builder.geometry.set_surface_geometry(
        6, cur_coor_surface, facet, facet_ptr );

    return true;
}

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Test geomodel build" );

        GeoModel3D in;
        build_geomodel( in );
        test_build_geomodel( in );
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
