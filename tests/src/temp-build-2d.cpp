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

#include <ringmesh/basic/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_api.h>
#include <ringmesh/geomodel/geomodel_builder.h>
#include <ringmesh/io/io.h>

int main()
{
    using namespace RINGMesh;

    try {

        // This line stands for the initialization
        // of Geogram and the factories of RINGMesh
        // IT IS MANDATORY
        default_configure();

        // Say Hello
        print_header_information();

        // We instantiate the class GeoModel
        GeoModel< 2 > geomodel;

        // To build the model, we have to use the class
        // GeoModelBuilder, which is a safety for the user
        // Indeed, the GeoModel can't be directly modified,
        // It has to be done using the GeoModelBuilder
        GeoModelBuilder< 2 > builder( geomodel );

        //#############################
        // Declaration of the Entities#
        //#############################

        index_t nb_interfaces = 7;
        index_t nb_layers = 3;

        index_t nb_corners = 12;
        index_t nb_lines = 17;
        index_t nb_surfaces = 6;

        // Create the Interfaces
        for( index_t interface_itr = 0; interface_itr < nb_interfaces;
            interface_itr++ ) {
            builder.geology.create_geological_entity(
                Interface< 2 >::type_name_static() );
        }

        // Create the Layers
        for( index_t layer = 0; layer < nb_layers; layer++ ) {
            builder.geology.create_geological_entity(
                Layer< 2 >::type_name_static() );
        }

        // Then we create the GeoModelMeshEntity
        // Create the Corners
        for( index_t corner = 0; corner < nb_corners; corner++ ) {
            builder.topology.create_mesh_entity< Corner >();
        }

        // Create the Lines
        for( index_t lines = 0; lines < nb_lines; lines++ ) {
            builder.topology.create_mesh_entity< Line >();
        }

        // Create the Surfaces
        for( index_t surface = 0; surface < nb_surfaces; surface++ ) {
            builder.topology.create_mesh_entity< Surface >();

        }

        //#############################
        // Setting the Geometry       #
        //#############################

        // We declare the coordinates of the corners. We arrange the corner in a
        // table
        std::vector< vec2 > corners_table( 12 );
        corners_table[0] = vec2( 0, 0 );
        corners_table[1] = vec2( 5, 0 );
        corners_table[2] = vec2( 10, 0 );
        corners_table[3] = vec2( 0, 3 );
        corners_table[4] = vec2( 5, 3 );
        corners_table[5] = vec2( 10, 3 );
        corners_table[6] = vec2( 0, 6 );
        corners_table[7] = vec2( 5, 6 );
        corners_table[8] = vec2( 10, 6 );
        corners_table[9] = vec2( 0, 10 );
        corners_table[10] = vec2( 5, 10 );
        corners_table[11] = vec2( 10, 10 );

        // We associate the coordinates with the corners
        for( index_t corner = 0; corner < nb_corners; corner++ ) {
            builder.geometry.set_corner( corner, corners_table[corner] );
        }

        // We associate the coordinates with the lines
        // We create a vector cur_coor_line containing the 2 vertices
        // for each line. Of course, you can have more vertices in a Line
        std::vector< vec2 > cur_coor_line( 2 );
        cur_coor_line[0] = corners_table[0];
        cur_coor_line[1] = corners_table[1];
        builder.geometry.set_line( 0, cur_coor_line );

        cur_coor_line[0] = corners_table[1];
        cur_coor_line[1] = corners_table[2];
        builder.geometry.set_line( 1, cur_coor_line );

        cur_coor_line[0] = corners_table[0];
        cur_coor_line[1] = corners_table[3];
        builder.geometry.set_line( 2, cur_coor_line );

        cur_coor_line[0] = corners_table[1];
        cur_coor_line[1] = corners_table[4];
        builder.geometry.set_line( 3, cur_coor_line );

        cur_coor_line[0] = corners_table[2];
        cur_coor_line[1] = corners_table[5];
        builder.geometry.set_line( 4, cur_coor_line );

        cur_coor_line[0] = corners_table[3];
        cur_coor_line[1] = corners_table[4];
        builder.geometry.set_line( 5, cur_coor_line );

        cur_coor_line[0] = corners_table[4];
        cur_coor_line[1] = corners_table[5];
        builder.geometry.set_line( 6, cur_coor_line );

        cur_coor_line[0] = corners_table[3];
        cur_coor_line[1] = corners_table[6];
        builder.geometry.set_line( 7, cur_coor_line );

        cur_coor_line[0] = corners_table[4];
        cur_coor_line[1] = corners_table[7];
        builder.geometry.set_line( 8, cur_coor_line );

        cur_coor_line[0] = corners_table[5];
        cur_coor_line[1] = corners_table[8];
        builder.geometry.set_line( 9, cur_coor_line );

        cur_coor_line[0] = corners_table[6];
        cur_coor_line[1] = corners_table[7];
        builder.geometry.set_line( 10, cur_coor_line );

        cur_coor_line[0] = corners_table[7];
        cur_coor_line[1] = corners_table[8];
        builder.geometry.set_line( 11, cur_coor_line );

        cur_coor_line[0] = corners_table[6];
        cur_coor_line[1] = corners_table[9];
        builder.geometry.set_line( 12, cur_coor_line );

        cur_coor_line[0] = corners_table[7];
        cur_coor_line[1] = corners_table[10];
        builder.geometry.set_line( 13, cur_coor_line );

        cur_coor_line[0] = corners_table[8];
        cur_coor_line[1] = corners_table[11];
        builder.geometry.set_line( 14, cur_coor_line );

        cur_coor_line[0] = corners_table[9];
        cur_coor_line[1] = corners_table[10];
        builder.geometry.set_line( 15, cur_coor_line );

        cur_coor_line[0] = corners_table[10];
        cur_coor_line[1] = corners_table[11];
        builder.geometry.set_line( 16, cur_coor_line );

        // We associate the coordinates with the Surfaces
        // We create a vector cur_coor_surface containing 4 vertices.
        // These 4 vertices delimits each surface so each surface
        // will contain one unique quad as a facet.
        // You can defined a more complicated mesh (for example a
        // triangular mesh) with these methods.
        std::vector< index_t > facet( 6 );
        facet[0] = 0;
        facet[1] = 1;
        facet[2] = 2;
        facet[3] = 0;
        facet[4] = 2;
        facet[5] = 3;
        std::vector< index_t > facet_ptr( 3 );
        facet_ptr[0] = 0;
        facet_ptr[1] = 2;
        facet_ptr[2] = 4;
        std::vector< vec2 > cur_coor_surface( 4 );
        cur_coor_surface[0] = corners_table[0];
        cur_coor_surface[1] = corners_table[1];
        cur_coor_surface[2] = corners_table[4];
        cur_coor_surface[3] = corners_table[3];
        builder.geometry.set_surface_geometry( 0, cur_coor_surface, facet,
            facet_ptr );

        cur_coor_surface[0] = corners_table[1];
        cur_coor_surface[1] = corners_table[2];
        cur_coor_surface[2] = corners_table[5];
        cur_coor_surface[3] = corners_table[4];
        builder.geometry.set_surface_geometry( 1, cur_coor_surface, facet,
            facet_ptr );

        cur_coor_surface[0] = corners_table[3];
        cur_coor_surface[1] = corners_table[4];
        cur_coor_surface[2] = corners_table[7];
        cur_coor_surface[3] = corners_table[6];
        builder.geometry.set_surface_geometry( 2, cur_coor_surface, facet,
            facet_ptr );

        cur_coor_surface[0] = corners_table[4];
        cur_coor_surface[1] = corners_table[5];
        cur_coor_surface[2] = corners_table[8];
        cur_coor_surface[3] = corners_table[7];
        builder.geometry.set_surface_geometry( 3, cur_coor_surface, facet,
            facet_ptr );

        cur_coor_surface[0] = corners_table[6];
        cur_coor_surface[1] = corners_table[7];
        cur_coor_surface[2] = corners_table[10];
        cur_coor_surface[3] = corners_table[9];
        builder.geometry.set_surface_geometry( 4, cur_coor_surface, facet,
            facet_ptr );

        cur_coor_surface[0] = corners_table[7];
        cur_coor_surface[1] = corners_table[8];
        cur_coor_surface[2] = corners_table[11];
        cur_coor_surface[3] = corners_table[10];
        builder.geometry.set_surface_geometry( 5, cur_coor_surface, facet,
            facet_ptr );

        //###################################
        // Setting the Boundaries relations #
        //###################################

        //We set the Corners which are incident entities of the lines
        // The add_mesh_entity_boundary_relation method take as first argument the
        // gme_t of the boundary and in second argument
        // the id of the GeoModelMeshentity bounded by the boundary
        // Remember :
        // Lines are bounded by Corners
        // Surfaces are bounded by Lines
        // Region are bounded by Surfaces

        // For corner 0
        gmme_id corner0( Corner< 2 >::type_name_static(), 0 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 0 ), corner0 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 2 ), corner0 );

        // For corner 1
        gmme_id corner1( Corner< 2 >::type_name_static(), 1 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 0 ), corner1 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 3 ), corner1 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 1 ), corner1 );

        // For corner 2
        gmme_id corner2( Corner< 2 >::type_name_static(), 2 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 1 ), corner2 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 4 ), corner2 );

        // For corner 3
        gmme_id corner3( Corner< 2 >::type_name_static(), 3 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 2 ), corner3 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 5 ), corner3 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 7 ), corner3 );

        // For corner 4
        gmme_id corner4( Corner< 2 >::type_name_static(), 4 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 3 ), corner4 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 5 ), corner4 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 6 ), corner4 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 8 ), corner4 );

        // For corner 5
        gmme_id corner5( Corner< 2 >::type_name_static(), 5 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 4 ), corner5 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 6 ), corner5 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 9 ), corner5 );

        // For corner 6
        gmme_id corner6( Corner< 2 >::type_name_static(), 6 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 10 ), corner6 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 12 ), corner6 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 7 ), corner6 );

        // For corner 7
        gmme_id corner7( Corner< 2 >::type_name_static(), 7 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 8 ), corner7 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 10 ), corner7 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 11 ), corner7 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 13 ), corner7 );

        // For corner 8
        gmme_id corner8( Corner< 2 >::type_name_static(), 8 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 11 ), corner8 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 14 ), corner8 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 9 ), corner8 );

        // For corner 9
        gmme_id corner9( Corner< 2 >::type_name_static(), 9 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 12 ), corner9 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 15 ), corner9 );

        // For corner 10
        gmme_id corner10( Corner< 2 >::type_name_static(), 10 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 16 ), corner10 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 13 ), corner10 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 15 ), corner10 );

        // For corner 11
        gmme_id corner11( Corner< 2 >::type_name_static(), 11 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 14 ), corner11 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Line< 2 >::type_name_static(), 16 ), corner11 );

        /////////////////////////////////////////////////////////

        // For line 0
        gmme_id line0( Line< 2 >::type_name_static(), 0 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 0 ), line0, false );

        // For line 1
        gmme_id line1( Line< 2 >::type_name_static(), 1 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 1 ), line1, false );

        // For line 2
        gmme_id line2( Line< 2 >::type_name_static(), 2 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 0 ), line2, true );

        // For line 3
        gmme_id line3( Line< 2 >::type_name_static(), 3 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 0 ), line3, false );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 1 ), line3, true );

        // For line 4
        gmme_id line4( Line< 2 >::type_name_static(), 4 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 1 ), line4, false );

        // For line 5
        gmme_id line5( Line< 2 >::type_name_static(), 5 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 0 ), line5, true );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 2 ), line5, false );

        // For line 6
        gmme_id line6( Line< 2 >::type_name_static(), 6 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 1 ), line6, true );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 3 ), line6, false );

        // For line 7
        gmme_id line7( Line< 2 >::type_name_static(), 7 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 2 ), line7, true );

        // For line 8
        gmme_id line8( Line< 2 >::type_name_static(), 8 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 2 ), line8, false );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 3 ), line8, true );

        // For line 9
        gmme_id line9( Line< 2 >::type_name_static(), 9 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 3 ), line9, false );

        // For line 10
        gmme_id line10( Line< 2 >::type_name_static(), 10 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 2 ), line10, true );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 4 ), line10, false );

        // For line 11
        gmme_id line11( Line< 2 >::type_name_static(), 11 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 3 ), line11, true );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 5 ), line11, false );

        // For line 12
        gmme_id line12( Line< 2 >::type_name_static(), 12 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 4 ), line12, true );

        // For line 13
        gmme_id line13( Line< 2 >::type_name_static(), 13 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 5 ), line13, false );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 4 ), line13, true );

        // For line 14
        gmme_id line14( Line< 2 >::type_name_static(), 14 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 5 ), line14, false );

        // For line 15
        gmme_id line15( Line< 2 >::type_name_static(), 15 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 4 ), line15, true );

        // For line 16
        gmme_id line16( Line< 2 >::type_name_static(), 16 );
        builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Surface< 2 >::type_name_static(), 5 ), line16, true );

        // For the Universe Boundary
        builder.topology.add_universe_boundary( 0, true );
        builder.topology.add_universe_boundary( 1, true );
        builder.topology.add_universe_boundary( 4, true );
        builder.topology.add_universe_boundary( 9, true );
        builder.topology.add_universe_boundary( 14, true );
        builder.topology.add_universe_boundary( 16, false );
        builder.topology.add_universe_boundary( 15, false );
        builder.topology.add_universe_boundary( 12, false );
        builder.topology.add_universe_boundary( 7, false );
        builder.topology.add_universe_boundary( 2, false );

        //#####################################
        // Setting the parent/child relations #
        //#####################################

        // For Interface 0
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 0 ),
            gmme_id( Line< 2 >::type_name_static(), 0 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 0 ),
            gmme_id( Line< 2 >::type_name_static(), 1 ) );

        // For Interface 1
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 1 ),
            gmme_id( Line< 2 >::type_name_static(), 5 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 1 ),
            gmme_id( Line< 2 >::type_name_static(), 6 ) );

        // For Interface 2
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 2 ),
            gmme_id( Line< 2 >::type_name_static(), 10 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 2 ),
            gmme_id( Line< 2 >::type_name_static(), 11 ) );

        // For Interface 3
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 3 ),
            gmme_id( Line< 2 >::type_name_static(), 15 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 3 ),
            gmme_id( Line< 2 >::type_name_static(), 16 ) );

        // For Interface 4
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 4 ),
            gmme_id( Line< 2 >::type_name_static(), 2 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 4 ),
            gmme_id( Line< 2 >::type_name_static(), 7 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 4 ),
            gmme_id( Line< 2 >::type_name_static(), 12 ) );

        // For Interface 5
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 5 ),
            gmme_id( Line< 2 >::type_name_static(), 13 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 5 ),
            gmme_id( Line< 2 >::type_name_static(), 8 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 5 ),
            gmme_id( Line< 2 >::type_name_static(), 3 ) );

        // For Interface 6
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 6 ),
            gmme_id( Line< 2 >::type_name_static(), 14 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 6 ),
            gmme_id( Line< 2 >::type_name_static(), 9 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Interface< 2 >::type_name_static(), 6 ),
            gmme_id( Line< 2 >::type_name_static(), 4 ) );

        ///////////////////////////////////////////////////

        // For Layer 0
        builder.geology.add_parent_children_relation(
            gmge_id( Layer< 2 >::type_name_static(), 0 ),
            gmme_id( Surface< 2 >::type_name_static(), 0 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Layer< 2 >::type_name_static(), 0 ),
            gmme_id( Surface< 2 >::type_name_static(), 1 ) );

        // For Layer 1
        builder.geology.add_parent_children_relation(
            gmge_id( Layer< 2 >::type_name_static(), 1 ),
            gmme_id( Surface< 2 >::type_name_static(), 2 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Layer< 2 >::type_name_static(), 1 ),
            gmme_id( Surface< 2 >::type_name_static(), 3 ) );

        // For Layer 2
        builder.geology.add_parent_children_relation(
            gmge_id( Layer< 2 >::type_name_static(), 2 ),
            gmme_id( Surface< 2 >::type_name_static(), 4 ) );
        builder.geology.add_parent_children_relation(
            gmge_id( Layer< 2 >::type_name_static(), 2 ),
            gmme_id( Surface< 2 >::type_name_static(), 5 ) );

        // Then, we end the model building
        // This method will set the missing information for the boundaries
        // and parent/child relation. e. g., if you decide to use the
        // add_parent_children_relation (like above), the child has no information of who
        // is his parent. This method deal with that by filling the missing information
        builder.end_geomodel();

        // We save the builded model
//        geomodel_save( geomodel, "builded_model.gm" );

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category(), e.what() );
        return 1;
    } catch( const std::exception& e ) {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    return 0;
}
