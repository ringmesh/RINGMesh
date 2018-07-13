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

#pragma once

#include <ringmesh/io/common.h>

namespace {
    /*!
     * In VTU/VTK, edges, polygons and cells (in the RINGMesh vocbulary)
     * are considered as... cells
     */
    index_t total_number_of_elements( const GeoModel3D& geomodel) {
        
        return geomodel.nb_corners() +
            geomodel.mesh.edges.nb() +
            geomodel.mesh.polygons.nb() +
            geomodel.mesh.cells.nb();
    }
    struct RINGMesh2VTK
    {
        index_t entity_type;
        index_t vertices[8];
    };

    static RINGMesh2VTK vertex_descriptor_vtk = {
        1, // type
        { 0 } // vertices
    };

    static RINGMesh2VTK edge_descriptor_vtk = {
        3, // type
        { 0,1  } // vertices
    };

    static RINGMesh2VTK triangle_descriptor_vtk = {
        5, // type
        { 0, 1, 2 } // vertices
    };

    static RINGMesh2VTK quad_descriptor_vtk = {
        9, // type
        { 0, 1, 2, 3 } // vertices
    };

    static RINGMesh2VTK hex_descriptor_vtk = {
        12, // type
        { 0, 4, 5, 1, 2, 6, 7, 3 } // vertices
    };

    static RINGMesh2VTK prism_descriptor_vtk = {
        13, // type
        { 0, 2, 1, 3, 5, 4 } // vertices
    };

    static RINGMesh2VTK pyramid_descriptor_vtk = {
        14, // type
        { 0, 1, 2, 3, 4 } // vertices
    };

    static RINGMesh2VTK tet_descriptor_vtk = {
        10, // type
        { 0, 1, 2, 3 } // vertices
    };

    static RINGMesh2VTK no_cell = {
        0, // type
        { } // vertices
    };

    static RINGMesh2VTK* ringmesh2vtk_vertex[2] = {
        &no_cell, &vertex_descriptor_vtk
    };

    static RINGMesh2VTK* ringmesh2vtk_edges[3] = {
        &no_cell, &no_cell, &edge_descriptor_vtk
    };

    static RINGMesh2VTK* ringmesh2vtk_polygons[5] = {
        &no_cell, &no_cell, &no_cell, &triangle_descriptor_vtk, &quad_descriptor_vtk
    };

    static RINGMesh2VTK* ringmesh2vtk_cells[9] = {
        &no_cell, &no_cell, &no_cell, &no_cell,
        &tet_descriptor_vtk, &pyramid_descriptor_vtk, &prism_descriptor_vtk, &no_cell,
        &hex_descriptor_vtk
    };

} // namespace RINGMesh
