/*
*  Copyright (c) 2012-2014, Bruno Levy
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*
*  * Redistributions of source code must retain the above copyright notice,
*  this list of conditions and the following disclaimer.
*  * Redistributions in binary form must reproduce the above copyright notice,
*  this list of conditions and the following disclaimer in the documentation
*  and/or other materials provided with the distribution.
*  * Neither the name of the ALICE Project-Team nor the names of its
*  contributors may be used to endorse or promote products derived from this
*  software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
*  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
*  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
*  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
*  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
*  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
*  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
*  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
*  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
*  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: Bruno Levy
*
*     Bruno.Levy@inria.fr
*     http://www.loria.fr/~levy
*
*     ALICE Project
*     LORIA, INRIA Lorraine,
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX
*     FRANCE
*
*
*  MODIFIED BY
*
* Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
* All rights reserved.
*
*     http://www.ring-team.org
*
*     RING Project
*     Ecole Nationale Superieure de Geologie - GeoRessources
*     2 Rue du Doyen Marcel Roubault - TSA 70605
*     54518 VANDOEUVRE-LES-NANCY
*     FRANCE
*/

#include <ringmesh/geogram_mesh_repair.h>

#include <geogram/points/colocate.h>

/*!
* @file Implementation of high level repairing functions on GEO::Mesh 
* @note Code modified from geogram/mesh/mesh_repair.cpp
* @author Jeanne Pellerin
*/

namespace RINGMesh {
    index_t detect_mesh_colocated_vertices(
        const GEO::Mesh& M, double tolerance, GEO::vector< index_t >& old2new )
    {
        index_t nb_unique_vertices = 0 ;
        if( tolerance == 0.0 ) {
            nb_unique_vertices = GEO::Geom::colocate_by_lexico_sort(
                M.vertices.point_ptr( 0 ), 3, M.vertices.nb(), old2new,
                M.vertices.dimension() ) ;
        } else {
            nb_unique_vertices = GEO::Geom::colocate(
                M.vertices.point_ptr( 0 ), 3, M.vertices.nb(), old2new,
                tolerance, M.vertices.dimension() ) ;
        }
        index_t nb_colocated_vertices( M.vertices.nb() - nb_unique_vertices ) ;
        return nb_colocated_vertices ;
    }

    bool has_mesh_colocate_vertices( const GEO::Mesh& M, double tolerance )
    {
        GEO::vector< index_t > old2new ;
        index_t nb_colocated_vertices = detect_mesh_colocated_vertices( M, tolerance, old2new ) ;
        if( nb_colocated_vertices == 0 ) {
            return false ;
        } else {
            return true ;
        }
    }

    void update_mesh_edges_vertices( GEO::Mesh& M, const GEO::vector<index_t>& old2new )
    {
        for( index_t e = 0; e < M.edges.nb(); ++e ) {
            M.edges.set_vertex( e, 0, old2new[ M.edges.vertex( e, 0 ) ] );
            M.edges.set_vertex( e, 1, old2new[ M.edges.vertex( e, 1 ) ] );
        }
    }
    void update_mesh_facets_vertices( GEO::Mesh& M, const GEO::vector<index_t>& old2new )
    {
        for( index_t c = 0; c < M.facet_corners.nb(); ++c ) {
            M.facet_corners.set_vertex( c, old2new[ M.facet_corners.vertex( c ) ] );
        }
    }

    void update_mesh_cells_vertices( GEO::Mesh& M, const GEO::vector<index_t>& old2new )
    {
        for( index_t ce = 0; ce < M.cells.nb(); ++ce ) {
            for( index_t c = M.cells.corners_begin( ce );
                 c<M.cells.corners_end( ce ); ++c
                 ) {
                M.cell_corners.set_vertex( c, old2new[ M.cell_corners.vertex( c ) ] );
            }
        }
    }
    void delete_colocated_vertices( GEO::Mesh& M, GEO::vector< index_t >& old2new )
    {
        for( index_t i = 0; i < old2new.size(); i++ ) {
            if( old2new[ i ] == i ) {
                old2new[ i ] = 0;
            } else {
                old2new[ i ] = 1;
            }
        }
        M.vertices.delete_elements( old2new );
    }


    void repair_colocate_vertices( GEO::Mesh& M, double colocate_epsilon )
    {
        GEO::vector<index_t> old2new;
        index_t nb_colocated_vertices = detect_mesh_colocated_vertices( M, colocate_epsilon, old2new ) ;
        if( nb_colocated_vertices == 0 ) {
            return ;
        }

        GEO::Logger::out( "GeoModel" ) << "Removing "
            << nb_colocated_vertices
            << " duplicated vertices" << std::endl;

        update_mesh_edges_vertices( M, old2new ) ;
        update_mesh_facets_vertices( M, old2new ) ;
        update_mesh_cells_vertices( M, old2new ) ;

        delete_colocated_vertices( M, old2new ) ;
    }
}