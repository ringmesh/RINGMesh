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
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
*/


#include <ringmesh/tetra_gen.h>
#include <ringmesh/boundary_model_element.h>
#include <ringmesh/well.h>

#include <geogram/mesh/mesh_tetrahedralize.h>

#include <iomanip>
#include <stack>
#include <sstream>
#include <cstdio>

#ifdef WIN32
    #include <io.h>
#endif 


namespace RINGMesh {

    void start_redirect( fpos_t& pos, FILE* out, int& fd )
    {
#ifndef RINGMESH_DEBUG
        //Save position of current standard output
        fgetpos( out, &pos ) ;
#ifdef WIN32
        fd = _dup( fileno( out ) ) ;
        freopen( "nul", "w", out ) ;
#else
        fd = dup( fileno( out ) ) ;
        FILE* f = freopen( "/dev/null", "w", out ) ;
#endif
#endif
    }

    void stop_redirect( fpos_t& pos, FILE* out, int& fd )
    {
#ifndef RINGMESH_DEBUG
        //Flush stdout so any buffered messages are delivered
        fflush( out ) ;
        //Close file and restore standard output to stdout - which should be the terminal
#ifdef WIN32
        _dup2( fd, fileno( out ) ) ;
#else
        dup2( fd, fileno( out ) ) ;
#endif
        close( fd ) ;
        clearerr( out ) ;
        fsetpos( out, &pos ) ;
#endif
    }

    TetraGen_var TetraGen::instantiate(
        const TetraMethod& method,
        GEO::Mesh& tetmesh,
        const BoundaryModelElement* region,
        bool add_steiner_points,
        const std::vector< vec3 >& internal_vertices,
        const WellGroup* wells )
    {
        switch( method ) {
            case TetGen:
                return new TetraGen_TetGen( tetmesh, region, add_steiner_points,
                    internal_vertices, wells ) ;
#ifdef USE_MG_TETRA
            case MG_Tetra:
                return new TetraGen_MG_Tetra( tetmesh, region, add_steiner_points,
                    internal_vertices, wells ) ;
#endif
            default:
                return nil ;
        }
    }

    TetraGen::TetraGen(
        GEO::Mesh& tetmesh,
        const BoundaryModelElement* region,
        bool refine,
        const std::vector< vec3 >& internal_vertices,
        const WellGroup* wells )
        :
            tetmesh_( tetmesh ),
            region_( region ),
            wells_( wells ),
            refine_( refine )
    {


        index_t nb_surfaces = region->nb_boundaries() ;
        std::vector< const BoundaryModelElement* > unique_surfaces ;
        unique_surfaces.reserve( nb_surfaces ) ;
        std::vector< index_t > surface_id ;
        index_t nb_points = 0, nb_facets = 0 ;
        for( index_t s = 0; s < nb_surfaces; s++ ) {
            const BoundaryModelElement& surface = region->boundary( s ) ;
            if( Utils::contains( surface_id, surface.bme_id().index ) ) continue ;
            nb_points += surface.nb_vertices() ;
            nb_facets += surface.nb_cells() ;

            surface_id.push_back( surface.bme_id().index ) ;
            unique_surfaces.push_back( &surface ) ;
        }

        std::vector< std::vector< Edge > > well_edges ;
        if( wells_ ) {
            wells->get_region_edges( region->bme_id().index, well_edges ) ;
        }
        index_t nb_points_without_well = nb_points ;
        nb_points += well_edges.size() ;
        MakeUnique uniqueID( unique_surfaces, true ) ;
        for( index_t w = 0; w < well_edges.size(); w++ ) {
            uniqueID.add_edges( well_edges[w] ) ;
        }
        uniqueID.unique() ;

        std::vector< vec3 > unique_points ;
        uniqueID.unique_points( unique_points ) ;
        const std::vector< index_t >& unique_indices = uniqueID.indices() ;
        internal_vertices_ptr_ = unique_points.size() ;
        tetmesh_.vertices.create_vertices(
            internal_vertices_ptr_ + internal_vertices.size() ) ;
        for( index_t p = 0; p < unique_points.size(); p++ ) {
            tetmesh_.vertices.point( p ) = unique_points[p] ;
        }
        for( index_t p = 0; p < internal_vertices.size(); p++ ) {
            tetmesh_.vertices.point( internal_vertices_ptr_ + p ) =
                internal_vertices[p] ;
        }

        if( !well_edges.empty() ) {
            index_t nb_well_edges = 0 ;
            for( index_t w = 0; w < well_edges.size(); w++ ) {
                nb_well_edges += well_edges[w].size() ;
            }
//            tetmesh_.edges.create_edges( nb_well_edges ) ;
            GEO::Attribute< index_t > edge_region( tetmesh_.edges.attributes(),
                surface_att_name ) ;
            index_t cur_vertex_id = nb_points_without_well ;
            for( index_t w = 0; w < well_edges.size(); w++ ) {
                for( index_t e = 0; e < well_edges[w].size(); e++ ) {
                    const Edge& cur_edge = well_edges[w][e] ;
                    index_t cur_id = tetmesh_.edges.create_edge(
                        unique_indices[cur_vertex_id++ ],
                        unique_indices[cur_vertex_id++ ] ) ;
                    edge_region[cur_id] = w ;
                }
            }
        }

        index_t offset = 0 ;
        index_t cur_facet = 0 ;
        tetmesh_.facets.create_triangles( nb_facets ) ;
        GEO::Attribute< index_t > surface_region( tetmesh_.facets.attributes(),
            surface_att_name ) ;
        for( index_t s = 0; s < unique_surfaces.size(); s++ ) {
            const Surface& surface = dynamic_cast< const Surface& >( *unique_surfaces[s] ) ;
            for( index_t t = 0; t < surface.nb_cells(); t++ ) {
                ringmesh_debug_assert( surface.is_triangle( t ) ) ;
                for( index_t v = 0; v < 3; v++ ) {
                    tetmesh_.facets.set_vertex( cur_facet, v,
                        unique_indices[offset + surface.surf_vertex_id( t, v )] ) ;
                }
                surface_region[cur_facet++] = surface.bme_id().index ;

            }
            offset += surface.nb_vertices() ;
        }
        tetmesh_.facets.connect() ;
    }

    TetraGen::~TetraGen()
    {
    }

    void TetraGen::initialize_storage(
        index_t nb_points,
        index_t nb_tets )
    {
        tetmesh_.vertices.clear( true, false ) ;
        tetmesh_.vertices.create_vertices( nb_points ) ;
        tetmesh_.cells.create_tets( nb_tets ) ;
    }

    void TetraGen::set_point( index_t index, const double* point )
    {
        for( index_t i = 0; i < 3; i++ ) {
            tetmesh_.vertices.point_ptr( index )[i] = point[i] ;
        }
    }

    void TetraGen::set_tetra( index_t index, int* tet, index_t nb_lines, index_t nb_triangles )
    {
        index_t corner_begin = tetmesh_.cells.corners_begin( index ) ;
        for( index_t v = 0; v < 4; v++ ) {
            tetmesh_.cell_corners.set_vertex( corner_begin++, tet[v] - 1  ) ;
        }
    }

    TetraGen_TetGen::TetraGen_TetGen(
        GEO::Mesh& tetmesh,
        const BoundaryModelElement* region,
        bool add_steiner_points,
        const std::vector< vec3 >& internal_vertices,
        const WellGroup* wells )
        :
            TetraGen( tetmesh, region, add_steiner_points, internal_vertices, wells )
    {
    }

    bool TetraGen_TetGen::tetrahedralize()
    {
        GEO::mesh_tetrahedralize( tetmesh_, false, refine_, 1.0 ) ;
        Utils::check_and_repair_mesh_consistency( *region_, tetmesh_ ) ;
        return true ;
    }

#ifdef USE_MG_TETRA

    TetraGen_MG_Tetra::TetraGen_MG_Tetra(
        GEO::Mesh& tetmesh,
        const BoundaryModelElement* region,
        bool add_steiner_points,
        const std::vector< vec3 >& internal_vertices,
        const WellGroup* wells )
        :
            TetraGen( tetmesh, region, add_steiner_points, internal_vertices, wells ),
            mesh_output_( nil )
    {
        fpos_t pos ;
        int fd = 0 ;
        start_redirect( pos, stdout, fd ) ;
        fpos_t pos_err ;
        int fd_err = 0 ;
        start_redirect( pos_err, stderr, fd_err ) ;

        context_ = context_new() ;
        mesh_input_ = mesh_new_in_memory( context_ ) ;
        context_set_message_callback(context_, my_message_cb, 0);

        mesh_set_vertex_count( mesh_input_, nb_total_points() ) ;
        for( index_t p = 0; p < tetmesh_.vertices.nb(); p++ ) {
            mesh_set_vertex_coordinates( mesh_input_, p + 1,
                tetmesh_.vertices.point_ptr( p ) ) ;
        }

        mesh_set_edge_count( mesh_input_, tetmesh_.edges.nb() ) ;
        for( index_t e = 0; e < tetmesh_.edges.nb(); e++ ) {
            meshgems_integer edge_indices[2] ;
            edge_indices[0] = tetmesh_.edges.vertex( e, 0 ) ;
            edge_indices[1] = tetmesh_.edges.vertex( e, 1 ) ;
            mesh_set_edge_vertices( mesh_input_, e + 1, edge_indices ) ;
        }

        mesh_set_triangle_count( mesh_input_, tetmesh_.facets.nb() ) ;
        for( index_t t = 0; t < tetmesh_.facets.nb(); t++ ) {
            meshgems_integer triangle_indices[3] ;
            triangle_indices[0] = tetmesh_.facets.vertex( t, 0 ) ;
            triangle_indices[1] = tetmesh_.facets.vertex( t, 1 ) ;
            triangle_indices[2] = tetmesh_.facets.vertex( t, 2 ) ;
            mesh_set_triangle_vertices( mesh_input_, t + 1, triangle_indices ) ;
        }

        tms_ = tetra_session_new( context_ ) ;
        tetra_set_surface_mesh( tms_, mesh_input_ ) ;
        tetra_set_param( tms_, "verbose", "4" ) ;
        tetra_set_param( tms_, "components", "all" ) ;
        tetra_set_param( tms_, "optimisation_level", "standard" ) ;
        tetra_set_param( tms_, "gradation", "1.1" ) ;
        tetra_set_param( tms_, "pthreads_mode", "aggressive" ) ;
        tetra_set_param( tms_, "max_number_of_threads", "8" ) ;
        tetra_set_param( tms_, "max_error_count", "5" ) ;

        stop_redirect( pos, stdout, fd ) ;
        stop_redirect( pos_err, stderr, fd_err ) ;
    }

    TetraGen_MG_Tetra::~TetraGen_MG_Tetra()
    {
        fpos_t pos ;
        int fd = 0 ;
        start_redirect( pos, stdout, fd ) ;
        fpos_t pos_err ;
        int fd_err = 0 ;
        start_redirect( pos_err, stderr, fd_err ) ;

        tetra_regain_mesh( tms_, mesh_output_ ) ;
        tetra_session_delete( tms_ ) ;
        mesh_delete( mesh_input_ ) ;
        context_delete( context_ ) ;

        stop_redirect( pos, stdout, fd ) ;
        stop_redirect( pos_err, stderr, fd_err ) ;
    }

    bool TetraGen_MG_Tetra::tetrahedralize()
    {
        fpos_t pos ;
        int fd = 0 ;
        start_redirect( pos, stdout, fd ) ;
        fpos_t pos_err ;
        int fd_err = 0 ;
        start_redirect( pos_err, stderr, fd_err ) ;

        status_t ret = tetra_mesh_boundary( tms_ ) ;
        if( ret != STATUS_OK ) {
            std::cout << "Encountered a problem while meshing boundary..."
                << std::endl ;
            return false ;
        }
        if( refine_ ) {
            ret = tetra_insert_volume_vertices( tms_ ) ;
            if( ret != STATUS_OK ) {
                std::cout << "Encountered a problem while meshing inside..."
                    << std::endl ;
                return false ;
            }
            ret = tetra_optimise_volume_regular( tms_ ) ;
            if( ret != STATUS_OK ) {
                std::cout << "Encountered a problem while meshing inside..."
                    << std::endl ;
                return false ;
            }
        }
        tetra_get_mesh( tms_, &mesh_output_ ) ;
        signed_index_t nb_points = 0 ;
        mesh_get_vertex_count( mesh_output_, &nb_points ) ;
        signed_index_t nb_tets = 0 ;
        mesh_get_tetrahedron_count( mesh_output_, &nb_tets ) ;
        signed_index_t nb_triangles = 0 ;
        mesh_get_triangle_count( mesh_output_, &nb_triangles ) ;
        signed_index_t nb_lines = 0 ;
        mesh_get_edge_count( mesh_output_, &nb_lines ) ;

        initialize_storage( nb_points, nb_tets ) ;
#pragma omp parallel for
        for( index_t t = 0; t < nb_tets; t++ ) {
            signed_index_t tet[4] ;
            mesh_get_tetrahedron_vertices( mesh_output_, t+1, tet ) ;
            set_tetra( t, tet, nb_lines, nb_triangles ) ;
        }

#pragma omp parallel for
        for( index_t p = 0; p < nb_points; p++ ) {
            double point[3] ;
            mesh_get_vertex_coordinates( mesh_output_, p+1, point ) ;
            set_point( p, point ) ;
        }
        tetmesh_.cells.connect() ;
        Utils::check_and_repair_mesh_consistency( *region_, tetmesh_ ) ;

        stop_redirect( pos, stdout, fd ) ;
        stop_redirect( pos_err, stderr, fd_err ) ;

        return true ;
    }

    status_t TetraGen_MG_Tetra::my_message_cb( message_t * msg, void *user_data )
    {
        char *desc ;
        integer e, ibuff[6] ;
        real rbuff[3] ;

        message_get_description( msg, &desc ) ;
        message_get_number( msg, &e ) ;

        if( e == 0 ) {
            std::cerr << desc << std::endl ;
        } else if( e == MESHGEMS_TETRA_CODE( -5110 ) ) {
            message_get_integer_data( msg, 1, 4, ibuff ) ;
            std::cerr << "two surface edges are intersecting : " << ibuff[0] << " "
                << ibuff[1] << " intersects " << ibuff[2] << " " << ibuff[3]
                << std::endl ;
        } else if( e == MESHGEMS_TETRA_CODE( -5120 ) ) {
            message_get_integer_data( msg, 1, 5, ibuff ) ;
            std::cerr << "surface edge intersects a surface face : " << ibuff[0]
                << " " << ibuff[1] << " intersects " << ibuff[2] << " " << ibuff[3]
                << " " << ibuff[4] << std::endl ;
        } else if( e == MESHGEMS_TETRA_CODE( -5150 ) ) {
            message_get_integer_data( msg, 1, 4, ibuff ) ;
            std::cerr << "boundary point inside a surface face : " << ibuff[0]
                << " in " << ibuff[1] << " " << ibuff[2] << " " << ibuff[3]
                << std::endl ;
        } else if( e == MESHGEMS_TETRA_CODE( 5200 ) ) {
            message_get_integer_data( msg, 1, 3, ibuff ) ;
            std::cerr << "duplicated face : " << ibuff[0] << " " << ibuff[1] << " "
                << ibuff[2] << " " << ibuff[3] << std::endl ;
        } else if( e == MESHGEMS_TETRA_CODE( -5621 ) ) {
            message_get_integer_data( msg, 1, 4, ibuff ) ;
            message_get_real_data( msg, 1, 1, rbuff ) ;
            std::cerr << "degenerated face : face " << ibuff[0] << " (" << ibuff[1]
                << ", " << ibuff[2] << ", " << ibuff[3] << ") with small inradius = "
                << rbuff[0] << std::endl ;
        } else if( e == MESHGEMS_TETRA_CODE( -5820 ) ) {
            message_get_integer_data( msg, 1, 2, ibuff ) ;
            std::cerr << "edge bounding a hole : " << ibuff[0] << " " << ibuff[1]
                << std::endl ;
        } else if( e == MESHGEMS_TETRA_CODE( 8423 ) ) {
            message_get_integer_data( msg, 1, 3, ibuff ) ;
            std::cerr << "constrained face cannot be enforced : " << ibuff[0] << " "
                << ibuff[1] << " " << ibuff[2] << std::endl ;
        } else if( e == MESHGEMS_TETRA_CODE( 8441 ) ) {
            message_get_integer_data( msg, 1, 2, ibuff ) ;
            std::cerr << "constrained edge cannot be enforced : " << ibuff[0] << " "
                << ibuff[1] << std::endl ;
        } else {
            std::cerr << "Error message not directly handle" << std::endl ;
            std::cerr << "Error(" << e << ") : " << desc << std::endl ;
        }
        return STATUS_OK ;
    }
#endif

}
