/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/tetra_gen.h>

#include <cstdio>
#include <iomanip>
#include <stack>
#include <sstream>

#ifdef WIN32
#   include <io.h>
#endif 


#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>

#include <ringmesh/algorithm.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_builder.h>
#include <ringmesh/geo_model_entity.h>
#include <ringmesh/geometry.h>
#include <ringmesh/tetgen_mesher.h>
#include <ringmesh/well.h>

/*!
 * @file Implementation of tetrahedral meshing region per region of a GeoModel
 * @author Arnaud Botella
 */

namespace RINGMesh {

#ifdef RINGMESH_WITH_TETGEN
    class RINGMESH_API TetraGen_TetGen: public TetraGen {
    public:
        TetraGen_TetGen()
            : TetraGen()
        {
        }
        virtual ~TetraGen_TetGen()
        {
        }

        virtual bool tetrahedralize( bool refine )
        {
            tetrahedralize_mesh_tetgen( tetmesh_constrain, refine, 1.0 ) ;
            initialize_storage( tetmesh_constrain.vertices.nb(),
                tetmesh_constrain.cells.nb() ) ;
            write_vertices_in_ringmesh_data_structure() ;
            write_tet_in_ringmesh_data_structure() ;
            return true ;
        }

        virtual void write_vertices_in_ringmesh_data_structure()
        {
            for( index_t v = 0; v < tetmesh_constrain.vertices.nb(); v++ ) {
                set_point( v, tetmesh_constrain.vertices.point_ptr( v ) ) ;
            }
        }
        virtual void write_tet_in_ringmesh_data_structure()
        {
            for( index_t tet = 0; tet < tetmesh_constrain.cells.nb(); tet++ ) {
                int vertex_indices[4] ;
                for( index_t v = 0; v < 4; v++ ) {
                    vertex_indices[v] = static_cast< int >( tetmesh_constrain.cells.vertex( tet, v ) ) ;
                }
                set_tetra( tet, vertex_indices ) ;
            }
            builder_->compute_region_adjacencies( output_region_ ) ;
        }
    } ;
#endif

#ifdef USE_MG_TETRA

    void start_redirect( fpos_t& pos, FILE* out, int& fd )
    {
#ifndef RINGMESH_DEBUG
        //Save position of current standard output
        fgetpos( out, &pos ) ;
#   ifdef WIN32
        fd = _dup( fileno( out ) ) ;
        freopen( "nul", "w", out ) ;
#   else
        fd = dup( fileno( out ) ) ;
        FILE* f = freopen( "/dev/null", "w", out ) ;
        ringmesh_unused(f) ;
#   endif
#endif
    }

    void stop_redirect( fpos_t& pos, FILE* out, int& fd )
    {
#ifndef RINGMESH_DEBUG
        //Flush stdout so any buffered messages are delivered
        fflush( out ) ;
        //Close file and restore standard output to stdout - which should be the terminal
#   ifdef WIN32
        _dup2( fd, fileno( out ) ) ;
#   else
        dup2( fd, fileno( out ) ) ;
#   endif
        close( fd ) ;
        clearerr( out ) ;
        fsetpos( out, &pos ) ;
#endif
    }

    class RINGMESH_API TetraGen_MG_Tetra: public TetraGen {
    public:
        TetraGen_MG_Tetra()
            :
                TetraGen(),
                context_( nil ),
                mesh_input_( nil ),
                mesh_output_( nil ),
                tms_( nil )
        {

        }

        virtual ~TetraGen_MG_Tetra()
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

        bool tetrahedralize( bool refine )
        {
            fpos_t pos ;
            int fd = 0 ;
            start_redirect( pos, stdout, fd ) ;
            fpos_t pos_err ;
            int fd_err = 0 ;
            start_redirect( pos_err, stderr, fd_err ) ;

            initialize_mgtetra_variables() ;

            set_mesh_in_mgtetra() ;

            set_meshing_parameters() ;

            generate_mesh( refine ) ;

            write_vertices_in_ringmesh_data_structure() ;
            write_tet_in_ringmesh_data_structure() ;

            stop_redirect( pos, stdout, fd ) ;
            stop_redirect( pos_err, stderr, fd_err ) ;

            return true ;
        }

        static status_t my_message_cb( message_t * msg, void *user_data )
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
                std::cerr << "two surface edges are intersecting : " << ibuff[0]
                    << " " << ibuff[1] << " intersects " << ibuff[2] << " "
                    << ibuff[3] << std::endl ;
            } else if( e == MESHGEMS_TETRA_CODE( -5120 ) ) {
                message_get_integer_data( msg, 1, 5, ibuff ) ;
                std::cerr << "surface edge intersects a surface face : " << ibuff[0]
                    << " " << ibuff[1] << " intersects " << ibuff[2] << " "
                    << ibuff[3] << " " << ibuff[4] << std::endl ;
            } else if( e == MESHGEMS_TETRA_CODE( -5150 ) ) {
                message_get_integer_data( msg, 1, 4, ibuff ) ;
                std::cerr << "boundary point inside a surface face : " << ibuff[0]
                    << " in " << ibuff[1] << " " << ibuff[2] << " " << ibuff[3]
                    << std::endl ;
            } else if( e == MESHGEMS_TETRA_CODE( 5200 ) ) {
                message_get_integer_data( msg, 1, 3, ibuff ) ;
                std::cerr << "duplicated face : " << ibuff[0] << " " << ibuff[1]
                    << " " << ibuff[2] << " " << ibuff[3] << std::endl ;
            } else if( e == MESHGEMS_TETRA_CODE( -5621 ) ) {
                message_get_integer_data( msg, 1, 4, ibuff ) ;
                message_get_real_data( msg, 1, 1, rbuff ) ;
                std::cerr << "degenerated face : face " << ibuff[0] << " ("
                    << ibuff[1] << ", " << ibuff[2] << ", " << ibuff[3]
                    << ") with small inradius = " << rbuff[0] << std::endl ;
            } else if( e == MESHGEMS_TETRA_CODE( -5820 ) ) {
                message_get_integer_data( msg, 1, 2, ibuff ) ;
                std::cerr << "edge bounding a hole : " << ibuff[0] << " " << ibuff[1]
                    << std::endl ;
            } else if( e == MESHGEMS_TETRA_CODE( 8423 ) ) {
                message_get_integer_data( msg, 1, 3, ibuff ) ;
                std::cerr << "constrained face cannot be enforced : " << ibuff[0]
                    << " " << ibuff[1] << " " << ibuff[2] << std::endl ;
            } else if( e == MESHGEMS_TETRA_CODE( 8441 ) ) {
                message_get_integer_data( msg, 1, 2, ibuff ) ;
                std::cerr << "constrained edge cannot be enforced : " << ibuff[0]
                    << " " << ibuff[1] << std::endl ;
            } else {
                std::cerr << "Error message not directly handle" << std::endl ;
                std::cerr << "Error(" << e << ") : " << desc << std::endl ;
            }
            return STATUS_OK ;
        }

    private:
        context_t* context_ ;
        mesh_t* mesh_input_ ;
        mesh_t* mesh_output_ ;
        tetra_session_t* tms_ ;

    private:

        void initialize_ringmesh_storage()
        {
            tetra_get_mesh( tms_, &mesh_output_ ) ;
            signed_index_t nb_points = 0 ;
            mesh_get_vertex_count( mesh_output_, &nb_points ) ;
            signed_index_t nb_tets = 0 ;
            mesh_get_tetrahedron_count( mesh_output_, &nb_tets ) ;
            initialize_storage( nb_points, nb_tets ) ;
        }

        void initialize_mgtetra_variables() {
            context_ = context_new() ;
            mesh_input_ = mesh_new_in_memory( context_ ) ;
        }

        void set_mesh_in_mgtetra()
        {
            set_vertices() ;
            set_edges() ;
            set_triangles() ;
        }

        bool generate_mesh( bool refine )
        {
            if( !create_boundary_mesh() ) {
                return false ;
            }
            if( refine ) {
                return refine_mesh() ;
            }
            return true ;
        }

        void set_vertices()
        {
            mesh_set_vertex_count( mesh_input_, tetmesh_->vertices.nb() ) ;
            for( index_t p = 0; p < tetmesh_->vertices.nb(); p++ ) {
                mesh_set_vertex_coordinates( mesh_input_, p + 1,
                    tetmesh_->vertices.point_ptr( p ) ) ;
            }
        }

        void set_edges()
        {
            mesh_set_edge_count( mesh_input_, tetmesh_->edges.nb() ) ;
            for( index_t e = 0; e < tetmesh_->edges.nb(); e++ ) {
                meshgems_integer edge_indices[2] ;
                edge_indices[0] = tetmesh_->edges.vertex( e, 0 ) + 1 ;
                edge_indices[1] = tetmesh_->edges.vertex( e, 1 ) + 1 ;
                mesh_set_edge_vertices( mesh_input_, e + 1, edge_indices ) ;
            }

        }

        void set_triangles()
        {
            mesh_set_triangle_count( mesh_input_, tetmesh_->facets.nb() ) ;
            for( index_t t = 0; t < tetmesh_->facets.nb(); t++ ) {
                meshgems_integer triangle_indices[3] ;
                triangle_indices[0] = tetmesh_->facets.vertex( t, 0 ) + 1 ;
                triangle_indices[1] = tetmesh_->facets.vertex( t, 1 ) + 1 ;
                triangle_indices[2] = tetmesh_->facets.vertex( t, 2 ) + 1 ;
                mesh_set_triangle_vertices( mesh_input_, t + 1, triangle_indices ) ;
            }
        }

        void set_meshing_parameters()
        {
            tms_ = tetra_session_new( context_ ) ;
            tetra_set_surface_mesh( tms_, mesh_input_ ) ;
            tetra_set_param( tms_, "verbose", "4" ) ;
            tetra_set_param( tms_, "components", "all" ) ;
            tetra_set_param( tms_, "optimisation_level", "standard" ) ;
            tetra_set_param( tms_, "gradation", "1.1" ) ;
            tetra_set_param( tms_, "pthreads_mode", "aggressive" ) ;
            tetra_set_param( tms_, "max_number_of_threads", "8" ) ;
            tetra_set_param( tms_, "max_error_count", "5" ) ;
        }

        bool create_boundary_mesh()
        {
            status_t ret = tetra_mesh_boundary( tms_ ) ;
            if( ret != STATUS_OK ) {
                Logger::err( "TetraGen" )
                    << "Encountered a problem while meshing boundary..."
                    << std::endl ;
                return false ;
            }
            return true ;
        }

        bool refine_mesh()
        {
            status_t ret = tetra_insert_volume_vertices( tms_ ) ;
            if( ret != STATUS_OK ) {
                Logger::err( "TetraGen" )
                    << "Encountered a problem while meshing inside..." << std::endl ;
                return false ;
            }
            ret = tetra_optimise_volume_regular( tms_ ) ;
            if( ret != STATUS_OK ) {
                Logger::err( "TetraGen" )
                    << "Encountered a problem while meshing inside..." << std::endl ;
                return false ;
            }
            return true ;

        }

        void write_vertices_in_ringmesh_data_structure()
        {
            RINGMESH_PARALLEL_LOOP
            for( index_t p = 0; p < tetmesh_->vertices.nb(); p++ ) {
                double point[3] ;
                mesh_get_vertex_coordinates( mesh_output_, p + 1, point ) ;
                set_point( p, point ) ;
            }
        }

        void write_tet_in_ringmesh_data_structure()
        {
            RINGMESH_PARALLEL_LOOP
            for( index_t t = 0; t < tetmesh_->cells.nb(); t++ ) {
                signed_index_t tet[4] ;
                mesh_get_tetrahedron_vertices( mesh_output_, t + 1, tet ) ;
                set_tetra( t, tet, tetmesh_->edges.nb(), tetmesh_->facets.nb() ) ;
            }
            tetmesh_->cells.connect() ;
        }

    } ;
#endif

    TetraGen* TetraGen::create(
        GeoModel& M,
        index_t region_id,
        const std::string& algo_name )
    {
        TetraGen* mesher = TetraGenFactory::create_object( algo_name ) ;
        if( !mesher ) {
#ifdef RINGMESH_WITH_TETGEN
            Logger::warn( "TetraGen" ) << "Could not create TetraGen mesher: "
                << algo_name << std::endl ;
            Logger::warn( "TetraGen" ) << "Falling back to TetGen mode"
                << std::endl ;
            mesher = new TetraGen_TetGen() ;
#else
            std::vector< std::string > names ;
            TetraGenFactory::list_creators( names ) ;
            Logger::err( "I/O" ) << "Currently supported mesher are: " ;
            for( index_t i = 0; i < names.size(); i++ ) {
                Logger::out( "I/O" ) << " " << names[ i ] ;
            }
            Logger::out( "I/O" ) << std::endl ;
            throw RINGMeshException( "TetraGen", "Could not create TetraGen mesher: "
                + algo_name ) ;
#endif
        }

        mesher->builder_ = new GeoModelBuilder( M ) ;
        mesher->output_region_ = region_id ;
        return mesher ;
    }

    TetraGen::TetraGen()
        : builder_( nil ), output_region_( NO_ID ), region_( nil ), wells_( nil )
    {
    }

    /*!
     * Sets the boundaries of the domain
     * @param[in] region The Region of the GeoModel to mesh
     * @param[in] wells the wells to be conformal to
     */
    void TetraGen::set_boundaries(
        const Region& region,
        const WellGroup* wells )
    {
        region_ = &region ;
        index_t nb_surfaces = region_->nb_boundaries() ;
        std::vector< const GeoModelMeshEntity* > unique_surfaces ;
        unique_surfaces.reserve( nb_surfaces ) ;
        std::vector< index_t > surface_id ;
        surface_id.reserve( nb_surfaces ) ;
        index_t nb_surface_points = 0, nb_facets = 0 ;
        for( index_t s = 0; s < nb_surfaces; s++ ) {
            const Surface& surface =
                dynamic_cast< const Surface& >( region_->boundary( s ) ) ;
            if( contains( surface_id, surface.index() ) ) continue ;
            nb_surface_points += surface.nb_vertices() ;
            nb_facets += surface.nb_mesh_elements() ;

            surface_id.push_back( surface.index() ) ;
            unique_surfaces.push_back( &surface ) ;
        }

        MakeUnique uniqueID( unique_surfaces, true ) ;
        std::vector< std::vector< Edge > > well_edges ;
        if( wells ) {
            wells->get_region_edges( region.index(), well_edges ) ;
            // Copy result of porting. Stupid, I know, but because of the interface
            // of MakeUnique. This Edge class is a pain [JP]
            for( index_t w = 0; w < well_edges.size(); w++ ) {
                std::vector< std::pair< vec3, vec3 > > wells_copy(
                    well_edges[w].size() ) ;
                for( index_t i = 0; i < wells_copy.size(); ++i ) {
                    wells_copy[i] = std::pair< vec3, vec3 >(
                        well_edges[w][i].value( 0 ), well_edges[w][i].value( 1 ) ) ;
                }
                uniqueID.add_edges( wells_copy ) ;
            }
        }

        std::vector< vec3 > region_vertices( region.nb_vertices() ) ;
        for( index_t v = 0; v < region.nb_vertices(); v++ ) {
            region_vertices[v] = region.vertex( v ) ;
        }
        uniqueID.add_points( region_vertices ) ;

        uniqueID.unique() ;
        const std::vector< index_t >& unique_indices = uniqueID.indices() ;
        std::vector< vec3 > unique_points ;
        uniqueID.unique_points( unique_points ) ;

        index_t starting_index = tetmesh_constrain.vertices.create_vertices(
            unique_points.size() ) ;
        GEO::Memory::copy( tetmesh_constrain.vertices.point_ptr( starting_index ),
            unique_points.data()->data(),
            3 * sizeof(double) * unique_points.size() ) ;

        if( !well_edges.empty() ) {
            index_t nb_well_edges = 0 ;
            for( index_t w = 0; w < well_edges.size(); w++ ) {
                nb_well_edges += well_edges[w].size() ;
            }
            tetmesh_constrain.edges.create_edges( nb_well_edges ) ;
            GEO::Attribute< index_t > edge_region( tetmesh_constrain.edges.attributes(),
                surface_att_name ) ;
            index_t cur_vertex_id = nb_surface_points ;
            index_t cur_edge = 0 ;
            for( index_t w = 0; w < well_edges.size(); w++ ) {
                for( index_t e = 0; e < well_edges[w].size(); e++ ) {
                    tetmesh_constrain.edges.set_vertex( cur_edge, 0,
                        starting_index + unique_indices[cur_vertex_id++ ] ) ;
                    tetmesh_constrain.edges.set_vertex( cur_edge, 1,
                        starting_index + unique_indices[cur_vertex_id++ ] ) ;
                    edge_region[cur_edge++ ] = w ;
                }
            }
        }

        index_t offset_vertices = 0 ;
        index_t offset_facets = 0 ;
        tetmesh_constrain.facets.create_triangles( nb_facets ) ;
        GEO::Attribute< index_t > surface_region( tetmesh_constrain.facets.attributes(),
            surface_att_name ) ;
        for( index_t s = 0; s < unique_surfaces.size(); s++ ) {
            const Surface& surface =
                dynamic_cast< const Surface& >( *unique_surfaces[s] ) ;
            RINGMESH_PARALLEL_LOOP
            for( index_t t = 0; t < surface.nb_mesh_elements(); t++ ) {
                ringmesh_assert( surface.facet_is_triangle( t ) ) ;
                for( index_t v = 0; v < 3; v++ ) {
                    tetmesh_constrain.facets.set_vertex( offset_facets + t, v,
                        starting_index
                            + unique_indices[offset_vertices
                                + surface.mesh_element_vertex_index( t, v )] ) ;
                }
                surface_region[offset_facets + t] = surface.index() ;

            }
            offset_vertices += surface.nb_vertices() ;
            offset_facets += surface.nb_mesh_elements() ;
        }
        tetmesh_constrain.facets.connect() ;
    }

    /*!
     * Set additional points to be in the output tetrahedral mesh
     * @param[in] points the points to add
     */
    void TetraGen::set_internal_points( const std::vector< vec3 >& points )
    {
        if( points.empty() ) return ;
        index_t start = tetmesh_constrain.vertices.create_vertices( points.size() ) ;
        GEO::Memory::copy( tetmesh_constrain.vertices.point_ptr( start ),
            points.front().data(), points.size() * 3 * sizeof(double) ) ;
    }

    TetraGen::~TetraGen()
    {
        if( builder_ ) delete builder_ ;
    }

    void TetraGen::set_point( index_t index, const double* point )
    {
        bool update = false ;
        vec3 vertex( point ) ;
        builder_->set_mesh_entity_vertex( GME::gme_t( Region::type_name_, output_region_ ),
            index, vertex, update ) ;
    }

    void TetraGen::set_tetra(
        index_t tetra_index,
        int* vertex_indices )
    {
        std::vector< index_t > corners( 4 ) ;
        for( index_t v = 0; v < 4; v++ ) {
            index_t vertex_id = static_cast< index_t >( vertex_indices[v] ) ;
            corners[v] = vertex_id ;
        }
        builder_->set_region_element_geometry( output_region_, tetra_index, corners ) ;
    }

    void TetraGen::initialize_storage( index_t nb_points, index_t nb_tets )
    {
        GME::gme_t region_id( Region::type_name_, output_region_ ) ;
        builder_->delete_mesh_entity_mesh( region_id ) ;
        builder_->create_mesh_entity_vertices( region_id, nb_points ) ;
        builder_->create_region_cells( output_region_, GEO::MESH_TET,
            nb_tets ) ;
    }

    void TetraGen::initialize()
    {
#ifdef RINGMESH_WITH_TETGEN
        ringmesh_register_tetragen( TetraGen_TetGen, "TetGen" ) ;
#endif

#ifdef USE_MG_TETRA
        ringmesh_register_tetragen( TetraGen_MG_Tetra, "MG_Tetra" );
#endif
    }
}
