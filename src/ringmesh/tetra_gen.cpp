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

#include <ringmesh/tetra_gen.h>
#include <ringmesh/geo_model_element.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/well.h>
#include <ringmesh/algorithm.h>
#include <ringmesh/geometry.h>

#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/mesh/mesh_AABB.h>

#include <iomanip>
#include <stack>
#include <sstream>
#include <cstdio>

#ifdef WIN32
#include <io.h>
#endif 

namespace RINGMesh {
    /*!
     * Tests if two adjacent facets have the same orientation
     * @param[in] mesh the mesh
     * @param[in] f1 the first facet index
     * @param[in] c11 the corner index in the first facet
     * @param[in] f2 the second facet index
     * @return the result of the test
     *
     * @todo Check that this code is not a duplicate of what is used to check validity
     *       of Geogram functions [JP]
     */
    bool facets_have_same_orientation(
        const GEO::Mesh& mesh,
        index_t f1,
        index_t c11,
        index_t f2 )
    {
        index_t c12 = mesh.facets.next_corner_around_facet( f1, c11 ) ;
        index_t v11 = mesh.facet_corners.vertex( c11 ) ;
        index_t v12 = mesh.facet_corners.vertex( c12 ) ;
        for( index_t c21 = mesh.facets.corners_begin( f2 );
            c21 < mesh.facets.corners_end( f2 ); c21++ ) {
            index_t c22 = mesh.facets.next_corner_around_facet( f2, c21 ) ;
            index_t v21 = mesh.facet_corners.vertex( c21 ) ;
            index_t v22 = mesh.facet_corners.vertex( c22 ) ;
            if( v11 == v21 && v12 == v22 ) {
                return false ;
            }
            if( v11 == v22 && v12 == v21 ) {
                return true ;
            }
        }
        return true ;
    }

    /*
     * @todo Check that this code is not a duplicate of what is used to check validity
     *       of of Geogram functions[ JP ]
     */
    void mesh_facet_connect( GEO::Mesh& mesh )
    {
        std::vector< index_t > temp ;
        temp.reserve( 7 ) ;
        std::vector< std::vector< index_t > > stars( mesh.vertices.nb(), temp ) ;
        for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
            for( index_t c = mesh.facets.corners_begin( f );
                c < mesh.facets.corners_end( f ); c++ ) {
                stars[mesh.facet_corners.vertex( c )].push_back( f ) ;
            }
        }
        for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
            for( index_t c = mesh.facets.corners_begin( f );
                c < mesh.facets.corners_end( f ); c++ ) {
                index_t f_adj = mesh.facet_corners.adjacent_facet( c ) ;
                if( f_adj != GEO::NO_FACET ) continue ;
                const std::vector< index_t >& star0 =
                    stars[mesh.facet_corners.vertex( c )] ;
                const std::vector< index_t >& star1 =
                    stars[mesh.facet_corners.vertex(
                        mesh.facets.next_corner_around_facet( f, c ) )] ;
                std::vector< index_t > intersect(
                    std::min( star0.size(), star1.size() ) ) ;
                intersect.erase(
                    std::set_intersection( star0.begin(), star0.end(), star1.begin(),
                        star1.end(), intersect.begin() ), intersect.end() ) ;
                if( intersect.size() > 1 ) {
                    for( index_t i = 0; i < intersect.size(); i++ ) {
                        index_t cur_f = intersect[i] ;
                        if( cur_f != f ) {
                            f_adj = cur_f ;
                        }
                    }
                    mesh.facet_corners.set_adjacent_facet( c, f_adj ) ;
                }
            }
        }
    }

    /*!
     * Repair the consistency between a GeoModel region
     * and its volumetric Mesh. It repairs duplicated facets and facet orientation
     * @param[in] region the GeoModel region
     * @param[in] mesh the mesh to repair
     * @param[in] check_duplicated_facet the test of duplicated facets is optional
     *
     * @todo Why not use mesh_repair functions of Geogram ?? Please comment. [JP]
     *       Check that this code is not a duplicate of what is used to check validity [JP]
     */
    void check_and_repair_mesh_consistency(
        const GeoModelElement& region,
        GEO::Mesh& mesh,
        bool check_duplicated_facet = false )
    {
        if( mesh.facets.nb() == 0 ) return ;

        GEO::Attribute< index_t > attribute( mesh.facets.attributes(),
            surface_att_name ) ;

        /// 0 - Remove duplicated facets (optionnal)
        if( check_duplicated_facet ) {
            std::vector< vec3 > barycenters( mesh.facets.nb(), vec3( 0, 0, 0 ) ) ;
            for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                barycenters[f] = GEO::Geom::mesh_facet_center( mesh, f ) ;
            }

            MakeUnique unique( barycenters ) ;
            unique.unique() ;
            const std::vector< index_t > indices = unique.indices() ;
            GEO::vector< index_t > facet_to_remove( mesh.facets.nb(), 0 ) ;
            signed_index_t cur_id = 0 ;
            for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                if( cur_id == indices[f] ) {
                    cur_id++ ;
                } else {
                    facet_to_remove[f] = 1 ;
                }
            }
            mesh.facets.delete_elements( facet_to_remove ) ;

            // I am pretty sure this Attribute resizing is done by Geogram [JP] 
            if( GEO::Attribute< index_t >::is_defined( mesh.facets.attributes(),
                surface_att_name ) ) {
                // Review : already defined above
                GEO::Attribute< index_t > attribute( mesh.facets.attributes(),
                    surface_att_name ) ;
                index_t offset = 0 ;
                cur_id = 0 ;
                for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                    if( cur_id == indices[f] ) {
                        cur_id++ ;
                        attribute[f - offset] = attribute[f] ;
                    } else {
                        offset++ ;
                    }
                }
                attribute.redim( attribute.size() - offset ) ;
            }
            mesh.facets.connect() ;
        }

        /// 1 - Check facet adjacencies for non-manifold surfaces
        std::vector< index_t > temp ;
        temp.reserve( 6 ) ;
        std::vector< std::vector< index_t > > stars( mesh.vertices.nb(), temp ) ;
        for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
            for( index_t v = 0; v < mesh.facets.nb_vertices( f ); v++ ) {
                stars[mesh.facets.vertex( f, v )].push_back( f ) ;
            }
        }

        for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
            index_t surface_id = attribute[f] ;
            for( index_t c = mesh.facets.corners_begin( f );
                c < mesh.facets.corners_end( f ); c++ ) {
                index_t f_adj = mesh.facet_corners.adjacent_facet( c ) ;
                if( f_adj != GEO::NO_FACET && attribute[f_adj] != surface_id ) {
                    f_adj = GEO::NO_FACET ;
                }
                if( f_adj == GEO::NO_FACET ) {
                    const std::vector< index_t >& star0 =
                        stars[mesh.facet_corners.vertex( c )] ;
                    const std::vector< index_t >& star1 =
                        stars[mesh.facet_corners.vertex(
                            mesh.facets.next_corner_around_facet( f, c ) )] ;
                    std::vector< index_t > intersect(
                        std::min( star0.size(), star1.size() ) ) ;
                    intersect.erase(
                        std::set_intersection( star0.begin(), star0.end(),
                            star1.begin(), star1.end(), intersect.begin() ),
                        intersect.end() ) ;
                    if( intersect.size() > 1 ) {
                        for( index_t i = 0; i < intersect.size(); i++ ) {
                            index_t cur_f = intersect[i] ;
                            if( cur_f != f && attribute[cur_f] == surface_id ) {
                                f_adj = cur_f ;
                            }
                        }
                    }
                }
                mesh.facet_corners.set_adjacent_facet( c, f_adj ) ;
            }
        }

        /// 2 - Reorient in the same direction using propagation
        std::vector< bool > facet_visited( mesh.facets.nb(), false ) ;
        for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
            if( facet_visited[f] ) continue ;
            index_t surface_id = attribute[f] ;
            std::stack< index_t > S ;
            S.push( f ) ;
            do {
                index_t cur_f = S.top() ;
                S.pop() ;
                if( facet_visited[cur_f] ) continue ;
                facet_visited[cur_f] = true ;
                for( index_t c = mesh.facets.corners_begin( cur_f );
                    c < mesh.facets.corners_end( cur_f ); c++ ) {
                    index_t f_adj = mesh.facet_corners.adjacent_facet( c ) ;
                    if( f_adj == GEO::NO_FACET || attribute[f_adj] != surface_id
                        || facet_visited[f_adj] ) continue ;
                    if( !facets_have_same_orientation( mesh, cur_f, c, f_adj ) ) {
                        mesh.facets.flip( f_adj ) ;
                    }
                    S.push( f_adj ) ;
                }
            } while( !S.empty() ) ;
        }

        /// 3 - Check for consistent orientation with GeoModel
        GEO::MeshFacetsAABB aabb( mesh ) ;
        std::vector< bool > flip_surface( region.model().nb_surfaces(), false ) ;
        bool flip_sthg = false ;
        for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
            const Surface& surface = dynamic_cast< const Surface& >( region.boundary(
                s ) ) ;
            vec3 barycenter = model_element_cell_center( surface, 0 ) ;
            vec3 nearest_point ;
            float64 distance ;
            index_t f = aabb.nearest_facet( barycenter, nearest_point, distance ) ;
            ringmesh_debug_assert( surface.index() == attribute[f] ) ;

            vec3 ori_normal = surface.normal( 0 ) ;
            vec3 test_normal = GEO::Geom::mesh_facet_normal( mesh, f ) ;
            if( dot( ori_normal, test_normal ) < 0 ) {
                flip_surface[surface.index()] = true ;
                flip_sthg = true ;
            }
        }
        if( flip_sthg ) {
            for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                index_t surface_id = attribute[f] ;
                if( flip_surface[surface_id] ) {
                    mesh.facets.flip( f ) ;
                }
            }
        }
    }

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
            GEO::mesh_tetrahedralize( *tetmesh_, false, refine, 1.0 ) ;
            check_and_repair_mesh_consistency( *region_, *tetmesh_ ) ;
            return true ;
        }
    } ;

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

            set_mesh_in_ringmesh() ;


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

        void set_mesh_in_ringmesh()
        {
            initialize_ringmesh_storage() ;
            write_vertices_in_ringmesh_data_structure() ;
            write_tet_in_ringmesh_data_structure() ;
            check_and_repair_mesh_consistency( *region_, *tetmesh_ ) ;
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
                GEO::Logger::err( "TetraGen" )
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
                GEO::Logger::err( "TetraGen" )
                    << "Encountered a problem while meshing inside..." << std::endl ;
                return false ;
            }
            ret = tetra_optimise_volume_regular( tms_ ) ;
            if( ret != STATUS_OK ) {
                GEO::Logger::err( "TetraGen" )
                    << "Encountered a problem while meshing inside..." << std::endl ;
                return false ;
            }
            return true ;

        }

        void initialize_ringmesh_storage()
        {
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

    /*!
     * Creates an instance of the tetrahedral mesher
     * @param[in,out] tetmesh this mesh will be filled with
     * the generated tetrahedral mesh
     * @param[in] algo_name the name of the algorithm to use
     * @return the corresponding instance
     */
    TetraGen* TetraGen::create( GEO::Mesh& tetmesh, const std::string& algo_name )
    {
        TetraGen* mesher = TetraGenFactory::create_object( algo_name ) ;
        if( !mesher ) {
            GEO::Logger::warn( "TetraGen" ) << "Could not create TetraGen mesher: "
                << algo_name << std::endl ;
            GEO::Logger::warn( "TetraGen" ) << "Falling back to TetGen mode"
                << std::endl ;
            mesher = new TetraGen_TetGen() ;
        }

        mesher->set_mesh( tetmesh ) ;
        return mesher ;
    }

    TetraGen::TetraGen()
        : tetmesh_( nil ), region_( nil ), wells_( nil )
    {
    }

    /*!
     * Sets the output tetrahedral mesh to fill
     * @param[out] tetmesh the mesh to tetrahedralize
     */
    void TetraGen::set_mesh( GEO::Mesh& tetmesh )
    {
        tetmesh_ = &tetmesh ;
    }

    /*!
     * Sets the boundaries of the domain
     * @param[in] region The Region of the GeoModel to mesh
     * @param[in] wells the wells to be conformal to
     */
    void TetraGen::set_boundaries(
        const GeoModelElement& region,
        const WellGroup* wells )
    {
        region_ = &region ;
        index_t nb_surfaces = region_->nb_boundaries() ;
        std::vector< const GeoModelMeshElement* > unique_surfaces ;
        unique_surfaces.reserve( nb_surfaces ) ;
        std::vector< index_t > surface_id ;
        surface_id.reserve( nb_surfaces ) ;
        index_t nb_surface_points = 0, nb_facets = 0 ;
        for( index_t s = 0; s < nb_surfaces; s++ ) {
            const Surface& surface =
                dynamic_cast< const Surface& >( region_->boundary( s ) ) ;
            if( contains( surface_id, surface.index() ) ) continue ;
            nb_surface_points += surface.nb_vertices() ;
            nb_facets += surface.nb_cells() ;

            surface_id.push_back( surface.index() ) ;
            unique_surfaces.push_back( &surface ) ;
        }

        MakeUnique uniqueID( unique_surfaces, true ) ;
        std::vector< std::vector< Edge > > well_edges ;
        if( wells ) {
            wells->get_region_edges( region.index(), well_edges ) ;
            // Copy result of porting. Stupid, I know, but because of the interface
            // of MakeUnique. This Edge class is a pain [JP]
            std::vector< std::pair< vec3, vec3 > > wells_copy ;
            for( index_t w = 0; w < well_edges.size(); w++ ) {
                wells_copy.resize( well_edges.size() ) ;
                for( index_t i = 0; i < wells_copy.size(); ++i ) {
                    wells_copy[i] = std::pair< vec3, vec3 >(
                        well_edges[w][i].value( 0 ), well_edges[w][i].value( 1 ) ) ;
                }
                uniqueID.add_edges( wells_copy ) ;
            }
        }
        uniqueID.unique() ;
        const std::vector< index_t >& unique_indices = uniqueID.indices() ;
        std::vector< vec3 > unique_points ;
        uniqueID.unique_points( unique_points ) ;

        index_t starting_index = tetmesh_->vertices.create_vertices(
            unique_points.size() ) ;
        GEO::Memory::copy( tetmesh_->vertices.point_ptr( starting_index ),
            unique_points.data()->data(),
            3 * sizeof(double) * unique_points.size() ) ;

        if( !well_edges.empty() ) {
            index_t nb_well_edges = 0 ;
            for( index_t w = 0; w < well_edges.size(); w++ ) {
                nb_well_edges += well_edges[w].size() ;
            }
            tetmesh_->edges.create_edges( nb_well_edges ) ;
            GEO::Attribute< index_t > edge_region( tetmesh_->edges.attributes(),
                surface_att_name ) ;
            index_t cur_vertex_id = nb_surface_points ;
            index_t cur_edge = 0 ;
            for( index_t w = 0; w < well_edges.size(); w++ ) {
                for( index_t e = 0; e < well_edges[w].size(); e++ ) {
                    tetmesh_->edges.set_vertex( cur_edge, 0,
                        starting_index + unique_indices[cur_vertex_id++ ] ) ;
                    tetmesh_->edges.set_vertex( cur_edge, 1,
                        starting_index + unique_indices[cur_vertex_id++ ] ) ;
                    edge_region[cur_edge++ ] = w ;
                }
            }
        }

        index_t offset_vertices = 0 ;
        index_t offset_facets = 0 ;
        tetmesh_->facets.create_triangles( nb_facets ) ;
        GEO::Attribute< index_t > surface_region( tetmesh_->facets.attributes(),
            surface_att_name ) ;
        for( index_t s = 0; s < unique_surfaces.size(); s++ ) {
            const Surface& surface =
                dynamic_cast< const Surface& >( *unique_surfaces[s] ) ;
            RINGMESH_PARALLEL_LOOP
            for( index_t t = 0; t < surface.nb_cells(); t++ ) {
                ringmesh_debug_assert( surface.is_triangle( t ) ) ;
                for( index_t v = 0; v < 3; v++ ) {
                    tetmesh_->facets.set_vertex( offset_facets + t, v,
                        starting_index
                            + unique_indices[offset_vertices
                                + surface.surf_vertex_id( t, v )] ) ;
                }
                surface_region[offset_facets + t] = surface.index() ;

            }
            offset_vertices += surface.nb_vertices() ;
            offset_facets += surface.nb_cells() ;
        }
        tetmesh_->facets.connect() ;
    }

    /*!
     * Set additional points to be in the output tetrahedral mesh
     * @param[in] points the points to add
     */
    void TetraGen::set_internal_points( const std::vector< vec3 >& points )
    {
        if( points.empty() ) return ;
        index_t start = tetmesh_->vertices.create_vertices( points.size() ) ;
        GEO::Memory::copy( tetmesh_->vertices.point_ptr( start ),
            points.front().data(), points.size() * 3 * sizeof(double) ) ;
    }

    TetraGen::~TetraGen()
    {
    }

    void TetraGen::initialize_storage( index_t nb_points, index_t nb_tets )
    {
        tetmesh_->vertices.clear( true, false ) ;
        tetmesh_->vertices.create_vertices( nb_points ) ;
        tetmesh_->cells.create_tets( nb_tets ) ;
    }

    void TetraGen::set_point( index_t index, const double* point )
    {
        for( index_t i = 0; i < 3; i++ ) {
            tetmesh_->vertices.point_ptr( index )[i] = point[i] ;
        }
    }

    void TetraGen::set_tetra(
        index_t index,
        int* tet,
        index_t nb_lines,
        index_t nb_triangles )
    {
        index_t corner_begin = tetmesh_->cells.corners_begin( index ) ;
        for( index_t v = 0; v < 4; v++ ) {
            tetmesh_->cell_corners.set_vertex( corner_begin++, tet[v] - 1 ) ;
        }
    }

    void TetraGen::initialize()
    {
        ringmesh_register_tetragen( TetraGen_TetGen, "TetGen" ) ;

#ifdef USE_MG_TETRA
        ringmesh_register_tetragen( TetraGen_MG_Tetra, "MG_Tetra" );
#endif
    }
}
