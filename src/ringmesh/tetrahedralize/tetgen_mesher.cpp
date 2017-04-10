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

#include <ringmesh/tetrahedralize/tetgen_mesher.h>

#include <cstring>

#include <geogram/mesh/mesh.h>
#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_builder.h>

#ifdef RINGMESH_WITH_TETGEN
/*!
 * @file Implementation of the interface of GEO::Mesh with Tetgen
 * @author Jeanne Pellerin
 */

namespace {
    using namespace RINGMesh ;

    bool is_mesh_tetrahedralizable( const GEO::Mesh& M )
    {
        if( M.facets.nb() == 0 ) {
            Logger::err( "RING", "Mesh to tetrahedralize has no facets " ) ;
            return false ;
        }
        if( !M.facets.are_simplices() ) {
            Logger::err( "RING", "Mesh to tetrahedralize is not triangulated" ) ;
            return false ;
        }
        if( M.cells.nb() != 0 ) {
            Logger::warn( "RING", "Mesh to tetrahedralize already have cells" ) ;
        }
        return true ;
    }
}

namespace RINGMesh {

    TetgenMesher::~TetgenMesher()
    {
        // Take over facet deletion of tetgen that does not set to 
        // nullptr pointers to polygonlist or holelist in facet
        delete[] tetgen_in_.facetlist ;
        tetgen_in_.facetlist = nullptr ;
        tetgen_in_.numberoffacets = 0 ;

        delete[] polygon_corners_ ;
        polygon_corners_ = nullptr ;

        delete[] polygons_ ;
        polygons_ = nullptr ;
    }

    void TetgenMesher::tetrahedralize(
        const GEO::Mesh& input_mesh,
        Mesh3DBuilder& output_mesh_builder )
    {
        initialize() ;
        copy_mesh_to_tetgen_input( input_mesh ) ;
        tetrahedralize() ;
        assign_result_tetmesh_to_mesh( output_mesh_builder ) ;
    }

    void TetgenMesher::tetrahedralize(
        const Mesh0D& input_mesh,
        Mesh3DBuilder& output_mesh_builder )
    {
        initialize() ;
        copy_vertices_to_tetgen_input( input_mesh ) ;
        tetrahedralize() ;
        assign_result_tetmesh_to_mesh( output_mesh_builder ) ;
    }

    void TetgenMesher::initialize()
    {
        initialize_tetgen_args() ;
        tetgen_in_.initialize() ;
        tetgen_out_.initialize() ;
    }

    void TetgenMesher::tetrahedralize()
    {
        try {
            GEO_3rdParty::tetrahedralize( &tetgen_args_, &tetgen_in_,
                &tetgen_out_ ) ;
        } catch( int code ) {
            Logger::err( "Tetgen", "Encountered a problem: " ) ;
            switch( code ) {
                case 1:
                    Logger::err( "Tetgen", "Out of memory" ) ;
                    break ;
                case 2:
                    Logger::err( "Tetgen",
                        "Please report this bug to Hang.Si@wias-berlin.de. Include" ) ;
                    Logger::err( "Tetgen",
                        "  the message above, your input data set, and the exact" ) ;
                    Logger::err( "Tetgen",
                        "  command line you used to run this program, thank you" ) ;
                    ;
                    break ;
                case 3:
                    Logger::err( "Tetgen",
                        "A self-intersection was detected. Program stopped" ) ;
                    Logger::err( "Tetgen",
                        "Hint: use -d option to detect all self-intersections" ) ;
                    break ;
                case 4:
                    Logger::err( "Tetgen",
                        "A very small input feature size was detected. Program stopped." ) ;
                    Logger::err( "Tetgen",
                        "Hint: use -T option to set a smaller tolerance." ) ;
                    break ;
                case 5:
                    Logger::err( "Tetgen",
                        "Two very close input facets were detected. Program stopped." ) ;
                    Logger::err( "Tetgen",
                        "Hint: use -Y option to avoid adding Steiner points in boundary." ) ;
                    break ;
                case 10:
                    Logger::err( "Tetgen",
                        "An input error was detected. Program stopped." ) ;
                    break ;
            }
        }
    }

    void TetgenMesher::copy_mesh_to_tetgen_input( const GEO::Mesh& M )
    {
        if( M.vertices.nb() != 0 ) {
            copy_vertices_to_tetgen_input( M ) ;
        }
        if( M.edges.nb() != 0 ) {
            copy_edges_to_tetgen_input( M ) ;
        }
        if( M.facets.nb() != 0 ) {
            copy_facets_to_tetgen_input( M ) ;
        }
    }

    void TetgenMesher::copy_vertices_to_tetgen_input( const GEO::Mesh& M )
    {
        tetgen_in_.numberofpoints = static_cast< int >( M.vertices.nb() ) ;
        tetgen_in_.pointlist = new double[3 * tetgen_in_.numberofpoints] ;
        GEO::Memory::copy( tetgen_in_.pointlist, M.vertices.point_ptr( 0 ),
            M.vertices.nb() * 3 * sizeof(double) ) ;
    }

    void TetgenMesher::copy_vertices_to_tetgen_input( const Mesh0D& M )
    {
        if( M.nb_vertices() != 0 ) {
            tetgen_in_.numberofpoints = static_cast< int >( M.nb_vertices() ) ;
            tetgen_in_.pointlist = new double[M.nb_vertices() * 3] ;
            for( index_t v = 0; v < M.nb_vertices(); v++ ) {
                for( index_t i = 0; i < 3; i++ ) {
                    tetgen_in_.pointlist[3 * v + 1] = M.vertex( v )[i] ;
                }
            }
        }
    }

    void TetgenMesher::copy_edges_to_tetgen_input( const GEO::Mesh& M )
    {
        tetgen_in_.numberofedges = static_cast< int >( M.edges.nb() ) ;
        tetgen_in_.edgelist = new int[2 * tetgen_in_.numberofedges] ;
        GEO::Memory::copy( tetgen_in_.edgelist, M.edges.vertex_index_ptr( 0 ),
            M.edges.nb() * 2 * sizeof(int) ) ;
    }

    void TetgenMesher::copy_facets_to_tetgen_input( const GEO::Mesh& M )
    {
        polygons_ = new GEO_3rdParty::tetgenio::polygon[M.facets.nb()] ;

        tetgen_in_.numberoffacets = static_cast< int >( M.facets.nb() ) ;
        tetgen_in_.facetlist =
            new GEO_3rdParty::tetgenio::facet[tetgen_in_.numberoffacets] ;

        polygon_corners_ = new int[M.facet_corners.nb()] ;
        GEO::Memory::copy( polygon_corners_, M.facet_corners.vertex_index_ptr( 0 ),
            M.facet_corners.nb() * sizeof(int) ) ;

        for( index_t f = 0; f < M.facets.nb(); ++f ) {
            GEO_3rdParty::tetgenio::facet& F = tetgen_in_.facetlist[f] ;
            GEO_3rdParty::tetgenio::init( &F ) ;
            F.numberofpolygons = 1 ;
            F.polygonlist = &polygons_[f] ;

            GEO_3rdParty::tetgenio::polygon& P = F.polygonlist[0] ;
            GEO_3rdParty::tetgenio::init( &P ) ;
            P.numberofvertices = static_cast< int >( M.facets.nb_corners( f ) ) ;
            P.vertexlist = &polygon_corners_[M.facets.corners_begin( f )] ;
        }
    }

    void TetgenMesher::set_regions(
        const std::vector< vec3 >& one_point_in_each_region )
    {
        index_t nb_regions = one_point_in_each_region.size() ;
        tetgen_in_.numberofregions = static_cast< int >( nb_regions ) ;
        tetgen_in_.regionlist = new double[5 * nb_regions] ;

        for( index_t i = 0; i != nb_regions; ++i ) {
            tetgen_in_.regionlist[5 * i] = one_point_in_each_region[i].x ;
            tetgen_in_.regionlist[5 * i + 1] = one_point_in_each_region[i].y ;
            tetgen_in_.regionlist[5 * i + 2] = one_point_in_each_region[i].z ;
            tetgen_in_.regionlist[5 * i + 3] = i ;
            tetgen_in_.regionlist[5 * i + 4] = DBL_MAX ; // Used only with the a switch
        }
    }

    void TetgenMesher::initialize_tetgen_args()
    {
        char* copy = new char[tetgen_command_line_.length() + 1] ;
        std::strcpy( copy, tetgen_command_line_.c_str() ) ;
        tetgen_args_.parse_commandline( copy ) ;
        delete[] copy ;
    }

    void TetgenMesher::assign_result_tetmesh_to_mesh(
        Mesh3DBuilder& output_mesh_builder ) const
    {
        output_mesh_builder.assign_vertices( get_result_tetmesh_points() ) ;
        output_mesh_builder.assign_cell_tet_mesh( get_result_tetmesh_tets() ) ;
        output_mesh_builder.remove_isolated_vertices() ;
        output_mesh_builder.connect_cells() ;
    }

    std::vector< double > TetgenMesher::get_result_tetmesh_points() const
    {
        index_t nb_points = static_cast< index_t >( tetgen_out_.numberofpoints ) ;
        std::vector< double > points( 3 * nb_points ) ;
        double* points_ptr = tetgen_out_.pointlist ;
        RINGMESH_PARALLEL_LOOP
        for( index_t i = 0; i < 3 * nb_points; ++i ) {
            points[i] = points_ptr[i] ;
        }
        return points ;
    }

    std::vector< index_t > TetgenMesher::get_result_tetmesh_tets() const
    {
        std::vector< index_t > tets_to_keep = determine_tets_to_keep() ;

        index_t nb_tets = static_cast< index_t >( tets_to_keep.size() ) ;
        std::vector< index_t > tets( 4 * nb_tets ) ;
        int* tets_ptr = tetgen_out_.tetrahedronlist ;
        RINGMESH_PARALLEL_LOOP
        for( index_t i = 0; i < nb_tets; ++i ) {
            index_t tetra = tets_to_keep[i] ;
            for( index_t v = 0; v < 4; v++ ) {
                tets[4 * i + v] = static_cast< index_t >( tets_ptr[4 * tetra + v] ) ;
            }
        }
        return tets ;
    }

    std::set< double > TetgenMesher::determine_tet_regions_to_keep() const
    {
        // Determine which regions are incident to
        // the 'exterior' (neighbor = -1).
        // The region Id of tet t is determined by:
        //  tetgen_out_.tetrahedronattributelist[t]
        std::set< double > regions_to_keep ;
        index_t nb_tets = static_cast< index_t >( tetgen_out_.numberoftetrahedra ) ;
        for( index_t t = 0; t < nb_tets; ++t ) {
            for( index_t f = 0; f < 4; ++f ) {
                signed_index_t n = tetgen_out_.neighborlist[t * 4 + f] ;
                if( n == -1 ) {
                    regions_to_keep.insert(
                        tetgen_out_.tetrahedronattributelist[t] ) ;
                    break ;
                }
            }
        }
        return regions_to_keep ;
    }

    std::vector< index_t > TetgenMesher::determine_tets_to_keep() const
    {
        std::vector< index_t > tets_to_keep ;
        std::set< double > regions_to_keep = determine_tet_regions_to_keep() ;

        index_t nb_tets = static_cast< index_t >( tetgen_out_.numberoftetrahedra ) ;
        tets_to_keep.reserve( nb_tets ) ;
        for( index_t t = 0; t < nb_tets; ++t ) {
            if( regions_to_keep.find( tetgen_out_.tetrahedronattributelist[t] )
                != regions_to_keep.end() ) {
                tets_to_keep.push_back( t ) ;
            }
        }
        return tets_to_keep ;
    }

    void tetrahedralize_mesh_tetgen(
        Mesh3DBuilder& out_tet_mesh,
        const GEO::Mesh& in_mesh,
        bool refine,
        double quality )
    {
        if( !is_mesh_tetrahedralizable( in_mesh ) ) {
            throw RINGMeshException( "TetGen", "Mesh cannot be tetrahedralized" ) ;
        }
        TetgenMesher mesher ;
        if( refine ) {
            mesher.add_points_to_match_quality( quality ) ;
        }
        mesher.tetrahedralize( in_mesh, out_tet_mesh ) ;
    }

    void tetrahedralize_mesh_tetgen(
        Mesh3DBuilder& out_tet_mesh,
        const Mesh0D& in_point_cloud,
        bool refine,
        double quality )
    {
        TetgenMesher mesher ( ) ;
        if( refine ) {
            mesher.add_points_to_match_quality( quality ) ;
        }
        mesher.tetrahedralize( in_point_cloud, out_tet_mesh ) ;

    }

    void TetgenMesher::add_points_to_match_quality( double quality )
    {
        tetgen_command_line_ += "q" + GEO::String::to_string( quality ) ;
    }
}

#endif

