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

#include <ringmesh/geogram_extension.h>

#include <geogram/basic/line_stream.h>
#include <geogram/basic/logger.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>

#include <geogram/points/colocate.h>

/*!
 * @todo Re-orgarnize this mess
 */

namespace RINGMesh {

    using GEO::vec3 ;
    using GEO::Mesh ;
        

    /***********************************************************************/
    /* Loading and saving a GEO::Mesh                                      */

    /*! 
     * @brief TSurfMeshIOHandler for importing .ts files into a mesh.
     */
    class TSurfMeshIOHandler : public GEO::MeshIOHandler {
    public:   
        TSurfMeshIOHandler() :
            mesh_dimension_(3),
            nb_vertices_(0), 
            nb_triangles_(0),
            z_sign_(1)
        {
        }
        
        /*!
         * @brief Load a TSurf saved in .ts format
         * @warning Assumes there is only one TSurf in the file.
         *          Will undoubtedly crash if it is not the case.
         * @todo Prevent crashing if the file is not as expected
         *
         * @param filename the name of the .ts file to be processed.
         * @param mesh to which vertices and facets will be assigned
         * @param flags Not used for now.
         */
        virtual bool load(
            const std::string& filename,
            GEO::Mesh& mesh,
            const GEO::MeshIOFlags& flag = GEO::MeshIOFlags()
            )
        {
            filename_ = filename ;
            if( !is_file_valid()) {
                return false ;
            } else {
                read_number_of_vertices_and_triangles() ;
                allocate_vertices() ;
                allocate_triangles() ;
                read_vertices_and_triangles() ;
                assign_and_repair_mesh( mesh ) ;
                return true ;
            }
        }

        /*!
         * @brief Save a Mesh in .ts format
         * @todo To be implemented.
         */
        virtual bool save( const GEO::Mesh&, const std::string&, const GEO::MeshIOFlags& )
        {
            GEO::Logger::err( "I/O" )
                << "Saving a Mesh into TSurf format not implemented yet"
                << std::endl ;
            return false ;
        }

    private:
        // This function read the z_sign too [PA]
        void read_number_of_vertices_and_triangles()
        {
            GEO::LineInput in( filename_ ) ;
            while( !in.eof() && in.get_line() ) {
                in.get_fields() ;
                if( in.nb_fields() > 0 ) {
                    if( in.field_matches( 0, "ZPOSITIVE" ) ) {
                        if( in.field_matches( 1, "Elevation" ) ) {
                            z_sign_ = 1 ;
                        } else if( in.field_matches( 1, "Depth" ) ) {
                            z_sign_ = -1 ;
                        }
                    }
                    else if( in.field_matches( 0, "VRTX" ) || in.field_matches( 0, "PVRTX" ) ) {
                        nb_vertices_++ ;
                    } else if( in.field_matches( 0, "PATOM" ) || in.field_matches( 0, "ATOM" ) ) {
                        nb_vertices_++ ;
                    } else if( in.field_matches( 0, "TRGL" ) ) {
                        nb_triangles_++ ;
                    }
                }
            }
        }
        
        void read_vertices_and_triangles()
        {
            GEO::LineInput in( filename_ ) ;
            index_t v = 0 ;
            index_t t = 0 ;
            while( !in.eof() && in.get_line() ) {
                in.get_fields() ;
                if( in.nb_fields() > 0 ) {
                    if( in.field_matches( 0, "VRTX" ) || in.field_matches( 0, "PVRTX" ) ) {
                        vertices_[ mesh_dimension_*v ]     =  in.field_as_double(2);
                        vertices_[ mesh_dimension_*v + 1 ] = in.field_as_double(3);
                        vertices_[ mesh_dimension_*v + 2 ] = in.field_as_double(4) * z_sign_ ;
                        ++v ;
                    } else if( in.field_matches( 0, "PATOM" ) || in.field_matches( 0, "ATOM" ) ) {
                        index_t v0 = in.field_as_uint( 2 ) - 1 ;
                        vertices_[ mesh_dimension_*v ]     = vertices_[ mesh_dimension_*v0 ] ;
                        vertices_[ mesh_dimension_*v + 1 ] = vertices_[ mesh_dimension_*v0 + 1 ] ;
                        vertices_[ mesh_dimension_*v + 2 ] = vertices_[ mesh_dimension_*v0 + 2 ] ;
                        ++v ;
                    } else if( in.field_matches( 0, "TRGL" ) ) {
                        triangles_[ 3*t ]     = index_t( in.field_as_uint(1)-1 ) ;
                        triangles_[ 3*t + 1 ] = index_t( in.field_as_uint(2)-1 ) ;
                        triangles_[ 3*t + 2 ] = index_t( in.field_as_uint(3)-1 ) ;
                        t++ ;
                    }
                }
            }
        }

        void assign_and_repair_mesh( GEO::Mesh& mesh )
        {
            mesh.facets.assign_triangle_mesh( mesh_dimension_, vertices_, triangles_, true ) ;
            // Do not use GEO::MESH_REPAIR_DEFAULT because it glues the 
            // disconnected edges along internal boundaries
            GEO::mesh_repair( mesh, GEO::MESH_REPAIR_DUP_F ) ;
        }
        
        bool is_file_valid()
        {
            GEO::LineInput in( filename_ ) ;
            if( !in.OK() ) {
                return false ;
            } else {
                return true ;
            }
        }

        void allocate_vertices()
        {
            vertices_.resize( mesh_dimension_ * nb_vertices_ ) ;
        }
        
        void allocate_triangles()
        {
            triangles_.resize( 3 * nb_triangles_ ) ;
        }

    private:
        index_t mesh_dimension_ ;
        index_t nb_vertices_ ;
        index_t nb_triangles_ ;
        int z_sign_ ;
        std::string filename_ ;
        GEO::vector< double > vertices_ ;
        GEO::vector< index_t > triangles_ ;
    } ;

    class LINMeshIOHandler: public GEO::MeshIOHandler {
    public:
        virtual bool load(
            const std::string& filename,
            GEO::Mesh& mesh,
            const GEO::MeshIOFlags& flag = GEO::MeshIOFlags() )
        {
            GEO::LineInput file( filename ) ;

            while( !file.eof() && file.get_line() ) {
                file.get_fields() ;
                if( file.nb_fields() > 0 ) {
                    if( file.field_matches( 0, "v" ) ) {
                        vec3 vertex = load_vertex( file, 1 ) ;
                        mesh.vertices.create_vertex( vertex.data() ) ;
                    } else if( file.field_matches( 0, "s" ) ) {
                        mesh.edges.create_edge(
                            file.field_as_uint( 1 ) - 1,
                            file.field_as_uint( 2 ) - 1 ) ;
                    }
                }
            }
            return true ;

        }
        virtual bool save(
            const GEO::Mesh& M,
            const std::string& filename,
            const GEO::MeshIOFlags& ioflags = GEO::MeshIOFlags() )
        {
            throw RINGMeshException( "I/O",
                "Saving a Mesh into .lin format not implemented yet" ) ;
            return false ;
        }

    private:
        vec3 load_vertex( GEO::LineInput& file, index_t field ) const
        {
            double x = file.field_as_double( field++ ) ;
            double y = file.field_as_double( field++ ) ;
            double z = file.field_as_double( field++ ) ;
            return vec3( x, y, z ) ;
        }
    } ;

    void ringmesh_mesh_io_initialize()
    {
        geo_register_MeshIOHandler_creator( TSurfMeshIOHandler, "ts" ) ;
        geo_register_MeshIOHandler_creator( LINMeshIOHandler, "lin" ) ;
    }


    /***********************************************************************/
#ifdef RINGMESH_WITH_TETGEN
#ifdef __GNUC__
#   pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

    bool is_mesh_tetrahedralizable( const GEO::Mesh& M ) 
    {         
        if( M.facets.nb() == 0 ) {
            GEO::Logger::err( "RING" ) << "Mesh to tetrahedralize has no facets "
                << std::endl ;
            return false ;
        }
        if( !M.facets.are_simplices() ) {
            GEO::Logger::err( "RING" ) << "Mesh to tetrahedralize is not triangulated"
                << std::endl ;
            return false ;
        }
        if( M.cells.nb() != 0 ) {
            GEO::Logger::warn( "RING" ) << "Mesh to tetrahedralize already have cells"
                << std::endl ;         
        }
        return true ;
    }


    TetgenMesher::~TetgenMesher()
    {
        // Take over facet deletion of tetgen that does not set to 
        // nil pointers to polygonlist or holelist in facet
        delete[] tetgen_in_.facetlist ;
        tetgen_in_.facetlist = nil ;
        tetgen_in_.numberoffacets = 0 ;

        delete[] polygon_corners_ ;
        polygon_corners_ = nil ;

        delete[] polygons_ ;
        polygons_ = nil ;
    }

    void TetgenMesher::tetrahedralize( const Mesh& input_mesh, 
                                       const std::string& command_line,
                                       Mesh& output_mesh )
    {
        set_command_line( command_line ) ;
        initialize() ;
        copy_mesh_to_tetgen_input( input_mesh ) ;
        tetrahedralize() ;
        assign_result_tetmesh_to_mesh( output_mesh ) ;
    }

    void  TetgenMesher::tetrahedralize( const GEO::Mesh& input_mesh,
                                        const std::vector< vec3 >& one_point_per_region,
                                        const std::string& command_line,
                                        GEO::Mesh& output_mesh ) 
    {
        set_command_line( command_line ) ;
        initialize() ;
        copy_mesh_to_tetgen_input( input_mesh ) ;
        set_regions( one_point_per_region ) ;
        tetrahedralize() ;
        assign_result_tetmesh_to_mesh( output_mesh ) ;
        fill_region_attribute_on_mesh_cells( output_mesh, "region" ) ;
    }



    void TetgenMesher::initialize()
    {
        initialize_tetgen_args() ;
        tetgen_in_.initialize() ;
        tetgen_out_.initialize() ;
    }

    void TetgenMesher::set_command_line( const std::string& command_line )
    {
        tetgen_command_line_ = command_line.c_str() ;
    }

    void TetgenMesher::tetrahedralize()
    {
        try {
            GEO_3rdParty::tetrahedralize( &tetgen_args_, &tetgen_in_,
                &tetgen_out_ ) ;
        } catch( int code ) {
            GEO::Logger::err( "Tetgen" ) << "Encountered a problem: " ;
            switch( code ) {
                case 1:
                    GEO::Logger::err( "Tetgen" ) << "Out of memory" ;
                    break ;
                case 2:
                    GEO::Logger::err( "Tetgen" )
                        << "Please report this bug to Hang.Si@wias-berlin.de. Include\n" ;
                    GEO::Logger::err( "Tetgen" )
                        << "  the message above, your input data set, and the exact\n" ;
                    GEO::Logger::err( "Tetgen" )
                        << "  command line you used to run this program, thank you" ;
                    break ;
                case 3:
                    GEO::Logger::err( "Tetgen" )
                        << "A self-intersection was detected. Program stopped\n" ;
                    GEO::Logger::err( "Tetgen" )
                        << "Hint: use -d option to detect all self-intersections" ;
                    break ;
                case 4:
                    GEO::Logger::err( "Tetgen" )
                        << "A very small input feature size was detected. Program stopped.\n" ;
                    GEO::Logger::err( "Tetgen" )
                        << "Hint: use -T option to set a smaller tolerance." ;
                    break ;
                case 5:
                    GEO::Logger::err( "Tetgen" )
                        << "Two very close input facets were detected. Program stopped.\n" ;
                    GEO::Logger::err( "Tetgen" )
                        << "Hint: use -Y option to avoid adding Steiner points in boundary." ;
                    break ;
                case 10:
                    GEO::Logger::err( "Tetgen" )
                        << "An input error was detected. Program stopped." ;
                    break ;
            }
            GEO::Logger::err( "Tetgen" ) << std::endl ;
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
        tetgen_in_.numberofpoints = static_cast<int>(M.vertices.nb()) ;
        tetgen_in_.pointlist = new double[ 3 * tetgen_in_.numberofpoints ] ;
        GEO::Memory::copy(
            tetgen_in_.pointlist, M.vertices.point_ptr( 0 ),
            M.vertices.nb() * 3 * sizeof( double )
            ) ;
    }

    void TetgenMesher::copy_edges_to_tetgen_input( const GEO::Mesh& M )
    {
        tetgen_in_.numberofedges = static_cast<int>(M.edges.nb()) ;
        tetgen_in_.edgelist = new int[ 2 * tetgen_in_.numberofedges ] ;
        GEO::Memory::copy(
            tetgen_in_.edgelist, M.edges.vertex_index_ptr( 0 ),
            M.edges.nb() * 2 * sizeof( int )
            ) ;
    }

    void TetgenMesher::copy_facets_to_tetgen_input( const GEO::Mesh& M )
    {
        polygons_ = new GEO_3rdParty::tetgenio::polygon[ M.facets.nb() ] ;

        tetgen_in_.numberoffacets = static_cast<int>(M.facets.nb()) ;
        tetgen_in_.facetlist = new GEO_3rdParty::tetgenio::facet[ tetgen_in_.numberoffacets ] ;

        polygon_corners_ = new int[ M.facet_corners.nb() ] ;
        GEO::Memory::copy(
            polygon_corners_, M.facet_corners.vertex_index_ptr( 0 ),
            M.facet_corners.nb()*sizeof( int )
            ) ;

        for( index_t f = 0; f < M.facets.nb(); ++f ) {
            GEO_3rdParty::tetgenio::facet& F = tetgen_in_.facetlist[ f ] ;
            GEO_3rdParty::tetgenio::init( &F ) ;
            F.numberofpolygons = 1 ;
            F.polygonlist = &polygons_[ f ] ;

            GEO_3rdParty::tetgenio::polygon& P = F.polygonlist[ 0 ] ;
            GEO_3rdParty::tetgenio::init( &P ) ;
            P.numberofvertices = static_cast< int >( M.facets.nb_corners( f ) ) ;
            P.vertexlist = &polygon_corners_[ M.facets.corners_begin( f ) ] ;
        }
    }

    void TetgenMesher::set_regions( const std::vector< vec3 >& one_point_in_each_region ) 
    {
        index_t nb_regions = one_point_in_each_region.size() ;
        tetgen_in_.numberofregions = static_cast< int >( nb_regions ) ;
        tetgen_in_.regionlist = new double[5*nb_regions] ;

        for( index_t i = 0; i != nb_regions; ++i ){
            tetgen_in_.regionlist[5*i] = one_point_in_each_region[i].x ;
            tetgen_in_.regionlist[5*i+1] = one_point_in_each_region[i].y ;
            tetgen_in_.regionlist[5*i+2] = one_point_in_each_region[i].z ;
            tetgen_in_.regionlist[5*i+3] = i ;
            tetgen_in_.regionlist[5*i+4] = DBL_MAX ; // Used only with the a switch
        }
    }

    void TetgenMesher::fill_region_attribute_on_mesh_cells( GEO::Mesh& M, const std::string& attribute_name ) const
    {
        double* tet_attributes = tetgen_out_.tetrahedronattributelist ;
        int one_tet_attribute_size = tetgen_out_.numberoftetrahedronattributes ;
        GEO::Attribute< index_t > region_id( M.cells.attributes(), attribute_name ) ;
        for( index_t i = 0; i < M.cells.nb(); ++i ) {
            // Nothing says where it is, so we hope that the shell id is the first 
            // attribute stored in tetgen [JP]
            region_id[ i ] = tet_attributes[ one_tet_attribute_size*i ] ;
        }
        region_id.unbind() ;
    }

    void TetgenMesher::initialize_tetgen_args()
    {
        char* copy = new char[ tetgen_command_line_.length() + 1 ] ;
        std::strcpy( copy, tetgen_command_line_.c_str() ) ;
        tetgen_args_.parse_commandline( copy ) ;
    }

    void TetgenMesher::assign_result_tetmesh_to_mesh( GEO::Mesh& M ) const
    {
        GEO::vector<double> points ;
        get_result_tetmesh_points( points ) ;

        GEO::vector<index_t> tets ;
        get_result_tetmesh_tets( tets ) ;

        M.cells.assign_tet_mesh( 3, points, tets, true ) ;
        M.cells.connect() ;
    }

    void TetgenMesher::get_result_tetmesh_points( GEO::vector< double >& points ) const
    {
        index_t nb_points = static_cast< index_t >( tetgen_out_.numberofpoints ) ;
        points.resize( 3 * nb_points ) ;
        double* points_ptr = tetgen_out_.pointlist ;
        RINGMESH_PARALLEL_LOOP
        for( index_t i = 0; i < 3 * nb_points; ++i ) {
            points[ i ] = points_ptr[ i ] ;
        }
    }

    void TetgenMesher::get_result_tetmesh_tets( GEO::vector< index_t>& tets ) const
    {
        index_t nb_tets = static_cast< index_t >( tetgen_out_.numberoftetrahedra );
        tets.resize( 4 * nb_tets );

        int* tets_ptr = tetgen_out_.tetrahedronlist ;
        int one_tet_size = tetgen_out_.numberofcorners ;
        RINGMESH_PARALLEL_LOOP
        for( index_t i = 0; i < nb_tets; ++i ) {
            tets[ 4 * i + 0 ] = tets_ptr[ one_tet_size*i + 0 ] ;
            tets[ 4 * i + 1 ] = tets_ptr[ one_tet_size*i + 1 ] ;
            tets[ 4 * i + 2 ] = tets_ptr[ one_tet_size*i + 2 ] ;
            tets[ 4 * i + 3 ] = tets_ptr[ one_tet_size*i + 3 ] ;
        }
    }

    void tetrahedralize_mesh_tetgen( GEO::Mesh& M, bool refine, double quality )
    {
        if( !is_mesh_tetrahedralizable( M ) ) {
            throw RINGMeshException( "TetGen", "Mesh cannot be tetrahedralized" ) ;
        }               
        TetgenMesher mesher ;
        if( refine ) {
            mesher.tetrahedralize( M, "QpYAq" + GEO::String::to_string( quality ),
                M ) ;
        } else {
            mesher.tetrahedralize( M, "QpYYA", M ) ;
        }
    }

#ifdef __GNUC__
#   pragma GCC diagnostic warning "-Wsign-conversion"
#endif
#endif
    
    /***********************************************************************/

    /*!
    * Computes the volume of a Mesh cell
    * @param[in] M the mesh
    * @param[in] c the cell index
    * @return the volume of the cell
    */
    double mesh_cell_volume( const GEO::Mesh& M, index_t c )
    {
        switch( M.cells.type( c ) ) {
            case GEO::MESH_TET:
                return GEO::Geom::tetra_volume(
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 1 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 3 ) ) ) ;
            case GEO::MESH_PYRAMID:
                return GEO::Geom::tetra_volume(
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 1 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 4 ) ) )
                    + GEO::Geom::tetra_volume(
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 3 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 4 ) ) ) ;
            case GEO::MESH_PRISM:
            case GEO::MESH_HEX: {
                vec3 ori( 0, 0, 0 ) ;
                double volume = 0 ;
                for( index_t f = 0; f < M.cells.nb_facets( c ); f++ ) {
                    switch( M.cells.facet_nb_vertices( c, f ) ) {
                        case 3:
                            volume += GEO::Geom::tetra_signed_volume(
                                GEO::Geom::mesh_vertex( M,
                                M.cells.facet_vertex( c, f, 0 ) ),
                                GEO::Geom::mesh_vertex( M,
                                M.cells.facet_vertex( c, f, 1 ) ),
                                GEO::Geom::mesh_vertex( M,
                                M.cells.facet_vertex( c, f, 2 ) ), ori ) ;
                            break ;
                        case 4:
                            volume += GEO::Geom::tetra_signed_volume(
                                GEO::Geom::mesh_vertex( M,
                                M.cells.facet_vertex( c, f, 0 ) ),
                                GEO::Geom::mesh_vertex( M,
                                M.cells.facet_vertex( c, f, 1 ) ),
                                GEO::Geom::mesh_vertex( M,
                                M.cells.facet_vertex( c, f, 2 ) ), ori ) ;
                            volume += GEO::Geom::tetra_signed_volume(
                                GEO::Geom::mesh_vertex( M,
                                M.cells.facet_vertex( c, f, 0 ) ),
                                GEO::Geom::mesh_vertex( M,
                                M.cells.facet_vertex( c, f, 2 ) ),
                                GEO::Geom::mesh_vertex( M,
                                M.cells.facet_vertex( c, f, 3 ) ), ori ) ;
                            break ;
                        default:
                            ringmesh_assert_not_reached;
                            return 0 ;
                    }
                }
                ringmesh_assert( volume > 0 ) ;
                return volume ;
            }
            default:
                return 0 ;
        }
    }

    /*!
    * Computes the Mesh cell facet barycenter
    * @param[in] M the mesh
    * @param[in] cell the cell index
    * @param[in] f the facet index in the cell
    * @return the cell facet center
    */
    vec3 mesh_cell_facet_center( const GEO::Mesh& M, index_t cell, index_t f )
    {
        vec3 result( 0., 0., 0. ) ;
        index_t nb_vertices = M.cells.facet_nb_vertices( cell, f ) ;        
        for( index_t v = 0; v < nb_vertices; ++v ) {
            result += GEO::Geom::mesh_vertex( M, M.cells.facet_vertex( cell, f, v ) ) ;
        }
        ringmesh_assert( nb_vertices > 0 );

        return result/ nb_vertices ;
    }

    /*!
    * Computes the non weighted barycenter of a volumetric 
    * cell of a Mesh
    * @param[in] M the mesh
    * @param[in] cell the cell index
    * @return the cell center
    */
    vec3 mesh_cell_center( const GEO::Mesh& M, index_t cell )
    {
        vec3 result( 0.0, 0.0, 0.0 ) ;
        double count = 0.0 ;
        for( index_t v = 0; v < M.cells.nb_vertices( cell ); ++v ) {
            result += GEO::Geom::mesh_vertex( M, M.cells.vertex( cell, v ) ) ;
            count += 1.0 ;
        }
        return ( 1.0 / count ) * result ;
    }


    /*!
    * Tests if a tetrahedron has an egde between two given points
    * @param[in] mesh the mesh
    * @param[in] t Tetrahedron index
    * @param[in] p0 First vertex index
    * @param[in] p1 Second vertex index
    * @param[out] edge Output edge index
    * @return The result of the test
    */
    bool has_edge(
        const GEO::Mesh& mesh,
        index_t t,
        index_t p0,
        index_t p1,
        index_t& edge )
    {
        using GEO::Numeric::uint8 ;
        for( uint8 e = 0; e < 6; e++ ) {
            index_t v0 = mesh.cells.edge_vertex( t, e, 0 ) ;
            index_t v1 = mesh.cells.edge_vertex( t, e, 1 ) ;
            if( ( p0 == v0 && p1 == v1 ) || ( p0 == v1 && p1 == v0 ) ) {
                edge = e ;
                return true ;
            }
        }
        return false ;
    }

    /*!
    * Gets all the next adjacent tetrahedra sharing an edge
    * @param[in] mesh the mesh
    * @param[in] t Starting tetrahedron index to test, should contain the edge
    * @param[in] prev Previous tetrahedron index
    * (if propagation around the edge, prevent to go back were we came from)
    * @param[in] p0 First vertex index of the edge
    * @param[in] p1 Second vertex index of the edge
    * @return The edge index
    * \pre the mesh needs to be tetrahedralized
    */
    index_t next_around_edge(
        const GEO::Mesh& mesh,
        index_t t,
        index_t prev,
        index_t p0,
        index_t p1 )
    {
        for( index_t adj = 0; adj < mesh.cells.nb_facets( t ); adj++ ) {
            index_t t_adj = mesh.cells.adjacent( t, adj ) ;
            if( t_adj == GEO::NO_CELL || t_adj == prev ) continue ;
            index_t edge ;
            if( has_edge( mesh, t_adj, p0, p1, edge ) ) {
                //todo handles any cell type
                return 6 * t_adj + edge ;
            }
        }
        return GEO::NO_CELL ;
    }

    /*!
    * Gets all the edge indices around one edge
    * @param[in] mesh the mesh
    * @param[in] t First tetrahedron index to test, should include the edge
    * @param[in] p0 First vertex index of the edge
    * @param[in] p1 Second vertex index of the edge
    * @param[out] result Output list of edge indices
    */
    void edges_around_edge(
        const GEO::Mesh& mesh,
        index_t t,
        index_t p0,
        index_t p1,
        std::vector< index_t >& result )
    {
        index_t prev = t ;
        index_t cur = t ;
        do {
            index_t info = next_around_edge( mesh, cur, prev, p0, p1 ) ;
            if( info == GEO::NO_CELL ) return ;
            result.push_back( info ) ;
            prev = cur ;
            cur = info / 6 ;
        } while( cur != t ) ;
    }

    /*!
    * Get vertices when an edge is divide into \p nb_parts parts
    * @param[in] mesh the mesh
    * @param[in] edge the edge id in \p mesh
    * @param[in] nb_parts the number of edge division
    * @param[out] points the points which divide the edge
    */
    void divide_edge_in_parts(
        const GEO::Mesh& mesh,
        index_t edge,
        index_t nb_parts,
        std::vector< vec3 >& points )
    {
        points.resize( nb_parts - 1 ) ;
        double pond = 1. / nb_parts ;
        vec3 node0 = GEO::Geom::mesh_vertex( mesh, mesh.edges.vertex( edge, 0 ) ) ;
        vec3 node1 = GEO::Geom::mesh_vertex( mesh, mesh.edges.vertex( edge, 1 ) ) ;
        for( index_t i = 0; i < nb_parts - 1; i++ ) {
            for( index_t j = 0; j < 3; j++ ) {
                points[ i ][ j ] = ( i + 1 ) * pond * node1[ j ]
                    + ( 1. - ( i + 1 ) * pond ) * node0[ j ] ;
            }
        }

    }

    void divide_edge_in_parts(
        vec3& node0,
        vec3& node1,
        index_t nb_parts,
        std::vector< vec3 >& points )
    {
        if( nb_parts > 0 ) {
            points.resize( nb_parts - 1 ) ;
            double pond = 1. / nb_parts ;
            for( index_t i = 0; i < nb_parts - 1; i++ ) {
                for( index_t j = 0; j < 3; j++ ) {
                    points[ i ][ j ] = ( i + 1 ) * pond * node1[ j ]
                        + ( 1. - ( i + 1 ) * pond ) * node0[ j ] ;
                }
            }
        }
    }
    /*!
    * Gets the closest local vertex index in a mesh cell of a point
    * @param[in] mesh the mesh
    * @param[in] p the point to test
    * @param[in] t the cell index
    * @return the local vertex index
    */
    index_t get_nearest_vertex_index(
        const GEO::Mesh& mesh,
        const vec3& p,
        index_t t )
    {
        float64 dist = GEO::Numeric::max_float64() ;
        index_t result = NO_ID ;
        for( index_t v = 0; v < mesh.cells.nb_vertices( t ); v++ ) {
            float64 distance = length2(
                GEO::Geom::mesh_vertex( mesh, mesh.cells.vertex( t, v ) ) - p ) ;
            if( distance < dist ) {
                result = v ;
            }
        }
        return result ;
    }




    /*!
    * Rotation of all the vertices of a mesh following
    * a defined rotational matrix.
    *
    * @param mesh[in,out] the mesh to rotate.
    *
    * @param[in] rot_mat matrix which defines the rotation.
    */
    void rotate_mesh(
        GEO::Mesh& mesh,
        const GEO::Matrix< float64, 4 >& rot_mat )
    {
        for( index_t v = 0; v < mesh.vertices.nb(); v++ ) {
            float64 old_coords[ 4 ] ;
            for( index_t i = 0; i < 3; i++ ) {
                old_coords[ i ] = mesh.vertices.point_ptr( v )[ i ] ;
            }
            old_coords[ 3 ] = 1. ;
            float64 new_coords[ 4 ] ;
            GEO::mult( rot_mat, old_coords, new_coords ) ;

            for( index_t i = 0; i < 3; i++ ) {
                mesh.vertices.point_ptr( v )[ i ] = new_coords[ i ] ;
            }
            ringmesh_assert( new_coords[ 3 ] == 1. ) ;
        }
    }



    void print_bounded_attributes( const GEO::Mesh& M )
    {
        {
            GEO::vector< std::string > names ;
            M.vertices.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size(), false ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    if( names[ a ] == "point" ) continue ;
                    is_bounded[ a ] = M.vertices.attributes().find_attribute_store(
                        names[ a ] )->has_observers() ;
                    if( is_bounded[ a ] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on vertices:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[ a ] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[ a ] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.edges.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[ a ] = M.edges.attributes().find_attribute_store(
                        names[ a ] )->has_observers() ;
                    if( is_bounded[ a ] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on edges:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[ a ] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[ a ] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.facets.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[ a ] = M.facets.attributes().find_attribute_store(
                        names[ a ] )->has_observers() ;
                    if( is_bounded[ a ] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on facets:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[ a ] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[ a ] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.facet_corners.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[ a ] =
                        M.facet_corners.attributes().find_attribute_store( names[ a ] )->has_observers() ;
                    if( is_bounded[ a ] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on facet_corners:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[ a ] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[ a ] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.cells.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[ a ] = M.cells.attributes().find_attribute_store(
                        names[ a ] )->has_observers() ;
                    if( is_bounded[ a ] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on cells:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[ a ] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[ a ] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.cell_corners.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[ a ] = M.cell_corners.attributes().find_attribute_store(
                        names[ a ] )->has_observers() ;
                    if( is_bounded[ a ] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on cell_corners:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[ a ] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[ a ] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.cell_facets.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[ a ] = M.cell_facets.attributes().find_attribute_store(
                        names[ a ] )->has_observers() ;
                    if( is_bounded[ a ] ) failed = true ;
                }
                if( failed ) {
                    GEO::Logger::err( "Attributes" )
                        << "Attributes still bounded on cell_facets:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[ a ] ) {
                            GEO::Logger::err( "Attributes" ) << " " << names[ a ] ;
                        }
                    }
                    GEO::Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
    }



    /*
    * @note Code modified from geogram/mesh/mesh_repair.cpp
    */
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
        index_t nb_colocated_vertices( M.vertices.nb() - nb_unique_vertices) ;
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

    void update_mesh_edges_vertices( GEO::Mesh& M, const GEO::vector<index_t>& old2new)
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
