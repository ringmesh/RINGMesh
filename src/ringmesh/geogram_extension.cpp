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

#include <ringmesh/geogram_extension.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/line_stream.h>


namespace RINGMesh {

    using GEO::vec3 ;
        

    // Copied from geo_model_builder.cpp. Still needed ?? check Geogram
    static double read_double( GEO::LineInput& in, index_t field )
    {
        double result ;
        std::istringstream iss( in.field( field ) ) ;
        iss >> result >> std::ws ;
        return result ;
    }
	
    /***********************************************************************/
    /* Loading and saving a GEO::Mesh                                      */

	/** 
	 * @brief TSurfMeshIOHandler for importing .ts files into a mesh.
	 */
	class RINGMESH_API TSurfMeshIOHandler : public GEO::MeshIOHandler {
	public:
		/** @brief creates a TSurfMeshIOHandler
		 */
		TSurfMeshIOHandler(){
		}

		/** @brief TSurfMeshIOHandler destructor
		 */
		~TSurfMeshIOHandler(){
		}

		/** @brief for loading a TSurf saved in .ts format
		 * @warning Assumes there is only one TSurf in the file.
		 * @param filename the name of the file to be processed.
		 * @param mesh the mesh where the surface will be created.
		 * @param flags Some flags, not used for now.
		 */
		virtual bool load( 
                        const std::string& filename,
                        GEO::Mesh& mesh,
                        const GEO::MeshIOFlags& flag = GEO::MeshIOFlags()
                        ) ;
		
		/** @brief for saving a Mesh in .ts format
		 * @todo to be implemented.
		 */
        virtual bool save( const GEO::Mesh& , const std::string& , const GEO::MeshIOFlags& )
        {
            GEO::Logger::err( "I/O" )
                << "Saving a Mesh into TSurf format not implemented yet"
                << std::endl ;
            return false ;
		}

	protected:
		/** @brief loads a .ts
		 * @note we assume there is only one TSURF saved in the file
		 */
		bool load_ts_file( GEO::Mesh& mesh, const std::string& filename ) ;
		
	} ;

	bool TSurfMeshIOHandler::load( const std::string& filename, GEO::Mesh& mesh, const GEO::MeshIOFlags& )
    {
		return load_ts_file( mesh, filename ) ;
    }

    bool TSurfMeshIOHandler::load_ts_file( GEO::Mesh& mesh, const std::string& filename )
    {
        // Count the number of triangles and vertices
        index_t nb_points = 0 ;
        index_t nb_triangles = 0 ;
        int z_sign = 1 ;
        {
            GEO::LineInput in( filename ) ;
            if( !in.OK() ) {
                return false ;
            }
            while( !in.eof() && in.get_line() ) {
                in.get_fields() ;
                if( in.nb_fields() > 0 ) {
                    if( in.field_matches( 0, "ZPOSITIVE" ) ) {
                        if( in.field_matches( 1, "Elevation" ) ) {
                            z_sign = 1 ;
                        } else if( in.field_matches( 1, "Depth" ) ) {
                            z_sign = -1 ;
                        } else {
                            ringmesh_assert_not_reached;
                        }
                    }
                    /// 2.1 Read the surface vertices and facets (only triangles in Gocad Model3d files)
                    else if( in.field_matches( 0,
                        "VRTX" ) || in.field_matches( 0, "PVRTX" ) ) {
                        nb_points++ ;
                    } else if( in.field_matches( 0,
                        "PATOM" ) | in.field_matches( 0, "ATOM" ) ) {
                        nb_points++ ;
                    } else if( in.field_matches( 0, "TRGL" ) ) {
                        nb_triangles++ ;
                    }
                }
            }
        }
        index_t dim = 3 ;
        GEO::vector< double > vertices( dim * nb_points ) ;
        GEO::vector< index_t > triangles( 3 * nb_triangles ) ;
        {
            GEO::LineInput in( filename ) ;
            if( !in.OK() ) {
                return false ;
            }
            index_t v = 0 ;
            index_t t = 0 ;
            while( !in.eof() && in.get_line() ) {
                in.get_fields() ;
                if( in.nb_fields() > 0 ) {
                    /// 2.1 Read the surface vertices and facets (only triangles in Gocad Model3d files)
                    if( in.field_matches( 0, "VRTX" )
                        || in.field_matches( 0, "PVRTX" ) ) {
                        vertices[ dim * v ] = read_double( in, 2 ) ;
                        vertices[ dim * v + 1 ] = read_double( in, 3 ) ;
                        vertices[ dim * v + 2 ] = z_sign * read_double( in, 4 ) ;
                        ++v ;
                    } else if( in.field_matches( 0, "PATOM" )
                               || in.field_matches( 0, "ATOM" ) ) {
                        index_t v0 = in.field_as_uint( 2 ) - 1 ;
                        vertices[ dim * v ] = vertices[ dim * v0 ] ;
                        vertices[ dim * v + 1 ] = vertices[ dim * v0 + 1 ] ;
                        vertices[ dim * v + 2 ] = vertices[ dim * v0 + 2 ] ;
                        ++v ;
                    } else if( in.field_matches( 0, "TRGL" ) ) {
                        triangles[ 3 * t ] = static_cast< index_t >( in.field_as_uint(
                            1 ) - 1 ) ;
                        triangles[ 3 * t + 1 ] =
                            static_cast< index_t >( in.field_as_uint( 2 ) - 1 ) ;
                        triangles[ 3 * t + 2 ] =
                            static_cast< index_t >( in.field_as_uint( 3 ) - 1 ) ;
                        t++ ;
                    }
                }
            }
        }

        mesh.facets.assign_triangle_mesh( dim, vertices, triangles, true ) ;

        GEO::mesh_repair( mesh, GEO::MESH_REPAIR_DEFAULT ) ;
        return true ;
    }
	
	void RINGMESH_API ringmesh_mesh_io_initialize() 
	{
		geo_register_MeshIOHandler_creator( TSurfMeshIOHandler, "ts" );
	}
	
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
                ringmesh_debug_assert( volume > 0 ) ;
                return volume ;
            }
            default:
                ringmesh_assert_not_reached;
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
        vec3 result( 0.0, 0.0, 0.0 ) ;
        double count = 0.0 ;
        for( index_t v = 0; v < M.cells.facet_nb_vertices( cell, f ); ++v ) {
            result += GEO::Geom::mesh_vertex( M,
                                              M.cells.facet_vertex( cell, f, v ) ) ;
            count += 1.0 ;
        }
        return ( 1.0 / count ) * result ;
    }

    /*!
    * Computes the Mesh cell facet normal
    * @param[in] M the mesh
    * @param[in] c the cell index
    * @param[in] f the facet index in the cell
    * @return the cell facet normal
    */
    vec3 mesh_cell_facet_normal( const GEO::Mesh& M, index_t c, index_t f )
    {
        const vec3& p1 = GEO::Geom::mesh_vertex( M,
                                                 M.cells.facet_vertex( c, f, 0 ) ) ;
        const vec3& p2 = GEO::Geom::mesh_vertex( M,
                                                 M.cells.facet_vertex( c, f, 1 ) ) ;
        const vec3& p3 = GEO::Geom::mesh_vertex( M,
                                                 M.cells.facet_vertex( c, f, 2 ) ) ;
        return cross( p2 - p1, p3 - p1 ) ;
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
        int cur = t ;
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
        index_t result = -1 ;
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
            ringmesh_debug_assert( new_coords[ 3 ] == 1. ) ;
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


}