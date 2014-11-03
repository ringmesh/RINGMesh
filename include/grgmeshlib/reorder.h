/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#ifndef __GRGMESH_REORDER__
#define __GRGMESH_REORDER__

#include <grgmeshlib/common.h>
#include <grgmeshlib/boundary_model_element.h>

namespace GRGMesh {

    template< class F > inline void parallel_for(
        F& f,
        uint32 from,
        uint32 to )
    {
#pragma omp parallel for
        for( uint32 i = from; i < to; i++ ) {
            f( i ) ;
        }
    }

    template< class IT, class CMP > inline IT split( IT begin, IT end, CMP cmp )
    {
        if( begin >= end ) return begin ;
        IT middle = begin + ( end - begin ) / 2 ;
        std::nth_element( begin, middle, end, cmp ) ;
        return middle ;
    }

    template< int32 COORD, bool UP, class MESH > struct Morton_vertex_cmp {
        Morton_vertex_cmp( const MESH& mesh )
            : mesh_( mesh )
        {
        }
        bool operator()( int32 i1, int32 i2 )
        {
            return mesh_.vertex( i1 )[COORD] < mesh_.vertex( i2 )[COORD] ;
        }
        const MESH& mesh_ ;
    } ;

    template< int32 COORD, bool UP > struct Morton_facet_vertex_cmp {
        Morton_facet_vertex_cmp( const SurfacePart& mesh )
            : mesh_( mesh )
        {
        }
        bool operator()( int32 i1, int32 i2 )
        {
            return mesh_.vertex( i1 )[COORD] < mesh_.vertex( i2 )[COORD] ;
        }
        const SurfacePart& mesh_ ;
    } ;

    template< int32 COORD, bool UP, class MESH > struct Morton_cell_cmp {
    public:
        Morton_cell_cmp( const MESH& mesh )
            : mesh_( mesh )
        {
        }
        float64 center( int32 t ) const
        {
            float64 result = 0.0 ;
            for( uint32 p = 0; p < mesh_.nb_vertices_in_cell( t ); p++ ) {
                result += mesh_.cell_vertex( t, p )[COORD] ;
            }
            return result ;
        }
        bool operator()( int32 t1, int32 t2 )
        {
            return ( center( t1 ) < center( t2 ) ) ;
        }
        const MESH& mesh_ ;
    } ;

    template< int32 COORD, bool UP, class MESH > struct Morton_facet_cmp {
    public:
        Morton_facet_cmp( const MESH& mesh )
            : mesh_( mesh )
        {
        }
        float64 center( int32 t ) const
        {
            float64 result = 0.0 ;
            for( uint32 p = 0; p < 3; p++ ) {
                result += mesh_.point( t, p )[COORD] ;
            }
            return result ;
        }
        bool operator()( int32 t1, int32 t2 )
        {
            return ( center( t1 ) < center( t2 ) ) ;
        }
        const MESH& mesh_ ;
    } ;

    template< class MESH, template< int32 COORD, bool UP, class MESH > class CMP > struct HilbertSort {
        template< int32 COORDX, bool UPX, bool UPY, bool UPZ, class IT >
        static void sort( const MESH& M, IT begin, IT end, uint32 limit = 1 )
        {
            const int32 COORDY = ( COORDX + 1 ) % 3, COORDZ = ( COORDY + 1 ) % 3 ;
            if( end - begin <= int32( limit ) ) return ;
            IT m0 = begin, m8 = end ;
            IT m4 = split( m0, m8, CMP< COORDX, UPX, MESH >( M ) ) ;
            IT m2 = split( m0, m4, CMP< COORDY, UPY, MESH >( M ) ) ;
            IT m1 = split( m0, m2, CMP< COORDZ, UPZ, MESH >( M ) ) ;
            IT m3 = split( m2, m4, CMP< COORDZ, !UPZ, MESH >( M ) ) ;
            IT m6 = split( m4, m8, CMP< COORDY, !UPY, MESH >( M ) ) ;
            IT m5 = split( m4, m6, CMP< COORDZ, UPZ, MESH >( M ) ) ;
            IT m7 = split( m6, m8, CMP< COORDZ, !UPZ, MESH >( M ) ) ;
            sort< COORDZ, UPZ, UPX, UPY >( M, m0, m1 ) ;
            sort< COORDY, UPY, UPZ, UPX >( M, m1, m2 ) ;
            sort< COORDY, UPY, UPZ, UPX >( M, m2, m3 ) ;
            sort< COORDX, UPX, !UPY, !UPZ >( M, m3, m4 ) ;
            sort< COORDX, UPX, !UPY, !UPZ >( M, m4, m5 ) ;
            sort< COORDY, !UPY, UPZ, !UPX >( M, m5, m6 ) ;
            sort< COORDY, !UPY, UPZ, !UPX >( M, m6, m7 ) ;
            sort< COORDZ, !UPZ, !UPX, UPY >( M, m7, m8 ) ;
        }

        HilbertSort(
            const MESH& M,
            std::vector< int32 >& sorted_indices,
            uint32 limit = 1 )
            : M_( M )
        {
            if( sorted_indices.size() <= limit ) return ;
            m0_ = sorted_indices.begin() ;
            m8_ = sorted_indices.end() ;
            m4_ = split( m0_, m8_, CMP< 0, false, MESH >( M ) ) ;
            parallel_for( *this, 0, 2 ) ; // computes m2,m6 in parallel
            parallel_for( *this, 10, 14 ) ; // computes m1,m3,m5,m7 in parallel
            parallel_for( *this, 20, 28 ) ; // sorts the 8 subsets in parallel
        }

        void operator()( int32 i )
        {
            const int32 COORDX = 0, COORDY = 1, COORDZ = 2 ;
            const bool UPX = false, UPY = false, UPZ = false ;
            switch( i ) {
                case 0:
                    m2_ = split( m0_, m4_, CMP< COORDY, UPY, MESH >( M_ ) ) ;
                    break ;
                case 1:
                    m6_ = split( m4_, m8_, CMP< COORDY, !UPY, MESH >( M_ ) ) ;
                    break ;
                case 10:
                    m1_ = split( m0_, m2_, CMP< COORDZ, UPZ, MESH >( M_ ) ) ;
                    break ;
                case 11:
                    m3_ = split( m2_, m4_, CMP< COORDZ, !UPZ, MESH >( M_ ) ) ;
                    break ;
                case 12:
                    m5_ = split( m4_, m6_, CMP< COORDZ, UPZ, MESH >( M_ ) ) ;
                    break ;
                case 13:
                    m7_ = split( m6_, m8_, CMP< COORDZ, !UPZ, MESH >( M_ ) ) ;
                    break ;
                case 20:
                    sort< COORDZ, UPZ, UPX, UPY >( M_, m0_, m1_ ) ;
                    break ;
                case 21:
                    sort< COORDY, UPY, UPZ, UPX >( M_, m1_, m2_ ) ;
                    break ;
                case 22:
                    sort< COORDY, UPY, UPZ, UPX >( M_, m2_, m3_ ) ;
                    break ;
                case 23:
                    sort< COORDX, UPX, !UPY, !UPZ >( M_, m3_, m4_ ) ;
                    break ;
                case 24:
                    sort< COORDX, UPX, !UPY, !UPZ >( M_, m4_, m5_ ) ;
                    break ;
                case 25:
                    sort< COORDY, !UPY, UPZ, !UPX >( M_, m5_, m6_ ) ;
                    break ;
                case 26:
                    sort< COORDY, !UPY, UPZ, !UPX >( M_, m6_, m7_ ) ;
                    break ;
                case 27:
                    sort< COORDZ, !UPZ, !UPX, UPY >( M_, m7_, m8_ ) ;
                    break ;
                default:
                    grgmesh_assert_not_reached ;
                    break ;
            }
        }

    private:
        const MESH& M_ ;
        std::vector< int32 >::iterator m0_, m1_, m2_, m3_, m4_, m5_, m6_, m7_, m8_ ;
    } ;

    inline void morton_vertex_sort(
        const MixedMesh& M,
        std::vector< int32 >& sorted_indices,
        bool use_line,
        bool use_trgl,
        bool use_quad,
        bool use_tetra,
        bool use_pyramid,
        bool use_prism,
        bool use_hexa )
    {        
        std::vector< int32 > used_vertices( M.nb_vertices(), 0 ) ;
        if( use_line ) {
            for( uint32 l = 0; l < M.nb_lines(); l++ ) {
                used_vertices[M.line_vertex_index( l, 0 )]++ ;
                used_vertices[M.line_vertex_index( l, 1 )]++ ;
            }
        }
        if( use_trgl ) {
            for( uint32 t = 0 ; t < M.nb_triangles(); t++ ) {
                used_vertices[M.triangle_vertex_index( t, 0 )]++ ;
                used_vertices[M.triangle_vertex_index( t, 1 )]++ ;
                used_vertices[M.triangle_vertex_index( t, 2 )]++ ;
                used_vertices[M.triangle_vertex_index( t, 3 )]++ ;
            }
        }
        if( use_quad ) {
            for( uint32 q = 0; q < M.nb_quad(); q++ ) {
                used_vertices[M.quad_vertex_index( q, 0 )]++ ;
                used_vertices[M.quad_vertex_index( q, 1 )]++ ;
                used_vertices[M.quad_vertex_index( q, 2 )]++ ;
                used_vertices[M.quad_vertex_index( q, 3 )]++ ;
            }
        }
        if( use_tetra ) {
            for( uint32 t = 0 ; t < M.nb_tetra(); t++ ) {
                used_vertices[M.tetra_vertex_index( t, 0 )]++ ;
                used_vertices[M.tetra_vertex_index( t, 1 )]++ ;
                used_vertices[M.tetra_vertex_index( t, 2 )]++ ;
                used_vertices[M.tetra_vertex_index( t, 3 )]++ ;
            }
        }
        if( use_pyramid ) {
            for( uint32 p = 0; p < M.nb_pyramids(); p++ ) {
                used_vertices[M.pyramid_vertex_index( p, 0 )]++ ;
                used_vertices[M.pyramid_vertex_index( p, 1 )]++ ;
                used_vertices[M.pyramid_vertex_index( p, 2 )]++ ;
                used_vertices[M.pyramid_vertex_index( p, 3 )]++ ;
                used_vertices[M.pyramid_vertex_index( p, 4 )]++ ;
            }
        }
        if( use_prism ) {
            for( uint32 p = 0; p < M.nb_prisms(); p++ ) {
                used_vertices[M.prism_vertex_index( p, 0 )]++ ;
                used_vertices[M.prism_vertex_index( p, 1 )]++ ;
                used_vertices[M.prism_vertex_index( p, 2 )]++ ;
                used_vertices[M.prism_vertex_index( p, 3 )]++ ;
                used_vertices[M.prism_vertex_index( p, 4 )]++ ;
                used_vertices[M.prism_vertex_index( p, 5 )]++ ;
            }
        }
        if( use_hexa ) {
            for( uint32 h = 0; h < M.nb_hexa(); h++ ) {
                used_vertices[M.hexa_vertex_index( h, 0 )]++ ;
                used_vertices[M.hexa_vertex_index( h, 1 )]++ ;
                used_vertices[M.hexa_vertex_index( h, 2 )]++ ;
                used_vertices[M.hexa_vertex_index( h, 3 )]++ ;
                used_vertices[M.hexa_vertex_index( h, 4 )]++ ;
                used_vertices[M.hexa_vertex_index( h, 5 )]++ ;
                used_vertices[M.hexa_vertex_index( h, 6 )]++ ;
                used_vertices[M.hexa_vertex_index( h, 7 )]++ ;
            }
        }

        sorted_indices.clear() ;
        sorted_indices.reserve( M.nb_vertices() ) ;
        for( uint32 i = 0; i < used_vertices.size(); i++ ) {
            if( used_vertices[i] > 0 ) sorted_indices.push_back( i ) ; 
        }
        HilbertSort< MixedMesh, Morton_vertex_cmp >( M, sorted_indices ) ;
    }

    inline void morton_cell_sort(
        const MixedMesh& M,
        std::vector< int32 >& sorted_indices,
        bool use_line,
        bool use_trgl,
        bool use_quad,
        bool use_tetra,
        bool use_pyramid,
        bool use_prism,
        bool use_hexa )
    {
        uint32 nb_points = 0 ;
        if( use_line )    nb_points +=  M.nb_lines() ;
        if( use_trgl )    nb_points +=  M.nb_triangles() ;
        if( use_quad )    nb_points +=  M.nb_quad() ;
        if( use_tetra )   nb_points +=  M.nb_tetra() ;
        if( use_pyramid ) nb_points +=  M.nb_pyramids() ;
        if( use_prism )   nb_points +=  M.nb_prisms() ;
        if( use_hexa )    nb_points +=  M.nb_hexa() ;
        sorted_indices.resize( nb_points ) ;

        uint32 i = 0 ;
        uint32 offset = 0 ;
        if( use_line ) {
            for( uint32 l = 0; l < M.nb_lines(); l++ ) {
                sorted_indices[i++] = offset + l ;
            }
        }
        offset += M.nb_lines() ;
        if( use_trgl ) {
            for( uint32 t = 0 ; t < M.nb_triangles(); t++ ) {
                sorted_indices[i++] = offset + t ;
            }
        }
        offset += M.nb_triangles() ;
        if( use_quad ) {
            for( uint32 q = 0; q < M.nb_quad(); q++ ) {
                sorted_indices[i++] = offset + q ;
            }
        }
        offset += M.nb_quad() ;
        if( use_tetra ) {
            for( uint32 t = 0 ; t < M.nb_tetra(); t++ ) {
                sorted_indices[i++] = offset + t ;
            }
        }
        offset += M.nb_tetra() ;
        if( use_pyramid ) {
            for( uint32 p = 0; p < M.nb_pyramids(); p++ ) {
                sorted_indices[i++] = offset + p ;
            }
        }
        offset += M.nb_pyramids() ;
        if( use_prism ) {
            for( uint32 p = 0; p < M.nb_prisms(); p++ ) {
                sorted_indices[i++] = offset + p ;
            }
        }
        offset += M.nb_prisms() ;
        if( use_hexa ) {
            for( uint32 h = 0; h < M.nb_hexa(); h++ ) {
                sorted_indices[i++] = offset + h ;
            }
        }

        HilbertSort< MixedMesh, Morton_cell_cmp >( M, sorted_indices ) ;
    }

    inline void morton_vertex_sort(
        const SurfacePart& M,
        std::vector< int32 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_vertices() ) ;
        for( uint32 i = 0; i < M.nb_vertices(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< SurfacePart, Morton_vertex_cmp >( M,
            sorted_indices ) ;
    }

    inline void morton_cell_sort(
        const SurfacePart& M,
        std::vector< int32 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_simplices() ) ;
        for( uint32 i = 0; i < M.nb_simplices(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< SurfacePart, Morton_facet_cmp >( M, sorted_indices ) ;
    }
}

#endif
