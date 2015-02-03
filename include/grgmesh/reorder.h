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

#include <grgmesh/common.h>
#include <grgmesh/boundary_model_element.h>

#include <geogram/mesh/mesh.h>

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

    template< int32 COORD, bool UP > struct Morton_vertex_cmp {
        Morton_vertex_cmp( const Surface& mesh )
            : mesh_( mesh )
        {
        }
        bool operator()( int32 i1, int32 i2 ) const
        {
            return mesh_.vertex( i1 )[COORD] < mesh_.vertex( i2 )[COORD] ;
        }
        const Surface& mesh_ ;
    } ;


    template< int32 COORD, bool UP > struct Morton_mesh_vertex_cmp {
        Morton_mesh_vertex_cmp( const GEO::Mesh& mesh )
            : mesh_( mesh )
        {
        }
        bool operator()( int32 i1, int32 i2 ) const
        {
            return mesh_.vertex_ptr( i1 )[COORD] < mesh_.vertex_ptr( i2 )[COORD] ;
        }
        const GEO::Mesh& mesh_ ;
    } ;

    template< int32 COORD, bool UP > struct Morton_facet_vertex_cmp {
        Morton_facet_vertex_cmp( const Surface& mesh )
            : mesh_( mesh )
        {
        }
        bool operator()( int32 i1, int32 i2 ) const
        {
            return mesh_.vertex( i1 )[COORD] < mesh_.vertex( i2 )[COORD] ;
        }
        const Surface& mesh_ ;
    } ;

    template< int32 COORD, bool UP > struct Morton_cell_cmp {
    public:
        Morton_cell_cmp( const GEO::Mesh& mesh )
            : mesh_( mesh )
        {
        }
        float64 center( int32 t ) const
        {
            float64 result = 0.0 ;
            for( uint32 p = 0; p < mesh_.cell_nb_vertices( t ); p++ ) {
                result += mesh_.vertex_ptr( mesh_.cell_vertex_index( t, p ) )[COORD] ;
            }
            return result ;
        }
        bool operator()( int32 t1, int32 t2 ) const
        {
            return ( center( t1 ) < center( t2 ) ) ;
        }
        const GEO::Mesh& mesh_ ;
    } ;

    template< int32 COORD, bool UP > struct Morton_facet_cmp {
    public:
        Morton_facet_cmp( const Surface& mesh )
            : mesh_( mesh )
        {
        }
        float64 center( int32 t ) const
        {
            float64 result = 0.0 ;
            for( uint32 p = 0; p < 3; p++ ) {
                result += mesh_.vertex( t, p )[COORD] ;
            }
            return result ;
        }
        bool operator()( int32 t1, int32 t2 ) const
        {
            return ( center( t1 ) < center( t2 ) ) ;
        }
        const Surface& mesh_ ;
    } ;

    template< class MESH, template< int32 COORD, bool UP > class CMP > struct HilbertSort {
        template< int32 COORDX, bool UPX, bool UPY, bool UPZ, class IT >
        static void sort( const MESH& M, IT begin, IT end, uint32 limit = 1 )
        {
            const int32 COORDY = ( COORDX + 1 ) % 3, COORDZ = ( COORDY + 1 ) % 3 ;
            if( end - begin <= int32( limit ) ) return ;
            IT m0 = begin, m8 = end ;
            IT m4 = split( m0, m8, CMP< COORDX, UPX >( M ) ) ;
            IT m2 = split( m0, m4, CMP< COORDY, UPY >( M ) ) ;
            IT m1 = split( m0, m2, CMP< COORDZ, UPZ >( M ) ) ;
            IT m3 = split( m2, m4, CMP< COORDZ, !UPZ >( M ) ) ;
            IT m6 = split( m4, m8, CMP< COORDY, !UPY >( M ) ) ;
            IT m5 = split( m4, m6, CMP< COORDZ, UPZ >( M ) ) ;
            IT m7 = split( m6, m8, CMP< COORDZ, !UPZ >( M ) ) ;
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
            m4_ = split( m0_, m8_, CMP< 0, false >( M ) ) ;
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
                    m2_ = split( m0_, m4_, CMP< COORDY, UPY >( M_ ) ) ;
                    break ;
                case 1:
                    m6_ = split( m4_, m8_, CMP< COORDY, !UPY >( M_ ) ) ;
                    break ;
                case 10:
                    m1_ = split( m0_, m2_, CMP< COORDZ, UPZ >( M_ ) ) ;
                    break ;
                case 11:
                    m3_ = split( m2_, m4_, CMP< COORDZ, !UPZ >( M_ ) ) ;
                    break ;
                case 12:
                    m5_ = split( m4_, m6_, CMP< COORDZ, UPZ >( M_ ) ) ;
                    break ;
                case 13:
                    m7_ = split( m6_, m8_, CMP< COORDZ, !UPZ >( M_ ) ) ;
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
        const Surface& M,
        std::vector< int32 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_vertices() ) ;
        for( uint32 i = 0; i < M.nb_vertices(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< Surface, Morton_vertex_cmp >( M,
            sorted_indices ) ;
    }

    inline void morton_cell_sort(
        const Surface& M,
        std::vector< int32 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_cells() ) ;
        for( uint32 i = 0; i < M.nb_cells(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< Surface, Morton_facet_cmp >( M, sorted_indices ) ;
    }


    inline void morton_vertex_sort(
        const GEO::Mesh& M,
        std::vector< int32 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_vertices() ) ;
        for( uint32 i = 0; i < M.nb_vertices(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< GEO::Mesh, Morton_mesh_vertex_cmp >( M,
            sorted_indices ) ;
    }

    inline void morton_cell_sort(
        const GEO::Mesh& M,
        std::vector< int32 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_cells() ) ;
        for( uint32 i = 0; i < M.nb_cells(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< GEO::Mesh, Morton_cell_cmp >( M, sorted_indices ) ;
    }
}

#endif
