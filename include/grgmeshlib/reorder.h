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
        unsigned int from,
        unsigned int to )
    {
#pragma omp parallel for
        for( uint64 i = from; i < to; i++ ) {
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

    /*
    template< int COORD, bool UP > struct Morton_vertex_cmp {
        Morton_vertex_cmp( const TetraMesh& mesh )
            : mesh_( mesh )
        {
        }
        bool operator()( int i1, int i2 )
        {
            return mesh_.vertex( i1 )[COORD] < mesh_.vertex( i2 )[COORD] ;
        }
        const TetraMesh& mesh_ ;
    } ;
    */
    template< int COORD, bool UP > struct Morton_facet_vertex_cmp {
        Morton_facet_vertex_cmp( const GRGMesh::SurfacePart& mesh )
            : mesh_( mesh )
        {
        }
        bool operator()( int i1, int i2 )
        {
            return mesh_.point( i1 )[COORD] < mesh_.point( i2 )[COORD] ;
        }
        const GRGMesh::SurfacePart& mesh_ ;
    } ;
/*
    template< int COORD, bool UP > struct Morton_tet_cmp {
    public:
        Morton_tet_cmp( const TetraMesh& mesh )
            : mesh_( mesh )
        {
        }
        double center( int t ) const
        {
            double result = 0.0 ;
            for( unsigned int p = 0; p < 4; p++ ) {
                result += mesh_.vertex( t, p )[COORD] ;
            }
            return result ;
        }
        bool operator()( int t1, int t2 )
        {
            return ( center( t1 ) < center( t2 ) ) ;
        }
        const TetraMesh& mesh_ ;
    } ;
*/
    template< int COORD, bool UP > struct Morton_facet_cmp {
    public:
        Morton_facet_cmp( const GRGMesh::SurfacePart& mesh )
            : mesh_( mesh )
        {
        }
        double center( int t ) const
        {
            double result = 0.0 ;
            for( unsigned int p = 0; p < 3; p++ ) {
                result += mesh_.point( t, p )[COORD] ;
            }
            return result ;
        }
        bool operator()( int t1, int t2 )
        {
            return ( center( t1 ) < center( t2 ) ) ;
        }
        const GRGMesh::SurfacePart& mesh_ ;
    } ;

    template< class MESH, template< int COORD, bool UP > class CMP > struct HilbertSort {
        template< int COORDX, bool UPX, bool UPY, bool UPZ, class IT >
        static void sort( const MESH& M, IT begin, IT end, unsigned int limit = 1 )
        {
            const int COORDY = ( COORDX + 1 ) % 3, COORDZ = ( COORDY + 1 ) % 3 ;
            if( end - begin <= int( limit ) ) return ;
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
            std::vector< int64 >& sorted_indices,
            uint64 limit = 1 )
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

        void operator()( int i )
        {
            const int COORDX = 0, COORDY = 1, COORDZ = 2 ;
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
        std::vector< int64 >::iterator m0_, m1_, m2_, m3_, m4_, m5_, m6_, m7_, m8_ ;
    } ;
/*
    inline void morton_vertex_sort(
        const TetraMesh& M,
        std::vector< int64 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_points() ) ;
        for( uint64 i = 0; i < M.nb_points(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< TetraMesh, Morton_vertex_cmp >( M, sorted_indices ) ;
    }

    inline void morton_tet_sort(
        const TetraMesh& M,
        std::vector< int64 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_tetra() ) ;
        for( uint64 i = 0; i < M.nb_tetra(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< TetraMesh, Morton_tet_cmp >( M, sorted_indices ) ;
    }
*/
    inline void morton_facet_vertex_sort(
        const GRGMesh::SurfacePart& M,
        std::vector< int64 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_points() ) ;
        for( uint64 i = 0; i < M.nb_points(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< SurfacePart, Morton_facet_vertex_cmp >( M,
            sorted_indices ) ;
    }

    inline void morton_facet_sort(
        const GRGMesh::SurfacePart& M,
        std::vector< int64 >& sorted_indices )
    {
        sorted_indices.resize( M.nb_simplices() ) ;
        for( uint64 i = 0; i < M.nb_simplices(); i++ ) {
            sorted_indices[i] = i ;
        }
        HilbertSort< SurfacePart, Morton_facet_cmp >( M, sorted_indices ) ;
    }
}

#endif
