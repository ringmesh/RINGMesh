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

#ifndef __GRGMESH_UTILS__
#define __GRGMESH_UTILS__

#include <grgmeshlib/common.h>
#include <grgmeshlib/vecn.h>
#include <ANN/ANN.h>

#include <iostream>
#include <sstream>
#include <vector>

namespace GRGMesh {

    class Edge ;
    class SurfacePart ;
    class ContactPart ;
    class MixedMesh ;
}


namespace GRGMesh {

    class Box3d {
    public:
        Box3d()
            :
                initialized_( false ),
                x_min_( 1e30 ),
                y_min_( 1e30 ),
                z_min_( 1e30 ),
                x_max_( -1e30 ),
                y_max_( -1e30 ),
                z_max_( -1e30 )
        {
        }
        bool initialized() const
        {
            return initialized_ ;
        }
        void clear()
        {
            initialized_ = false ;
        }
        float64 x_min() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return x_min_ ;
        }
        float64 y_min() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return y_min_ ;
        }
        float64 z_min() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return z_min_ ;
        }
        float64 x_max() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return x_max_ ;
        }
        float64 y_max() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return y_max_ ;
        }
        float64 z_max() const
        {
            grgmesh_debug_assert( initialized_ ) ;
            return z_max_ ;
        }
        float64 min( unsigned axis ) const
        {
            return ( axis == 0 ) ? x_min_ : ( ( axis == 1 ) ? y_min_ : z_min_ ) ;
        }
        float64 max( unsigned axis ) const
        {
            return ( axis == 0 ) ? x_max_ : ( ( axis == 1 ) ? y_max_ : z_max_ ) ;
        }
        float64 width() const
        {
            return x_max() - x_min() ;
        }
        float64 height() const
        {
            return y_max() - y_min() ;
        }
        float64 depth() const
        {
            return z_max() - z_min() ;
        }
        vec3 center() const
        {
            return vec3( 0.5 * ( x_max() + x_min() ),
                0.5 * ( y_max() + y_min() ), 0.5 * ( z_max() + z_min() ) ) ;
        }
        double radius() const
        {
            return 0.5
                * ::sqrt(
                    ::sqrt( x_max() - x_min() ) + ::sqrt( y_max() - y_min() )
                        + ::sqrt( z_max() - z_min() ) ) ;
        }
        void add_point( const vec3& p )
        {
            if( !initialized_ ) {
                x_min_ = p[0] ;
                y_min_ = p[1] ;
                z_min_ = p[2] ;
                x_max_ = p[0] ;
                y_max_ = p[1] ;
                z_max_ = p[2] ;
                initialized_ = true ;
            } else {
                x_min_ = std::min( x_min_, p[0] ) ;
                y_min_ = std::min( y_min_, p[1] ) ;
                z_min_ = std::min( z_min_, p[2] ) ;
                x_max_ = std::min( x_max_, p[0] ) ;
                y_max_ = std::min( y_max_, p[1] ) ;
                z_max_ = std::min( z_max_, p[2] ) ;
            }
        }
        void add_box( const Box3d& b )
        {
            if( b.initialized() ) {
                add_point( vec3( b.x_min(), b.y_min(), b.z_min() ) ) ;
                add_point( vec3( b.x_max(), b.y_max(), b.z_max() ) ) ;
            }
        }

    private:
        bool initialized_ ;
        float64 x_min_ ;
        float64 y_min_ ;
        float64 z_min_ ;
        float64 x_max_ ;
        float64 y_max_ ;
        float64 z_max_ ;
    } ;

    namespace Utils {
        inline double triangle_area( const vec3& p1, const vec3& p2, const vec3& p3 )
        {
            return 0.5 * length( cross( p2 - p1, p3 - p1 ) ) ;
        }
        template< class T > static bool contains(
            const std::vector< T >& v,
            const T& t )
        {
            return find( v, t ) != -1 ;
        }
        template< class T > static int find( const std::vector< T >& v, const T& t )
        {
            for( unsigned int i = 0; i < v.size(); i++ ) {
                if( v[i] == t ) return i ;
            }
            return -1 ;
        }
        template< class T1, class T2 > static bool inexact_equal( const T1& v1, const T2& v2 ) {
            for( unsigned int i = 0; i < 3; i++ ) {
                double diff( v1[i]-v2[i] ) ;
                if( diff > epsilon || diff < -epsilon ) {
                    return false ;
                }
            }
            return true ;
        }
        template< class T1, class T2 > static bool triple_equal(
            const T1& rhs1,
            const T1& rhs2,
            const T1& rhs3,
            const T2& lhs1,
            const T2& lhs2,
            const T2& lhs3 )
        {
            if( rhs1 == lhs1 ) {
                if( rhs2 == lhs2 && rhs3 == lhs3 ) {
                    return true ;
                } else if( rhs2 == lhs3 && rhs3 == lhs2 ) {
                    return true ;
                }
            } else if( rhs1 == lhs2 ) {
                if( rhs2 == lhs1 && rhs3 == lhs3 ) {
                    return true ;
                } else if( rhs2 == lhs3 && rhs3 == lhs1 ) {
                    return true ;
                }
            } else if( rhs1 == lhs3 ) {
                if( rhs2 == lhs1 && rhs3 == lhs2 ) {
                    return true ;
                } else if( rhs2 == lhs2 && rhs3 == lhs1 ) {
                    return true ;
                }
            }
            return false ;
        }
    }

    class InputStream {
    public:
        InputStream(
            std::istream& in
        ) : in_(in), line_in_(nil) {   }
        ~InputStream() {
            delete line_in_ ; line_in_ = nil ;
        }
        bool eof() const { return in_.eof() ; }
        bool eol() const { return line_in_ == nil || line_in_->eof() ; }
        bool ok() const { return in_ != 0; }

        void get_line() {
            in_.getline(buffer_, 65536) ;
            bool check_multiline = true ;
            int total_length = 65536 ;
            char* ptr = buffer_ ;

            // If the line ends with a backslash, append
            // the next line to the current line.
            while(check_multiline) {
                int L = (int)strlen(ptr) ;
                total_length -= L ;
                ptr = ptr + L - 2;
                if(*ptr == '\\' && total_length > 0) {
                    *ptr = ' ' ;
                    ptr++ ;
                    in_.getline(ptr, total_length) ;
                } else {
                    check_multiline = false ;
                }
            }

            if(total_length < 0) {
                std::cerr << "MultiLine longer than 65536 bytes" << std::endl ;
            }

            delete line_in_ ; line_in_ = new std::istringstream(buffer_) ;
        }

        std::istream& line() {
            grgmesh_assert(line_in_ != nil) ;
            return *line_in_ ;
        }

        const char *current_line() const {
            return buffer_;
        }

        template <class T> InputStream& operator>>(T& param) {
            param = T() ; // reset do default value, in case *line_in_ is EOF.
            *line_in_ >> param;
            return *this;
        }

    private:
        std::istream& in_ ;
        std::istringstream* line_in_ ;
        char buffer_[65536] ;
    } ;

    class GRGMESH_API MakeUnique {
    public:
        MakeUnique( const std::vector< vec3 >& data ) ;
        template< class T > MakeUnique( const std::vector< T >& data )
        {
            int nb_points = 0 ;
            for( unsigned int i = 0; i < data.size(); i++ ) {
                nb_points += data[i].points().size() ;
            }
            points_.resize( nb_points ) ;
            indices_.resize( nb_points ) ;
            int cur_id = 0 ;
            for( unsigned int i = 0; i < data.size(); i++ ) {
                for( unsigned int p = 0; p < data[i].points().size(); p++, cur_id++ ) {
                    points_[cur_id] = data[i].points()[p] ;
                    indices_[cur_id] = cur_id ;
                }
            }
        }
        template< class T > MakeUnique( const std::vector< T >& data, bool T_is_a_pointer )
        {
            int nb_points = 0 ;
            for( unsigned int i = 0; i < data.size(); i++ ) {
                nb_points += data[i]->nb_points() ;
            }
            points_.resize( nb_points ) ;
            indices_.resize( nb_points ) ;
            int cur_id = 0 ;
            for( unsigned int i = 0; i < data.size(); i++ ) {
                for( unsigned int p = 0; p < data[i]->nb_points(); p++, cur_id++ ) {
                    points_[cur_id] = data[i]->point( p ) ;
                    indices_[cur_id] = cur_id ;
                }
            }
        }

        void add_edges( const std::vector< Edge >& points ) ;
        void add_points( const std::vector< vec3 >& points ) ;

        void unique( int nb_neighbors = 5 ) ;

        const std::vector< vec3 >& points() const { return points_ ; }
        void unique_points( std::vector< vec3 >& results ) const ;
        const std::vector< int >& indices() const { return indices_ ; }

    private:
        std::vector< vec3 > points_ ;
        std::vector< int > indices_ ;
    } ;

    class GRGMESH_API ColocaterANN {
    public:
        ColocaterANN( const SurfacePart& mesh ) ;
        ColocaterANN( const ContactPart& mesh ) ;
        ColocaterANN( const MixedMesh& mesh ) ;
        ColocaterANN( const MixedMesh& mesh, bool use_surface_id ) ;
        ColocaterANN( const std::vector< vec3 >& vertices ) ;
        ColocaterANN( const std::vector< Edge >& edges ) ;

        ~ColocaterANN() {
            annDeallocPts( ann_points_ ) ;
            delete ann_tree_;
            annClose() ;
        }

        void get_mapped_colocated(
            vec3& v,
            std::vector< unsigned int >& result,
            int nb_neighbors = 2 ) const ;

        bool get_colocated(
            const vec3& v,
            std::vector< unsigned int >& result,
            int nb_neighbors = 2 ) const {
            return get_colocated( const_cast< vec3& >( v ), result, nb_neighbors ) ;
        }
        bool get_colocated(
            vec3& v,
            std::vector< unsigned int >& result,
            int nb_neighbors = 2 ) const ;
        void get_neighbors(
            const vec3& v,
            std::vector< int >& result,
            int nb_neighbors = 2 ) const {
            return get_neighbors( const_cast< vec3& >( v ), result, nb_neighbors ) ;
        }
        void get_neighbors(
            vec3& v,
            std::vector< int >& result,
            int nb_neighbors = 2 ) const ;

        vec3 point( int i ) const {
            return vec3( ann_points_[i] ) ;
        }

    private:
        std::vector< int > mapped_indices_ ;
        ANNpointArray ann_points_ ;
        ANNkd_tree* ann_tree_ ;
    } ;


    template< class T, int n >
    class GRGMESH_API Array {
    public:
        void assign( const std::vector< T >& values ) {
            grgmesh_debug_assert( values.size() < n+1 ) ;
            for( unsigned int i = 0; i < values.size(); i++ ) {
                values_[i] = values[i] ;
            }
        }

        T value( int i ) const { return values_[i] ; }
        T& value( int i ) { return values_[i] ; }

        double normalized_value( int i, double max, double min, double scale ) const {
            double s = values_[i] ;
            s = (s - min) / (max - min) ;
            s = std::min(s, 1.0) ;
            s = std::max(s, 0.0) ;
            s *= scale ;
        return s ;
        }

    protected:
        T values_[n] ;
    } ;

    template< int n >
    class GRGMESH_API intArrayTmpl: public Array< int, n > {
    public:
        intArrayTmpl() {
            for( unsigned int i = 0; i < n; i++ ) {
                this->values_[i] = -1 ;
            }
        }
    } ;
    typedef intArrayTmpl< 6 > intArray ;
    typedef intArrayTmpl< 12 > edgeArray ;

    template< int n >
    class GRGMESH_API boolArrayTmpl: public Array< bool, n > {
    public:
        boolArrayTmpl() {
            for( unsigned int i = 0; i < n; i++ ) {
                this->values_[i] = false ;
            }
        }
    } ;
    typedef boolArrayTmpl< 6 > boolArray ;

    class GRGMESH_API Edge: public Array< vec3, 2 > {
    public:
        Edge( const vec3& v0, const vec3& v1 ) {
            values_[0] = v0 ;
            values_[1] = v1 ;
        }
        vec3 barycenter() const {
            return (values_[0]+values_[1]) / static_cast< double >( 2 ) ;
        }
    } ;

}

#endif
