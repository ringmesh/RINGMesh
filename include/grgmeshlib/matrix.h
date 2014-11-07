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

#ifndef __GRGMESH_MATRIX__
#define __GRGMESH_MATRIX__

#include <grgmeshlib/common.h>

#include <vector>
#include <algorithm>

namespace GRGMesh {

    template< class T >
    struct ElementImpl {
        ElementImpl()
            : index( -1 )
        {
        }
        T value ;
        int32 index ;

    } ;

    template< class T >
    class GRGMESH_API RowImpl {
    public:
        typedef ElementImpl< T > Element ;

        RowImpl()
            : nb_elements_( 0 ), capacity_( 4 )
        {
            reallocate( capacity_ ) ;
        }
        ~RowImpl() {
            delete[] elements_ ;
        }

        void set_element( uint32 j, const T& value ) {
            int32 index = find( j ) ;
            if( index == -1 ) {
                if( nb_elements_ == capacity_ ) grow() ;
                Element& elt = elements_[nb_elements_++] ;
                elt.index = j ;
                elt.value = value ;
            } else {
                elements_[index].value = value ;
            }
        }

        int32 find( uint32 j ) const {
            for( uint32 i = 0; i < nb_elements_; i++ ) {
                if( elements_[i].index == j ) return i ;
            }
            return -1 ;
        }

        bool get_element( uint32 j, T& value ) const {
            int32 index = find( j ) ;
            if( index == -1 ) return false ;
            value = elements_[index].value ;
            return true ;
        }

    private:
        void reallocate( uint32 new_capacity ) {
            Element* new_elements = new Element[new_capacity] ;
            std::copy( elements_, elements_ + nb_elements_, new_elements ) ;
            delete[] elements_ ;
            elements_ = new_elements ;
        }
        void grow() {
            capacity_ = capacity_ * 2 ;
            reallocate( capacity_ ) ;
        }

    private:
        Element* elements_ ;
        uint32 nb_elements_ ;
        uint32 capacity_ ;
    } ;

    template< class T >
    class GRGMESH_API SparseMatrix {
    public:
        typedef RowImpl< T > Row ;
        typedef SparseMatrix< T > thisclass ;

        SparseMatrix()
        {
        }
        void build_matrix( uint32 n )
        {
             rows_.resize( n ) ;
        }

        uint32 n() const { return rows_.size() ; }
        const Row& row( uint32 i ) const { return rows_[i] ; }
        void set_element( uint32 i, uint32 j, const T& value ) {
            rows_[i].set_element( j, value ) ;
        }
        bool get_element( uint32 i, uint32 j, T& value ) const {
            return rows_[i].get_element( j, value ) ;
        }
    private:
        SparseMatrix( const thisclass &rhs ) ;
        thisclass& operator=( const thisclass &rhs ) ;

    private:
        std::vector< Row > rows_ ;
    } ;
}

#endif
