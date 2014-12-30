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

#include <grgmesh/common.h>

namespace GRGMesh {

    template< class T >
    struct ElementImpl {
        ElementImpl()
            : index( -1 )
        {
        }
        T value ;
        signed_index_t index ;

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

        void set_element( index_t j, const T& value ) {
            signed_index_t index = find( j ) ;
            if( index == -1 ) {
                if( nb_elements_ == capacity_ ) grow() ;
                Element& elt = elements_[nb_elements_++] ;
                elt.index = j ;
                elt.value = value ;
            } else {
                elements_[index].value = value ;
            }
        }

        signed_index_t find( index_t j ) const {
            for( index_t i = 0; i < nb_elements_; i++ ) {
                if( elements_[i].index == j ) return i ;
            }
            return -1 ;
        }

        bool get_element( index_t j, T& value ) const {
            signed_index_t index = find( j ) ;
            if( index == -1 ) return false ;
            value = elements_[index].value ;
            return true ;
        }
        void element( index_t i, T& value ) const {
            grgmesh_debug_assert( i < nb_elements_ ) ;
            value = elements_[i] ;
        }


    private:
        void reallocate( index_t new_capacity ) {
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
        index_t nb_elements_ ;
        index_t capacity_ ;
    } ;

    template< class T >
    class GRGMESH_API SparseMatrix {
    public:
        typedef RowImpl< T > Row ;
        typedef SparseMatrix< T > thisclass ;

        SparseMatrix() : rows_( nil ), nb_rows_( 0 )
        {
        }
        ~SparseMatrix() {
            if( rows_ ) delete[] rows_ ;
        }
        void build_matrix( index_t n )
        {
            nb_rows_= n ;
            rows_ = new Row[n] ;
        }

        index_t n() const { return nb_rows_ ; }
        const Row& row( index_t i ) const { return rows_[i] ; }
        void set_element( index_t i, index_t j, const T& value ) {
            rows_[i].set_element( j, value ) ;
        }
        bool get_element( index_t i, index_t j, T& value ) const {
            return rows_[i].get_element( j, value ) ;
        }
    private:
        SparseMatrix( const thisclass &rhs ) ;
        thisclass& operator=( const thisclass &rhs ) ;

    private:
        Row* rows_ ;
        index_t nb_rows_ ;
    } ;
}

#endif
