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
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Sup�rieure de G�ologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
 */

#ifndef __GRGMESH_MATRIX__
#define __GRGMESH_MATRIX__

#include <grgmesh/common.h>
#include <deque>

namespace GRGMesh {
    /**
     * @brief enum of MatrixType, This is useful to further specialize the template in the future
     * */
    enum MatrixType {
        heavy = 0, light = 1, other = 2
    } ;

    /**
     * @Brief Basic container for the sparse matrix, i.e. the "elements".
     * */
    template< class T >
    struct ElementImpl {
        const static index_t NOT_USED = index_t( -1 ) ;
        ElementImpl()
            : index( NOT_USED )
        {
        }
        T value ;
        index_t index ;

    } ;

    /**
     * @brief Basic "Row" of the matrix, this stores the elements of the matrix in a line-oriented way
     * */
    template< class T >
    class RowImpl {
    public:
        typedef ElementImpl< T > Element ;

        RowImpl()
            : nb_elements_( 0 ), capacity_( 0 ), elements_( nullptr )
        {
        }
        ~RowImpl()
        {
            if( elements_ != nullptr ) {
                delete[] elements_ ; // Constructor does not allocate elements_
            }
        }

        void set_element( index_t j, const T& value )
        {
            index_t index ;
            if( !find( j, index ) ) {
                if( nb_elements_ == capacity_ ) grow() ;
                Element& elt = elements_[nb_elements_++ ] ;
                elt.index = j ;
                elt.value = value ;
            } else {

                elements_[index].value = value ;

            }
        }

        bool find( index_t j, index_t& index = dummy_index_t ) const
        {
            for( index_t e = 0; e < nb_elements_; e++ ) {
                if( elements_[e].index == j ) {
                    index = e ;
                    return true ;
                }
            }
            return false ;
        }
        bool exist( index_t j )
        {
            for( index_t e = 0; e < nb_elements_; e++ ) {
                if( elements_[e].index == j ) {
                    return true ;
                }
            }
            return false ;
        }

        bool get_element( index_t j, T& value ) const
        {
            index_t index ;
            if( !find( j, index ) ) return false ;
            value = elements_[index].value ;
            return true ;
        }
        void element( index_t e, T& value ) const
        {
            grgmesh_debug_assert( e < nb_elements_ ) ;
            value = elements_[e].value ;
        }
        index_t index( index_t e ) const
        {
            grgmesh_debug_assert( e < nb_elements_ ) ;
            return elements_[e].index ;
        }
        T& operator[]( index_t i ) const
        {
            grgmesh_debug_assert( i < nb_elements_ ) ;
            return elements_[i].value ;
        }
        index_t nb_elements() const
        {
            return nb_elements_ ;
        }

        void reallocate( index_t new_capacity )
        {
            Element* new_elements = new Element[new_capacity] ;

            if( elements_ != nullptr ) {
                std::copy( elements_, elements_ + nb_elements_, new_elements ) ;
                delete[] elements_ ; // Constructor does not allocate elements_
            }
            elements_ = new_elements ;
        }

    private:
        void grow()
        {
            capacity_ = capacity_ * 2 ;
            reallocate( capacity_ ) ;
        }

    private:
        Element* elements_ ;
        index_t nb_elements_ ;
        index_t capacity_ ;
    } ;

    /*!
     *  @brief This is the parent class for sparse matrices, the main difference between light and heavy type matrices
     * depend on the contents of rows elements: Light will contain type T objects, while heavy an index to access a std::deque.
     * */
    template< class T, typename RowType > class SparseMatrixImpl {
    public:
        typedef RowImpl< RowType > Row ;
        SparseMatrixImpl( bool is_symmetrical = false )
            : rows_( NULL ), ni_( 0 ), nj_( 0 ), is_symmetrical_( is_symmetrical )
        {
        }
        ~SparseMatrixImpl()
        {
            if( rows_ ) delete[] rows_ ;
        }
        /*!
         * test the existence of a given i-j element
         * @param[in] index_t i the given row
         * @param[in] index_t j the given column
         * @return bool true if it exists, false if it does not exist
         */
        bool exist( index_t i, index_t j ) const
        { // test existence of the i-j element
            grgmesh_debug_assert( i < ni_ && j < nj_ && i >= 0 && j >= 0 ) ;
            return rows_[i].exist( j ) ;
        }
        /*!
         * @brief gets number of elements within a row
         * @param[in] i row index
         * @retrun index_t number of elements
         */
        index_t get_nb_elements_in_line( index_t i ) const
        {
            grgmesh_debug_assert( i < ni_ ) ;
            return rows_[i].nb_elements() ;
        }
        /*!
         * @brief gets the j that correspond to the given index within the row
         * @param[in] i row index
         * @param[in] e the index within the row
         * @retrun index_t of the corresponding j column in the matrix
         */
        index_t get_column_in_line( index_t i, index_t e ) const
        {
            return rows_[i].index( e ) ;
        }
        /*!
         * @brief gets the rows_ index corresponding to a given i-j couple
         * @param[in] i row index
         * @param[in] j column index
         * @param[out] index the index within a row
         * @retrun true if success, false if the i-j couple is empty
         */
        bool get_index_in_line( index_t i, index_t j, index_t& index ) const
        {
            grgmesh_debug_assert( i < ni_ && j < nj_ && i >= 0 && j >= 0 ) ;
            return rows_[i].find( j, index ) ;
        }
        /*!
         * returns the number of lines
         * @return index_t ni number of lines of the matrix
         */
        index_t ni() const
        {
            return ni_ ;
        }
        /*!
         * returns the number of columns
         * @return index_t nj number of columns of the matrix
         */
        index_t nj() const
        {
            return nj_ ;
        }
        /*!
         * returns true if the matrix is symmetrical
         * @return bool is_symmetrical (true if yes)
         */
        bool is_symmetrical() const
        {
            return is_symmetrical_ ;
        }
        /*!
         * build the matrix, in allocate the ni number of lines and sets the matrix dimensions
         * @param[in] index_t ni number of lines of the matrix
         * @param[in] index_t nj number of columns of the matrix
         */
        void build_matrix( index_t ni, index_t nj )
        {
            ni_ = ni ;
            nj_ = nj ;
            rows_ = new Row[ni] ;
            for( index_t i = 0; i < ni; ++i ) {
                rows_[i].reallocate( nj ) ;
            }
        }

    protected:
        Row* rows_ ;
        index_t ni_, nj_ ; // matrix dimensions
        bool is_symmetrical_ ;
    } ;

    /*!
     * declaration of a template class SparseMatrix, it will be specialazed for the different MatrixType
     * */
    template< class T, MatrixType Light = MatrixType(
        2 * sizeof(T) <= 2 * sizeof(index_t) + sizeof(T) ) >
    class SparseMatrix: public SparseMatrixImpl< T, T > {
    } ;

    /*!
     * specialization of SparseMatrix for MatrixType "light"
     * */
    template< class T > class SparseMatrix< T, light > : public SparseMatrixImpl< T,
        T > {
    public:
        typedef SparseMatrix< T, light > thisclass ;
        SparseMatrix( bool is_symetrical = false )
            : SparseMatrixImpl< T, T >( is_symetrical )
        {
        }
        /*!
         * set the value of element i-j in the matrix
         * @param[in] i row index
         * @param[in] j column index
         * @param[in] value to store
         * @return bool true (for instance no checks for errors...)
         * */
        bool set_element( index_t i, index_t j, const T& value )
        {
            this->rows_[i].set_element( j, value ) ;
            if( this->is_symmetrical_ ) {
                this->rows_[j].set_element( i, value ) ;
            }
            return true ;
        }
        /*
         * get the value of element i-j in the matrix, returns false if no element is found
         * @param[in] i row index
         * @param[in] j column index
         * @param[out] value to retrieve
         * @return bool true if success and false if the element is inexistent
         * */
        bool get_element( index_t i, index_t j, T& value ) const
        { // valeur du i_j_ieme element stocke sur la ligne
            return this->rows_[i].get_element( j, value ) ;
        }
        /*!
         *  get the value of e-element (index within the row, not in the matrix) on line i
         * @param[in] i row index
         * @param[in] e index within the row
         * @param[out] value to retrieve
         * */
        void get_element_in_line( index_t i, index_t e, T& value ) const
        { // valeur du e_ieme element stocke sur la ligne
            value = this->rows_[i][e] ;
            return ;
        }

    private:
        SparseMatrix( const thisclass &rhs ) ;
        thisclass& operator=( const thisclass &rhs ) ;
    } ;

    /*!
     * specialization of SparseMatrix for MatrixType "heavy" ,
     * the main difference with light is that here we store only
     * one copy of the data in the case of a symmetrical matrix.
     * The data are strored in a std::deque and the rows contains the
     * ids of the values within the deque.
     * */
    template< class T > class SparseMatrix< T, heavy > : public SparseMatrixImpl< T,
        index_t > {
    public:
        typedef SparseMatrix< T, heavy > thisclass ;
        SparseMatrix( bool is_symetrical = false )
            : SparseMatrixImpl< T, index_t >( is_symetrical )
        {
        }
        /*!
         * set the value of element i-j in the matrix
         * @param[in] i row index
         * @param[in] j column index
         * @param[in] value to store
         * @return bool true (for instance no checks for errors...)
         * */
        bool set_element( index_t i, index_t j, const T& value )
        {
            // small difference with light type: we fill a deque
            index_t value_id ;
            if( this->exist( i, j ) ) {
                value_id = get_value_id( i, j ) ;
                values_[value_id] = value ;
            } else {
                values_.push_back( value ) ;
                value_id = values_.size() - 1 ;
            }
            this->rows_[i].set_element( j, value_id ) ;
            if( this->is_symmetrical_ ) {
                this->rows_[j].set_element( i, value_id ) ;
            }
            return true ;
        }
        /*
         * get the value of element i-j in the matrix, returns false if no element is found
         * @param[in] i row index
         * @param[in] j column index
         * @param[out] value to retrieve
         * @return bool true if success and false if the element is inexistent
         * */
        bool get_element( index_t i, index_t j, T& value ) const
        {
            index_t value_id ;
            if( !this->rows_[i].get_element( j, value_id ) ) {
                return false ;
            } ;
            value = values_[value_id] ;
            return true ;
        }
        /*!
         *  get the value of e-element (index within the row, not in the matrix) on line i
         * @param[in] i row index
         * @param[in] e index within the row
         * @param[out] value to retrieve
         * */
        void get_element_in_line( index_t i, index_t e, T& value ) const
        {
            index_t value_id = this->rows_[i][e] ;
            value = values_[value_id] ;
            return ;
        }
        void print_matrix( void )
        {
            std::cout << "deque size = " << values_.size() << std::endl ;
        }

    private:
        index_t get_value_id( index_t i, index_t j )
        {
            return this->rows_[i].find( j ) ;
        }
        SparseMatrix( const thisclass &rhs ) ;
        thisclass& operator=( const thisclass &rhs ) ;
    private:
        std::deque< T > values_ ;
    } ;

    template< class T >
    void product_matrix_by_vector(
        const SparseMatrix< T >& mat1,
        const std::vector< T >& mat2,
        std::vector< T >& result )
    {
        grgmesh_assert( mat1.nj() == mat2.size() ) ;
        result.clear() ;
        result.resize( mat1.ni(), 0 ) ;

        for( index_t i = 0; i < mat1.ni(); ++i ) {
            grgmesh_assert( i >= 0 && i < result.size() ) ;
            for( index_t e = 0; e < mat1.get_nb_elements_in_line( i ); ++e ) {
                index_t j = mat1.get_column_in_line( i, e ) ;
                T i_j_result ;
                mat1.get_element_in_line( i, e, i_j_result ) ;
                i_j_result *= mat2[j] ;
                result[i] += i_j_result ;
            }
        }
    }
}
#endif
