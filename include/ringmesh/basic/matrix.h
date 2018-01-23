/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#pragma once

#include <ringmesh/basic/common.h>
#include <ringmesh/basic/task_handler.h>

#include <deque>
#include <memory>

/*!
 * @file Template matix declarations and definitions
 * @author Arnaud Botella and Andrea Borghi
 */

namespace RINGMesh
{
    /*!
     * @brief enum of MatrixType, This is useful to further specialize the
     * template in the future
     */
    enum MatrixType
    {
        heavy = 0,
        light = 1,
        other = 2
    };

    /*!
     * @brief Basic container for the sparse matrix, i.e. the "elements".
     */
    template < typename T >
    struct ElementImpl
    {
        const static index_t NOT_USED = index_t( -1 );
        ElementImpl() : index( NOT_USED ) {}

        T value;
        index_t index;
    };

    /*!
     * @brief Basic "Row" of the matrix, this stores the elements of the matrix
     * in a line-oriented way
     */
    template < typename T >
    class RowImpl
    {
    public:
        using Element = ElementImpl< T >;
        RowImpl() : elements_( new Element[4] ) {}

        void set_element( index_t j, const T& value )
        {
            index_t index = find( j );
            if( index == NO_ID )
            {
                push_element( j, value );
            }
            else
            {
                elements_[index].value = value;
            }
        }

        void push_element( index_t j, const T& value )
        {
            if( nb_elements_ == capacity_ )
            {
                grow();
            }
            Element& elt = elements_[nb_elements_++];
            elt.index = j;
            elt.value = value;
        }

        index_t find( index_t j ) const
        {
            for( auto e : range( nb_elements_ ) )
            {
                if( elements_[e].index == j )
                {
                    return e;
                }
            }
            return NO_ID;
        }

        bool exist( index_t j )
        {
            for( auto e : range( nb_elements_ ) )
            {
                if( elements_[e].index == j )
                {
                    return true;
                }
            }
            return false;
        }

        std::tuple< bool, T > get_element( index_t j ) const
        {
            index_t index = find( j );
            if( index == NO_ID )
            {
                return std::make_tuple( false, T() );
            }
            T value = elements_[index].value;
            return std::make_tuple( true, value );
        }

        T element( index_t e ) const
        {
            ringmesh_assert( e < nb_elements_ );
            return elements_[e].value;
        }

        index_t index( index_t e ) const
        {
            ringmesh_assert( e < nb_elements_ );
            return elements_[e].index;
        }

        T& operator[]( index_t i ) const
        {
            ringmesh_assert( i < nb_elements_ );
            return elements_[i].value;
        }

        index_t nb_elements() const
        {
            return nb_elements_;
        }

    private:
        void reallocate( index_t new_capacity )
        {
            auto new_elements = new Element[new_capacity];
            std::copy(
                elements_.get(), elements_.get() + nb_elements_, new_elements );
            elements_.reset( new_elements );
        }

        void grow()
        {
            ringmesh_assert( capacity_ != 0 );
            capacity_ = capacity_ * 2;
            reallocate( capacity_ );
        }

    private:
        std::unique_ptr< Element[] > elements_{};
        index_t nb_elements_{ 0 };
        index_t capacity_{ 4 };
    };

    /*!
     *  @brief This is the parent class for sparse matrices, the main difference
     * between light and heavy type matrices
     * depend on the contents of rows elements: Light will contain type T
     * objects, while heavy an index to access a std::deque.
     */
    template < typename T, typename RowType >
    class SparseMatrixImpl
    {
        ringmesh_disable_copy_and_move( SparseMatrixImpl );

    public:
        using Row = RowImpl< RowType >;
        explicit SparseMatrixImpl( bool is_symmetrical = false )
            : is_symmetrical_( is_symmetrical )
        {
        }

        ~SparseMatrixImpl() = default;

        /*!
         * Test the existence of a given i-j element
         * @param[in] i the given row
         * @param[in] j the given column
         * @return bool true if it exists, false if it does not exist
         */
        bool exist( index_t i, index_t j ) const
        { // test existence of the i-j element
            ringmesh_assert( i < ni_ && j < nj_ && i >= 0 && j >= 0 );
            return rows_[i].exist( j );
        }

        /*!
         * @brief gets number of elements within a row
         * @param[in] i row index
         * @return index_t number of elements
         */
        index_t get_nb_elements_in_line( index_t i ) const
        {
            ringmesh_assert( i < ni_ );
            return rows_[i].nb_elements();
        }

        /*!
         * @brief gets the j that correspond to the given index within the row
         * @param[in] i row index
         * @param[in] e the index within the row
         * @return index_t of the corresponding j column in the matrix
         */
        index_t get_column_in_line( index_t i, index_t e ) const
        {
            return rows_[i].index( e );
        }

        /*!
         * @brief gets the rows_ index corresponding to a given i-j couple
         * @param[in] i row index
         * @param[in] j column index
         * @return a tuple containing:
         * - a boolean: true if success, else false if the i-j couple is empty.
         * - the index the index within a row if any, else NO_ID.
         */
        std::tuple< bool, index_t > get_index_in_line(
            index_t i, index_t j ) const
        {
            ringmesh_assert( i < ni_ && j < nj_ && i >= 0 && j >= 0 );
            index_t index = rows_[i].find( j );
            return std::make_tuple( index != NO_ID, index );
        }

        /*!
         * returns the number of lines
         * @return index_t ni number of lines of the matrix
         */
        index_t ni() const
        {
            return ni_;
        }

        /*!
         * returns the number of columns
         * @return index_t nj number of columns of the matrix
         */
        index_t nj() const
        {
            return nj_;
        }

        /*!
         * returns true if the matrix is symmetrical
         * @return bool is_symmetrical (true if yes)
         */
        bool is_symmetrical() const
        {
            return is_symmetrical_;
        }

        /*!
         * build the matrix, in allocate the ni number of lines and sets the
         * matrix dimensions
         * @param[in] ni number of lines of the matrix
         * @param[in] nj number of columns of the matrix
         */
        void build_matrix( index_t ni, index_t nj )
        {
            ni_ = ni;
            nj_ = nj;
            rows_.reset( new Row[ni] );
        }

        /*!
         *  get the value of e-element (index within the row, not in the matrix)
         * on line i.
         *  this code should never be reached.
         * @param[in] i row index
         * @param[in] e index within the row
         * @return value to retrieve
         */
        T get_element_in_line( index_t i, index_t e ) const
        {
            ringmesh_unused( i );
            ringmesh_unused( e );
            ringmesh_assert_not_reached;
            return T();
        }

    protected:
        std::unique_ptr< Row[] > rows_{};
        index_t ni_{ 0 }, nj_{ 0 }; // matrix dimensions
        bool is_symmetrical_;
    };

    /*!
     * declaration of a template class SparseMatrix, it will be specialazed for
     * the different MatrixType
     */
    template < typename T,
        MatrixType Light = MatrixType(
            2 * sizeof( T ) <= 2 * sizeof( index_t ) + sizeof( T ) ) >
    class SparseMatrix : public SparseMatrixImpl< T, T >
    {
    };

    /*!
     * specialization of SparseMatrix for MatrixType "light"
     */
    template < typename T >
    class basic_api SparseMatrix< T, light > : public SparseMatrixImpl< T, T >
    {
    public:
        using thisclass = SparseMatrix< T, light >;
        explicit SparseMatrix( bool is_symetrical = false )
            : SparseMatrixImpl< T, T >( is_symetrical )
        {
        }

        /*!
         * set the value of element i-j in the matrix
         * @param[in] i row index
         * @param[in] j column index
         * @param[in] value to store
         */
        void set_element( index_t i, index_t j, const T& value )
        {
            this->rows_[i].set_element( j, value );
            if( this->is_symmetrical_ )
            {
                this->rows_[j].set_element( i, value );
            }
        }

        /*!
         * set the value of element i-j in the matrix without verifying
         * if the element i-j already exists !!! BE CAREFULL
         * @param[in] i row index
         * @param[in] j column index
         * @param[in] value to store
         */
        void push_element( index_t i, index_t j, const T& value )
        {
            this->rows_[i].push_element( j, value );
            if( this->is_symmetrical_ )
            {
                this->rows_[j].push_element( i, value );
            }
        }

        /*
         * get the value of element i-j in the matrix, returns false if no
         * element is found
         * @param[in] i row index
         * @param[in] j column index
         * @return a tuple containing:
         * - a boolean: true if success and false if the element is inexistent
         * - value to retrieve if any.
         */
        std::tuple< bool, T > get_element( index_t i, index_t j ) const
        { // valeur du i_j_ieme element stocke sur la ligne
            return this->rows_[i].get_element( j );
        }

        /*!
         *  get the value of e-element (index within the row, not in the matrix)
         * on line i
         * @param[in] i row index
         * @param[in] e index within the row
         * @return value to retrieve
         */
        T get_element_in_line( index_t i, index_t e ) const
        { // valeur du e_ieme element stocke sur la ligne
            return this->rows_[i][e];
        }
    };

    /*!
     * specialization of SparseMatrix for MatrixType "heavy" ,
     * the main difference with light is that here we store only
     * one copy of the data in the case of a symmetrical matrix.
     * The data are stored in a std::deque and the rows contains the
     * ids of the values within the deque.
     */
    template < typename T >
    class basic_api SparseMatrix< T, heavy > : public SparseMatrixImpl< T, index_t >
    {
    public:
        using thisclass = SparseMatrix< T, heavy >;
        explicit SparseMatrix( bool is_symetrical = false )
            : SparseMatrixImpl< T, index_t >( is_symetrical )
        {
        }

        /*!
         * set the value of element i-j in the matrix
         * @param[in] i row index
         * @param[in] j column index
         * @param[in] value to store
         */
        void set_element( index_t i, index_t j, const T& value )
        {
            // small difference with light type: we fill a deque
            index_t value_id;
            if( this->exist( i, j ) )
            {
                value_id = get_value_id( i, j );
                values_[value_id] = value;
            }
            else
            {
                values_.push_back( value );
                value_id = values_.size() - 1;
            }
            this->rows_[i].set_element( j, value_id );
            if( this->is_symmetrical_ )
            {
                this->rows_[j].set_element( i, value_id );
            }
        }

        /*
         * get the value of element i-j in the matrix, returns false if no
         * element is found
         * @param[in] i row index
         * @param[in] j column index
         * @return the value to retrieve
         */
        T get_element( index_t i, index_t j ) const
        {
            bool value_exists{ false };
            index_t value_id;
            std::tie( value_exists, value_id ) =
                this->rows_[i].get_element( j );
            if( !value_exists )
            {
                return T();
            }
            return values_[value_id];
        }

        /*!
         *  get the value of e-element (index within the row, not in the matrix)
         * on line i
         * @param[in] i row index
         * @param[in] e index within the row
         * @return value to retrieve
         */
        T get_element_in_line( index_t i, index_t e ) const
        {
            index_t value_id = this->rows_[i][e];
            return values_[value_id];
        }

    private:
        index_t get_value_id( index_t i, index_t j )
        {
            return this->rows_[i].find( j );
        }

    private:
        std::deque< T > values_{};
    };

    // Note: without light or heavy, it does not compile on Windows.
    // Error C2770. BC
    template < typename T >
    std::vector< T > product_matrix_by_vector(
        const SparseMatrix< T, light >& mat1, const std::vector< T >& mat2 )
    {
        ringmesh_assert( mat1.nj() == mat2.size() );
        std::vector< T > result( mat1.ni(), 0 );
        parallel_for( mat1.ni(), [&mat1, &mat2, &result]( index_t i ) {
            for( auto e : range( mat1.get_nb_elements_in_line( i ) ) )
            {
                index_t j = mat1.get_column_in_line( i, e );
                T i_j_result = mat1.get_element_in_line( i, e );
                i_j_result *= mat2[j];
                result[i] += i_j_result;
            }
        } );
        return result;
    }

} // namespace RINGMesh
