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
#include <deque>
#include <map>
#include <list>


namespace GRGMesh {

	template< class T > class SparseMatrix {
	public:
		SparseMatrix(bool is_symmetrical = false){}
		~SparseMatrix(){}
		index_t get_nb_element_line( index_t i ) ; // nb d'elements stockes sur la ligne
		 void get_value_line(index_t i, index_t e, const T& value ){}// valeur du e_ieme element stocke sur la ligne
		 index_t get_index_line(index_t i, index_t e ){} // j du e_ieme element stocke sur la ligne
		 bool get_element(index_t i, index_t j, const T& value){} // valeur du i_j_ieme element stocke sur la ligne
		 bool exist(index_t i, index_t j){}
		 void print_matrix(){}
		 bool set_element(index_t i, index_t j, const T& value){}
		 void build_matrix(index_t ni, index_t nj = ni){}
	private:
		SparseMatrixImpl* sparse_matrix_impl_;
	};

	template < class T > class SparseMatrixImpl{
	public:
		SparseMatrixImpl(bool is_symmetrical = false){}
		~SparseMatrixImpl(){}
		virtual void get_line(index_t i, std::list<T>&)=0;
		virtual bool get_element(index_t i, index_t j, const T& value)=0;
		virtual bool find_element(index_t i, index_t j)=0;
		virtual void print_matrix()=0;
		virtual bool set_element(index_t i, index_t j, const T& value)=0;
		virtual void set_size(index_t ni, index_t nj)=0;
	};

	template< class T > class SparseMatrixHeavyType;
	template< class T > class SparseMatrixLightType;

	template< class T > class LightMatrix;

	template< class T > struct ElementImpl;
	template< class T > class GRGMESH_API RowOrColumnImpl;



	template < class T > class MatrixFactory
	{
	  public:
	    static SparseMatrix < T > *new_matrix(bool is_symmetrical = false)
	    {
	      if (2*sizeof(T) <= 2*sizeof(index_t)+sizeof(T))
	        return new SparseMatrixLightType<T>(is_symmetrical);
	      else
	        return new SparseMatrixHeavyType<T>(is_symmetrical);
	      return 0;
	    }
	};


	typedef std::map<index_t, index_t>  Rows ;
	typedef std::map<index_t, index_t>::iterator  RowsIterator ;

	class RowKeyIterator : public RowsIterator{
	public:
		RowKeyIterator() : RowsIterator() {};
		RowKeyIterator(RowsIterator s) : RowsIterator(s) {};
		index_t* operator->() { return (index_t* const)&(RowsIterator::operator->()->first); }
		index_t operator*() { return RowsIterator::operator*().first; }
	};


	template< class T > class GRGMESH_API SparseMatrixHeavyType : public SparseMatrix < T > {
	public:
		typedef SparseMatrixHeavyType< T > thisclass ;

		SparseMatrixHeavyType(bool is_symmetrical = false ) :
			is_symmetrical_( is_symmetrical )
		{
			//values_.push_back(-1); // undefined values == -1 ...
		}
		~SparseMatrixHeavyType(){};

		void set_size(index_t ni, index_t nj){
			ni_ = ni;
			nj_ = nj;
		}

		bool get_element(index_t i, index_t j, T& value){
			grgmesh_debug_assert( i < ni_ ) ;
			grgmesh_debug_assert( j < nj_ ) ;
			grgmesh_debug_assert( i >= 0 ) ;
			grgmesh_debug_assert( j >= 0 ) ;

			value = values_[value_id_matrix_[i][j]];
			return true;
		}
		void get_line(index_t i, std::list<T>& output){
			RowKeyIterator it;
			//std::cout << "output.size() = " << output.size() << std::endl;
			int c=0;
			for (it=value_id_matrix_[i].begin(); it != value_id_matrix_[i].end(); ++it){
				//std::cout << " --------- ITERATION nb " << c << "---------" << std::endl;
				output.push_back(values_[value_id_matrix_[i][*it]]);
				//std::cout << "output.size() = " << output.size() << std::endl;
				c++;

			}
		}

		void print_matrix(){
			std::cout << "values contained in DEQUE:"<< std::endl;
			print_deque();
			std::cout << "matrix representation:"<< std::endl;
			std::cout << "value_id_matrix_.size() = " << value_id_matrix_.size() << "\n " ;
			for (int i = 0; i<ni_ ; i++){
				for (int j = 0; j<ni_ ; j++){
					std::cout << values_[value_id_matrix_[i][j]] << " " ;
				}
				std::cout <<std::endl ;
			}
			std::cout << "value_id_matrix_.size() = " << value_id_matrix_.size() << "\n " ;

		}

		bool set_element( index_t i, index_t j, const T& value ){
			if ( ! insert_value( i, j, value ) ){ return false; }
			//insert_value( i, j, value );
			if (is_symmetrical_){
				//insert_value( j, i, value );
     			if ( ! insert_value( j, i, value ) ){ return false; }
			}
			return true;
		}
	private:
		SparseMatrixHeavyType( const thisclass &rhs ) ;
		thisclass& operator=( const thisclass &rhs ) ;
		std::deque< T > values_;
		std::map< index_t, Rows > value_id_matrix_; //  < I , < J, value> >
		index_t ni_,nj_; // matrix dimensions
		bool is_symmetrical_;
	private:
		void print_deque(){
			int c = 0;
			for (c=0 ; c<values_.size(); c++){
				std::cout << c << " "<< values_[c] << std::endl;

			}
		}
		bool insert_value( index_t i, index_t j, const T& value ){
			grgmesh_debug_assert( i < ni_ ) ;
			grgmesh_debug_assert( j < nj_ ) ;
			grgmesh_debug_assert( i >= 0 ) ;
			grgmesh_debug_assert( j >= 0 ) ;

			std::map < index_t, Rows >::iterator    it_i;
			RowsIterator it_j;

			// todo, smart way to erase the content of values_ if a new element replaces the old one
			values_.push_back(value);
			index_t value_id = values_.size() - 1;

			std::pair < index_t, index_t >  j_to_insert ;
			j_to_insert.first = j;
			j_to_insert.second= value_id;

			if ( ! bool(value_id_matrix_.count(i))){ // create empty line
				std::map < index_t, index_t >  map_j_to_insert ;
				map_j_to_insert.insert(j_to_insert);

				std::pair < index_t, Rows > ij_to_insert;
				ij_to_insert.first = i;
				ij_to_insert.second = map_j_to_insert;
									//j_to_insert;
				value_id_matrix_.insert(ij_to_insert);
			}
			else { // add column if line already exists
				//std::cout << "value_id_matrix_[i].count(j) : " << value_id_matrix_[i].count(j) << std::endl;
				if ( bool( value_id_matrix_[i].count(j))){value_id_matrix_[i].erase(j);}
				value_id_matrix_[i].insert(j_to_insert);
			}
			return true;
		}
	};






// Arnaud's implementation
	// incomplete class declarations, less compilation problems...
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

    template< class T >
    class RowImpl {
    public:
        typedef ElementImpl< T > Element ;

        RowImpl()
            : nb_elements_( 0 ), capacity_( 4 )
        {
            elements_ = new Element[capacity_] ;
        }
        ~RowImpl() {
            delete[] elements_ ;
        }

        void set_element( index_t j, const T& value ) {
            index_t index ;
            if( !find( j, index ) ) {
                if( nb_elements_ == capacity_ ) grow() ;
                Element& elt = elements_[nb_elements_++] ;
                elt.index = j ;
                elt.value = value ;
            } else {
                elements_[index].value = value ;
            }
        }

        bool find( index_t j, index_t& index = dummy_index_t ) const {
            for( index_t i = 0; i < nb_elements_; i++ ) {
                if( elements_[i].index == j ) {
                    index = i ;
                    return true ;
                }
            }
            return false ;
        }

        bool get_element( index_t j, T& value ) const {
            index_t index ;
            if( !find( j, index ) ) return false ;
            value = elements_[index].value ;
            return true ;
        }
        void element( index_t i, T& value ) const {
            grgmesh_debug_assert( i < nb_elements_ ) ;
            value = elements_[i].value ;
        }
        index_t index( index_t i ) const {
            grgmesh_debug_assert( i < nb_elements_ ) ;
            return elements_[i].index ;
        }
        T& operator[]( index_t i ) const {
            grgmesh_debug_assert( i < nb_elements_ ) ;
            return elements_[i].value ;
        }
        index_t nb_elements() const {
            return nb_elements_ ;
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
    class SparseMatrixLightType  : public SparseMatrix < T >{
        //todo need to handle symmetric matrix
    public:
        typedef RowImpl< T > Row ;
        typedef SparseMatrixLightType< T > thisclass ;

        SparseMatrixLightType(bool is_symmetrical = false) : rows_( nil ), nb_rows_( 0 )
        {
        	is_symmetrical_ = is_symmetrical;
        }
        ~SparseMatrixLightType() {
            if( rows_ ) delete[] rows_ ;
        }
        void build_matrix( index_t n )
        {
            nb_rows_= n ;
            rows_ = new Row[n] ;
        }

        index_t n() const { return nb_rows_ ; }
        const Row& row( index_t i ) const { return rows_[i] ; }
        bool set_element( index_t i, index_t j, const T& value ) {
            grgmesh_debug_assert( i < nb_rows_ ) ;
            rows_[i].set_element( j, value ) ;
            if (is_symmetrical_){rows_[j].set_element( i, value ) ; }
            return true;
        }
        bool get_element( index_t i, index_t j, T& value ) const {
            grgmesh_debug_assert( i < nb_rows_ ) ;
            return rows_[i].get_element( j, value ) ;
        }
        bool find_element( index_t i, index_t j ) const{
            grgmesh_debug_assert( i < nb_rows_ ) ;
            return rows_[i].find( j ) ;
        }
        // todo: print_matrix and get_line
        void get_line(index_t i, std::list<T>&){}
    private:
        SparseMatrixLightType( const thisclass &rhs ) ;
        thisclass& operator=( const thisclass &rhs ) ;

    private:
        Row* rows_ ;
        index_t nb_rows_ ;
        bool is_symmetrical_;
    } ;
}

#endif
