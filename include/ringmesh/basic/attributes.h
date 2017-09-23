/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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
#include <ringmesh/basic/logger.h>
#include <geogram/basic/geofile.h>
#include <map>
#include <typeinfo>
#include <set>

/**
 * \file ringmesh/basic/attributes.h
 * \brief Generic mechanism for attributes.
 * This file has been strongly inspired from the attributes of the geogram library.
 * Many thanks to Bruno Levy (Bruno.Levy@inria.fr) who did the first implementation in Geogram.
 */

namespace RINGMesh {

    class AttributeStore;

    /**
     * \brief Internal class for creating an AttributeStore
     *  from the type name of its elements.
     */
    class RINGMESH_API AttributeStoreCreator {
    public:

        /**
         * \brief AttributeStoreCreator destructor.
         */
        virtual ~AttributeStoreCreator() = default;

        /**
         * \brief Creates a new attribute store.
         * \param[in] dimension number of elements in each item
         * \return a pointer to the newly created AttributeStore
         */
        virtual AttributeStore* create_attribute_store( index_t dimension ) = 0;

    private:
        std::string element_type_name_;
        std::string element_typeid_name_;
    };

    /**
     * \brief Notifies a set of AttributeStoreObservers
     *  each time the stored array changes size and/or
     *  base address and/or dimension.
     */
    class RINGMESH_API AttributeStore {
    public:
        /**
         * \brief AttributeStore constructor.
         * \param[in] elemsize size of one element,
         *  in bytes.
         * \param[in] dim number of elements in
         *  each item. Default is 1 for standard
         *  attributes and can be greater for vector
         *  attributes.
         */
        AttributeStore( index_t elemsize, index_t dim = 1 );

        /**
         * \brief AttributeStore destructor.
         */
        virtual ~AttributeStore();

        /**
         * \brief Tests whether this AttributeStore stores
         *  elements of a given type.
         * \param[in] type_name the name of the type, as given by
         *     typeid(T).name()
         * \retval true if this AttributeStore stores elements
         *   of type \p type_name
         * \retval false otherwise
         */
        virtual bool elements_type_matches( const std::string& type_name ) const = 0;

        /**
         * \brief Gets the typeid name of the element type stored
         *  in this AttributeStore.
         * \return the typeid name, as a string.
         */
        virtual std::string element_typeid_name() const = 0;

        /**
         * \brief Gets the size.
         * \return the number of items
         */
        index_t size() const
        {
            return cached_size_;
        }

        /**
         * \brief Resizes this AttributeStore
         * \param[in] new_size new number of items
         */
        virtual void resize( index_t new_size ) = 0;

        /**
         * \brief Resizes this AttributeStore to 0.
         * \param[in] keep_memory if true, then memory
         *  is kept reserved for future use.
         */
        virtual void clear( bool keep_memory = false ) = 0;

        /**
         * \brief Gets the dimension.
         * \details The dimension is 1 for standard attributes and
         *  can be greater for vector attributes.
         */
        index_t dimension() const
        {
            return dimension_;
        }

        /**
         * \brief Sets the dimension.
         * \details The dimension is 1 for standard attributes and
         *  can be greater for vector attributes. The existing
         *  fields are kept. If the new dimension is greater than
         *  the old one, then new fields are initialized to the default
         *  value for the attribute type.
         * \param[in] dim the new dimension
         */
        virtual void redim( index_t dim ) = 0;

        /**
         * \brief Applies a permutation to the stored attributes.
         * \details Applying a permutation to the data is equivalent
         *  to:
         * \code
         * for(i=0; i<permutation.size(); i++) {
         *    data2[i] = data[permutation[i]]
         * }
         * data = data2 ;
         * \endcode
         * But it is done in-place.
         * \param[in] permutation the permutation.
         *  It is temporarily changed during execution of the
         *  function, but identical to the input on exit.
         * \note This function uses memcpy(). If required, it
         *  can be overloaded in derived classes.
         */
        virtual void apply_permutation( const std::vector< index_t >& permutation );

        /**
         * \brief Compresses the stored attributes, by
         *  applying an index mapping that fills-in the gaps.
         * \details This is equivalent to:
         * \code
         * for(i=0; i<size(); i++) {
         *    if(old2new[i] != index_t(-1)) {
         *       data2[old2new[i]] = data[i];
         *    }
         * }
         * data = data2 ;
         * \endcode
         * \param[in] old2new the index mapping to be applied.
         * \pre old2new[i] <= i || old2new[i] == index_t(-1)
         * \note This function uses memcpy(). If required, it
         *  can be overloaded in derived classes.
         */
        virtual void compress( const std::vector< index_t >& old2new );

        /**
         * \brief Zeroes all the memory associated with this
         *  AttributeStore.
         * \details Subclasses may overload this function for
         *  attributes that have non "plain ordinary datatypes"
         *  and that need a more elaborate initialization mechanism.
         */
        virtual void zero();

        /**
         * \brief Creates a new AttributeStore that is a carbon copy
         *  of this AttributeStore.
         * \details Only the data is copied, observers are not copied.
         */
        virtual AttributeStore* clone() const = 0;

        /**
         * \brief Copies an item
         * \param[in] to index of the destination item
         * \param[in] from index of the source item
         */
        void copy_item( index_t to, index_t from )
        {
            ringmesh_assert( from < cached_size_ );
            ringmesh_assert( to < cached_size_ );
            index_t item_size = element_size_ * dimension_;
            for( auto i : range( item_size ) ) {
                *( cached_base_addr_ + to * item_size + i ) = *( cached_base_addr_
                    + from * item_size + i );
            }
        }

        /**
         * \brief Gets a pointer to the stored data.
         * \return A pointer to the memory block
         */
        pointer data()
        {
            return cached_base_addr_;
        }

        /**
         * \brief Gets a pointer to the stored data.
         * \return A const pointer to the memory block
         */
        pointer data() const
        {
            return cached_base_addr_;
        }

        /**
         * \brief Gets the element size.
         * \return the size of an element, in bytes
         */
        size_t element_size() const
        {
            return size_t( element_size_ );
        }

        /**
         * \brief Gets the total number of elements.
         * \details This corresponds to one position past
         *  the last valid index.
         * \return the total number of elements.
         */
        index_t nb_elements() const
        {
            return cached_size_ * dimension_;
        }

        /**
         * \brief Tests whether a given element type is registered in
         *   the system.
         * \param[in] element_type_name a const reference to a string
         *   with the C++ type name
         * \retval true if the element type was registered
         * \retval false otherwise
         */
        static bool element_type_name_is_known(
            const std::string& element_type_name )
        {
            return ( type_name_to_creator_.find( element_type_name )
                != type_name_to_creator_.end() );
        }

        /**
         * \brief Tests whether a given element type is registered in
         *   the system.
         * \param[in] element_typeid_name a const reference to a string
         *   with the mangled type, as given by typeid(T).name()
         * \retval true if the element type was registered
         * \retval false otherwise
         */
        static bool element_typeid_name_is_known(
            const std::string& element_typeid_name )
        {
            return ( typeid_name_to_type_name_.find( element_typeid_name )
                != typeid_name_to_type_name_.end() );
        }

        /**
         * \brief Creates an attribute store of a given type
         * \param[in] element_type_name a const reference to a string with
         *  the C++ type of the elements to be stored in the attribute
         * \param[in] dimension number of elements in each item
         */
        static AttributeStore* create_attribute_store_by_element_type_name(
            const std::string& element_type_name,
            index_t dimension )
        {
            ringmesh_assert( element_type_name_is_known( element_type_name ) );
            return type_name_to_creator_[element_type_name]->create_attribute_store(
                dimension );
        }

        /**
         * \brief Gets an element type name from its mangled name
         * \param[in] element_typeid_name a const reference to a
         *  string with the mangled type name, as given by typeid(T).name()
         * \return a string with the C++ type name
         * \pre element_typeid_name_is_known(element_typeid_name)
         */
        static std::string element_type_name_by_element_typeid_name(
            const std::string& element_typeid_name )
        {
            ringmesh_assert( element_typeid_name_is_known( element_typeid_name ) );
            return typeid_name_to_type_name_[element_typeid_name];
        }

        /**
         * \brief Gets an element mangled type name from its C++ name.
         * \param[in] element_type_name a reference to a string with
         *  the C++ type name
         * \return a string with the mangled type name, as given by
         *  typeid(T).name()
         * \pre element_type_name_is_known(element_type_name)
         */
        static std::string element_typeid_name_by_element_type_name(
            const std::string& element_type_name )
        {
            ringmesh_assert( element_type_name_is_known( element_type_name ) );
            return type_name_to_typeid_name_[element_type_name];
        }

        /**
         * \brief Registers a new element type
         * \note Internal use function, one should use
         *  geo_register_attribute_type instead
         * \param[in] creator a pointer to the AttributeStoreCreator
         * \param[in] element_type_name a const reference to a string with the
         *  C++ type name of the elements
         * \param[in] element_typeid_name a const reference to a string with
         *  the mangled type name of the elements, as given by typeid(T).name()
         */
        static void register_attribute_creator(
            AttributeStoreCreator* creator,
            const std::string& element_type_name,
            const std::string& element_typeid_name )
        {
            if( element_type_name_is_known( element_type_name ) ) {
                Logger::warn( "Attributes", element_type_name,
                    " already registered" );
                if( element_typeid_name_is_known( element_typeid_name ) ) {
                    bool already_registered_attribute_has_same_type =
                        ( type_name_to_typeid_name_[element_type_name]
                            == element_typeid_name );
                    ringmesh_assert( already_registered_attribute_has_same_type );
                    ringmesh_unused( already_registered_attribute_has_same_type );
                }
            }
            type_name_to_creator_[element_type_name] = std::unique_ptr<
                AttributeStoreCreator >( creator );
            typeid_name_to_type_name_[element_typeid_name] = element_type_name;
            type_name_to_typeid_name_[element_type_name] = element_typeid_name;
        }

    protected:
        /**
         * \brief If size or base address differ from the
         *  cached values, notify all the observers,
         *  and update the cached base address and size.
         * \param[in] base_addr the new base address
         * \param[in] size the new size
         * \param[in] dim the new dimension
         */
        virtual void notify( pointer base_addr, index_t size, index_t dim );

    protected:
        index_t element_size_;
        index_t dimension_;
        pointer cached_base_addr_ { nullptr };
        index_t cached_size_ { 0 };
        std::mutex lock_;

        static std::map< std::string, std::unique_ptr< AttributeStoreCreator > > type_name_to_creator_;

        static std::map< std::string, std::string > typeid_name_to_type_name_;

        static std::map< std::string, std::string > type_name_to_typeid_name_;

        friend class AttributeStoreObserver;
    };

    /*********************************************************************/

    /**
     * \brief Stores an array of elements of a given type,
     *  and notifies a set of AttributeStoreObservers each time the
     *  storead array changes size and/or base address.
     */
    template< class T > class TypedAttributeStore: public AttributeStore {
    public:

        /**
         * \brief Creates a new empty attribute store.
         * \param[in] dim number of elements in each item,
         *  default value is 1, can be greater for vector
         *  attributes.
         */
        TypedAttributeStore( index_t dim = 1 )
            : AttributeStore( index_t( sizeof(T) ), dim )
        {
        }

        virtual bool elements_type_matches( const std::string& type_name ) const
        {
            return type_name == typeid(T).name();
        }

        virtual std::string element_typeid_name() const
        {
            return typeid(T).name();
        }

    };

    template< class T > class VectorStore: public TypedAttributeStore< T > {
    public:

        /**
         * \brief Creates a new empty attribute store.
         * \param[in] dim number of elements in each item,
         *  default value is 1, can be greater for vector
         *  attributes.
         */
        VectorStore( index_t dim = 1 )
            : TypedAttributeStore< T >( dim )
        {
        }

        virtual void resize( index_t new_size )
        {
            store_.resize( new_size * this->dimension() );
        }


        virtual void clear( bool keep_memory = false )
        {
            if( keep_memory ) {
                store_.resize( 0 );
            } else {
                store_.clear();
            }
        }

        virtual void redim( index_t dim )
        {
            if( dim == this->dimension() ) {
                return;
            }
            std::vector< T > new_store( this->size() * dim );
            index_t copy_dim = GEO::geo_min( dim, this->dimension() );
            for( index_t i = 0; i < this->size(); ++i ) {
                for( index_t c = 0; c < copy_dim; ++c ) {
                    new_store[dim * i + c] = store_[this->dimension() * i + c];
                }
            }
            store_.swap( new_store );
        }

        virtual AttributeStore* clone() const
        {
            VectorStore< T >* result = new VectorStore< T >(
                this->dimension() );
            result->resize( this->size() );
            result->store_ = store_;
            return result;
        }

    private:
        std::vector< T > store_;
    };


    /*********************************************************************/

    /**
     * \brief Implementation of AttributeStoreCreator for a specific type.
     * \tparam T type of the elements
     */
    template< class T >
    class TypedAttributeStoreCreator: public AttributeStoreCreator {
    public:
        /**
         * \copydoc AttributeStoreCreator::create_attribute_store()
         */
        virtual AttributeStore* create_attribute_store( index_t dim )
        {
            return new TypedAttributeStore< T >( dim );
        }
    };

    /*********************************************************************/

    /**
     * \brief Helper class to register new attribute types
     * \tparam T attribute element type
     */
    template< class T > class geo_register_attribute_type {
    public:
        /**
         * \brief geo_register_attribute_type constructor
         * \param[in] type_name a const reference to a string with
         *  the C++ type name.
         * \details If the attribute is already registered with the same
         *  \p type_name and same \p T, then a warning message is issued.
         *  If the attribute is already registered with the same \p type_name
         *  but a different \p T, then an assertion failure is triggered.
         */
        geo_register_attribute_type( const std::string& type_name )
        {
//            AttributeStore::register_attribute_creator(
//                new TypedAttributeStoreCreator<T>, type_name, typeid( T ).name()
//                );
//            if( type_name == "bool" ) {
//                GEO::GeoFile::register_ascii_attribute_serializer(
//                    type_name,
//                    read_ascii_attribute<bool>,
//                    write_ascii_attribute<bool>
//                    );
//            } else {
//                GEO::GeoFile::register_ascii_attribute_serializer(
//                    type_name,
//                    read_ascii_attribute<T>,
//                    write_ascii_attribute<T>
//                    );
//            }
        }
    };

    /*********************************************************************/

    /**
     * \brief Managers a set of attributes attached to
     *  an object.
     */
    class RINGMESH_API AttributesManager {
    public:
        /**
         * \brief Constructs a new empty AttributesManager.
         */
        AttributesManager() = default;

        /**
         * \brief AttributesManager destructor.
         */
        ~AttributesManager();

        /**
         * \brief Gets the number of attributes.
         * \return The number of attributes managed by this
         *   AttributesManager.
         */
        index_t nb() const
        {
            return index_t( attributes_.size() );
        }

        /**
         * \brief Gets the names of all the attributes in this
         *   AttributeStore.
         * \param[out] names a vector of all attribute names
         */
        void list_attribute_names( std::vector< std::string >& names ) const;

        /**
         * \brief Gets the size.
         * \details All attributes stored in an AttributesManager have
         *  the same number of items.
         * \return the number of items of each attribute.
         */
        index_t size() const
        {
            return size_;
        }

        /**
         * \brief Resizes all the attributes managed by this
         *  AttributesManager.
         * \param[in] new_size the new number of items for
         *  all attributes.
         */
        void resize( index_t new_size );

        /**
         * \brief Clears this AttributesManager
         * \param[in] keep_attributes if true, then all
         *  attributes are resized to 0 but their names are
         *  kept.
         * \param[in] keep_memory if true, allocated memory
         *  is kept reserved.
         */
        void clear( bool keep_attributes, bool keep_memory = false );

        /**
         * \brief Zeroes all the attributes.
         */
        void zero();

        /**
         * \brief Binds an AttributeStore with the specified name.
         *  Ownership of this AttributeStore is transfered to
         *  the AttributesManager.
         * \param[in] name the name
         * \param[in] as a pointer to the AttributeStore to be bound
         * \pre No AttributeStore is already bound to the same name
         */
        void bind_attribute_store( const std::string& name, AttributeStore* as );

        /**
         * \brief Finds an AttributeStore by name.
         * \param[in] name the name under which the AttributeStore was bound
         * \return a pointer to the attribute store or nullptr if is is undefined.
         */
        AttributeStore* find_attribute_store( const std::string& name );

        /**
         * \brief Finds an AttributeStore by name.
         * \param[in] name the name under which the AttributeStore was bound
         * \return a const pointer to the attribute store or nullptr if is is
         *  undefined.
         */
        const AttributeStore* find_attribute_store( const std::string& name ) const;

        /**
         * \brief Tests whether an attribute is defined.
         * \param[in] name name of the attribute
         * \retval true if an attribute with the specified name exists
         * \retval false otherwise
         */
        bool is_defined( const std::string& name )
        {
            return ( find_attribute_store( name ) != nullptr );
        }

        /**
         * \brief Deletes an AttributeStore by name.
         * \param[in] name the name of the attribute store
         *  to be deleted.
         */
        void delete_attribute_store( const std::string& name );

        /**
         * \brief Deletes an AttributeStore.
         * \param[in] as a pointer to the attribute store
         *  to be deleted.
         */
        void delete_attribute_store( AttributeStore* as );

        /**
         * \brief Applies a permutation to the stored attributes.
         * \details Applying a permutation to the data is equivalent
         *  to:
         * \code
         * for(i=0; i<permutation.size(); i++) {
         *    data2[i] = data[permutation[i]]
         * }
         * data = data2 ;
         * \endcode
         * But it is done in-place.
         * \param[in] permutation the permutation.
         *  It is temporarily changed during execution of the
         *  function, but identical to the input on exit.
         */
        void apply_permutation( const std::vector< index_t >& permutation );

        /**
         * \brief Compresses the stored attributes, by
         *  applying an index mapping that fills-in the gaps.
         * \details This is equivalent to:
         * \code
         * for(i=0; i<size(); i++) {
         *    if(old2new[i] != index_t(-1)) {
         *       data2[old2new[i]] = data[i];
         *    }
         * }
         * data = data2 ;
         * \endcode
         * \param[in] old2new the index mapping to be applied.
         * \pre old2new[i] <= i || old2new[i] == index_t(-1)
         */
        void compress( const std::vector< index_t >& old2new );

        /**
         * \brief Copies all the attributes from another AttributesManager.
         * \details Previous content of this AttributesManager is erased.
         */
        void copy( const AttributesManager& rhs );

        /**
         * \brief Copies all the attributes of an item into another one.
         * \param[in] to index of the destination item
         * \param[in] from index of the source item
         * \note This function is not efficient.
         */
        void copy_item( index_t to, index_t from );

    private:
        /**
         * \brief Forbids copy.
         * \details This is to make sure that client code does
         *   not unintentionlly copies an AttributesManager (for
         *   instance by passing it by-value to a function).
         *   Use copy() instead.
         */
        AttributesManager( const AttributesManager& rhs );

        /**
         * \brief Forbids copy.
         * \details This is to make sure that client code does
         *   not unintentionlly copies an AttributesManager (for
         *   instance by passing it by-value to a function).
         *   Use copy() instead.
         */
        const AttributesManager& operator=( const AttributesManager& rhs );

    private:
        index_t size_ { 0 };
        std::map< std::string, AttributeStore* > attributes_;
    };

    /*********************************************************************/

    /**
     * \brief Base class for Attributes, that manipulates an
     *  attribute stored in an AttributesManager.
     */
    template< class T > class AttributeBase {
    public:

        /**
         * \brief Creates an unitialized (unbound) Attribute.
         */
        AttributeBase()
            : manager_( nullptr ), store_( nullptr )
        {
        }

        /**
         * \brief Creates or retreives a persistent attribute attached to
         *  a given AttributesManager.
         * \details If the attribute already exists with the specified
         *  name in the AttributesManager then it is retreived, else
         *  it is created and bound to the name.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name name of the attribute
         */
        AttributeBase( AttributesManager& manager, const std::string& name )
            : manager_( nullptr ), store_( nullptr )
        {
            bind( manager, name );
        }

        /**
         * \brief Tests whether an Attribute is bound.
         * \retval true if this Attribute is bound
         * \retval false otherwise
         */
        bool is_bound() const
        {
            return ( store_ != nullptr );
        }

        /**
         * \brief Unbinds this Attribute.
         * \pre is_bound()
         */
        void unbind()
        {
            ringmesh_assert( is_bound() );
            manager_ = nullptr;
            store_ = nullptr;
        }

        /**
         * \brief Binds this Attribute to an AttributesManager.
         * \details If the attribute already exists with the specified
         *  name in the AttributesManager then it is retreived, else
         *  it is created and bound to the name.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name name of the attribute
         * \pre !is_bound()
         */
        void bind( AttributesManager& manager, const std::string& name )
        {
            ringmesh_assert( !is_bound() );
            manager_ = &manager;
            store_ = manager_->find_attribute_store( name );
            if( store_ == nullptr ) {
                store_ = new VectorStore< T >();
                manager_->bind_attribute_store( name, store_ );
            } else {
                ringmesh_assert( store_->elements_type_matches( typeid(T).name() ) );
            }
        }

        /**
         * \brief Binds this Attribute to an AttributesManager if it
         *  already exists in the AttributesManager.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name name of the attribute
         * \pre !is_bound()
         */
        void bind_if_is_defined(
            AttributesManager& manager,
            const std::string& name )
        {
            ringmesh_assert( !is_bound() );
            manager_ = &manager;
            store_ = manager_->find_attribute_store( name );
            if( store_ != nullptr ) {
                ringmesh_assert( store_->elements_type_matches( typeid(T).name() ) );
            }
        }

        /**
         * \brief Creates and binds a new vector attribute.
         * \param[in] manager the attribute manager
         * \param[in] name the name of the attribute
         * \param[in] dimension the number of elements per item
         */
        void create_vector_attribute(
            AttributesManager& manager,
            const std::string& name,
            index_t dimension )
        {
            ringmesh_assert( !is_bound() );
            manager_ = &manager;
            ringmesh_assert( manager_->find_attribute_store( name ) == nullptr );
            store_ = new VectorStore< T >( dimension );
            manager_->bind_attribute_store( name, store_ );
        }

        /**
         * \brief Destroys this attribute in the AttributesManager.
         * \details On exit, the attribute is no-longer accessible in
         *  the AttributesManager, its name is available again, and
         *  this attribute is in the unbound state.
         */
        void destroy()
        {
            ringmesh_assert( is_bound() );
            manager_->delete_attribute_store( store_ );
            store_ = nullptr;
            manager_ = nullptr;
        }

        index_t dimension() const
        {
            return store_->dimension();
        }

        /**
         * \brief Sets the dimension.
         * \details The dimension is 1 for standard attributes and
         *  can be greater for vector attributes. The existing
         *  fields are kept. If the new dimension is greater than
         *  the old one, then new fields are initialized to the default
         *  value for the attribute type.
         * \param[in] new_dim the new dimension
         */
        void redim( index_t new_dim )
        {
            ringmesh_assert( is_bound() );
            store_->redim( new_dim );
        }

        /**
         * \brief Attribute destructor
         * \details
         *  The attribute is not destroyed, it can be retreived later
         *  by binding with the same name. To destroy the attribute,
         *  use detroy() instead.
         */
        ~AttributeBase()
        {
            if( is_bound() ) {
                unbind();
            }
        }

        /**
         * \brief Tests whether an attribute with the specified name and with
         *  corresponding type exists in an AttributesManager.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name the name of the attribute
         * \param[in] dim dimension, or 0 if any dimension can match
         */
        static bool is_defined(
            AttributesManager& manager,
            const std::string& name,
            index_t dim = 0 )
        {
            AttributeStore* store = manager.find_attribute_store( name );
            return ( store != nullptr
                && store->elements_type_matches( typeid(T).name() )
                && ( ( dim == 0 ) || ( store->dimension() == dim ) ) );
        }

        /**
         * \brief Gets the size.
         * \return The number of items in this attribute.
         */
        index_t size() const
        {
            return store_->size();
        }

        /**
         * \brief Sets all the elements of this Attribute to
         *   zero.
         */
        void zero()
        {
            ringmesh_assert( is_bound() );
            store_->zero();
        }

        /**
         * \brief Gets the AttributeManager this Attribute is bound to.
         * \return a pointer to the attributes manager.
         */
        AttributesManager* manager() const
        {
            return manager_;
        }

    protected:
        AttributesManager* manager_;
        AttributeStore* store_;
    };

    /*********************************************************************/

    /**
     * \brief Manages an attribute attached to a set of object.
     * \tparam T type of the attributes. Needs to be a basic type
     *  or a plain ordinary datatype (classes that do dynamic
     *  memory allocation are not allowed here).
     */
    template< class T > class Attribute: public AttributeBase< T > {
    public:
        typedef AttributeBase< T > superclass;

        /**
         * \brief Creates an unitialized (unbound) Attribute.
         */
        Attribute()
            : superclass()
        {
        }

        /**
         * \brief Creates or retreives a persistent attribute attached to
         *  a given AttributesManager.
         * \details If the attribute already exists with the specified
         *  name in the AttributesManager then it is retreived, else
         *  it is created and bound to the name.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name name of the attribute
         */
        Attribute( AttributesManager& manager, const std::string& name )
            : superclass( manager, name )
        {
        }

        /**
         * \brief Gets a modifiable element by index
         * \param [in] i index of the element
         * \return a modifiable reference to the \p i%th element
         */
        T& operator[]( unsigned int i )
        {
            ringmesh_assert( i < superclass::nb_elements() );
            return ( (T*) superclass::store_->data() )[i];
        }

        /**
         * \brief Gets an element by index
         * \param [in] i index of the element
         * \return a const reference to the \p i%th element
         */
        const T& operator[]( unsigned int i ) const
        {
            ringmesh_assert( i < superclass::nb_elements() );
            return ( (const T*) superclass::store_->data() )[i];
        }

        /**
         * \brief Sets all the elements in this attribute
         *   to a specified value.
         * \param[in] val the value
         */
        void fill( const T& val )
        {
            for( index_t i = 0; i < superclass::nb_elements(); ++i ) {
                ( *this )[i] = val;
            }
        }

    private:
        /**
         * \brief Forbids copy.
         */
        Attribute( const Attribute< T >& rhs );
        /**
         * \brief Forbids copy.
         */
        Attribute< T >& operator=( const Attribute< T >& rhs );
    };
    /**
     * \brief Specialization of Attribute for booleans
     * \details Attribute needs a specialization for bool, since
     *   vector<bool> uses compressed storage (1 bit per boolean),
     *   that is not compatible with the attribute management
     *   mechanism. This wrapper class uses an Attribute<Numeric::uint8>
     *   and does the appropriate conversions, using an accessor class.
     */
    template <> class Attribute<bool> : public AttributeBase<Byte> {
    public:
        typedef AttributeBase<Byte> superclass;

        Attribute() : superclass() {
        }

        Attribute(AttributesManager& manager, const std::string& name) :
            superclass(manager,name) {
        }

        class BoolAttributeAccessor;


        /**
         * \brief Accessor class for adapting Attribute<bool>
         *  indexing.
         */
        class ConstBoolAttributeAccessor {
        public:
            /**
             * \brief ConstBoolAttributeAccessor constructor.
             */
            ConstBoolAttributeAccessor(
                const Attribute<bool>& attribute,
                index_t index
            ) :
                attribute_(&attribute),
                index_(index) {
            }

            /**
             * \brief Converts a BoolAttributeAccessor to a bool.
             * \details Performs the actual lookup.
             */
            operator bool() const {
                return (attribute_->element(index_) != 0);
            }

        private:
            const Attribute<bool>* attribute_;
            index_t index_;

            friend class BoolAttributeAccessor;
        };

        /**
         * \brief Accessor class for adapting Attribute<bool>
         *  indexing.
         */
        class BoolAttributeAccessor {
        public:
            /**
             * \brief BoolAttributeAccessor constructor.
             */
            BoolAttributeAccessor(
                Attribute<bool>& attribute,
                index_t index
            ) :
                attribute_(&attribute),
                index_(index) {
            }

            /**
             * \brief Converts a BoolAttributeAccessor to a bool.
             * \details Performs the actual lookup.
             */
            operator bool() const {
                return (attribute_->element(index_) != 0);
            }

            /**
             * \brief Copy-constructor.
             * \param[in] rhs a const reference to the
             *  BoolAttributeAccessor to be copied.
             */
            BoolAttributeAccessor(const BoolAttributeAccessor& rhs) {
                attribute_ = rhs.attribute_;
                index_ = rhs.index_;
            }

            /**
             * \brief Assigns a bool to a BoolAttributeAccessor.
             * \details Stores the boolean into the Attribute.
             */
            BoolAttributeAccessor& operator=(bool x) {
                attribute_->element(index_) = Byte(x);
                return *this;
            }

            /**
             * \brief Copies a bool from another attribute.
             * \param[in] rhs a const reference to the BoolAttributeAccessor
             *  to be copied.
             */
            BoolAttributeAccessor& operator=(
                const BoolAttributeAccessor& rhs
            ) {
                if(&rhs != this) {
                    attribute_->element(index_) =
                        rhs.attribute_->element(rhs.index_);
                }
                return *this;
            }

            /**
             * \brief Copies a bool from another attribute.
             * \param[in] rhs a const reference to the
             *  ConstBoolAttributeAccessor to be copied.
             */
            BoolAttributeAccessor& operator=(
                const ConstBoolAttributeAccessor& rhs
            ) {
                attribute_->element(index_) =
                    rhs.attribute_->element(rhs.index_);
                return *this;
            }

        private:
            Attribute<bool>* attribute_;
            index_t index_;
        };


        BoolAttributeAccessor operator[](index_t i) {
            return BoolAttributeAccessor(*this,i);
        }

        ConstBoolAttributeAccessor operator[](index_t i) const {
            return ConstBoolAttributeAccessor(*this,i);
        }

        /**
         * \brief Sets all the elements in this attribute
         *   to a specified value.
         * \param[in] val the value
         */
        void fill(bool val) {
            for(index_t i=0; i<superclass::store_->nb_elements(); ++i) {
                element(i) = Byte(val);
            }
        }

    protected:

        friend class BoolAttributeAccessor;
        friend class ConstBoolAttributeAccessor;

        /**
         * \brief Gets a modifiable element by index
         * \param [in] i index of the element
         * \return a modifiable reference to the \p i%th element
         */
        Byte& element(unsigned int i) {
            geo_debug_assert(i < superclass::store_->nb_elements());
            return superclass::store_->data()[i];
        }

        /**
         * \brief Gets an element by index
         * \param [in] i index of the element
         * \return a const reference to the \p i%th element
         */
        const Byte& element(unsigned int i) const {
            geo_debug_assert(i < superclass::store_->nb_elements());
            return superclass::store_->data()[i];
        }

    private:
        /**
         * \brief Forbids copy.
         */
        Attribute(const Attribute<bool>& rhs);
        /**
         * \brief Forbids copy.
         */
        Attribute<bool>& operator=(const Attribute<bool>& rhs);
    } ;

    /***********************************************************/

    /**
     * \brief Access to an attribute as a double regardless its type.
     * \details The attribute can be an element of a vector attribute.
     */
    class RINGMESH_API ReadOnlyScalarAttributeAdapter {

    public:
        /**
         * \brief Internal representation of the attribute.
         */
        enum ElementType {
            ET_NONE = 0,
            ET_UINT8 = 1,
            ET_INT8 = 2,
            ET_UINT32 = 3,
            ET_INT32 = 4,
            ET_FLOAT32 = 5,
            ET_FLOAT64 = 6,
            ET_VEC2 = 7,
            ET_VEC3 = 8
        };

        /**
         * \brief ReadOnlyScalarAttributeAdapter constructor.
         */
        ReadOnlyScalarAttributeAdapter()
            :
                manager_( nullptr ),
                store_( nullptr ),
                element_type_( ET_NONE ),
                element_index_( index_t( -1 ) )
        {
        }

        /**
         * \brief ReadOnlyScalarAttributeAdapter constructor.
         * \details Retreives a persistent attribute attached to
         *  a given AttributesManager.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name name of the attribute with an optional index,
         *   for instance, "foobar[5]" refers to the 5th coordinate of
         *   the "foobar" vector attribute.
         */
        ReadOnlyScalarAttributeAdapter(
            const AttributesManager& manager,
            const std::string& name )
            : manager_( nullptr ), store_( nullptr )
        {
            bind_if_is_defined( manager, name );
        }

        /**
         * \brief Tests whether an Attribute is bound.
         * \retval true if this Attribute is bound
         * \retval false otherwise
         */
        bool is_bound() const
        {
            return ( store_ != nullptr );
        }

        /**
         * \brief Unbinds this Attribute.
         * \pre is_bound()
         */
        void unbind()
        {
            ringmesh_assert( is_bound() );
            manager_ = nullptr;
            store_ = nullptr;
            element_type_ = ET_NONE;
            element_index_ = index_t( -1 );
        }

        /**
         * \brief Binds this Attribute to an AttributesManager if it
         *  already exists in the AttributesManager.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name name of the attribute with an optional index,
         *   for instance, "foobar[5]" refers to the 5th coordinate of
         *   the "foobar" vector attribute.
         * \pre !is_bound()
         */
        void bind_if_is_defined(
            const AttributesManager& manager,
            const std::string& name );

        /**
         * \brief ReadonlyScalarAttributeAdapter destructor
         * \details
         *  The attribute is not destroyed, it can be retreived later
         *  by binding with the same name. To destroy the attribute,
         *  use detroy() instead.
         */
        ~ReadOnlyScalarAttributeAdapter()
        {
            if( is_bound() ) {
                unbind();
            }
        }

        /**
         * \brief Tests whether an attribute with the specified name and with
         *  a type that can be converted to double exists in an
         *  AttributesManager.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name the name of the attribute with an optional index,
         *   for instance, "foobar[5]" refers to the 5th coordinate of
         *   the "foobar" vector attribute.
         */
        static bool is_defined(
            const AttributesManager& manager,
            const std::string& name );

        /**
         * \brief Gets the size.
         * \return The number of items in this attribute.
         */
        index_t size() const
        {
            return ( store_ == nullptr ) ? 0 : store_->size();
        }

        /**
         * \brief Gets the internal representation of the
         *  elements.
         * \return one of ET_NONE (if unbound), ET_UINT8,
         *  ET_INT8, ET_UINT32, ET_INT32, ET_FLOAT32,
         *  ET_FLOAT64, ET_VEC2, ET_VEC3
         */
        ElementType element_type() const
        {
            return element_type_;
        }

        /**
         * \brief Gets the element index.
         * \return the index of the elements accessed in
         *  the bound attribute, or 0 if the bound attribute
         *  is scalar.
         */
        index_t element_index() const
        {
            return element_index_;
        }

        /**
         * \brief Gets the AttributeStore.
         * \return a const pointer to the AttributeStore.
         */
        const AttributeStore* attribute_store() const
        {
            return store_;
        }

        /**
         * \brief Gets a property value
         * \param[in] i the index of the item
         * \return the value of the property, convertex
         *  into a double
         * \pre is_bound() && i < size()
         */
        double operator[]( index_t i )
        {
            double result = 0.0;
            switch( element_type_ ) {
                case ET_UINT8:
                    result = get_element< uint8 >( i );
                    break;
                case ET_INT8:
                    result = get_element< int8 >( i );
                    break;
                case ET_UINT32:
                    result = get_element< uint32 >( i );
                    break;
                case ET_INT32:
                    result = get_element< int32 >( i );
                    break;
                case ET_FLOAT32:
                    result = get_element< float32 >( i );
                    break;
                case ET_FLOAT64:
                    result = get_element< float64 >( i );
                    break;
                case ET_VEC2:
                    result = get_element< float64 >( i, 2 );
                    break;
                case ET_VEC3:
                    result = get_element< float64 >( i, 3 );
                    break;
                case ET_NONE:
                    ringmesh_assert_not_reached;
            }
            return result;
        }

        /**
         * \brief Tests whether a ReadOnlyScalarAttributeAdapter can
         *  be bound to a given attribute store.
         * \param[in] store a pointer to the attribute store.
         * \retval true if it can be bound
         * \retval false otherwise
         */
        static bool can_be_bound_to( const AttributeStore* store )
        {
            return element_type( store ) != ET_NONE;
        }

        /**
         * \brief Gets the number of scalar components per item in an
         *  AttributeStore.
         * \param[in] store a pointer to the attribute store.
         * \return the number of scalar components per item in an
         *  AttributeStore.
         */
        static index_t nb_scalar_elements_per_item( const AttributeStore* store );

    protected:
        /**
         * \brief Gets the base attribute name from a compound name.
         * \param[in] name the string with the attribute name and optional
         *  index. For instance, "foobar[5]" refers to the 5th coordinate of
         *  the "foobar" vector attribute.
         * \return the attribute name. For instance, for "foobar[5]", it
         *  returns "foobar".
         */
        static std::string attribute_base_name( const std::string& name );

        /**
         * \brief Gets the base attribute name from a compound name.
         * \param[in] name the string with the attribute name and optional
         *  index. For instance, "foobar[5]" refers to the 5th coordinate of
         *  the "foobar" vector attribute.
         * \return the index or zero if no index was specified. For instance,
         *  for "foobar[5]" it returns 5, and for "foobar" it returns 0
         */
        static index_t attribute_element_index( const std::string& name );

        /**
         * \brief Gets the element type stored in an AttributeStore.
         * \param[in] store a const pointer to the AttributeStore
         * \return one of ET_UINT8, ET_INT8, ET_UINT32, ET_INT32,
         *  ET_FLOAT32, ET_FLOAT64 if the type of the attribute is
         *  compatible with those types, or ET_NONE if it is incompatible.
         */
        static ElementType element_type( const AttributeStore* store );

        /**
         * \brief Gets an element.
         * \param[in] i index of the element
         * \param[in] multiplier multiplier applied to the index before fetching
         *  the raw pointer.
         */
        template< class T > double get_element( index_t i, index_t multiplier = 1 ) const
        {
            ringmesh_assert( is_bound() );
            ringmesh_assert( i < size() );
            return double(
                reinterpret_cast< const T* >( store_->data() )[( i * store_->dimension()
                    * multiplier ) + element_index_] );
        }

    private:
        const AttributesManager* manager_;
        const AttributeStore* store_;
        ElementType element_type_;
        index_t element_index_;
    };

/***********************************************************/
}

