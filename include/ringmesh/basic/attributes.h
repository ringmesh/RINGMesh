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

#include <geogram/basic/geofile.h>
#include <map>
#include <ringmesh/basic/common.h>
#include <ringmesh/basic/factory.h>
#include <ringmesh/basic/logger.h>
#include <ringmesh/basic/logger.h>
#include <set>
#include <typeinfo>

/**
 * \file ringmesh/basic/attributes.h
 * \brief Generic mechanism for attributes.
 * This file has been strongly inspired from the attributes of the geogram
 * library.
 * Many thanks to Bruno Levy (Bruno.Levy@inria.fr) who did the first
 * implementation in Geogram.
 */

namespace RINGMesh
{
    class AttributeStore;

    class RINGMESH_API Store
    {
    public:
        Store( index_t elemsize ) : element_size_( elemsize )
        {
        }

        virtual ~Store() = default;
        /**
        * \brief Gets the size.
        * \return the number of items
        */
        index_t size() const
        {
            return size_;
        }

        /**
        * \brief Gets the total number of elements.
        * \details This corresponds to one position past
        *  the last valid index.
        * \return the total number of elements.
        */
        index_t nb_elements() const
        {
            return size_;
        }
        virtual const void* data() const = 0;
        virtual void* data() = 0;

        virtual const void* element( index_t element ) const = 0;
        virtual void* element( index_t element ) = 0;

        virtual void compress( const std::vector< index_t >& old2new ) = 0;
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
        virtual void apply_permutation(
            const std::vector< index_t >& permutation ) = 0;
        /**
        * \brief Copies an item
        * \param[in] to index of the destination item
        * \param[in] from index of the source item
        */
        virtual void copy_item( index_t to, index_t from ) = 0;

        virtual bool elements_type_matches(
            const std::string& type_name ) const = 0;
        virtual std::string element_typeid_name() const = 0;

        virtual void resize( index_t new_size ) = 0;

        virtual void clear() = 0;

        virtual Store* clone() const = 0;

        index_t element_size() const
        {
            return element_size_;
        }

    protected:
        index_t size_{ 0 };
        index_t element_size_; // size of element in byte
    };
    /**
    * \brief Stores an array of elements of a given type.
    */
    template < class T >
    class TypedStore : public Store
    {
    public:
        /**
        * \brief Creates a new empty attribute store.
        */
        TypedStore() : Store( index_t( sizeof( T ) ) )
        {
        }

        bool elements_type_matches( const std::string& type_name ) const final
        {
            return type_name == element_typeid_name();
        }

        std::string element_typeid_name() const override
        {
            return typeid( T ).name();
        }
    };

    template < class T >
    class VectorStoreBase : public TypedStore< T >
    {
    public:
        /**
        * \brief Creates a new empty attribute store.
        */
        VectorStoreBase() = default;

        void resize( index_t new_size ) final
        {
            this->size_ = new_size;
            vector_.resize( new_size );
        }

        const void* data() const final
        {
            return vector_.data();
        }

        void* data() final
        {
            return vector_.data();
        }

        const void* element( index_t element ) const final
        {
            return &this->vector_[element];
        }
        void* element( index_t element ) final
        {
            return &this->vector_[element];
        }
        void clear() final
        {
            vector_.clear();
        }

        void compress( const std::vector< index_t >& /*old2new*/ ) final
        {
            return;
        }
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
        void apply_permutation(
            const std::vector< index_t >& permutation ) final
        {
            std::vector< T > tmp_vector;
            tmp_vector.resize( vector_.size() );
            for( auto p : range( permutation.size() ) )
            {
                tmp_vector[p] = vector_[permutation[p]];
            }
            vector_ = tmp_vector;
        }
        /**
        * \brief Copies an item
        * \param[in] to index of the destination item
        * \param[in] from index of the source item
        */
        void copy_item( index_t to, index_t from ) final
        {
            ringmesh_assert( to < vector_.size() && from < vector_.size() );
            vector_[from] = vector_[to];
        }

    protected:
        std::vector< T > vector_;
    };

    template < class T >
    class VectorStore : public VectorStoreBase< T >
    {
    public:
        VectorStore() = default;

        Store* clone() const final
        {
            VectorStore< T >* result = new VectorStore< T >();
            result->resize( this->size() );
            result->vector_ = this->vector_;
            return result;
        }
    };

    template <>
    class VectorStore< bool > : public VectorStoreBase< Byte >
    {
    public:
        /**
        * \brief Creates a new empty attribute store.
        */
        VectorStore() = default;

        Store* clone() const final
        {
            VectorStore< bool >* result = new VectorStore< bool >();
            result->resize( this->size() );
            result->vector_ = vector_;
            return result;
        }
        std::string element_typeid_name() const final
        {
            return std::string( "bool" );
        }
    };

    template < class T >
    class ConstantStore : public TypedStore< T >
    {
    public:
        /**
        * \brief Creates a new empty attribute store.
        */
        ConstantStore()
        {
            this->size_ = 1;
        }

        void resize( index_t ) override
        {
            Logger::warn( "Trying to resize constant attribute" );
        }

        const void* data() const final
        {
            return &constant_;
        }
        void* data() final
        {
            return &constant_;
        }
        const void* element( index_t element ) const final
        {
            ringmesh_unused( element );
            return &this->constant_;
        }
        void* element( index_t element ) final
        {
            ringmesh_unused( element );
            return &this->constant_;
        }
        void clear() override
        {
        }

        Store* clone() const override
        {
            ConstantStore< T >* result = new ConstantStore< T >();
            result->constant_ = constant_;
            return result;
        }
        void compress( const std::vector< index_t >& old2new ) final
        {
            return;
        }
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
        void apply_permutation(
            const std::vector< index_t >& permutation ) final
        {
            return;
        }
        /**
        * \brief Copies an item
        * \param[in] to index of the destination item
        * \param[in] from index of the source item
        */
        void copy_item( index_t to, index_t from ) final
        {
            return;
        }

    protected:
        T constant_;
    };

    /**
     * \brief Internal class for creating an AttributeStore
     *  from the type name of its elements.
     */
    class RINGMESH_API AttributeStoreCreator
    {
    public:
        /**
         * \brief AttributeStoreCreator destructor.
         */
        virtual ~AttributeStoreCreator() = default;

        /**
         * \brief Creates a new attribute store.
         * \return a pointer to the newly created AttributeStore
         */
        virtual AttributeStore* create_attribute_store() = 0;

    private:
        std::string element_type_name_;
        std::string element_typeid_name_;
    };

    class RINGMESH_API AttributeStore
    {
    public:
        /**
         * \brief AttributeStore constructor.
         * \param[in] elemsize size of one element,
         *  in bytes.
         */
        AttributeStore() = default;

        /**
         * \brief AttributeStore destructor.
         */
        ~AttributeStore();

        void set_store( Store* store )
        {
            store_ = std::unique_ptr< Store >( store );
        }

        const Store& get_store() const
        {
            return *store_;
        }

        /**
         * \brief Resizes this AttributeStore
         * \param[in] new_size new number of items
         */
        void resize( index_t new_size )
        {
            store_->resize( new_size );
        }

        /**
         * \brief Resizes this AttributeStore to 0.
         */
        void clear()
        {
            store_->clear();
        }

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
        void apply_permutation( const std::vector< index_t >& permutation )
        {
            store_->apply_permutation( permutation );
        }

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
        void compress( const std::vector< index_t >& old2new )
        {
            store_->compress( old2new );
        }

        /**
         * \brief Creates a new AttributeStore that is a carbon copy
         *  of this AttributeStore.
         * \details Only the data is copied.
         */
        AttributeStore* clone() const
        {
            AttributeStore* new_attstore =
                create_attribute_store_by_element_type_name(
                    store_->element_typeid_name() );

            new_attstore->set_store( store_->clone() );

            return new_attstore;
        }

        /**
         * \brief Gets a pointer to the stored data.
         * \return A pointer to the memory block
         */
        void* data()
        {
            return store_->data();
        }

        /**
         * \brief Gets a pointer to the stored data.
         * \return A const pointer to the memory block
         */
        void* data() const
        {
            return store_->data();
        }

        const void* element( index_t element ) const
        {
            return store_->element( element );
        }

        void* element( index_t element )
        {
            return store_->element( element );
        }

        bool elements_type_matches( const std::string& type_name ) const
        {
            return store_->elements_type_matches( type_name );
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
         */
        static AttributeStore* create_attribute_store_by_element_type_name(
            const std::string& element_type_name )
        {
            ringmesh_assert( element_type_name_is_known( element_type_name ) );
            return type_name_to_creator_[element_type_name]
                ->create_attribute_store();
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
            ringmesh_assert(
                element_typeid_name_is_known( element_typeid_name ) );
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
        static void register_attribute_creator( AttributeStoreCreator* creator,
            const std::string& element_type_name,
            const std::string& element_typeid_name )
        {
            if( element_type_name_is_known( element_type_name ) )
            {
                Logger::warn(
                    "Attributes", element_type_name, " already registered" );
                if( element_typeid_name_is_known( element_typeid_name ) )
                {
                    bool already_registered_attribute_has_same_type =
                        ( type_name_to_typeid_name_[element_type_name]
                            == element_typeid_name );
                    ringmesh_assert(
                        already_registered_attribute_has_same_type );
                    ringmesh_unused(
                        already_registered_attribute_has_same_type );
                }
            }
            type_name_to_creator_[element_type_name] =
                std::unique_ptr< AttributeStoreCreator >( creator );
            typeid_name_to_type_name_[element_typeid_name] = element_type_name;
            type_name_to_typeid_name_[element_type_name] = element_typeid_name;
        }

    protected:
        std::unique_ptr< Store > store_{ nullptr };

        static std::map< std::string, std::unique_ptr< AttributeStoreCreator > >
            type_name_to_creator_;

        static std::map< std::string, std::string > typeid_name_to_type_name_;

        static std::map< std::string, std::string > type_name_to_typeid_name_;
    };

    /*********************************************************************/

    /**
     * \brief Managers a set of attributes attached to
     *  an object.
     */
    class RINGMESH_API AttributesManager
    {
        ringmesh_disable_copy_and_move( AttributesManager );

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
        index_t nb_attributes() const
        {
            return static_cast< index_t >( attributes_.size() );
        }

        /**
         * \brief Gets the names of all the attributes in this
         *   AttributeStore.
         * \return a vector of all attribute names
         */
        std::vector< std::string > attribute_names() const;

        /**
         * \brief Gets the number of items in an Attribute.
         * \details All attributes stored in an AttributesManager have
         *  the same number of items.
         * \return the number of items of each attribute.
         */
        index_t nb_items() const
        {
            return nb_items_;
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
         */
        void clear( bool keep_attributes );

        /**
         * \brief Binds an AttributeStore with the specified name.
         *  Ownership of this AttributeStore is transfered to
         *  the AttributesManager.
         * \param[in] name the name
         * \param[in] as a pointer to the AttributeStore to be bound
         * \pre No AttributeStore is already bound to the same name
         */
        void bind_attribute_store(
            const std::string& name, AttributeStore* as );

        /**
         * \brief Finds an AttributeStore by name.
         * \param[in] name the name under which the AttributeStore was bound
         * \return a pointer to the attribute store or nullptr if is is
         * undefined.
         */
        AttributeStore* find_attribute_store( const std::string& name );

        /**
         * \brief Finds an AttributeStore by name.
         * \param[in] name the name under which the AttributeStore was bound
         * \return a const pointer to the attribute store or nullptr if is is
         *  undefined.
         */
        const AttributeStore* find_attribute_store(
            const std::string& name ) const;

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
        index_t nb_items_{ 0 };
        std::map< std::string, AttributeStore* > attributes_;
    };

    /*********************************************************************/

    /**
     * \brief Base class for Attributes, that manipulates an
     *  attribute stored in an AttributesManager.
     */
    template < class T >
    class Attribute
    {
        ringmesh_disable_copy_and_move( Attribute );

    public:
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
            if( store_ == nullptr )
            {
                store_ = new AttributeStore;
                store_->set_store( new VectorStore< T > );
                DEBUG( manager_->nb_items() );
                store_->resize( manager_->nb_items() );
                manager_->bind_attribute_store( name, store_ );
            }
            else
            {
                ringmesh_assert( store_->get_store().elements_type_matches(
                    typeid( T ).name() ) );
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
            AttributesManager& manager, const std::string& name )
        {
            ringmesh_assert( !is_bound() );
            manager_ = &manager;
            store_ = manager_->find_attribute_store( name );
            if( store_ != nullptr )
            {
                ringmesh_assert(
                    store_->elements_type_matches( typeid( T ).name() ) );
            }
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

        void set_constant_value( T value )
        {
            store_->set_store( new ConstantStore< T > );
            set_value( 0, value );
        }

        /**
         * \brief Attribute destructor
         * \details
         *  The attribute is not destroyed, it can be retreived later
         *  by binding with the same name. To destroy the attribute,
         *  use detroy() instead.
         */
        ~Attribute()
        {
            if( is_bound() )
            {
                unbind();
            }
        }

        /**
         * \brief Tests whether an attribute with the specified name and with
         *  corresponding type exists in an AttributesManager.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name the name of the attribute
         */
        static bool is_defined(
            AttributesManager& manager, const std::string& name )
        {
            AttributeStore* store = manager.find_attribute_store( name );
            return ( store != nullptr
                     && store->elements_type_matches( typeid( T ).name() ) );
        }

        /**
         * \brief Gets the size.
         * \return The number of items in this attribute.
         */
        index_t size() const
        {
            return store_->get_store().size();
        }

        /**
         * \brief Gets the AttributeManager this Attribute is bound to.
         * \return a pointer to the attributes manager.
         */
        AttributesManager* manager() const
        {
            return manager_;
        }

        /**
          * \brief Creates an unitialized (unbound) Attribute.
          */
        Attribute() = default;

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
        {
            bind( manager, name );
        }

        /**
         * \brief Gets a modifiable element by index
         * \param [in] i index of the element
         * \param [in] value to set at the \p i%th element
         */
        void set_value( index_t i, const T& value )
        {
            ringmesh_assert( i < size() );
            *reinterpret_cast< T* >( store_->element( i ) ) = value;
        }

        /**
         * \brief Gets an element by index
         * \param [in] i index of the element
         * \return a const reference to the \p i%th element
         */
        const T& operator[]( unsigned int i ) const
        {
            return value( i );
        }

        /**
        * \brief Gets an element by index
        * \param [in] i index of the element
        * \return a const reference to the \p i%th element
        */
        const T& value( unsigned int i ) const
        {
            ringmesh_assert( i < size() );
            return *reinterpret_cast< const T* >( store_->element( i ) );
        }

        /**
         * \brief Sets all the elements in this attribute
         *   to a specified value.
         * \param[in] val the value
         */
        void fill( const T& val )
        {
            for( index_t i = 0; i < store_->get_store().nb_elements(); ++i )
            {
                ( *this )[i] = val;
            }
        }

    protected:
        AttributesManager* manager_{ nullptr };
        AttributeStore* store_{ nullptr };
    };

    /***********************************************************/

    /**
     * \brief Access to an attribute as a double regardless its type.
     * \details The attribute can be an element of a vector attribute.
     */
    class RINGMESH_API ReadOnlyScalarAttributeAdapter
    {
    public:
        /**
         * \brief ReadOnlyScalarAttributeAdapter constructor.
         */
        ReadOnlyScalarAttributeAdapter() = delete;

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
            const AttributesManager& manager, const std::string& name )
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
            const AttributesManager& manager, const std::string& name );

        /**
         * \brief ReadonlyScalarAttributeAdapter destructor
         * \details
         *  The attribute is not destroyed, it can be retreived later
         *  by binding with the same name. To destroy the attribute,
         *  use detroy() instead.
         */
        virtual ~ReadOnlyScalarAttributeAdapter()
        {
            if( is_bound() )
            {
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
            const AttributesManager& manager, const std::string& name );

        /**
         * \brief Gets the size.
         * \return The number of items in this attribute.
         */
        index_t size() const
        {
            return ( store_ == nullptr ) ? 0 : store_->get_store().size();
        }

        /**
         * \brief Gets the internal representation of the
         *  elements.
         */
        // virtual AttributeElementType element_type() const = 0 ;

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
        virtual double operator[]( index_t i ) = 0;

        /**
         * \brief Gets the number of scalar components per item in an
         *  AttributeStore.
         * \param[in] store a pointer to the attribute store.
         * \return the number of scalar components per item in an
         *  AttributeStore.
         */
        virtual index_t nb_scalar_elements_per_item() const = 0;
        virtual bool is_integer_like_attribute() const = 0;

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
        // static ElementType element_type( const AttributeStore* store );

        /**
         * \brief Gets an element.
         * \param[in] i index of the element
         * \param[in] multiplier multiplier applied to the index before fetching
         *  the raw pointer.
         */
        template < class T >
        double get_element( index_t i, index_t multiplier = 1 ) const
        {
            ringmesh_assert( is_bound() );
            ringmesh_assert( i < size() );
            return double( reinterpret_cast< const T* >(
                store_->data() )[( i * multiplier ) + element_index_] );
        }

    protected:
        const AttributesManager* manager_;
        const AttributeStore* store_;
        index_t element_index_;
    };
    using ReadOnlyScalarAttributeAdapterFactory = Factory< std::string,
        ReadOnlyScalarAttributeAdapter,
        const AttributesManager&,
        const std::string& >;
    void register_read_only_scalar_attribute();
    /***********************************************************/
}
