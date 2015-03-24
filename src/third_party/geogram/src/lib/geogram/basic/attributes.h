/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef __GEOGRAM_BASIC_ATTRIBUTES__
#define __GEOGRAM_BASIC_ATTRIBUTES__


#include <geogram/basic/common.h>
#include <geogram/basic/memory.h>
#include <geogram/basic/numeric.h>
#include <geogram/basic/process.h>
#include <map>
#include <typeinfo>
#include <set>

/**
 * \file geogram/basic/attributes.h
 * \brief Generic mechanism for attributes.
 */

namespace GEO {

    class AttributeStore;
    
    /**
     * \brief Base class for attributes. They are notified
     *  whenever the AttributeStore is modified.
     */
    class GEOGRAM_API AttributeStoreObserver {
    public:

        /**
         * \brief Creates a new uninitialied AttributeStore.
         */
        AttributeStoreObserver() : base_addr_(nil), size_(0), dimension_(0) {
        }

        /**
         * \brief Callback function, called by the AttributeStore
         *  whenever it is modified.
         * \param[in] base_addr new base address of the AttributeStore
         * \param[in] size new number of items in the AttributeStore
         * \param[in] dim  new dimension, i.e. number of elements per item
         */
        void notify(
            Memory::pointer base_addr, index_t size, index_t dim
        ) {
            base_addr_ = base_addr;
            size_ = size;
            dimension_ = dim;
        }

        /**
         * \brief Gets the size.
         * \return the number of items
         */
        index_t size() const {
            return size_;
        }

        /**
         * \brief Gets the dimension.
         * \return the number of elements per item
         */
        index_t dimension() const {
            return dimension_;
        }

        /**
         * \brief Gets the total number of elements.
         * \details This corresponds to one position past
         *  the last valid index.
         * \return the total number of elements.
         */
        index_t nb_elements() const {
            return size_ * dimension_;
        }

        void register_me(AttributeStore* store);
        
        void unregister_me(AttributeStore* store);

        
    protected:
        Memory::pointer base_addr_;
        index_t size_;
        index_t dimension_;
    };


    /*********************************************************************/
    
    /**
     * \brief Notifies a set of AttributeStoreObservers 
     *  each time the stored array changes size and/or 
     *  base address and/or dimension.
     */
    class GEOGRAM_API AttributeStore {
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
        AttributeStore(index_t elemsize, index_t dim=1);

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
        virtual bool elements_type_matches(
            const std::string& type_name
        ) const = 0;

        /**
         * \brief Gets the size.
         * \return the number of items
         */
        index_t size() const {
            return cached_size_;
        }

        /**
         * \brief Resizes this AttributeStore
         * \param[in] new_size new number of items
         */
        virtual void resize(index_t new_size) = 0;


        /**
         * \brief Resizes this AttributeStore to 0.
         * \param[in] keep_memory if true, then memory
         *  is kept reserved for future use.
         */
        virtual void clear(bool keep_memory = false) = 0;
        
        /**
         * \brief Tests whether observers listen to this AttributeStore.
         * \retval true if at least one observer is bound to this AttributeStore
         * \retval false otherwise
         */
        bool has_observers() const {
            return !observers_.empty();
        }

        /**
         * \brief Gets the dimension.
         * \details The dimension is 1 for standard attributes and
         *  can be greater for vector attributes.
         */
        index_t dimension() const {
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
        virtual void redim(index_t dim) = 0;


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
         * \param[in] permutation_in the permutation.
         *  It is temporarily changed during execution of the
         *  function, but identical to the input on exit.
         * \note This function uses memcpy(). If required, it
         *  can be overloaded in derived classes.
         */
        virtual void apply_permutation(
            const vector<index_t>& permutation
        );

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
        virtual void compress(const vector<index_t>& old2new);

        /**
         * \brief Zeroes all the memory associated with this 
         *  AttributeStore.
         * \details Subclasses may overload this function for
         *  attributes that have non "plain ordinary datatypes"
         *  and that need a more elaborate initialization mechanism.
         */
        virtual void zero();
        
    protected:
        /**
         * \brief If size or base address differ from the
         *  cached values, notify all the observers, 
         *  and update the cached base address and size.
         * \param[in] base_addr the new base address
         * \param[in] size the new size
         * \param[in] dim the new dimension
         */
        void notify(Memory::pointer base_addr, index_t size, index_t dim);

        /**
         * \brief Registers an observer.
         * \details All the registered observers are notified whenever
         *  the size or base pointer in this AttributeStore change.
         *  The function is thread-safe.
         * \param[in] observer the AttributeStoreObserver to be
         *  registered.
         */
        void register_observer(AttributeStoreObserver* observer);

        /**
         * \brief Unregisters an observer.
         * \param[in] observer the AttributeStoreObserver to be
         *  unregistered.
         *  The function is thread-safe.
         * \pre \p observer is registered.
         */
        void unregister_observer(AttributeStoreObserver* observer);

        
    protected:
        index_t element_size_;
        index_t dimension_;        
        Memory::pointer cached_base_addr_;
        index_t cached_size_;
        std::set<AttributeStoreObserver*> observers_;
        Process::spinlock lock_;
        
        friend class AttributeStoreObserver;
    };

    /*********************************************************************/    
    
    /**
     * \brief Stores an array of elements of a given type, 
     *  and notifies a set of AttributeStoreObservers each time the
     *  storead array changes size and/or base address.
     */
    template <class T> class TypedAttributeStore : public AttributeStore {
    public:

        /**
         * \brief Creates a new empty attribute store.
         * \param[in] dim number of elements in each item,
         *  default value is 1, can be greater for vector
         *  attributes.
         */
        TypedAttributeStore(index_t dim=1) :
            AttributeStore(index_t(sizeof(T)),dim) {
        }

        virtual void resize(index_t new_size) {
            store_.resize(new_size*dimension_);
            notify(
                store_.empty() ? nil : Memory::pointer(store_.data()),
                new_size,
                dimension_
            );
        }

        virtual void clear(bool keep_memory=false) {
            if(keep_memory) {
                store_.resize(0);
            } else {
                store_.clear();
            }
            notify(nil, 0, dimension_);
        }

        
        virtual void redim(index_t dim) {
            if(dim == dimension()) {
                return;
            }
            vector<T> new_store(size()*dim);
            index_t copy_dim = GEO::geo_min(dim, dimension());
            for(index_t i = 0; i < size(); ++i) {
                for(index_t c = 0; c < copy_dim; ++c) {
                    new_store[dim * i + c] = store_[dimension_ * i + c];
                }
            }
            store_.swap(new_store);
            notify(
                store_.empty() ? nil : Memory::pointer(store_.data()),
                store_.size(),
                dim
            );
        }
        
        virtual bool elements_type_matches(const std::string& type_name) const {
            return type_name == typeid(T).name();
        }

    private:
        vector<T> store_;
    };


    /*********************************************************************/    

    class GEOGRAM_API AttributesManager {
    public:
        /**
         * \brief Constructs a new empty AttributesManager.
         */
        AttributesManager();


        /**
         * \brief AttributesManager destructor.
         */
        ~AttributesManager();

        /**
         * \brief Gets the number of attributes.
         * \return The number of attributes managed by this
         *   AttributeManager.
         */
        index_t nb() const {
            return index_t(attributes_.size());
        }

        /**
         * \brief Gets the names of all the attributes in this
         *   AttributeStore.
         * \param[out] names a vector of all attribute names
         */
        void list_attribute_names(vector<std::string>& names) const;
        
        /**
         * \brief Gets the size.
         * \details All attributes stored in an AttributeManager have
         *  the same number of items.
         * \return the number of items of each attribute.
         */
        index_t size() const {
            return size_;
        }
        
        /**
         * \brief Resizes all the attributes managed by this
         *  AttributeManager.
         * \param(in] new_size the new number of items for
         *  all attributes.
         */
        void resize(index_t new_size);

        /**
         * \brief Clears this AttributeManager
         * \param[in] keep_attributes if true, then all
         *  attributes are resized to 0 but their names are
         *  kept.
         * \param[in] keep_memory if true, allocated memory
         *  is kept reserved.
         */
        void clear(bool keep_attributes, bool keep_memory = false);
        

        /**
         * \brief Zeroes all the attributes.
         */
        void zero();
        
        /**
         * \brief Binds an AttributeStore with the specified name.
         *  Ownership of this AttributeStore is transfered to
         *  the AttributesManager.
         * \param[in] the name 
         * \param[in] as a pointer to the AttributeStore to be bound
         * \pre No AttributeStore is already bound to the same name
         */
        void bind_attribute_store(const std::string& name, AttributeStore* as);

        /**
         * \brief Finds an AttributeStore by name.
         * \param[in] name the name under which the AttributeStore was bound
         * \return a pointer to the attribute store or nil if is is undefined.
         */
        AttributeStore* find_attribute_store(const std::string& name);

        /**
         * \brief Tests whether an attribute is defined.
         * \param[in] name name of the attribute
         * \retval true if an attribute with the specified name exists
         * \retval false otherwise
         */
        bool is_defined(const std::string& name) {
            return (find_attribute_store(name) != nil);
        }
        
        /**
         * \brief Deletes an AttributeStore by name.
         * \param[in] name the name of the attribute store
         *  to be deleted.
         */
        void delete_attribute_store(const std::string& name);

        /**
         * \brief Deletes an AttributeStore.
         * \param[in] as a pointer to the attribute store
         *  to be deleted.
         */
        void delete_attribute_store(AttributeStore* as);

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
         * \param[in] permutation_in the permutation.
         *  It is temporarily changed during execution of the
         *  function, but identical to the input on exit.
         */
        void apply_permutation(
            const vector<index_t>& permutation
        );

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
        void compress(const vector<index_t>& old2new);

        
    private:
        /**
         * \brief Forbids copy.
         */
        AttributesManager(const AttributesManager* rhs);

        /**
         * \brief Forbids copy.
         */
        const AttributesManager& operator=(const AttributesManager& rhs);
        
    private:
        index_t size_;
        std::map<std::string, AttributeStore*> attributes_;
    } ;


    /*********************************************************************/    


    /**
     * \brief Base class for Attributes, that manipulates an 
     *  attribute stored in an AttributeManager.
     */
    template <class T> class AttributeBase : public AttributeStoreObserver {
    public:

        /**
         * \brief Creates an unitialized (unbound) Attribute.
         */
        AttributeBase() :
            manager_(nil),
            store_(nil) {
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
        AttributeBase(AttributesManager& manager, const std::string& name) :
            manager_(nil),
            store_(nil) {
            bind(manager, name);
        }

        /**
         * \brief Tests whether an Attribute is bound.
         * \retval true if this Attribute is bound
         * \retval false otherwise
         */
        bool is_bound() const {
            return (store_ != nil);
        }

        /**
         * \brief Unbinds this Attribute.
         * \pre is_bound()
         */
        void unbind() {
            geo_assert(is_bound());
            unregister_me(store_);
            manager_ = nil;
            store_ = nil;
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
        void bind(AttributesManager& manager, const std::string& name) {
            geo_assert(!is_bound());
            manager_ = &manager;
            store_ = manager_->find_attribute_store(name);
            if(store_ == nil) {
                store_ = new TypedAttributeStore<T>();
                manager_->bind_attribute_store(name,store_);
            } else {
                geo_assert(store_->elements_type_matches(typeid(T).name()));
            }
            register_me(store_);
        }


        /**
         * \brief Binds this Attribute to an AttributesManager if it
         *  already exists in the AttributeManager.
         * \param[in] manager a reference to the AttributesManager
         * \param[in] name name of the attribute
         * \pre !is_bound()
         */
        void bind_if_is_defined(
            AttributesManager& manager, const std::string& name
        ) {
            geo_assert(!is_bound());
            manager_ = &manager;
            store_ = manager_->find_attribute_store(name);
            if(store_ != nil) {
                geo_assert(store_->elements_type_matches(typeid(T).name()));
                register_me(store_);                
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
            index_t dimension
        ) {
            geo_assert(!is_bound());
            manager_ = &manager;            
            geo_assert(manager_->find_attribute_store(name) == nil);
            store_ = new TypedAttributeStore<T>(dimension);
            manager_->bind_attribute_store(name,store_);
            register_me(store_);
        }
        
        /**
         * \brief Destroys this attribute in the AttributeManager.
         * \details On exit, the attribute is no-longer accessible in
         *  the AttributeManager, its name is available again, and
         *  this attribute is in the unbound state.
         */
        void destroy() {
            geo_assert(is_bound());
            unregister_me(store_);
            manager_->delete_attribute_store(store_);
            store_ = nil;
            manager_ = nil;
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
        void redim(index_t new_dim) {
            geo_assert(is_bound());
            store_->redim(new_dim);
        }
        
        /**
         * \brief Attribute destructor
         * \details 
         *  The attribute is not destroyed, it can be retreived later 
         *  by binding with the same name. To destroy the attribute,
         *  use detroy() instead.
         */
        ~AttributeBase() {
            if(is_bound()) {
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
            AttributesManager& manager, const std::string& name,
            index_t dim = 0
        ) {
            AttributeStore* store = manager.find_attribute_store(name);
            return (
                store != nil &&
                store->elements_type_matches(typeid(T).name()) &&
                ((dim == 0) || (store->dimension() == dim))
            );
        }

        /**
         * \brief Gets the size.
         * \return The number of items in this attribute.
         */
        index_t size() const {
            return size_;
        }

        /**
         * \brief Sets all the elements of this Attribute to
         *   zero.
         */
        void zero() {
            geo_debug_assert(is_bound());
            store_->zero();
        }

    private:
        AttributesManager* manager_;
        AttributeStore* store_;
    } ;
    
    /*********************************************************************/

    template <class T> class Attribute : public AttributeBase<T> {
    public:
        typedef AttributeBase<T> superclass;

        /**
         * \brief Creates an unitialized (unbound) Attribute.
         */
        Attribute() : superclass() {
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
        Attribute(AttributesManager& manager, const std::string& name) :
            superclass(manager, name) {
        }

        /**
         * \brief Gets a modifiable element by index
         * \param [in] i index of the element
         * \return a modifiable reference to the \p i%th element
         */
        T& operator[](unsigned int i) {
            geo_debug_assert(i < superclass::nb_elements());
            return ((T*)superclass::base_addr_)[i];
        }

        /**
         * \brief Gets an element by index
         * \param [in] i index of the element
         * \return a const reference to the \p i%th element
         */
        const T& operator[](unsigned int i) const {
            geo_debug_assert(i < superclass::nb_elements());
            return ((const T*)superclass::base_addr_)[i];
        }

        /**
         * \brief Sets all the elements in this attribute
         *   to a specified value.
         * \param[in] val the value
         */
        void fill(const T& val) {
            for(index_t i=0; i<superclass::nb_elements(); ++i) {
                (*this)[i] = val;
            }
        }
    };
    
    /*********************************************************************/
    
    /**
     * \brief Specialization of Attribute for booleans
     * \details Attribute needs a specialization for bool, since
     *   vector<bool> uses compressed storage (1 bit per boolean),
     *   that is not compatible with the attribute management 
     *   mechanism. This wrapper class uses an Attribute<Numeric::uint8>
     *   and does the appropriate conversions, using an accessor class.
     */
    template <> class Attribute<bool> : public AttributeBase<Numeric::uint8> {
    public:
        typedef AttributeBase<Numeric::uint8> superclass;
        
        Attribute() : superclass() {
        }
        
        Attribute(AttributesManager& manager, const std::string& name) :
            superclass(manager,name) {
        }

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
                attribute_(attribute),
                index_(index) {
            }

            /**
             * \brief Converts a BoolAttributeAccessor to a bool.
             * \details Performs the actual lookup.
             */
            operator bool() const {
                return (attribute_.element(index_) != 0);
            }

            /**
             * \brief Assigns a bool to a BoolAttributeAccessor.
             * \details Stores the boolean into the Attribute.
             */
            BoolAttributeAccessor& operator=(bool x) {
                attribute_.element(index_) = Numeric::uint8(x);
                return *this;
            }
            
        private:
            Attribute<bool>& attribute_;
            index_t index_;
        };


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
                attribute_(attribute),
                index_(index) {
            }

            /**
             * \brief Converts a BoolAttributeAccessor to a bool.
             * \details Performs the actual lookup.
             */
            operator bool() const {
                return (attribute_.element(index_) != 0);
            }

        private:
            const Attribute<bool>& attribute_;
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
            for(index_t i=0; i<superclass::nb_elements(); ++i) {
                element(i) = Numeric::uint8(val);
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
        Numeric::uint8& element(unsigned int i) {
            geo_debug_assert(i < superclass::nb_elements());
            return ((Numeric::uint8*)superclass::base_addr_)[i];
        }

        /**
         * \brief Gets an element by index
         * \param [in] i index of the element
         * \return a const reference to the \p i%th element
         */
        const Numeric::uint8& element(unsigned int i) const {
            geo_debug_assert(i < superclass::nb_elements());
            return ((const Numeric::uint8*)superclass::base_addr_)[i];
        }
    } ;
 
    /*********************************************************************/
}

#endif

