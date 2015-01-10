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
*/

#ifndef __GRGMESH_ATTRIBUTE__
#define __GRGMESH_ATTRIBUTE__

#include <grgmesh/common.h>

#include <geogram/basic/counted.h>
#include <geogram/basic/smart_pointer.h>

#include <vector>
#include <map>
#include <typeinfo>

namespace GRGMesh {

     /** 
     * \brief Abstract base class to store an attribute as a vector of bytes.
     * 
     * Vector size depends on the size of the items stored.
     * Derives of Counted so that SmartPointers of this class may be used
     */
    class GRGMESH_API AttributeStore: public GEO::Counted {
    public:
        byte* data( index_t id)
        {
            return &data_[id * item_size_] ;
        }
        virtual const std::type_info& attribute_type_id() const = 0 ;
        index_t size() const { return data_.size() / item_size_ ; }
        index_t item_size() const { return item_size_ ; }
    protected:
        AttributeStore( index_t item_size, index_t size )
            : item_size_( item_size ), data_( item_size * size )
        {
        }
    protected:
        index_t item_size_ ;
        std::vector< byte > data_ ;
    } ;

    typedef GEO::SmartPointer< AttributeStore > AttributeStore_var ;

    template< class ATTRIBUTE >
    class AttributeStoreImpl: public AttributeStore {
    public:
        AttributeStoreImpl( index_t size )
            : AttributeStore( sizeof(ATTRIBUTE), size )
        {
        }
        virtual const std::type_info& attribute_type_id() const {
            return typeid( ATTRIBUTE ) ;
        }
    } ;


    /**
     * \brief Generic manager of the attributes stored for one location
     *
     * Template by an integer instead of an enum specific to a class
     * Each object on which attributes will be created object should 
     * define possible locations with an enum. Jeanne.
     */
    template< int32 LOCATION >
    class AttributeManager {
    public:
        enum Mode {
            FIND, CREATE, FIND_OR_CREATE
        } ;

    public :
        AttributeManager() {}
        virtual ~AttributeManager() {}

        void list_named_attributes( std::vector< std::string >& names )
        {
            names.clear() ;
            for( std::map< std::string, AttributeStore_var >::iterator it =
                attributes_.begin(); it != attributes_.end(); it++ ) {
                names.push_back( it->first ) ;
            }
        }
        bool named_attribute_is_bound( const std::string& name )
        {
            return ( attributes_.find( name ) != attributes_.end() ) ;
        }
        void delete_named_attribute( const std::string& name )
        {
            std::map< std::string, AttributeStore_var >::iterator it =
                attributes_.find( name ) ;
            grgmesh_debug_assert( it != attributes_.end() ) ;
            grgmesh_debug_assert( !it->second->is_shared() ) ;
            attributes_.erase( it ) ;
        }

        const ElementType& record_type_id() const  { return LOCATION ; }

        void bind_named_attribute_store(
            const std::string& name,
            AttributeStore* as )
        {
            grgmesh_debug_assert( !named_attribute_is_bound( name ) ) ;
            attributes_[name] = as ;
        }

        AttributeStore* resolve_named_attribute_store( const std::string& name )
        {
            std::map< std::string, AttributeStore_var >::iterator it =
                attributes_.find( name ) ;
            grgmesh_debug_assert( it != attributes_.end() ) ;
            return it->second ;
        }
    private:
        std::map< std::string, AttributeStore_var > attributes_ ;
    } ;


    /** 
     * \brief Generic attribute class - Storage on given elements of a given object 
     * The elements on which is defined the attribute are of the given ElementType.
     * Access to attribute value only using the index of the element in the object. Jeanne
     */ 
    template< int32 LOCATION, class ATTRIBUTE >
    class Attribute {
    public:
        typedef Attribute< LOCATION, ATTRIBUTE > thisclass ;
        typedef AttributeManager< LOCATION > Manager ;
        typedef AttributeStoreImpl< ATTRIBUTE > Store ;

        Attribute()
        {
        }
        Attribute( Manager* manager, index_t size ) {
            bind( manager, size ) ;
        }
        Attribute( Manager* manager, index_t size, const std::string& name )
        {
            bind( manager, size, name ) ;
        }
        Attribute( const thisclass& rhs )
            : store_( rhs.store_ )
        {
        }
        thisclass& operator=( const thisclass& rhs )
        {
            store_ = rhs.store_ ;
            return *this ;
        }

        index_t size() const { return store_->size() ; }
        bool is_bound() const
        {
            return !store_.is_nil() ;
        }

        void bind( Manager* manager, index_t size )
        {
            store_ = new Store( size ) ;
        }

        void bind( Manager* manager, index_t size, const std::string& name )
        {
            if( manager->named_attribute_is_bound( name ) ) {
                store_ = resolve_named_attribute_store( manager, name ) ;
                grgmesh_debug_assert( store_ != nil ) ;
                // Sanity check, checks the attribute type.
                AttributeStore* check_type = store_ ;
                grgmesh_debug_assert(
                    dynamic_cast< AttributeStore* >( check_type ) != nil ) ;
            } else {
                store_ = new Store( size ) ;
                bind_named_attribute_store( manager, name, store_ ) ;

            }
        }

        void unbind() {
            store_.reset() ;
        }

        ATTRIBUTE& operator[]( const index_t& id )
        {
            return *data( id ) ;
        }

        const ATTRIBUTE& operator[]( const index_t& id ) const
        {
            return *data( id ) ;
        }

        /**
         * Checks whether manager has an attribute of this type
         * bound with the specified name.
         */
        static bool is_defined( Manager* manager, const std::string& name )
        {
            return ( manager->named_attribute_is_bound( name )
                && dynamic_cast< Store* >( resolve_named_attribute_store(
                    manager, name ) ) != nil ) ;
        }

    protected:
        ATTRIBUTE* data( const index_t& id ) const
        {
            return reinterpret_cast< ATTRIBUTE* >( store_->data( id ) ) ;
        }

        // must be static because used in static function is_definedline 191
        static AttributeStore* resolve_named_attribute_store(
            Manager* manager,
            const std::string& name )
        {
            return manager->resolve_named_attribute_store( name ) ;
        }

        void bind_named_attribute_store(
            Manager* manager,
            const std::string& name,
            AttributeStore* as )
        {
            manager->bind_named_attribute_store( name, as ) ;
        }

    private:
        AttributeStore_var store_ ;
    } ;



}

#endif
