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

#ifndef __GRGMESH_ATTRIBUTE__
#define __GRGMESH_ATTRIBUTE__

#include <grgmeshlib/common.h>
#include <grgmeshlib/smart_pointer.h>

#include <vector>
#include <map>
#include <typeinfo>

namespace GRGMesh {


    class AttributeStore: public Counted {
    public:
        pointer data( uint32 id) const
        {
            return data_[id * item_size_] ;
        }
        virtual const std::type_info& attribute_type_id() const = 0 ;
        uint32 size() const { return data_.size() / item_size_ ; }
    protected:
        AttributeStore( uint32 item_size, uint32 size )
            :
                item_size_( item_size ), data_( item_size * size )
        {
        }
    protected:
        uint32 item_size_ ;
        std::vector< pointer > data_ ;
    } ;

    typedef SmartPointer< AttributeStore > AttributeStore_var ;

    template< class ATTRIBUTE >
    class AttributeStoreImpl: public AttributeStore {
    public:
        AttributeStoreImpl( uint32 size )
            : AttributeStore( sizeof(ATTRIBUTE), size )
        {
        }
        virtual const std::type_info& attribute_type_id() const {
            return typeid( ATTRIBUTE ) ;
        }
    } ;


    template< ElementType LOCATION >
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





    template< ElementType LOCATION, class ATTRIBUTE >
    class Attribute {
    public:
        typedef Attribute< LOCATION, ATTRIBUTE > thisclass ;
        typedef AttributeManager< LOCATION > Manager ;
        typedef AttributeStoreImpl< ATTRIBUTE > Store ;

        Attribute()
        {
        }
        Attribute( Manager* manager, uint32 size ) {
            bind( manager, size ) ;
        }
        Attribute( Manager* manager, uint32 size, const std::string& name )
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

        uint32 size() const { return store_->size() ; }
        bool is_bound() const
        {
            return !store_.is_nil() ;
        }

        void bind( Manager* manager, uint32 size )
        {
            store_ = new Store( size ) ;
        }

        void bind( Manager* manager, uint32 size, const std::string& name )
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
            store_.forget() ;
        }

        ATTRIBUTE& operator[]( const uint32& id )
        {
            return *data( id ) ;
        }

        const ATTRIBUTE& operator[]( const uint32& id ) const
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
        ATTRIBUTE* data( const uint32& id ) const
        {
            return reinterpret_cast< ATTRIBUTE* >( store_->data( id ) ) ;
        }

        AttributeStore* resolve_named_attribute_store(
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
