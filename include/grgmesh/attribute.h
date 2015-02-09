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


#ifndef __GRGMESH_ATTRIBUTE__
#define __GRGMESH_ATTRIBUTE__

#include <grgmesh/common.h>

#include <geogram/basic/counted.h>
#include <geogram/basic/smart_pointer.h>

#include <vector>
#include <map>
#include <typeinfo>
#include <fstream>
#include <iostream>
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
            grgmesh_debug_assert( id < nb_elements_ ) ;
            return &data_[id * item_size_] ;
        }
        virtual const std::type_info& attribute_type_id() const = 0 ;
        index_t size() const { return nb_elements_ ; }
        index_t item_size() const { return item_size_ ; }
    protected:
        AttributeStore( index_t item_size, index_t size )
            : item_size_( item_size ), data_( nil ), nb_elements_( size )
        {
            data_ = new byte[item_size * size] ;
        }
        virtual ~AttributeStore() {
            delete[] data_ ;
        }
    protected:
        index_t item_size_ ;
        byte* data_ ;
        index_t nb_elements_ ;
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


    class AttributeManager {
    public:
        AttributeManager() {}
        virtual ~AttributeManager() {}

        void list_named_attributes( std::vector< std::string >& names )
        {
            names.clear() ;
            for( std::map< std::string, AttributeStore_var >::iterator it =
                attributes_.begin(); it != attributes_.end(); ++it ) {
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

        virtual int32 record_type_id() const = 0 ;
    
    private:
        std::map< std::string, AttributeStore_var > attributes_ ;
    } ;

    /**
     * \brief Generic manager of the attributes stored for one location
     *
     * Template by an integer instead of an enum specific to a class
     * Each object on which attributes will be created object should 
     * define possible locations with an enum. Jeanne.
     */
    template< int32 LOCATION >
    class AttributeManagerImpl : public AttributeManager {
    public :
        virtual int32 record_type_id() const { return LOCATION ; }
    } ;


   class AttributeBase {
   protected:
       static AttributeStore* resolve_named_attribute_store(
           AttributeManager* manager, const std::string& name
       ) {
           return manager->resolve_named_attribute_store(name) ;
       }
   
       static void bind_named_attribute_store(
           AttributeManager* manager, 
           const std::string& name, AttributeStore* as
       ) {
           manager->bind_named_attribute_store(name,as) ;
       }
   } ;

    /** 
     * \brief Generic attribute class - Storage on given elements of a given object 
     * The elements on which is defined the attribute are of the given ElementType.
     * Access to attribute value only using the index of the element in the object. Jeanne
     */ 
    template< int32 LOCATION, class ATTRIBUTE >
    class Attribute : public AttributeBase {
    public:
        typedef Attribute< LOCATION, ATTRIBUTE > thisclass ;
        typedef AttributeManagerImpl< LOCATION > Manager ;
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
                grgmesh_assert( store_ != nil ) ;
                // Sanity check, checks the attribute type.
                AttributeStore* check_type = store_ ;
                grgmesh_assert(
                    dynamic_cast< Store* >( check_type ) != nil ) ;
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

       /* // must be static because used in static function is_definedline 191
        static AttributeStore* resolve_named_attribute_store(
            Manager* manager,
            const std::string& name )
        {
            return manager->resolve_named_attribute_store( name ) ;
        }

        static void bind_named_attribute_store(
            Manager* manager,
            const std::string& name,
            AttributeStore* as )
        {
            manager->bind_named_attribute_store( name, as ) ;
        }
        */

    private:
        AttributeStore_var store_ ;
    } ;


    
    // What follows come directly from Graphite Copyright Bruno L�vy
    /**
     * AttributeSerializer is used to save and load attributes attached
     * to an object. This is the base class for serializing the value of 
     * an attribute and creating an AttributeStore given its type name.
     */
    class AttributeSerializer : public GEO::Counted {
    public:
        static void initialize() ;
        static void terminate() ;

        static AttributeSerializer* resolve_by_type(const std::type_info& attribute_type) ;
        static AttributeSerializer* resolve_by_name(const std::string& attribute_type_name) ;
        static std::string find_name_by_type(const std::type_info& attribute_type) ;
        
        static void bind(
            const std::type_info& attribute_type,
            const std::string& attribute_type_name,
            AttributeSerializer* serializer
        ) ;

        virtual AttributeStore* create_attribute_store(AttributeManager* manager) = 0 ;
        virtual bool serialize_read(std::istream& in,   byte* addr) = 0 ;
        virtual bool serialize_write(std::ostream& out, byte* addr) = 0 ;

    private:
        typedef std::map<std::string, GEO::SmartPointer<AttributeSerializer> > SerializerMap ;
        typedef std::map<std::string, std::string> StringMap ;

        static SerializerMap* type_to_serializer_ ;
        static SerializerMap* name_to_serializer_ ;
        static StringMap*     type_to_name_ ;
    } ;

    typedef GEO::SmartPointer<AttributeSerializer> AttributeSerializer_var ;

    //_________________________________________________________________________________________________________

    /**
     * Default implementation of AttributeSerializer.
     */
    template< class ATTRIBUTE > class GenericAttributeSerializer : public AttributeSerializer {
    public:
        typedef AttributeStoreImpl< ATTRIBUTE > AttributeStore;
        
        virtual AttributeStore* create_attribute_store(AttributeManager* manager) {
            return new AttributeStore(manager) ;
        }
        virtual bool serialize_read(std::istream& in, byte* addr) {
            ATTRIBUTE& attr = *reinterpret_cast<ATTRIBUTE*>(addr) ;
            in >> attr ;
            return true ;
        }
        virtual bool serialize_write(std::ostream& out, byte* addr) {
            ATTRIBUTE& attr = *reinterpret_cast<ATTRIBUTE*>(addr) ;
            out << attr ;
            return true ;
        }

    } ;

    /**
     * Use this class to declare a new serializable attribute type.
     * In the common.cpp file of the library, add:
     * ogf_register_attribute_type<MyAttributeType>("MyAttributeType") ;
     */
    template <class T> class grgmesh_register_attribute_type {
    public:
        grgmesh_register_attribute_type(const std::string& type_name) {
            AttributeSerializer::bind(typeid(T), type_name, new GenericAttributeSerializer<T>()) ;
        }
    } ;

    /**
     * SerializedAttributeRef is what SerializedAttribute::operator[] returns.
     * It is just meant to overload operator<< and operator>>.
     */
    class SerializedAttributeRef {
    public:
        SerializedAttributeRef(
            AttributeSerializer* ser, byte* addr
        ) : serializer_(ser), addr_(addr) {
        }
        AttributeSerializer* serializer() const { return serializer_ ; }
        byte* addr() const { return addr_ ; }
    private:
        AttributeSerializer* serializer_ ;
        byte* addr_ ;
    } ;
    
    inline std::istream& operator>>(std::istream& in, const SerializedAttributeRef& r) {
        r.serializer()->serialize_read(in, r.addr()) ;
        return in ;
    }

    inline std::ostream& operator<<(std::ostream& out, const SerializedAttributeRef& r) {
        r.serializer()->serialize_write(out, r.addr()) ;
        return out ;
    }

    /**
     * SerializedAttribute allows writing attribute values into a stream,
     * reading attribute values from a stream, and creating an attribute
     * from its type name. 
     */
    template< int32 LOCATION > class SerializedAttribute : public AttributeBase {
    public:
        typedef AttributeManagerImpl< LOCATION > AttributeManagerT;

        void bind(AttributeManagerT* manager, const std::string& name) {
            attribute_manager_ = manager ;
            attribute_store_ = resolve_named_attribute_store(manager, name) ;
            if(attribute_store_ != nil) {
                serializer_ = AttributeSerializer::resolve_by_type(attribute_store_->attribute_type_id()) ;
            }
            name_ = name ;
        }

        void bind(AttributeManagerT* manager, const std::string& name, const std::string& attribute_type_name) {
            attribute_manager_ = manager ;
            serializer_ = AttributeSerializer::resolve_by_name(attribute_type_name) ;
            if(serializer_ != nil) {
                if(attribute_manager_->named_attribute_is_bound(name)) {
                    attribute_store_ = resolve_named_attribute_store(attribute_manager_,name) ;
                    grgmesh_assert(
                        AttributeSerializer::find_name_by_type(
                            attribute_store_->attribute_type_id()
                        ) == attribute_type_name
                    ) ;
                } else {
                    attribute_store_ = serializer_->create_attribute_store(attribute_manager_) ;
                    bind_named_attribute_store(attribute_manager_,name,attribute_store_) ;
                }
            }
            name_ = name ;
        }

        void unbind() {
            attribute_manager_ = nil ;
            attribute_store_ = nil ;
            serializer_ = nil ;
        }

        SerializedAttribute() {
            attribute_manager_ = nil ;
            attribute_store_ = nil ;
            serializer_ = nil ;
        }

        SerializedAttribute(AttributeManagerT* manager, const std::string& name) {
            bind(manager, name) ;
        }

        SerializedAttribute(
            AttributeManagerT* manager, const std::string& name, const std::string& attribute_type_name
        ) {
            bind(manager, name, attribute_type_name) ;
        }

        SerializedAttribute(const SerializedAttribute& rhs) {
            attribute_manager_ = rhs.attribute_manager_ ;
            attribute_store_   = rhs.attribute_store_ ;
            serializer_        = rhs.serializer_ ;
            name_              = rhs.name_ ;
        }

        bool is_bound() const {
            return (attribute_manager_ != nil) && (attribute_store_ != nil) && (serializer_ != nil) ;
        }

        const std::string& name() const { return name_ ; }

        std::string type_name() const {
            grgmesh_assert(attribute_store_ != nil) ;
            return AttributeSerializer::find_name_by_type(attribute_store_->attribute_type_id()) ;
        }

        SerializedAttributeRef operator[]( const index_t record ) {
            return SerializedAttributeRef(serializer_, attribute_store_->data(record)) ;
        }

    private:
        AttributeManagerT* attribute_manager_ ;
        AttributeStore* attribute_store_ ;
        AttributeSerializer* serializer_ ;
        std::string name_ ;
    } ;

       
    template< int32 T > inline bool get_serializable_attributes(
        AttributeManagerImpl<T>* manager, std::vector<SerializedAttribute<T> >& attributes,
        std::ostream& out
    ) {
        bool result = false ;
        std::vector<std::string> names ;
        manager->list_named_attributes(names) ;
        for(unsigned int i=0; i<names.size(); i++) {
            attributes.push_back(SerializedAttribute<T>()) ;
            attributes.rbegin()->bind(manager, names[i]) ;
            if(attributes.rbegin()->is_bound()) {
                std::cerr << "Attribute " << names[i] // << " on " << localisation << " : " 
                          << attributes.rbegin()->type_name() << std::endl ;
                out << /*attribute_kw << " " <<*/ names[i] << " " /*<< localisation << " " */
                    << attributes.rbegin()->type_name() << std::endl ;
                result = true ;
            } else {
                std::cerr << "Attribute " << names[i] /*<< " on " << localisation */
                          << " is not serializable" << std::endl ;
                attributes.pop_back() ;
            }
        }
        return result ;
    }
    
    template< int32 T > inline void serialize_read_attributes(
        std::istream& in, int32 item, std::vector< SerializedAttribute<T> >& attributes
    ) {
        for(unsigned int i=0; i<attributes.size(); i++) {
            in >> attributes[i][item] ;
        }
    }

    template< int32 T > inline void serialize_write_attributes(
        std::ostream& out, int32 item, std::vector< SerializedAttribute<T> >& attributes
    ) {
        for(unsigned int i=0; i<attributes.size(); i++) {
            out << attributes[i][item] << " " ;
        }
    }

}

#endif
