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
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_ATTRIBUTE__
#define __RINGMESH_ATTRIBUTE__

#include <ringmesh/common.h>

#include <geogram/basic/counted.h>
#include <geogram/basic/smart_pointer.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/attributes.h>

#include <vector>
#include <map>
#include <typeinfo>
#include <fstream>
#include <iostream>

namespace RINGMesh {


//    // What follows come directly from Graphite Copyright Bruno Levy
//    /**
//     * AttributeSerializer is used to save and load attributes attached
//     * to an object. This is the base class for serializing the value of
//     * an attribute and creating an AttributeStore given its type name.
//     */
//    class AttributeSerializer : public GEO::Counted {
//    public:
//        static void initialize() ;
//
//        static void terminate() ;
//
//        static AttributeSerializer* resolve_by_type(
//            const std::type_info& attribute_type ) ;
//
//        static AttributeSerializer* resolve_by_name(
//            const std::string& attribute_type_name ) ;
//
//        static std::string find_name_by_type( const std::type_info& attribute_type ) ;
//
//        static void bind(
//            const std::type_info& attribute_type,
//            const std::string& attribute_type_name,
//            AttributeSerializer* serializer ) ;
//
//        virtual AttributeStore* create_attribute_store( index_t record_size ) = 0 ;
//
//        virtual bool serialize_read(
//            std::istream& in,
//            byte* addr ) = 0 ;
//
//        virtual bool serialize_write(
//            std::ostream& out,
//            byte* addr ) = 0 ;
//
//    private:
//        typedef std::map< std::string,
//                          GEO::SmartPointer< AttributeSerializer > > SerializerMap ;
//        typedef std::map< std::string, std::string > StringMap ;
//
//        static SerializerMap* type_to_serializer_ ;
//        static SerializerMap* name_to_serializer_ ;
//        static StringMap*     type_to_name_ ;
//    } ;
//
//    typedef GEO::SmartPointer< AttributeSerializer > AttributeSerializer_var ;
//
//    // _________________________________________________________________________________________________________
//
//    /**
//     * Default implementation of AttributeSerializer.
//     */
//    template< class ATTRIBUTE >
//    class GenericAttributeSerializer : public
//                                       AttributeSerializer {
//    public:
//        virtual AttributeStoreImpl< ATTRIBUTE >* create_attribute_store(
//            index_t record_size )
//        {
//            return new AttributeStoreImpl< ATTRIBUTE >( record_size ) ;
//        }
//
//        virtual bool serialize_read(
//            std::istream& in,
//            byte* addr )
//        {
//            ATTRIBUTE& attr = *reinterpret_cast< ATTRIBUTE* >( addr ) ;
//            in >> attr ;
//            return true ;
//        }
//
//        virtual bool serialize_write(
//            std::ostream& out,
//            byte* addr )
//        {
//            ATTRIBUTE& attr = *reinterpret_cast< ATTRIBUTE* >( addr ) ;
//            out << attr ;
//            return true ;
//        }
//    } ;
//
//    /**
//     * Use this class to declare a new serializable attribute type.
//     * In the common.cpp file of the library, add:
//     * ringmesh_register_attribute_type<MyAttributeType>("MyAttributeType") ;
//     */
//    template< class T >
//    class ringmesh_register_attribute_type {
//    public:
//        ringmesh_register_attribute_type( const std::string& type_name )
//        {
//            AttributeSerializer::bind( typeid( T ), type_name,
//                new GenericAttributeSerializer< T >() ) ;
//        }
//    } ;
//
//    /**
//     * SerializedAttributeRef is what SerializedAttribute::operator[] returns.
//     * It is just meant to overload operator<< and operator>>.
//     */
//    class SerializedAttributeRef {
//    public:
//        SerializedAttributeRef(
//            AttributeSerializer* ser,
//            byte* addr ) : serializer_( ser ), addr_( addr )
//        {
//        }
//
//        AttributeSerializer* serializer() const { return serializer_ ;}
//        byte* addr() const { return addr_ ;}
//
//    private:
//        AttributeSerializer* serializer_ ;
//        byte* addr_ ;
//    } ;
//
//    inline std::istream& operator>>(
//        std::istream& in,
//        const SerializedAttributeRef& r )
//    {
//        r.serializer()->serialize_read( in, r.addr() ) ;
//        return in ;
//    }
//
//
//    inline std::ostream& operator<<(
//        std::ostream& out,
//        const SerializedAttributeRef& r )
//    {
//        r.serializer()->serialize_write( out, r.addr() ) ;
//        return out ;
//    }
//
//
//    /**
//     * SerializedAttribute allows writing attribute values into a stream,
//     * reading attribute values from a stream, and creating an attribute
//     * from its type name.
//     */
//    class SerializedAttribute : public AttributeBase {
//    public:
//        typedef GEO::AttributesManager AttributeManagerT ;
//
//        void bind(
//            AttributeManagerT* manager,
//            const std::string& name )
//        {
//            attribute_manager_ = manager ;
//            attribute_store_ = resolve_named_attribute_store( manager, name ) ;
//            if( attribute_store_ != nil ) {
//                serializer_ = AttributeSerializer::resolve_by_type(
//                    attribute_store_->attribute_type_id() ) ;
//            }
//            name_ = name ;
//        }
//
//        void bind(
//            AttributeManagerT* manager,
//            const std::string& name,
//            const std::string& attribute_type_name,
//            index_t record_size )
//        {
//            attribute_manager_ = manager ;
//            serializer_ = AttributeSerializer::resolve_by_name( attribute_type_name ) ;
//            if( serializer_ != nil ) {
//                if( attribute_manager_->named_attribute_is_bound( name ) ) {
//                    attribute_store_ = resolve_named_attribute_store(
//                        attribute_manager_, name ) ;
//                    ringmesh_assert(
//                        AttributeSerializer::find_name_by_type(
//                            attribute_store_->attribute_type_id()
//                            ) == attribute_type_name
//                        ) ;
//                } else {
//                    attribute_store_ = serializer_->create_attribute_store(
//                        record_size ) ;
//                    bind_named_attribute_store( attribute_manager_, name,
//                        attribute_store_ ) ;
//                }
//            }
//            name_ = name ;
//        }
//
//        void unbind()
//        {
//            attribute_manager_ = nil ;
//            attribute_store_ = nil ;
//            serializer_ = nil ;
//        }
//
//        SerializedAttribute()
//        {
//            attribute_manager_ = nil ;
//            attribute_store_ = nil ;
//            serializer_ = nil ;
//        }
//
//        SerializedAttribute(
//            AttributeManagerT* manager,
//            const std::string& name )
//        {
//            bind( manager, name ) ;
//        }
//
//        SerializedAttribute(
//            AttributeManagerT* manager,
//            const std::string& name,
//            const std::string& attribute_type_name )
//        {
//            bind( manager, name, attribute_type_name ) ;
//        }
//
//        SerializedAttribute( const SerializedAttribute& rhs )
//        {
//            attribute_manager_ = rhs.attribute_manager_ ;
//            attribute_store_   = rhs.attribute_store_ ;
//            serializer_        = rhs.serializer_ ;
//            name_              = rhs.name_ ;
//        }
//
//        bool is_bound() const
//        {
//            return ( attribute_manager_ != nil ) && ( attribute_store_ != nil ) &&
//                   ( serializer_ != nil ) ;
//        }
//
//        const std::string& name() const { return name_ ;}
//
//        std::string type_name() const
//        {
//            ringmesh_assert( attribute_store_ != nil ) ;
//            return AttributeSerializer::find_name_by_type(
//                attribute_store_->attribute_type_id() ) ;
//        }
//
//        SerializedAttributeRef operator[]( const index_t record )
//        {
//            return SerializedAttributeRef( serializer_,
//                attribute_store_->data( record ) ) ;
//        }
//
//    private:
//        AttributeManagerT* attribute_manager_ ;
//        AttributeStore* attribute_store_ ;
//        AttributeSerializer* serializer_ ;
//        std::string name_ ;
//    } ;
//
//    /*!
//     * \brief Get from the manager the attributes that can be saved in a file
//     * Those for which there is a type name has beed registered cf. ringmesh_register_attribute_type
//     * The names and types of the writable attributes are written in the output stream
//     * The corresponding serialized attributes are added to attributes
//     *
//     *  Remove this function and keep the next one ??
//     */
//    template< int32 T >
//    inline void get_serializable_attributes(
//        AttributeManagerImpl< T >* manager,
//        std::vector< SerializedAttribute< T > >& attributes,
//        std::ostream& out )
//    {
//        std::vector< std::string > names ;
//        manager->list_named_attributes( names ) ;
//
//        for( unsigned int i = 0; i < names.size(); i++ ) {
//            attributes.push_back( SerializedAttribute< T >() ) ;
//            attributes.back().bind( manager, names[ i ] ) ;
//            if( attributes.back().is_bound() ) {
//                out << names[ i ] << " "
//                    << attributes.back().type_name() ;
//            } else {
//                std::cerr << "Attribute " << names[ i ] << " is not serializable" <<
//                std::endl ;
//                attributes.pop_back() ;
//            }
//        }
//        out << std::endl ;
//    }
//
//
//    template< int32 T >
//    inline void get_serializable_attributes(
//        AttributeManagerImpl< T >* manager,
//        std::vector< SerializedAttribute< T > >& attributes )
//    {
//        std::vector< std::string > names ;
//        manager->list_named_attributes( names ) ;
//
//        for( unsigned int i = 0; i < names.size(); i++ ) {
//            attributes.push_back( SerializedAttribute< T >() ) ;
//            attributes.back().bind( manager, names[ i ] ) ;
//
//            // If the attribute is no bound remove it
//            if( !attributes.back().is_bound() ) {
//                attributes.pop_back() ;
//            }
//        }
//    }
//
//
//    /*!
//     * Generic reading of a attribute values in LineStream
//     *
//     * \param[in] in Input GEO::LineStream
//     * \param[in] start_field From which index should the attributes be read in in
//     * \param[in] item Index of the item to which the read values should be linked to
//     * \param[in] attributes The attributes to read
//     */
//    template< int32 T >
//    inline void serialize_read_attributes(
//        const GEO::LineInput& in,
//        index_t start_field,
//        int32 item,
//        std::vector< SerializedAttribute< T > >& attributes )
//    {
//        ringmesh_assert( start_field + attributes.size() - 1 < in.nb_fields()  ) ;
//        for( unsigned int i = 0; i < attributes.size(); i++ ) {
//            std::istringstream is( in.field( start_field + i ) ) ;
//            is >> attributes[ i ][ item ] ;
//        }
//    }
//
//
//    /*!
//     * Generic writing of attributes in a ostream
//     */
//    template< int32 T >
//    inline void serialize_write_attributes(
//        std::ostream& out,
//        int32 item,
//        std::vector< SerializedAttribute< T > >& attributes )
//    {
//        for( unsigned int i = 0; i < attributes.size(); i++ ) {
//            out << attributes[ i ][ item ] << " " ;
//        }
//    }
}

#endif
