/*
* This file has been strongly inspired from the attributes of the geogram
* library.
* Many thanks to Bruno Levy (Bruno.Levy@inria.fr) who did the first
* implementation in Geogram.
*
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

#include <ringmesh/geogram_extension/geogram_extension.h>

#include <geogram/basic/permutation.h>
#include <ringmesh/basic/attributes.h>

#include <algorithm>
#include <memory>

namespace RINGMesh
{
    /******************************************************************/

    std::map< std::string, std::unique_ptr< AttributeStoreCreator > >
        AttributeStore::type_name_to_creator_;

    std::map< std::string, std::string >
        AttributeStore::typeid_name_to_type_name_;

    std::map< std::string, std::string >
        AttributeStore::type_name_to_typeid_name_;

    AttributeStore::~AttributeStore()
    {
        // It is illegal to keep an Attribute<> active
        // when the object it is bound to is destroyed.
    }

    /*************************************************************************/

    AttributesManager::~AttributesManager()
    {
        clear( false );
    }

    void AttributesManager::resize( index_t new_size )
    {
        if( new_size == nb_items_ )
        {
            return;
        }
        for( std::map< std::string, AttributeStore* >::iterator it =
                 attributes_.begin();
             it != attributes_.end(); ++it )
        {
            it->second->resize( new_size );
        }
        nb_items_ = new_size;
    }

    void AttributesManager::apply_permutation(
        const std::vector< index_t >& permutation )
    {
        for( std::map< std::string, AttributeStore* >::iterator it =
                 attributes_.begin();
             it != attributes_.end(); ++it )
        {
            it->second->apply_permutation( permutation );
        }
    }

    void AttributesManager::compress( const std::vector< index_t >& old2new )
    {
        for( std::map< std::string, AttributeStore* >::iterator it =
                 attributes_.begin();
             it != attributes_.end(); ++it )
        {
            it->second->compress( old2new );
        }
    }

    void AttributesManager::bind_attribute_store(
        const std::string& name, AttributeStore* as )
    {
        ringmesh_assert( find_attribute_store( name ) == nullptr );
        attributes_[name] = as;
    }

    std::vector< std::string > AttributesManager::attribute_names() const
    {
        std::vector< std::string > names;
        names.reserve( attributes_.size() );
        for( auto it : attributes_ )
        {
            names.push_back( it.first );
        }
        return names;
    }

    AttributeStore* AttributesManager::find_attribute_store(
        const std::string& name )
    {
        std::map< std::string, AttributeStore* >::iterator it =
            attributes_.find( name );
        if( it == attributes_.end() )
        {
            return nil;
        }
        return it->second;
    }

    const AttributeStore* AttributesManager::find_attribute_store(
        const std::string& name ) const
    {
        std::map< std::string, AttributeStore* >::const_iterator it =
            attributes_.find( name );
        if( it == attributes_.end() )
        {
            return nil;
        }
        return it->second;
    }

    void AttributesManager::delete_attribute_store( const std::string& name )
    {
        std::map< std::string, AttributeStore* >::iterator it =
            attributes_.find( name );
        ringmesh_assert( it != attributes_.end() );
        delete it->second;
        attributes_.erase( it );
    }

    void AttributesManager::delete_attribute_store( AttributeStore* as )
    {
        for( std::map< std::string, AttributeStore* >::iterator it =
                 attributes_.begin();
             it != attributes_.end(); ++it )
        {
            if( it->second == as )
            {
                delete as;
                attributes_.erase( it );
                return;
            }
        }
        ringmesh_assert_not_reached;
    }

    void AttributesManager::clear( bool keep_attributes )
    {
        if( keep_attributes )
        {
            for( std::map< std::string, AttributeStore* >::iterator it =
                     attributes_.begin();
                 it != attributes_.end(); ++it )
            {
                it->second->clear();
            }
        }
        else
        {
            for( std::map< std::string, AttributeStore* >::iterator it =
                     attributes_.begin();
                 it != attributes_.end(); ++it )
            {
                delete it->second;
            }
            attributes_.clear();
        }
        nb_items_ = 0;
    }

    void AttributesManager::copy( const AttributesManager& rhs )
    {
        clear( false );
        resize( rhs.nb_items() );
        for( std::map< std::string, AttributeStore* >::const_iterator it =
                 rhs.attributes_.begin();
             it != rhs.attributes_.end(); ++it )
        {
            bind_attribute_store( it->first, it->second->clone() );
        }
    }

    /*   void AttributesManager::copy_item( index_t to, index_t from )
       {
           for(
               std::map<std::string, AttributeStore*>::iterator
               it = attributes_.begin();
           it != attributes_.end(); ++it
               ) {
               it->second->copy_item( to, from );
           }
       }*/

    /************************************************************************/

    std::string ReadOnlyScalarAttributeAdapter::attribute_base_name(
        const std::string& name )
    {
        size_t pos = name.find( '[' );
        if( pos == std::string::npos )
        {
            return name;
        }
        return name.substr( 0, pos );
    }

    index_t ReadOnlyScalarAttributeAdapter::attribute_element_index(
        const std::string& name )
    {
        index_t result = 0;
        size_t pos = name.find( '[' );
        if( pos != std::string::npos )
        {
            try
            {
                if( pos + 2 > name.length() )
                {
                    result = index_t( -1 );
                }
                else
                {
                    result = GEO::String::to_uint(
                        name.substr( pos + 1, name.length() - pos - 2 ) );
                }
            }
            catch( ... )
            {
                result = index_t( -1 );
            }
        }
        return result;
    }

    void ReadOnlyScalarAttributeAdapter::bind_if_is_defined(
        const AttributesManager& manager, const std::string& name )
    {
        ringmesh_assert( !is_bound() );
        manager_ = &manager;
        element_index_ = attribute_element_index( name );
        store_ = manager_->find_attribute_store( attribute_base_name( name ) );

        if( store_ == nil || element_index_ == index_t( -1 ) )
        {
            store_ = nil;
            element_index_ = index_t( -1 );
            return;
        }

        // Test element_index_ validity: should be smaller than
        // store's dimension (or 2*store dimension if a vec2,
        // or 3*store's dimension if a vec3)
        if( element_index_ >= nb_scalar_elements_per_item() )
        {
            store_ = nil;
            element_index_ = index_t( -1 );
            return;
        }
    }

    bool ReadOnlyScalarAttributeAdapter::is_defined(
        const AttributesManager& manager, const std::string& name )
    {
        std::string attribute_name = attribute_base_name( name );
        const AttributeStore* att_store =
            manager.find_attribute_store( attribute_name );

        if( att_store == nil )
        {
            return false;
        }
        std::unique_ptr< ReadOnlyScalarAttributeAdapter > adapter =
            ReadOnlyScalarAttributeAdapterFactory::create(
                att_store->get_store().element_typeid_name(), manager,
                attribute_name );

        index_t element_index = attribute_element_index( name );
        if( element_index == index_t( -1 ) )
        {
            return false;
        }

        if( element_index >= adapter->nb_scalar_elements_per_item() )
        {
            return false;
        }

        return true;
    }

    /************************************************************************/

    class ReadOnlyScalarAttributeAdapterET_UINT8
        : public ReadOnlyScalarAttributeAdapter
    {
    public:
        ReadOnlyScalarAttributeAdapterET_UINT8(
            const AttributesManager& manager, const std::string& name )
            : ReadOnlyScalarAttributeAdapter( manager, name )
        {
        }

        double operator[]( index_t i ) override
        {
            return get_element< uint8 >( i );
        }

        index_t nb_scalar_elements_per_item() const override
        {
            return 1;
        }
        /* AttributeElementType element_type() const override
         {
             return typeid( uint8 ).name();
         }*/
        bool is_integer_like_attribute() const override
        {
            return true;
        }
    };

    class ReadOnlyScalarAttributeAdapterET_INT8
        : public ReadOnlyScalarAttributeAdapter
    {
    public:
        ReadOnlyScalarAttributeAdapterET_INT8(
            const AttributesManager& manager, const std::string& name )
            : ReadOnlyScalarAttributeAdapter( manager, name )
        {
        }

        double operator[]( index_t i ) override
        {
            return get_element< int8 >( i );
        }

        index_t nb_scalar_elements_per_item() const override
        {
            return 1;
        }
        /*AttributeElementType element_type() const override
        {
            return typeid( int8 ).name();
        }*/
        bool is_integer_like_attribute() const override
        {
            return true;
        }
    };

    class ReadOnlyScalarAttributeAdapterET_UINT32
        : public ReadOnlyScalarAttributeAdapter
    {
    public:
        ReadOnlyScalarAttributeAdapterET_UINT32(
            const AttributesManager& manager, const std::string& name )
            : ReadOnlyScalarAttributeAdapter( manager, name )
        {
        }

        double operator[]( index_t i ) override
        {
            return get_element< uint32 >( i );
        }

        index_t nb_scalar_elements_per_item() const override
        {
            return 1;
        }
        /*AttributeElementType element_type() const override
        {
            return typeid( unsigned int ).name();
        }*/
        bool is_integer_like_attribute() const override
        {
            return true;
        }
    };

    class ReadOnlyScalarAttributeAdapterET_INT32
        : public ReadOnlyScalarAttributeAdapter
    {
    public:
        ReadOnlyScalarAttributeAdapterET_INT32(
            const AttributesManager& manager, const std::string& name )
            : ReadOnlyScalarAttributeAdapter( manager, name )
        {
        }

        double operator[]( index_t i ) override
        {
            return get_element< int32 >( i );
        }

        index_t nb_scalar_elements_per_item() const override
        {
            return 1;
        }
        /*AttributeElementType element_type() const override
        {
            return typeid( int ).name();
        }*/
        bool is_integer_like_attribute() const override
        {
            return true;
        }
    };

    class ReadOnlyScalarAttributeAdapterET_FLOAT32
        : public ReadOnlyScalarAttributeAdapter
    {
    public:
        ReadOnlyScalarAttributeAdapterET_FLOAT32(
            const AttributesManager& manager, const std::string& name )
            : ReadOnlyScalarAttributeAdapter( manager, name )
        {
        }

        double operator[]( index_t i ) override
        {
            return get_element< float32 >( i );
        }

        index_t nb_scalar_elements_per_item() const override
        {
            return 1;
        }
        /*AttributeElementType element_type() const override
        {
            return typeid( float ).name();
        }*/
        bool is_integer_like_attribute() const override
        {
            return false;
        }
    };

    class ReadOnlyScalarAttributeAdapterET_FLOAT64
        : public ReadOnlyScalarAttributeAdapter
    {
    public:
        ReadOnlyScalarAttributeAdapterET_FLOAT64(
            const AttributesManager& manager, const std::string& name )
            : ReadOnlyScalarAttributeAdapter( manager, name )
        {
        }

        double operator[]( index_t i ) override
        {
            return get_element< float64 >( i );
        }

        index_t nb_scalar_elements_per_item() const override
        {
            return 1;
        }
        /*virtual AttributeElementType element_type() const override
        {
            return  typeid( double ).name();
        }*/
        bool is_integer_like_attribute() const override
        {
            return false;
        }
    };

    class ReadOnlyScalarAttributeAdapterET_VEC2
        : public ReadOnlyScalarAttributeAdapter
    {
    public:
        ReadOnlyScalarAttributeAdapterET_VEC2(
            const AttributesManager& manager, const std::string& name )
            : ReadOnlyScalarAttributeAdapter( manager, name )
        {
        }

        double operator[]( index_t i ) override
        {
            return get_element< float64 >( i, 2 );
        }

        index_t nb_scalar_elements_per_item() const override
        {
            return 2;
        }
        /*AttributeElementType element_type() const override
        {
            return  typeid( vec2 ).name();
        }*/
        bool is_integer_like_attribute() const override
        {
            return false;
        }
    };

    class ReadOnlyScalarAttributeAdapterET_VEC3
        : public ReadOnlyScalarAttributeAdapter
    {
    public:
        ReadOnlyScalarAttributeAdapterET_VEC3(
            const AttributesManager& manager, const std::string& name )
            : ReadOnlyScalarAttributeAdapter( manager, name )
        {
        }

        double operator[]( index_t i ) override
        {
            return get_element< float64 >( i, 3 );
        }

        index_t nb_scalar_elements_per_item() const override
        {
            return 3;
        }
        /*AttributeElementType element_type() const override
        {
            return  typeid( vec3 ).name();
        }*/
        bool is_integer_like_attribute() const override
        {
            return false;
        }
    };

    void register_read_only_scalar_attribute()
    {
        ReadOnlyScalarAttributeAdapterFactory::
            register_creator< ReadOnlyScalarAttributeAdapterET_UINT8 >(
                typeid( uint8 ).name() );
        ReadOnlyScalarAttributeAdapterFactory::
            register_creator< ReadOnlyScalarAttributeAdapterET_INT8 >(
                typeid( char ).name() );
        ReadOnlyScalarAttributeAdapterFactory::
            register_creator< ReadOnlyScalarAttributeAdapterET_INT8 >(
                typeid( int8 ).name() );
        ReadOnlyScalarAttributeAdapterFactory::
            register_creator< ReadOnlyScalarAttributeAdapterET_UINT32 >(
                typeid( unsigned int ).name() );
        ReadOnlyScalarAttributeAdapterFactory::
            register_creator< ReadOnlyScalarAttributeAdapterET_INT32 >(
                typeid( int ).name() );
        ReadOnlyScalarAttributeAdapterFactory::
            register_creator< ReadOnlyScalarAttributeAdapterET_FLOAT32 >(
                typeid( float ).name() );
        ReadOnlyScalarAttributeAdapterFactory::
            register_creator< ReadOnlyScalarAttributeAdapterET_FLOAT64 >(
                typeid( double ).name() );
        ReadOnlyScalarAttributeAdapterFactory::
            register_creator< ReadOnlyScalarAttributeAdapterET_VEC2 >(
                typeid( vec2 ).name() );
        ReadOnlyScalarAttributeAdapterFactory::
            register_creator< ReadOnlyScalarAttributeAdapterET_VEC3 >(
                typeid( vec3 ).name() );
    }
    /************************************************************************/
}