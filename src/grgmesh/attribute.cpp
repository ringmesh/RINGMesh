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
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
*/

#include <grgmesh/attribute.h>
namespace GRGMesh {

    AttributeSerializer::SerializerMap* AttributeSerializer::type_to_serializer_ = nil ;
    AttributeSerializer::SerializerMap* AttributeSerializer::name_to_serializer_  = nil ;
    AttributeSerializer::StringMap*     AttributeSerializer::type_to_name_        = nil ;    

    void AttributeSerializer::initialize() {
        grgmesh_assert(type_to_serializer_ == nil) ;
        type_to_serializer_ = new SerializerMap ;
        grgmesh_assert(name_to_serializer_ == nil) ;
        name_to_serializer_ = new SerializerMap ;
        grgmesh_assert(type_to_name_ == nil) ;
        type_to_name_ = new StringMap ;
    }
    
    void AttributeSerializer::terminate() {
        delete type_to_serializer_ ;
        type_to_serializer_ = nil ;
        delete name_to_serializer_ ;
        name_to_serializer_ = nil ;
        delete type_to_name_ ;
        type_to_name_ = nil ;
    }
    
    AttributeSerializer* AttributeSerializer::resolve_by_type(const std::type_info& attribute_type) {
        grgmesh_assert(type_to_serializer_ != nil) ;
        SerializerMap::iterator it = type_to_serializer_->find(attribute_type.name()) ;
        if(it == type_to_serializer_->end()) {
            return nil ;
        }
        return it->second ;
    }

    AttributeSerializer* AttributeSerializer::resolve_by_name(const std::string& type_name) {
        grgmesh_assert(name_to_serializer_ != nil) ;       
        SerializerMap::iterator it = name_to_serializer_->find(type_name) ;
        if(it == name_to_serializer_->end()) {
            return nil ;
        }
        return it->second ;
    }
    

    std::string AttributeSerializer::find_name_by_type(const std::type_info& attribute_type) {
        grgmesh_assert(type_to_name_ != nil) ; 
        StringMap::iterator it = type_to_name_->find(attribute_type.name()) ;
        if(it == type_to_name_->end()) {
            return "unknown" ;
        }
        return it->second ;
    }

    void AttributeSerializer::bind(
        const std::type_info& attribute_type, const std::string& attribute_type_name, 
        AttributeSerializer* serializer
    ) {
        grgmesh_assert(resolve_by_type(attribute_type) == nil) ;
        grgmesh_assert(resolve_by_name(attribute_type_name) == nil) ;
        (*type_to_serializer_)[attribute_type.name()] = serializer ;
        (*name_to_serializer_)[attribute_type_name] = serializer ;
        (*type_to_name_)[attribute_type.name()] = attribute_type_name ;
    }
    
}