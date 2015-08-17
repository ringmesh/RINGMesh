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

#include <ringmesh/boundary_model.h>
#include <ringmesh/boundary_model_builder.h>
#include <ringmesh/io.h>
#include <ringmesh/utils.h>

namespace {
    using namespace RINGMesh ;

    bool check_with_expected_result(
        const BoundaryModelBuilder& builder,
        std::set< BME::bme_t >& initial_elts,
        const std::set< BME::bme_t >& result )
    {
        builder.get_dependent_elements( initial_elts ) ;
        if( initial_elts.size() != result.size() ) {
            return false ;
        }

        for( std::set< BME::bme_t >::const_iterator bme_itr = result.begin();
            bme_itr != result.end(); ++bme_itr ) {
            const BME::bme_t& cur_bme_t = *bme_itr ;
            if( initial_elts.find( cur_bme_t ) == initial_elts.end() ) {
                return false ;
            }
        }

        // If same size and all the elt in one set are in the other set, so we
        // are good!
        return true ;
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    GEO::Logger::out( "TEST" ) << "Test get_dependent_elements" << std::endl ;

    BoundaryModel in ;
    if( !RINGMeshIO::load( "../data/model1.ml", in ) ) {
        GEO::Logger::out( "TEST" ) << "FAILED" << std::endl ;
        return 1 ;
    }
    BoundaryModelBuilder builder( in ) ;

    // First test
    index_t reg_id = in.find_element( BoundaryModelElement::REGION, "Region_3" ) ;
    if( reg_id == NO_ID ) {
        GEO::Logger::out( "TEST" ) << "FAILED" << std::endl ;
        return 1 ;
    }
    std::set< BME::bme_t > elt_set ;
    elt_set.insert( in.region( reg_id ).bme_id() ) ;

    std::set< BME::bme_t > result ;
    result.insert( BME::bme_t( BoundaryModelElement::REGION, 1 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CORNER, 12 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CORNER, 13 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CORNER, 18 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CORNER, 19 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 12 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 13 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 14 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 25 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 26 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 27 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 34 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 35 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 3 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 8 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 12 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 16 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 20 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CONTACT, 14 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CONTACT, 20 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CONTACT, 22 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CONTACT, 23 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::INTERFACE, 8 ) ) ;

    if( !check_with_expected_result( builder, elt_set, result ) ) {
        GEO::Logger::out( "TEST" ) << "FAILED" << std::endl ;
        return 1 ;
    }

    elt_set.clear() ;
    result.clear() ;
    // Second test
    reg_id = in.find_element( BoundaryModelElement::REGION, "h1_model1_1" ) ;
    if( reg_id == NO_ID ) {
        GEO::Logger::out( "TEST" ) << "FAILED" << std::endl ;
        return 1 ;
    }
    elt_set.insert( in.region( reg_id ).bme_id() ) ;

    result.insert( BME::bme_t( BoundaryModelElement::REGION, 3 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CORNER, 14 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CORNER, 15 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CORNER, 16 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CORNER, 17 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 19 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 20 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 21 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 22 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 23 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 24 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 32 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::LINE, 33 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 6 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 7 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 11 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 15 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::SURFACE, 19 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CONTACT, 15 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CONTACT, 16 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CONTACT, 17 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::CONTACT, 18 ) ) ;
    result.insert( BME::bme_t( BoundaryModelElement::INTERFACE, 4 ) ) ;

    if( !check_with_expected_result( builder, elt_set, result ) ) {
        GEO::Logger::out( "TEST" ) << "FAILED" << std::endl ;
        return 1 ;
    }

    GEO::Logger::out( "TEST" ) << "SUCCESS" << std::endl ;
    return 0 ;
}
