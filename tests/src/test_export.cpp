/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/ringmesh_tests_config.h>

#include <vector>

#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/io/io.h>

/*!
 * @author Antoine MAzuyer
 */

using namespace RINGMesh ;

static const std::string model_to_test = "modelA6" ;

void test_file( const GeoModel& geomodel, const std::string& file_extension )
{
    std::string file_name_out( ringmesh_test_output_path ) ;
    file_name_out += model_to_test ;
    file_name_out += file_extension ;
    std::string file_name_ref( ringmesh_test_data_path ) ;
    file_name_ref += model_to_test ;
    file_name_ref += file_extension ;
    geomodel_save( geomodel, file_name_out ) ;

    if( !compare_files( file_name_out, file_name_ref ) ) {
        throw RINGMeshException( "TEST", "Aster export is broken" ) ;
    }
}

const static std::string aster_format = ".mail" ;
static const std::string file_format_to_test[1] = { aster_format } ;

int main()
{
    using namespace RINGMesh ;
    try {
        default_configure() ;

        Logger::out( "TEST" ) << "Test Export" << std::endl ;
        std::string file_name( ringmesh_test_data_path ) ;
        file_name += model_to_test ;
        file_name += ".ml" ;
        GeoModel geomodel ;
        geomodel_load( geomodel, file_name ) ;

        for( index_t i = 0; i < 1; i++ ) {
            test_file( geomodel, file_format_to_test[i] ) ;
        }

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    Logger::out( "TEST" ) << "SUCCESS" << std::endl ;
    return 0 ;
}
