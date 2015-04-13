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
#include <ringmesh/io.h>

#include <cstring>

bool compare_file( const std::string& f1, const std::string& f2 )
{
    const unsigned int MAX_LINE_LEN = 65535 ;

    std::ifstream lFile( f1.c_str() ) ;
    std::ifstream rFile( f2.c_str() ) ;

    char* lBuffer = new char[MAX_LINE_LEN]() ;
    char* rBuffer = new char[MAX_LINE_LEN]() ;

    do {
        lFile.read( lBuffer, MAX_LINE_LEN ) ;
        rFile.read( rBuffer, MAX_LINE_LEN ) ;
        unsigned int numberOfRead = lFile.gcount() ;

        if( std::memcmp( lBuffer, rBuffer, numberOfRead ) != 0 ) {
            delete[] lBuffer ;
            delete[] rBuffer ;
            return false ;
        }
    } while( lFile.good() || rFile.good() ) ;
    return true ;
}

int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    GEO::Logger::out("TEST") << "Test IO for a BoundaryModel in .ml" << std::endl ;

    BoundaryModel in ;
    RINGMeshIO::load( "../data/model1.ml", in ) ;
    RINGMeshIO::save( in, "out.ml" ) ;

    BoundaryModel in2 ;
    RINGMeshIO::load( "out.ml", in2 ) ;
    RINGMeshIO::save( in2, "out2.ml" ) ;

    bool res = compare_file( "out.ml", "out2.ml" ) ;
    if( res )
        GEO::Logger::out("TEST") << "SUCCES" << std::endl ;
    else
        GEO::Logger::out("TEST") << "FAILED" << std::endl ;
    return res ;
}
