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
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/io.h>

#include <fstream>
#include <cstring>

namespace RINGMesh {

    /*!
    * Compares the contains of two files
    * @param[in] f1 the first filename
    * @param[in] f2 the second filename
    * @return return True if the files are identical
    */
    bool compare_files( const std::string& f1, const std::string& f2 )
    {
        const unsigned int MAX_LINE_LEN = 65535 ;

        std::ifstream lFile( f1.c_str() ) ;
        std::ifstream rFile( f2.c_str() ) ;

        char* lBuffer = new char[ MAX_LINE_LEN ]() ;
        char* rBuffer = new char[ MAX_LINE_LEN ]() ;

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
        delete[] lBuffer ;
        delete[] rBuffer ;
        return true ;
    }
   
    void mesh_initialize()
    {
        GeoModelSurfaceIOHandler::initialize() ;
        GeoModelVolumeIOHandler::initialize() ;
        WellGroupIOHandler::initialize() ;
    }

    /*!
    * Read the coordinates system informations of files exported from Gocad.
    * @param[in] in The orientation of z-axis in Gocad. "Elevation" for
    * increasing z toward top and "Depth" for increasing z toward bottom.
    * @return Return 1 if Elevation direction, -1 if Depth direction.
    */
    int read_gocad_coordinates_system( const char* in )
    {
    	std::string str_in( in );
		if( str_in == "Elevation" ) {
			return 1 ;
		} else if( str_in == "Depth" ) {
			return -1 ;
		} else {
			ringmesh_assert_not_reached;
		}

    }
}
