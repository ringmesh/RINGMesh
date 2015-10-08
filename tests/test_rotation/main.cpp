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
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/geometry.h>
#include <geogram/basic/logger.h>

using namespace RINGMesh ;

namespace {
    bool are_equal_matrices(
        const GEO::Matrix< float64, 4 >& lhs,
        const GEO::Matrix< float64, 4 >& rhs )
    {
        for( index_t mat_i = 0; mat_i < 4; ++mat_i ) {
            for( index_t mat_j = 0; mat_j < 4; ++mat_j ) {
                float64 diff = lhs( mat_i, mat_j ) - rhs( mat_i, mat_j ) ;
                if( diff > epsilon || diff < -epsilon ) {
                    GEO::Logger::out( "TEST" ) << "Error at " << mat_i << " , "
                        << mat_j << std::endl ;
                    GEO::Logger::out( "TEST" ) << "One is " << lhs( mat_i, mat_j )
                        << std::endl ;
                    GEO::Logger::out( "TEST" ) << "Other is " << rhs( mat_i, mat_j )
                        << std::endl ;
                    return false ;
                }
            }
        }
        return true ;
    }
}

int main( int argc, char** argv )
{

    GEO::Logger::out( "TEST" )
        << "Test rotations of a boundary model and a macro mesh" << std::endl ;

    const vec3 origin( 0, 0, 0 ) ;
    const float64 pi = 3.14159265359 ;
    const float64 step = 0.1 ;

    GEO::Matrix< float64, 4 > rot_mat_degree ;
    GEO::Matrix< float64, 4 > rot_mat_radian ;
    GEO::Matrix< float64, 4 > result ;
    result( 0, 3 ) = 0 ;
    result( 1, 3 ) = 0 ;
    result( 2, 3 ) = 0 ;
    result( 3, 0 ) = 0 ;
    result( 3, 1 ) = 0 ;
    result( 3, 2 ) = 0 ;
    result( 3, 3 ) = 1 ;

    // Tests rotation along x axis
    vec3 axis( 1, 0, 0 ) ;
    for( float64 angle = 0.; angle <= 360.; angle += step ) {
        rotation_matrix_about_arbitrary_axis( origin, axis, angle, true,
            rot_mat_degree ) ;
        float64 angle_rad = angle * pi / 180. ;
        rotation_matrix_about_arbitrary_axis( origin, axis, angle_rad, false,
            rot_mat_radian ) ;
        result( 0, 0 ) = 1 ;
        result( 0, 1 ) = 0 ;
        result( 0, 2 ) = 0 ;

        result( 1, 0 ) = 0 ;
        result( 1, 1 ) = std::cos( angle_rad ) ;
        result( 1, 2 ) = -std::sin( angle_rad ) ;

        result( 2, 0 ) = 0 ;
        result( 2, 1 ) = std::sin( angle_rad ) ;
        result( 2, 2 ) = std::cos( angle_rad ) ;

        if( !are_equal_matrices( rot_mat_degree, result ) ) {
            GEO::Logger::out( "TEST" ) << "FAILED for axis x for angle " << angle
                << " degrees." << std::endl ;
            return 1 ;
        }

        if( !are_equal_matrices( rot_mat_radian, result ) ) {
            GEO::Logger::out( "TEST" ) << "FAILED for axis x for angle " << angle_rad
                << " radians." << std::endl ;
            return 1 ;
        }
    }

    // Tests rotation along y axis
    axis = vec3( 0, 1, 0 ) ;
    for( float64 angle = 0.; angle <= 360.; angle += step ) {
        rotation_matrix_about_arbitrary_axis( origin, axis, angle, true,
            rot_mat_degree ) ;
        float64 angle_rad = angle * pi / 180. ;
        rotation_matrix_about_arbitrary_axis( origin, axis, angle_rad, false,
            rot_mat_radian ) ;
        result( 0, 0 ) = std::cos( angle_rad ) ;
        result( 0, 1 ) = 0 ;
        result( 0, 2 ) = std::sin( angle_rad ) ;

        result( 1, 0 ) = 0 ;
        result( 1, 1 ) = 1 ;
        result( 1, 2 ) = 0 ;

        result( 2, 0 ) = -std::sin( angle_rad ) ;
        result( 2, 1 ) = 0 ;
        result( 2, 2 ) = std::cos( angle_rad ) ;

        if( !are_equal_matrices( rot_mat_degree, result ) ) {
            GEO::Logger::out( "TEST" ) << "FAILED for axis y for angle " << angle
                << " degrees." << std::endl ;
            return 1 ;
        }

        if( !are_equal_matrices( rot_mat_radian, result ) ) {
            GEO::Logger::out( "TEST" ) << "FAILED for axis y for angle " << angle_rad
                << " radians." << std::endl ;
            return 1 ;
        }
    }

    // Tests rotation along z axis
    axis = vec3( 0, 0, 1 ) ;
    for( float64 angle = 0.; angle <= 360.; angle += step ) {
        rotation_matrix_about_arbitrary_axis( origin, axis, angle, true,
            rot_mat_degree ) ;
        float64 angle_rad = angle * pi / 180. ;
        rotation_matrix_about_arbitrary_axis( origin, axis, angle_rad, false,
            rot_mat_radian ) ;
        result( 0, 0 ) = std::cos( angle_rad ) ;
        result( 0, 1 ) = -std::sin( angle_rad ) ;
        result( 0, 2 ) = 0 ;

        result( 1, 0 ) = std::sin( angle_rad ) ;
        result( 1, 1 ) = std::cos( angle_rad ) ;
        result( 1, 2 ) = 0 ;

        result( 2, 0 ) = 0 ;
        result( 2, 1 ) = 0 ;
        result( 2, 2 ) = 1 ;

        if( !are_equal_matrices( rot_mat_degree, result ) ) {
            GEO::Logger::out( "TEST" ) << "FAILED for axis z for angle " << angle
                << " degrees." << std::endl ;
            return 1 ;
        }

        if( !are_equal_matrices( rot_mat_radian, result ) ) {
            GEO::Logger::out( "TEST" ) << "FAILED for axis z for angle " << angle_rad
                << " radians." << std::endl ;
            return 1 ;
        }
    }

    GEO::Logger::out( "TEST" ) << "SUCCESS" << std::endl ;
    return 0 ;
}
