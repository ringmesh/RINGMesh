/*
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

#include <ringmesh/basic/frame.h>

#include <ringmesh/basic/geometry.h>

/*!
 * @file Basic classes for representing frame
 * @author Pierre Anquez
 */

namespace RINGMesh
{
    PlaneReferenceFrame3D::PlaneReferenceFrame3D( const Geometry::Plane& plane )
    {
        origin() = plane.origin;
        ( *this )[2] = plane.normal;

        // @todo A generic algorithm to find the first vector belonging to the
        // plane
        // can be designed using principal component of the plane normal.
        // However it is not a simple problem. The current version is not
        // generic
        // and is based on the idea that the plane is either a map section
        // or a cross-section. [PA]
        vec3 another_point_for_v_axis{ origin() };
        if( std::fabs( ( *this )[2].z )
            > ( std::fabs( ( *this )[2].x ) + std::fabs( ( *this )[2].y ) ) )
        {
            // Case where plane is sub-horizontal
            // (v axis is set towards 3D y direction)
            another_point_for_v_axis += vec3{ 0., 1., 0. };
        }
        else
        {
            // Case where plane is not sub-horizontal
            // (v axis is set towards 3D z direction)
            another_point_for_v_axis += vec3{ 0., 0., 1. };
        }
        vec3 v_axis_point;
        std::tie( std::ignore, v_axis_point ) = Distance::point_to_plane(
            another_point_for_v_axis, { ( *this )[2], origin() } );
        ringmesh_assert(
            ( origin() - v_axis_point ).length() > global_epsilon );
        ( *this )[1] = normalize( v_axis_point - origin() );
        ( *this )[0] = cross( ( *this )[1], ( *this )[2] );
    }

    template class basic_api FrameBase< 2 >;
    template class basic_api ReferenceFrame< 2 >;
    template class basic_api ReferenceFrameManipulator< 2 >;

    template class basic_api FrameBase< 3 >;
    template class basic_api ReferenceFrame< 3 >;
    template class basic_api ReferenceFrameManipulator< 3 >;
} // namespace RINGMesh
