/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

#pragma once

#include <array>

#include <ringmesh/basic/common.h>

/*!
 * @file Basic classes for representing frames
 * @author Pierre Anquez
 */

namespace RINGMesh
{
    namespace Geometry
    {
        struct Plane;
        FORWARD_DECLARATION_DIMENSION_STRUCT( Line );
        ALIAS_2D( Line );
    }
}

namespace RINGMesh
{
    template < index_t DIMENSION >
    class basic_api FrameBase
    {
    public:
        virtual ~FrameBase() = default;

        const vecn< DIMENSION >& operator[]( index_t coord ) const
        {
            ringmesh_assert( coord < DIMENSION );
            return axis_[coord];
        }

        vecn< DIMENSION >& operator[]( index_t coord )
        {
            ringmesh_assert( coord < DIMENSION );
            return axis_[coord];
        }

    protected:
        // Default Frame aligned on space axis
        FrameBase()
        {
            for( auto coord : RINGMesh::range( DIMENSION ) )
            {
                axis_[coord][coord] = 1.;
            }
        }

    private:
        std::array< vecn< DIMENSION >, DIMENSION > axis_;
    };
    ALIAS_2D_AND_3D( FrameBase );

    template < index_t DIMENSION >
    class basic_api Frame : public FrameBase< DIMENSION >
    {
    public:
        Frame() = default;
    };
    ALIAS_2D_AND_3D( Frame );

    template <>
    class basic_api Frame< 2 > : public FrameBase< 2 >
    {
    public:
        Frame() = default;

        Frame( vec2 x, vec2 y )
        {
            ( *this )[0] = std::move( x );
            ( *this )[1] = std::move( y );
        }
    };

    template <>
    class basic_api Frame< 3 > : public FrameBase< 3 >
    {
    public:
        Frame() = default;

        Frame( vec3 x, vec3 y, vec3 z )
        {
            ( *this )[0] = std::move( x );
            ( *this )[1] = std::move( y );
            ( *this )[2] = std::move( z );
        }
    };

    template < index_t DIMENSION >
    class basic_api ReferenceFrame : public Frame< DIMENSION >
    {
    public:
        ReferenceFrame() = default;

        ReferenceFrame(
            vecn< DIMENSION > frame_origin, Frame< DIMENSION > frame )
            : Frame< DIMENSION >( std::move( frame ) ),
              origin_( std::move( frame_origin ) )
        {
        }

        const vecn< DIMENSION >& origin() const
        {
            return origin_;
        }

        vecn< DIMENSION >& origin()
        {
            return origin_;
        }

        /*!
         * Two reference frames are considered equal if they have the same origin
         * and the same axis in the same order.
         */
        bool operator==(
            const ReferenceFrame< DIMENSION >& other_reference_frame ) const
        {
            if( other_reference_frame.origin() != origin_ )
            {
                return false;
            }
            for( auto i : range( DIMENSION ) )
            {
                if( other_reference_frame[i] != ( *this )[i] )
                {
                    return false;
                }
            }
            return true;
        }

        bool operator!=(
            const ReferenceFrame< DIMENSION >& other_reference_frame ) const
        {
            return !( other_reference_frame == *this );
        }

    private:
        vecn< DIMENSION > origin_{};
    };
    ALIAS_2D_AND_3D( ReferenceFrame );

    /*!
     * Verifies if the two input reference frames are almost the same:
     * @return true if the origin and each each of the vectors of the input
     * reference frames are close up to epsilon, false otherwise.
     */
    template < index_t DIMENSION >
    bool inexact_equal( const ReferenceFrame< DIMENSION >& v1,
        const ReferenceFrame< DIMENSION >& v2,
        double epsilon )
    {
        if( !inexact_equal( v1.origin(), v2.origin(), epsilon ) )
        {
            return false;
        }
        for( auto i : range( DIMENSION ) )
        {
            if( !inexact_equal( v1[i], v2[i], epsilon ) )
            {
                return false;
            }
        }
        return true;
    }

    /*!
     * @brief Class for manipulating reference frames coordinates and such.
     * All of its functions are static.
     */
    template < index_t DIMENSION >
    class basic_api ReferenceFrameManipulator
    {
    public:
    	/*!
    	 * @return the coordinates of the input point \global_coords in the
    	 * input \reference_frame.
    	 */
        static vecn< DIMENSION > coords_from_global_to_frame(
            const ReferenceFrame< DIMENSION >& reference_frame,
            const vecn< DIMENSION >& global_coords )
        {
            GEO::Matrix< DIMENSION, double > base_change_matrix =
                inverse_reference_matrix( reference_frame );
            vecn< DIMENSION > local_coords;
            for( auto i : range( DIMENSION ) )
            {
                for( auto x : range( DIMENSION ) )
                {
                    local_coords[i] +=
                        base_change_matrix( i, x )
                        * ( global_coords[x] - reference_frame.origin()[x] );
                }
            }
            return local_coords;
        }

    	/*!
    	 * @return the global coordinates of the input point from its coordinates
    	 * \local_coords in the input \reference_frame.
    	 */
        static vecn< DIMENSION > coords_from_frame_to_global(
            const ReferenceFrame< DIMENSION >& reference_frame,
            const vecn< DIMENSION >& local_coords )
        {
            vecn< DIMENSION > global_coords{ reference_frame.origin() };
            for( auto x : range( DIMENSION ) )
            {
                for( auto i : range( DIMENSION ) )
                {
                    global_coords[x] += local_coords[i] * reference_frame[i][x];
                }
            }
            return global_coords;
        }

        /*!
         * @return the inverse of the input reference_frame (equivalent to the
         * coordinates of the global frame in the input reference frame).
         */
        static ReferenceFrame< DIMENSION > reference_frame_from_global_to_local(
            const ReferenceFrame< DIMENSION >& reference_frame )
        {
            ReferenceFrame< DIMENSION > inverse_reference_frame;
            GEO::Matrix< DIMENSION, double > base_change_matrix{
                inverse_reference_matrix( reference_frame )
            };
            for( auto i : range( DIMENSION ) )
            {
                for( auto j : range( DIMENSION ) )
                {
                    inverse_reference_frame.origin()[i] -=
                        base_change_matrix( i, j )
                        * reference_frame.origin()[j];
                    inverse_reference_frame[i][j] = base_change_matrix( j, i );
                }
            }
            return inverse_reference_frame;
        }

        /*!
         * @return the inverse of the input reference_frame (equivalent to the
         * coordinates of the global frame in the input reference frame).
         * This method is optimized for orthogonal frames and works only on them.
         */
        static ReferenceFrame< DIMENSION >
            orthogonal_reference_frame_from_global_to_local(
                const ReferenceFrame< DIMENSION >& reference_frame )
        {
            ReferenceFrame< DIMENSION > inverse_reference_frame;
            for( auto i : range( DIMENSION ) )
            {
                double square_length_i{ reference_frame[i].length()
                                        * reference_frame[i].length() };
                for( auto j : range( DIMENSION ) )
                {
                    inverse_reference_frame.origin()[i] -=
                        reference_frame.origin()[j] * reference_frame[i][j]
                        / square_length_i;
                    inverse_reference_frame[j][i] =
                        reference_frame[i][j] / square_length_i;
                }
            }
            return inverse_reference_frame;
        }

        /*!
         * Verifies if the given reference frame is orthogonal.
         */
        static bool is_frame_orthogonal(
            const ReferenceFrame< DIMENSION >& reference_frame )
        {
            for( auto i : range( DIMENSION ) )
            {
                for( auto j : range( i + 1, DIMENSION ) )
                {
                    if( dot( reference_frame[i], reference_frame[j] ) != 0 )
                    {
                        return false;
                    }
                }
            }
            return true;
        }

    private:
        static GEO::Matrix< DIMENSION, double > inverse_reference_matrix(
            const ReferenceFrame< DIMENSION >& reference_frame )
        {
            GEO::Matrix< DIMENSION, double > matrix;
            for( auto x : range( DIMENSION ) )
            {
                for( auto i : range( DIMENSION ) )
                {
                    matrix( x, i ) = reference_frame[i][x];
                }
            }
            return matrix.inverse();
        }
    };
    ALIAS_2D_AND_3D( ReferenceFrameManipulator );

    /*!
     * @brief Reference frame aligned along the plane normal and whose u axis is
     * upward
     */
    class basic_api PlaneReferenceFrame3D : public ReferenceFrame< 3 >
    {
    public:
        PlaneReferenceFrame3D() = default;

        PlaneReferenceFrame3D( vec3 frame_origin, Frame3D frame )
            : ReferenceFrame3D( std::move( frame_origin ), std::move( frame ) )
        {
        }

        explicit PlaneReferenceFrame3D( const Geometry::Plane& plane );
    };

    /*!
     * @brief Reference frame aligned along the plane normal and whose u axis is
     * upward
     */
    class basic_api LineReferenceFrame2D : public ReferenceFrame< 2 >
    {
    public:
        LineReferenceFrame2D() = default;

        LineReferenceFrame2D( vec2 frame_origin, Frame2D frame )
            : ReferenceFrame2D( std::move( frame_origin ), std::move( frame ) )
        {
        }

        explicit LineReferenceFrame2D( const Geometry::Line2D& line );
    };
} // namespace RINGMesh
