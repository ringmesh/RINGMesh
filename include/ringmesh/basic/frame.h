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
    class RINGMESH_API FrameBase
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
    class RINGMESH_API Frame : public FrameBase< DIMENSION >
    {
    public:
        Frame() = default;
    };
    ALIAS_2D_AND_3D( Frame );

    template <>
    class RINGMESH_API Frame< 2 > : public FrameBase< 2 >
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
    class RINGMESH_API Frame< 3 > : public FrameBase< 3 >
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
    class RINGMESH_API ReferenceFrame : public Frame< DIMENSION >
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

        vecn< DIMENSION > coords_to_local(
            const vecn< DIMENSION > global_coords ) const;

        vecn< DIMENSION > coords_to_global(
            const vecn< DIMENSION > local_coords ) const
        {
            vecn< DIMENSION > global_coords;
            for( auto coord : RINGMesh::range( DIMENSION ) )
            {
                global_coords[coord] = ( *this ).origin()[coord];
                for( auto coor : RINGMesh::range( DIMENSION ) )
                {
                    global_coords[coord] +=
                        local_coords[coor] * ( *this )[coor][coord];
                }
            }
            return global_coords;
        }

    private:
        vecn< DIMENSION > origin_{};
    };
    ALIAS_2D_AND_3D( ReferenceFrame );

    class RINGMESH_API CGFrame3D : public ReferenceFrame< 3 >
    {
    public:
        CGFrame3D() = default;

        CGFrame3D( vec3 frame_origin, Frame3D frame )
            : ReferenceFrame3D( std::move( frame_origin ), std::move( frame ) )
        {
            base_change_ = ReferenceFrame3D();
            GEO::Matrix< 3, double > bchange;
            for( index_t i = 0; i < 3; i++ )
            {
                for( index_t j = 0; j < 3; j++ )
                {
                    bchange( i, j ) = frame[j][i];
                }
            }
            bchange = bchange.inverse();
            for( index_t i = 0; i < 3; i++ )
            {
                base_change_.origin()[i] = 0;
                for( index_t j = 0; j < 3; j++ )
                {
                    base_change_.origin()[i] -=
                        bchange( i, j ) * ( *this ).origin()[j];
                    base_change_[i][j] = bchange( j, i );
                }
            }
        }

        vec3 coords_to_local( const vec3 global_coords ) const
        {
            vec3 local_coords;
            local_coords[0] = base_change_.origin()[0]
                              + global_coords[0] * base_change_[0][0]
                              + global_coords[1] * base_change_[1][0]
                              + global_coords[2] * base_change_[2][0];
            local_coords[1] = base_change_.origin()[1]
                              + global_coords[0] * base_change_[0][1]
                              + global_coords[1] * base_change_[1][1]
                              + global_coords[2] * base_change_[2][1];
            local_coords[2] = base_change_.origin()[2]
                              + global_coords[0] * base_change_[0][2]
                              + global_coords[1] * base_change_[1][2]
                              + global_coords[2] * base_change_[2][2];
            return local_coords;
        }

    protected:
        ReferenceFrame3D base_change_;
    };

    /*!
     * @brief Reference frame aligned along the plane normal and whose u axis is
     * upward
     */
    class RINGMESH_API PlaneReferenceFrame3D : public ReferenceFrame< 3 >
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
    class RINGMESH_API LineReferenceFrame2D : public ReferenceFrame< 2 >
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
