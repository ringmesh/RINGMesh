/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#pragma once

#include <ringmesh/basic/common.h>

/*!
 * @file Box3D class declaration
 * @author Arnaud Botella
 */

namespace RINGMesh {

    class RINGMESH_API Box3d {
    public:
        Box3d()
            : initialized_( false )
        {
        }

        bool initialized() const
        {
            return initialized_;
        }

        void clear()
        {
            initialized_ = false;
        }

        double width() const
        {
            return max_[0] - min_[0];
        }

        double height() const
        {
            return max_[1] - min_[1];
        }

        double depth() const
        {
            return max_[2] - min_[2];
        }

        const vec3& min() const
        {
            return min_;
        }

        const vec3& max() const
        {
            return max_;
        }

        vec3 center() const
        {
            return 0.5 * ( min() + max() );
        }

        vec3 diagonal() const
        {
            return max() - min();
        }

        void add_point( const vec3& p );

        void add_box( const Box3d& b )
        {
            if( b.initialized() ) {
                add_point( b.min() );
                add_point( b.max() );
            }
        }

        bool bboxes_overlap( const Box3d& B ) const
        {
            for( index_t c = 0; c < 3; ++c ) {
                if( max()[c] < B.min()[c] ) {
                    return false;
                }
                if( min()[c] > B.max()[c] ) {
                    return false;
                }
            }
            return true;
        }

        Box3d bbox_union( const Box3d& B ) const
        {
            Box3d result = *this;
            result.add_box( B );
            return result;
        }

        bool contains( const vec3& b ) const
        {
            for( index_t c = 0; c < 3; ++c ) {
                if( b[c] < min()[c] ) {
                    return false;
                }
                if( b[c] > max()[c] ) {
                    return false;
                }
            }
            return true;
        }

        double distance_to_center( const vec3& p ) const;

        double signed_distance( const vec3& p ) const;

    private:
        bool initialized_;
        vec3 min_;
        vec3 max_;

    };

}
