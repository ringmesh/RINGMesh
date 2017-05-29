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
 * @file BoxND class declaration
 * @author Arnaud Botella
 */

namespace RINGMesh {

    template< index_t DIMENSION >
    class RINGMESH_API Box {
    public:
        Box()
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

        double depth() const //@todo not defined in 2D
        {
            return max_[2] - min_[2];
        }

        const vecn<DIMENSION>& min() const
        {
            return min_;
        }

        const vecn<DIMENSION>& max() const
        {
            return max_;
        }

        vecn<DIMENSION> center() const
        {
            return 0.5 * ( min() + max() );
        }

        vecn<DIMENSION> diagonal() const
        {
            return max() - min();
        }

        void add_point( const vecn<DIMENSION>& p );

        void add_box( const Box& b )
        {
            if( b.initialized() ) {
                add_point( b.min() );
                add_point( b.max() );
            }
        }

        inline bool bboxes_overlap( const Box& B ) const
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

        inline Box bbox_union( const Box& B ) const
        {
            Box result = *this;
            result.add_box( B );
            return result;
        }

        bool contains( const vecn<DIMENSION>& b ) const
        {
            for( index_t c = 0; c < DIMENSION; ++c ) {
                if( b[c] < min()[c] ) {
                    return false;
                }
                if( b[c] > max()[c] ) {
                    return false;
                }
            }
            return true;
        }

        double distance_to_center( const vecn<DIMENSION>& p ) const;

        double signed_distance( const vecn<DIMENSION>& p ) const;

    private:
        bool initialized_;
        vecn<DIMENSION> min_;
        vecn<DIMENSION> max_;

    };

}
