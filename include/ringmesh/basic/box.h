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

#include <ringmesh/basic/common.h>

/*!
 * @file Box class declaration
 * @author Arnaud Botella
 */

namespace RINGMesh
{
    template < index_t DIMENSION >
    class Box
    {
    public:
        bool initialized() const;

        void clear();

        const vecn< DIMENSION >& min() const;

        const vecn< DIMENSION >& max() const;

        vecn< DIMENSION > center() const;

        vecn< DIMENSION > diagonal() const;

        void add_point( const vecn< DIMENSION >& p );

        void add_box( const Box< DIMENSION >& b );

        bool bboxes_overlap( const Box< DIMENSION >& B ) const;

        Box< DIMENSION > bbox_union( const Box< DIMENSION >& B ) const;

        /*!
         * Computes the intersection between this box and another one
         * @param[in] B another box
         * @return a tuple containing:
         * - a boolean (true if the two boxes do intersect, false otherwise).
         * - the box corresponding to the intersection.
         */
        std::tuple< bool, Box< DIMENSION > > bbox_intersection(
            const Box< DIMENSION >& B ) const;

        bool contains( const vecn< DIMENSION >& b ) const;

        double distance_to_center( const vecn< DIMENSION >& p ) const;

        double signed_distance( const vecn< DIMENSION >& p ) const;

    private:
        bool initialized_{ false };
        vecn< DIMENSION > min_{};
        vecn< DIMENSION > max_{};
    };
    ALIAS_2D_AND_3D( Box );

} // namespace RINGMesh
