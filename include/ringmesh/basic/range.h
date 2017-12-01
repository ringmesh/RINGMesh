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

#include <ringmesh/basic/types.h>

namespace RINGMesh
{
    /*!
     * This class can be used to iterate over integer loop.
     * Example:
     *              = C++98 loop =
     *    for( index_t i = 0; i < n; i++ ) {
     *      // do something
     *    }
     *
     *            = C++11-like loop =
     *    for( index_t i : range( n ) ) {
     *      // do something
     *    }
     */
    class range
    {
    public:
        template < typename T1, typename T2 >
        range( T1 begin, T2 end )
            : iter_( static_cast< index_t >( begin ) ),
              last_( static_cast< index_t >( end ) )
        {
        }
        template < typename T >
        explicit range( T end ) : last_( static_cast< index_t >( end ) )
        {
        }
        // Iterable functions
        const range& begin() const
        {
            return *this;
        }
        const range& end() const
        {
            return *this;
        }
        // Iterator functions
        bool operator!=( const range& /*unused*/ ) const
        {
            return iter_ < last_;
        }
        void operator++()
        {
            ++iter_;
        }
        index_t operator*() const
        {
            return iter_;
        }

    protected:
        index_t iter_{ 0 };
        index_t last_{ 0 };
    };

} // namespace RINGMesh
