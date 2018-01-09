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

/*!
 * @file Implementation of visualization of GeoModelEntities
 * @author Benjamin Chauvin and Arnaud Botella
 */

#include <ringmesh/visualize/geomodel_gfx.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_entity.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

namespace RINGMesh
{
    template < index_t DIMENSION >
    GeoModelGfxBase< DIMENSION >::GeoModelGfxBase(
        GeoModelGfx< DIMENSION >& gfx )
        : corners( gfx ), lines( gfx ), surfaces( gfx ), attribute( gfx )
    {
    }

    template < index_t DIMENSION >
    void GeoModelGfxBase< DIMENSION >::set_geomodel(
        const GeoModel< DIMENSION >& geomodel )
    {
        geomodel_ = &geomodel;
        initialize();
    }

    template < index_t DIMENSION >
    const GeoModel< DIMENSION >* GeoModelGfxBase< DIMENSION >::geomodel() const
    {
        return geomodel_;
    }

    template < index_t DIMENSION >
    void GeoModelGfxBase< DIMENSION >::initialize()
    {
        ringmesh_assert( geomodel_ );
        corners.initialize();
        lines.initialize();
        surfaces.initialize();
    }

    template < index_t DIMENSION >
    GeoModelGfx< DIMENSION >::GeoModelGfx()
        : GeoModelGfxBase< DIMENSION >( *this )
    {
    }

    GeoModelGfx< 3 >::GeoModelGfx()
        : GeoModelGfxBase< 3 >( *this ), regions( *this )
    {
    }

    void GeoModelGfx< 3 >::initialize()
    {
        GeoModelGfxBase3D::initialize();
        regions.initialize();
    }

    template class visualize_api GeoModelGfxBase< 2 >;
    template class visualize_api GeoModelGfx< 2 >;

    template class visualize_api GeoModelGfxBase< 3 >;

} // namespace RINGMesh

#endif
