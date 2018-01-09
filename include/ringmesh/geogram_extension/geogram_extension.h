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

#include <ringmesh/geogram_extension/common.h>

#include <mutex>

#include <geogram/basic/attributes.h>
#include <geogram/basic/memory.h>

/*!
 * @file Helper functions on classes defined in Geogram
 * @author Various
 */

namespace GEO
{
    class Mesh;
} // namespace GEO

namespace RINGMesh
{
    /*!
     * Partial copy the content of a standrad library vector to a GEO::Vector.
     * A lot of copies, when we need to call Geogram functions.
     * @todo Could we set Geogram vector to be a std::vector ??
     */
    template < typename T, typename U = T >
    GEO::vector< U > copy_std_vector_to_geo_vector(
        const std::vector< T >& in, index_t from, index_t to )
    {
        ringmesh_assert( to < in.size() + 1 );
        ringmesh_assert( from < to );
        index_t nb_to_copy( to - from );
        GEO::vector< U > out( nb_to_copy );
        for( auto i : range( nb_to_copy ) )
        {
            out[i] = in[from + i];
        }
        return out;
    }

    /*!
     * Copy the content of a standard library vector to the memory aligned
     * GEO::Vector.
     * A lot of copies, when we need to call Geogram functions.
     * @todo Could we set Geogram vector to be a std::vector ??
     */
    template < typename T, typename U = T >
    GEO::vector< U > copy_std_vector_to_geo_vector( const std::vector< T >& in )
    {
        index_t size = static_cast< index_t >( in.size() );
        return copy_std_vector_to_geo_vector< T, U >( in, 0, size );
    }

    /***********************************************************************/
    /* Loading and saving a GEO::Mesh                                      */

    /*! @brief complement the available MeshIOHandler
     */
    void geogram_extension_api ringmesh_geogram_mesh_io_initialize();

    /******************************************************************/
    /* Operations on a GEO::Mesh                                      */

    /*!
     * Computes the signed volume of a Mesh cell
     * @param[in] M the mesh
     * @param[in] c the cell index
     * @return the signed volume of the cell
     */
    double geogram_extension_api mesh_cell_signed_volume(
        const GEO::Mesh& M, index_t c );
    /*!
     * Computes the volume of a Mesh cell
     * @param[in] M the mesh
     * @param[in] c the cell index
     * @return the volume of the cell
     */
    double geogram_extension_api mesh_cell_volume(
        const GEO::Mesh& M, index_t c );

    /*!
     * Computes the Mesh cell facet barycenter
     * @param[in] M the mesh
     * @param[in] cell the cell index
     * @param[in] f the facet index in the cell
     * @return the cell facet center
     */
    vec3 geogram_extension_api mesh_cell_facet_barycenter(
        const GEO::Mesh& M, index_t cell, index_t f );

    /*!
     * Computes the non weighted barycenter of a volumetric
     * cell of a Mesh
     * @param[in] M the mesh
     * @param[in] cell the cell index
     * @return the cell center
     */
    vec3 geogram_extension_api mesh_cell_barycenter(
        const GEO::Mesh& M, index_t cell );

    /*!
     * @brief Vector of pointers to Geogram attributes
     * @note Necessary since one cannot create, vectors of Geogram attributes
     * does
     * not compile, because @#$# (no idea) [JP]
     * @todo Probably extremely prone to bugs. Is it worth the risk?
     */
    template < typename T >
    class AttributeVector : public std::vector< GEO::Attribute< T >* >
    {
        ringmesh_disable_copy_and_move( AttributeVector );

    public:
        using base_class = std::vector< GEO::Attribute< T >* >;
        AttributeVector() = default;
        explicit AttributeVector( index_t size ) : base_class( size, nullptr )
        {
        }

        void bind_one_attribute( index_t i,
            GEO::AttributesManager& manager,
            const std::string& attribute_name )
        {
            base_class::operator[]( i ) =
                new GEO::Attribute< T >( manager, attribute_name );
        }

        GEO::Attribute< T >& operator[]( index_t i )
        {
            return *base_class::operator[]( i );
        }

        const GEO::Attribute< T >& operator[]( index_t i ) const
        {
            return *base_class::operator[]( i );
        }

        bool is_attribute_bound( index_t i ) const
        {
            return base_class::operator[]( i ) != nullptr;
        }

        void unbind( index_t i )
        {
            if( base_class::operator[]( i ) )
            {
                // I am not sure, but unbind should do the deallocation [JP]
                operator[]( i ).unbind();
                delete base_class::operator[]( i );
                base_class::operator[]( i ) = nullptr;
            }
        }

        ~AttributeVector()
        {
            for( auto i : range( base_class::size() ) )
            {
                unbind( i );
            }
        }
    };

} // namespace RINGMesh
