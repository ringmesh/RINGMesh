/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#ifndef __RINGMESH_GEOGRAM_EXTENSION__
#define __RINGMESH_GEOGRAM_EXTENSION__

#include <ringmesh/common.h>

#include <geogram/basic/memory.h>
#include <geogram/basic/attributes.h>
#include <geogram/mesh/mesh.h>

/*!
 * @file Helper functions on classes defined in Geogram
 * @author Various
 */

namespace RINGMesh {

    /*!
     * Copy the content of a standrad library vector to the memory aligned GEO::Vector. 
     * A lot of copies, when we need to call Geogram functions. 
     * @todo Could we set Geogram vector to be a std::vector ?? 
     */
    template< typename T >
    void copy_std_vector_to_geo_vector(
        const std::vector< T >& in,
        GEO::vector< T >& out )
    {
        out.resize( in.size() ) ;
        for( index_t i = 0; i < in.size(); ++i ) {
            out[i] = in[i] ;
        }
    }

    /*!
     * Partial copy the content of a standrad library vector to a GEO::Vector.
     * A lot of copies, when we need to call Geogram functions.
     * @todo Could we set Geogram vector to be a std::vector ??
     */
    template< typename T >
    void copy_std_vector_to_geo_vector(
        const std::vector< T >& in,
        index_t from,
        index_t to,
        GEO::vector< T >& out )
    {
        ringmesh_assert( to < in.size() + 1 ) ;
        ringmesh_assert( from < to ) ;
        index_t nb_to_copy( to - from ) ;
        out.resize( nb_to_copy ) ;
        for( index_t i = 0; i != nb_to_copy; ++i ) {
            out[i] = in[from + i] ;
        }
    }

    /***********************************************************************/
    /* Loading and saving a GEO::Mesh                                      */

    /*! @brief complement the available MeshIOHandler
     */
    void RINGMESH_API ringmesh_mesh_io_initialize() ;

    /******************************************************************/
    /* Operations on a GEO::Mesh                                      */


    void RINGMESH_API rotate_mesh(
        GEO::Mesh& mesh,
        const GEO::Matrix< double, 4 >& rot_mat ) ;

    double RINGMESH_API mesh_cell_signed_volume( const GEO::Mesh& M, index_t c ) ;
    double RINGMESH_API mesh_cell_volume( const GEO::Mesh& M, index_t c ) ;

    vec3 RINGMESH_API mesh_cell_facet_center(
        const GEO::Mesh& M,
        index_t cell,
        index_t f ) ;

    vec3 RINGMESH_API mesh_cell_center( const GEO::Mesh& M, index_t cell ) ;

    void RINGMESH_API divide_edge_in_parts(
        const GEO::Mesh& mesh,
        index_t edge,
        index_t nb_parts,
        std::vector< vec3 >& points ) ;

    void RINGMESH_API divide_edge_in_parts(
        vec3& node0,
        vec3& node1,
        index_t nb_parts,
        std::vector< vec3 >& points ) ;

    index_t RINGMESH_API get_nearest_vertex_index(
        const GEO::Mesh& mesh,
        const vec3& p,
        index_t t ) ;

    bool RINGMESH_API facets_have_same_orientation(
        const GEO::Mesh& mesh,
        index_t f1,
        index_t c11,
        index_t f2 ) ;

    void RINGMESH_API mesh_facet_connect( GEO::Mesh& mesh ) ;

    void RINGMESH_API barycentric_coordinates_point_inside_mesh_facet(
        const vec3& point_inside_facet,
        const GEO::Mesh& mesh,
        index_t facet,
        std::vector<double>& barycentric_coordinates ) ;
 

    /*!
     * @brief Vector of pointers to Geogram attributes
     * @note Necessary since one cannot create, vectors of Geogram attributes does
     * not compile, because @#$# (no idea) [JP]
     * @todo Probably extremely prone to bugs. Is it worth the risk?
     */
    template< typename T >
    class AttributeVector: public std::vector< GEO::Attribute< T >* > {
    ringmesh_disable_copy( AttributeVector ) ;
    public:
        typedef std::vector< GEO::Attribute< T >* > base_class ;
        AttributeVector()
            : base_class()
        {
        }
        AttributeVector( index_t size )
            : base_class( size, nil )
        {
        }

        void bind_one_attribute(
            index_t i,
            GEO::AttributesManager& manager,
            const std::string& attribute_name )
        {
            base_class::operator[]( i ) = new GEO::Attribute< T >( manager,
                attribute_name ) ;
        }

        GEO::Attribute< T >& operator[]( index_t i )
        {
            return *base_class::operator[]( i ) ;
        }

        const GEO::Attribute< T >& operator[]( index_t i ) const
        {
            return *base_class::operator[]( i ) ;
        }

        ~AttributeVector()
        {
            for( index_t i = 0; i < base_class::size(); i++ ) {
                if( base_class::operator[]( i ) ) {
                    // I am not sure, but unbind should do the deallocation [JP]
                    operator[]( i ).unbind() ;
                    delete base_class::operator[]( i ) ;
                }
            }
        }
    } ;

    /*!
     * @brief Typed attribute existence check
     */
    template< typename T >
    bool is_attribute_defined(
        GEO::AttributesManager& manager,
        const std::string& attribute_name )
    {
        GEO::AttributeStore* store = manager.find_attribute_store( attribute_name ) ;
        if( store == nil ) {
            return false ;
        } else {
            std::string T_type_name( typeid(T).name() ) ;
            return store->elements_type_matches( T_type_name ) ;
        }
    }

    /*!
     * @brief Type sensitive check of Attribute existence on a Mesh facets
     */
    template< typename T >
    bool is_facet_attribute_defined(
        const GEO::Mesh& mesh,
        const std::string& attribute_name )
    {
        GEO::AttributesManager& manager = mesh.facets.attributes() ;
        return is_attribute_defined< T >( manager, attribute_name ) ;
    }

    /*!
     * @brief Type sensitive check of Attribute existence on a Mesh cells
     */
    template< typename T >
    bool is_cell_attribute_defined(
        const GEO::Mesh& mesh,
        const std::string& attribute_name )
    {
        GEO::AttributesManager& manager = mesh.cells.attributes() ;
        return is_attribute_defined< T >( manager, attribute_name ) ;
    }

    void RINGMESH_API print_bounded_attributes( const GEO::Mesh& M ) ;

}

#endif
