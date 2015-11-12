/*
* Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*  Contacts:
*     Arnaud.Botella@univ-lorraine.fr
*     Antoine.Mazuyer@univ-lorraine.fr
*     Jeanne.Pellerin@wias-berlin.de
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

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {

    /*!
     * Copy the content of a standrad library vector to the
     * memory aligned GEO::Vector. 
     * A lot of copies, when we need to call Geogram functions. 
     * @todo Could we set Geogram vector to be a std::vector ?? 
     */
    template< class T >
    void copy_std_vector_to_geo_vector( const std::vector<T>& in, GEO::vector<T>& out )
    {
        out.resize( in.size() ) ;
        for( index_t i = 0; i < in.size(); ++i) {
            out[ i ] = in[ i ]  ; 
        }
    }

    template< class T >
    void copy_std_vector_to_geo_vector( 
        const std::vector<T>& in, index_t from, index_t to, GEO::vector<T>& out )
    {
        ringmesh_assert( to < in.size()+1 ) ;
        ringmesh_assert( from < to ) ;
        index_t nb_to_copy( to - from ) ;
        out.resize( nb_to_copy ) ;
        index_t count = 0 ;
        for( index_t i = 0; i != nb_to_copy; ++i) {
            out[ i ] = in[ from +i ] ;
        }
    }

	
    /***********************************************************************/
    /* Loading and saving a GEO::Mesh                                      */

    /** @brief complement the available MeshIOHandler
     */
    void RINGMESH_API ringmesh_mesh_io_initialize() ;
            
    /******************************************************************/
    /* Operations on a GEO::Mesh                                      */

    
    /*!
     * @brief Constrained tetrahedralize of the volumes defined by a triangulated surface mesh
     * @details Does not require this mesh to be a closed manifold
     * as the equivalent Geogram function does.
     */
    bool RINGMESH_API tetrahedralize_mesh_tetgen( GEO::Mesh& M ) ;

    
    void RINGMESH_API rotate_mesh(
        GEO::Mesh& mesh,
        const GEO::Matrix< float64, 4 >& rot_mat ) ;

  
    double RINGMESH_API mesh_cell_volume( const GEO::Mesh& M, index_t c ) ;

    vec3 RINGMESH_API mesh_cell_facet_normal(
        const GEO::Mesh& M,
        index_t c,
        index_t f ) ;

    vec3 RINGMESH_API mesh_cell_facet_center(
        const GEO::Mesh& M,
        index_t cell,
        index_t f ) ;

    vec3 RINGMESH_API mesh_cell_center( const GEO::Mesh& M, index_t cell ) ;

    bool RINGMESH_API has_edge(
        const GEO::Mesh& mesh,
        index_t t,
        index_t p0,
        index_t p1,
        index_t& edge ) ;

    index_t RINGMESH_API next_around_edge(
        const GEO::Mesh& mesh,
        index_t t,
        index_t prev,
        index_t p0,
        index_t p1 ) ;

    void RINGMESH_API edges_around_edge(
        const GEO::Mesh& mesh,
        index_t t,
        index_t p0,
        index_t p1,
        std::vector< index_t >& result ) ;

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
  

    /*!
    * @brief Returns true if there are colocated vertices in the Mesh
    * @details This is a wrapper around Geogram Colocate functions
    */
    bool RINGMESH_API has_mesh_colocate_vertices( const GEO::Mesh& M, double tolerance ) ;


    /*!
    * @brief Merges the vertices of a mesh that are at the same
    *  geometric location
    * @note Copied from geogram/mes/mesh_repair.cpp. No choice since BL
    * will not give access to it.
    */
    void RINGMESH_API repair_colocate_vertices( GEO::Mesh& M, double tolerance ) ;


    /*!
    * \brief Convenient class to manipulate vectors of Geogram attributes.
    * \details Used to ease the storage of a common attribute on several
    * meshes grouped in the same object, for example those stored by a MacroMesh.
    */
    template< class T >
    class AttributeHandler : public std::vector< GEO::Attribute< T >* > {
        ringmesh_disable_copy( AttributeHandler ) ;
    public:
        typedef std::vector< GEO::Attribute< T >* > base_class ;
        AttributeHandler()
            : base_class()
        {}
        AttributeHandler( index_t size )
            : base_class( size, nil )
        {}

        /*!
        * Allocate one attribute on one component of the vector
        * @param[in] m id of the GEO::Mesh
        * @param[in] name name of the attribute
        * @param[in] am attribute manager, saying where the attribute is (cells, facets...)
        */
        void allocate_attribute(
            const index_t m,
            GEO::AttributesManager& am,
            const std::string& name )
        {
            ringmesh_debug_assert( m < base_class::size() ) ;
            ringmesh_debug_assert( !base_class::operator[]( m ) ) ;
            base_class::operator[]( m ) = new GEO::Attribute< T >( am, name ) ;
        }

        /*!
        * Allocate one vector of attributes on one component of the vector
        * @param[in] m id of the GEO::Mesh
        * @param[in] name name of the attribute
        * @param[in] am attribute manager, saying where the attribute is (cells, facets...)
        * @param[in] size size of the vector of attributes
        */
        void allocate_attribute(
            const index_t m,
            GEO::AttributesManager& am,
            const std::string& name,
            index_t size )
        {
            ringmesh_debug_assert( m < base_class::size() ) ;
            ringmesh_debug_assert( !base_class::operator[]( m ) ) ;
            base_class::operator[]( m ) = new GEO::Attribute< T >() ;
            base_class::operator[]( m )->create_vector_attribute( am, name, size ) ;
        }

        GEO::Attribute< T >& operator[]( index_t i )
        {
            return *base_class::operator[]( i ) ;
        }

        const GEO::Attribute< T >& operator[]( index_t i ) const
        {
            return *base_class::operator[]( i ) ;
        }

        ~AttributeHandler()
        {
            for( index_t i = 0; i < base_class::size(); i++ ) {
                if( base_class::operator[]( i ) )
                    delete base_class::operator[]( i ) ;
            }
        }
    } ;


    void RINGMESH_API print_bounded_attributes( const GEO::Mesh& M ) ;

}


#endif
