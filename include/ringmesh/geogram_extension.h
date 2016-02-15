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
*
*
*
*
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

#ifdef RINGMESH_WITH_TETGEN
#   include <geogram/third_party/tetgen/tetgen.h>
#endif 

namespace RINGMesh {

    /*!
     * Copy the content of a standrad library vector to the memory aligned GEO::Vector. 
     * A lot of copies, when we need to call Geogram functions. 
     * @todo Could we set Geogram vector to be a std::vector ?? 
     */
    template< class T >
    void copy_std_vector_to_geo_vector( const std::vector<T>& in, GEO::vector<T>& out )
    {
        out.resize( in.size() ) ;
        for( index_t i = 0; i < in.size(); ++i ) {
            out[ i ] = in[ i ]  ; 
        }
    }

    /*!
     * Partial copy the content of a standrad library vector to a GEO::Vector.
     * A lot of copies, when we need to call Geogram functions.
     * @todo Could we set Geogram vector to be a std::vector ??
     */
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

    /*! @brief complement the available MeshIOHandler
     */
    void RINGMESH_API ringmesh_mesh_io_initialize() ;
            
    /******************************************************************/
    /* Operations on a GEO::Mesh                                      */

#ifdef RINGMESH_WITH_TETGEN
    /// @todo Move all tetgen related stuff in one or two files

    /*! 
     * @brief Utility class to set Tetgen switches and check their consistency
     * @details Tetgen arguments are a mess and this class helps set the basic options
     * @todo To implement!
     *
     * Q: quiet
     * p: input data is surfacic
     * q: desired quality
     * O0: do not optimize mesh at all -> a lot of flat tets
     * V: verbose - A LOT of information
     * Y: prohibit steiner points on boundaries
     * A: generate region tags for each shell.      
     *
     * Meshing with incomplete quality value "Qpq%fYA"
     */
    class TetgenCommandLine {
    public:
        const std::string command_line() const {
            return command_line_ ;
        }

    private:
        std::string command_line_ ;
    };


    /*!
     * @brief Tetgen wrapper
     * 
     */
    class TetgenMesher {
        ringmesh_disable_copy( TetgenMesher ) ;
    public:
        TetgenMesher()
            : polygons_( nil ), polygon_corners_( nil )
        {
        }
        ~TetgenMesher() ;

        void tetrahedralize( const GEO::Mesh& input_mesh, 
                             const std::string& command_line, 
                             GEO::Mesh& output_mesh ) ; 
        
        void tetrahedralize( const GEO::Mesh& input_mesh, 
                             const std::vector< vec3 >& one_point_per_region,
                             const std::string& command_line, 
                             GEO::Mesh& output_mesh ) ; 

    private:
        void initialize() ;
        void initialize_tetgen_args() ;         
        void set_command_line( const std::string& command_line ) ;
        void tetrahedralize() ;

        void copy_mesh_to_tetgen_input( const GEO::Mesh& M ) ;
        void copy_vertices_to_tetgen_input( const GEO::Mesh& M ) ;
        void copy_edges_to_tetgen_input( const GEO::Mesh& M ) ;
        void copy_facets_to_tetgen_input( const GEO::Mesh& M ) ;
        void set_regions( const std::vector< vec3 >& one_point_per_region ) ;

        void fill_region_attribute_on_mesh_cells( GEO::Mesh& M, const std::string& attribute_name ) const ;
        void assign_result_tetmesh_to_mesh( GEO::Mesh& M ) const;
        void get_result_tetmesh_points( GEO::vector< double >& points ) const ;
        void get_result_tetmesh_tets( GEO::vector< index_t>& tets ) const ;

    private:
        GEO_3rdParty::tetgenio tetgen_in_ ;
        GEO_3rdParty::tetgenio tetgen_out_ ;
        std::string tetgen_command_line_ ;
        GEO_3rdParty::tetgenbehavior tetgen_args_ ;

        GEO_3rdParty::tetgenio::polygon* polygons_ ;
        int* polygon_corners_ ;        
    };
   


    /*!
     * @brief Constrained tetrahedralize of the volumes defined by a triangulated surface mesh
     * @details Does not require this mesh to be a closed manifold
     */
    void RINGMESH_API tetrahedralize_mesh_tetgen( GEO::Mesh& M, bool refine, double quality ) ;

#endif

    
    void RINGMESH_API rotate_mesh( GEO::Mesh& mesh, const GEO::Matrix< float64, 4 >& rot_mat ) ;

  
    double RINGMESH_API mesh_cell_volume( const GEO::Mesh& M, index_t c ) ;

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
    * @details This is a wrapper around Geogram colocate functions
    */
    bool RINGMESH_API has_mesh_colocate_vertices( const GEO::Mesh& M, double tolerance ) ;


    /*!
    * @brief Merges the vertices of a mesh that are at the same geometric location
    * @note Copied from geogram/mes/mesh_repair.cpp. No choice since BL will not give access to it.
    */
    void RINGMESH_API repair_colocate_vertices( GEO::Mesh& M, double tolerance ) ;


    /*!
    * @brief Vector of pointers to Geogram attributes
    * @note Necessary since one cannot create, vectors of Geogram attributes does 
    * not compile, because @#$# (no idea) [JP]
    * @todo Probably extremely prone to bugs. Is it worth the risk? 
    */
    template< class T >
    class AttributeVector : public std::vector< GEO::Attribute< T >* > {
        ringmesh_disable_copy( AttributeVector ) ;
    public:
        typedef std::vector< GEO::Attribute< T >* > base_class ;
        AttributeVector()
            : base_class()
        {}
        AttributeVector( index_t size )
            : base_class( size, nil )
        {}

        void bind_one_attribute( index_t i,
                                 GEO::AttributesManager& manager,
                                 const std::string& attribute_name )
        {
            base_class::operator[]( i ) = new GEO::Attribute< T >( manager, attribute_name ) ;
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
    template< class T > 
    bool is_attribute_defined( GEO::AttributesManager& manager, 
                               const std::string& attribute_name )
    {
        GEO::AttributeStore* store = manager.find_attribute_store( attribute_name ) ;
        if( store == nil ) {
            return false ;
        } else {
            std::string T_type_name( typeid( T ).name() );
            return store->elements_type_matches( T_type_name ) ;
        }
    }

    /*!
     * @brief Type sensitive check of Attribute existence on a Mesh facets
     */
    template< class T >
    bool is_facet_attribute_defined( const GEO::Mesh& mesh,
                                     const std::string& attribute_name )
    {
        GEO::AttributesManager& manager = mesh.facets.attributes() ;
        return is_attribute_defined< T >( manager, attribute_name ) ;
    }

    /*!
    * @brief Type sensitive check of Attribute existence on a Mesh cells
    */
    template< class T >
    bool is_cell_attribute_defined( const GEO::Mesh& mesh,
                                     const std::string& attribute_name )
    {
        GEO::AttributesManager& manager = mesh.cells.attributes() ;      
        return is_attribute_defined< T >( manager, attribute_name ) ;
    }


    void RINGMESH_API print_bounded_attributes( const GEO::Mesh& M ) ;

}


#endif
