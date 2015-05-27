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
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*! \author Jeanne Pellerin */

#ifndef __RINGMESH_BOUNDARY_MODEL_FROM_SURFACE__
#define __RINGMESH_BOUNDARY_MODEL_FROM_SURFACE__

#include <ringmesh/common.h>
#include <ringmesh/boundary_model.h>

#include <vector>
#include <string>
#include <stack>

namespace RINGMesh {
    /*!
     * @brief Base class for all classes building a BoundaryModel.
     * @details Derive from this class to build or modify a BoundaryModel
     */
    class RINGMESH_API BoundaryModelBuilder {
    protected:
        BoundaryModelBuilder( BoundaryModel& model )
              : model_( model ) {}
        virtual ~BoundaryModelBuilder() {}

        
        void set_model_name( const std::string& name )
        {
            model_.name_ = name ;
        }

        /*!
         * @brief Generic accessor to a modifiable BoundaryModelElement
         */
        BoundaryModelElement& element(
            const BME::bme_t& t )
        {
            return const_cast< BoundaryModelElement& >( model_.element( t ) ) ;
        }

        /*! @}
         * \name Filling BoundaryModelElement attributes.
         * @{
         */
        void set_model(
            const BME::bme_t& t,
            BoundaryModel* m )
        {
            element( t ).set_model( m ) ;
        }

        void set_element_index(
            const BME::bme_t& t )
        {
            element( t ).set_id( t.index ) ;
        }

        void set_element_name(
            const BME::bme_t& t,
            const std::string& name )
        {
            element( t ).set_name( name ) ;
        }

        void set_element_geol_feature(
            const BME::bme_t& t,
            BME::GEOL_FEATURE geol )
        {
            element( t ).set_geological_feature( geol ) ;
        }

        void add_element_boundary(
            const BME::bme_t& t,
            const BME::bme_t& boundary,
            bool side = false )
        {
            if( t.type == BoundaryModelElement::REGION || t.type ==
                BoundaryModelElement::LAYER )
            {
                element( t ).add_boundary( boundary, side ) ;
            } else { element( t ).add_boundary( boundary ) ;}
        }

        void add_element_in_boundary(
            const BME::bme_t& t,
            const BME::bme_t& in_boundary )
        {
            element( t ).add_in_boundary( in_boundary ) ;
        }

        void set_parent(
            const BME::bme_t& t,
            const BME::bme_t& parent_index )
        {
            element( t ).set_parent( parent_index ) ;
        }

        void add_child(
            const BME::bme_t& t,
            const BME::bme_t& child_index )
        {
            element( t ).add_child( child_index ) ;
        }

        void set_element_vertex(
            BME::bme_t t,
            index_t v,
            const vec3& point )
        {
            element( t ).set_vertex( v, point, false ) ;
        }

        /*! @}
         * \name Find and/or create one BoundaryModelElement.
         * @{
         */
        BME::bme_t create_element( BME::TYPE e_type ) ;

        void erase_element( const BME::bme_t& t ) ;

        void resize_elements(
            BME::TYPE e_type,
            index_t nb ) ;

        // Corner
        BME::bme_t find_corner( const vec3& point) const ;
        BME::bme_t find_corner( index_t model_point_id ) const ;
        BME::bme_t create_corner( const vec3& point ) ;
        BME::bme_t find_or_create_corner( const vec3& point) ;

        // Line
        BME::bme_t find_line( const std::vector< vec3 >& vertices ) const ;

        BME::bme_t create_line( const std::vector< vec3 >& vertices ) ;

        BME::bme_t find_or_create_line( const std::vector< vec3 >& vertices ) ;

        // Surface
        BME::bme_t create_surface() ;

        // Contact
        BME::bme_t find_contact( const std::vector< index_t >& interfaces ) const ;

        BME::bme_t create_contact( const std::vector< index_t >& interfaces ) ;

        BME::bme_t find_or_create_contact( const std::vector< index_t >& interfaces ) ;

        // Interface
        BME::bme_t find_interface( const std::string& name ) const ;

        BME::bme_t create_interface(
            const std::string& name,
            BME::GEOL_FEATURE type = BME::NO_GEOL ) ;

        // Region
        BME::bme_t create_region() ;

        BME::bme_t create_region(
            const std::string& name,
            const std::vector< std::pair< index_t, bool > >& boundaries ) ;

        // Layers
        BME::bme_t create_layer( const std::string& name ) ;

        // Universe
        void set_universe( const std::vector< std::pair< index_t,
                                                         bool > >& boundaries ) ;

        /*! @}
         * \name Set element geometry
         * @{
         */
        void set_corner(
            const BME::bme_t& corner_id,
            const vec3& point ) ;

        void set_line(
            const BME::bme_t& id,
            const std::vector< vec3 >& vertices ) ;
                
        void set_surface_geometry(
            const BME::bme_t& surface_id,
            const std::vector< vec3 >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr,
            const std::vector< index_t >& surface_adjacencies = std::vector< index_t >() ) ;

        // Same functions but to call in the case where the vertices of the 
        // model are filled first 
        index_t add_unique_vertex( const vec3& p ) ;

        void set_corner(
            const BME::bme_t& corner_id,
            index_t unique_vertex ) ;

        void set_line(
            const BME::bme_t& id,
            const std::vector< index_t >& unique_vertices ) ;

        void set_surface_geometry(
            const BME::bme_t& surface_id,
            const std::vector< index_t >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr,
            const std::vector< index_t >& surface_adjacencies = std::vector< index_t >() ) ;

        void set_surface_geometry_bis(
            const BME::bme_t& surface_id,
            const std::vector< index_t >& corners,
            const std::vector< index_t >& facet_ptr,
            const std::vector< index_t >& corner_adjacent_facets = std::vector< index_t >() ) ;

        void set_surface_adjacencies( const BME::bme_t& surface_id ) ;

        /*! @}
         * \name Fix model - Check validity and fill missing stuff
         * @{
         */
        bool end_model() ;

        void update_all_ids() ;

        void init_global_model_element_access() ;

        bool complete_element_connectivity() ;  

        void fill_elements_boundaries( BME::TYPE type ) ;

        void fill_elements_in_boundaries( BME::TYPE type ) ;

        void fill_elements_parent( BME::TYPE type ) ;

        void fill_elements_children( BME::TYPE type ) ;

        /*!
         * @}
         */

    protected:
        BoundaryModel& model_ ;
    } ;

    /*!
     * @brief Build a BoundaryModel from a Gocad Model3D (file_model.ml)
     */
    class RINGMESH_API BoundaryModelBuilderGocad : public BoundaryModelBuilder {
    public:
        BoundaryModelBuilderGocad( BoundaryModel& model )
              : BoundaryModelBuilder( model ) {}
        virtual ~BoundaryModelBuilderGocad() {}

        void load_ml_file( const std::string& ml_file_name ) ;

        BME::bme_t determine_line_vertices(
            const Surface& S,
            index_t id0,
            index_t id1,
            std::vector< index_t >& border_vertex_model_ids ) const ;

        BME::bme_t determine_line_vertices(
            const Surface& S,
            index_t first_vertex,
            index_t second_vertex,
            std::vector< vec3 >& border_vertex_model_vertices ) const ;

    private:
        /*!
         * @brief Triangle that set the orientation of a TFACE
         */
        struct KeyFacet {
            KeyFacet(
                const vec3& p0,
                const vec3& p1,
                const vec3& p2 ) :
                p0_( p0 ), p1_( p1 ), p2_( p2 ) {}

        public:
            vec3 p0_ ;
            vec3 p1_ ;
            vec3 p2_ ;
        } ;

        void build_contacts() ;

        void create_surface(
            const std::string& interface_name,
            const std::string& type,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2 ) ;

        /*!
         * @brief Check if the surface triangle orientations match the one of the key facet
         */
        bool check_key_facet_orientation( index_t surface ) const;

        index_t find_key_facet(
            index_t surface_id,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            bool& same_orientation ) const ;

    private:
        std::vector< KeyFacet > key_facets_ ;
    } ;

    /*!
     * @brief Build a BoundaryModel from a file_model.bm
     */
    class RINGMESH_API BoundaryModelBuilderBM : public BoundaryModelBuilder {
    public:
        BoundaryModelBuilderBM( BoundaryModel& model )
              : BoundaryModelBuilder( model ) {}
        virtual ~BoundaryModelBuilderBM() {}

        bool load_file( const std::string& bm_file_name ) ;

    private:
        static BME::TYPE match_nb_elements( const char* s ) ;

        static BME::TYPE match_type( const char* s ) ;

        static bool match_high_level_type( const char* s )
        {
            return BME::child_allowed( match_type( s ) ) ;
        }
    } ;

    /*!
     * @brief Builder of a BoundaryModel from a conformal surface meshes
     *        in which the manifold connected components are disjoints
     * @todo A TESTER 
     */
    class RINGMESH_API BoundaryModelBuilderSurface : public BoundaryModelBuilder {
    public:
        BoundaryModelBuilderSurface( BoundaryModel& model  ) :
            BoundaryModelBuilder( model ) {}
        virtual ~BoundaryModelBuilderSurface() {}

        template< class MESH >
        void set_surfaces( const MESH& mesh ) ;

        void build_model() ;
    } ;


    /*!
     * @brief Create the model surfaces from the connected components of the input surfacic mesh
     * @details The class MESH should implement the following functions
     *  - vertices.nb()
     *  - const vec3& point( index_t i )
     *  - index_t facet_corners.nb()
     *  - index_t facets.nb()
     *  - index_t facets.corners_begin( index_t f )
     *  - index_t facets.corners_end( index_t f )
     *  - index_t facet_corners.vertex( index_t c )
     *  - signed_index_t facet_corners.adjacent_facet( index_t c )   -1 if no neighbor
     */
    template< class MESH >
    void BoundaryModelBuilderSurface::set_surfaces( const MESH& mesh )
    {
        /// 1. Copy the vertices of the input mesh to the model
        // reserve_vertices( mesh.nb_vertices() ) ;
        for( index_t i = 0; i < mesh.vertices.nb(); i++ ) {
            add_unique_vertex( mesh.vertices.point( i ) ) ;
        }

        /// 2. Propagate on the input mesh facet to determine its surface connected components
        // Vectors used in the loops
        std::vector< index_t > corners ;
        std::vector< index_t > facets_ptr ;

        corners.reserve( mesh.facet_corners.nb() ) ;
        facets_ptr.reserve( mesh.facets.nb() ) ;

        std::vector< bool > visited( mesh.facets.nb(), false ) ;
        for( index_t i = 0; i < mesh.facets.nb(); i++ ) {
            if( !visited[ i ] ) {
                // Index of the Surface to create form this facet
                index_t cc_index = model_.nb_surfaces() ;

                // Get the facets that are in the same connected component than the current facet
                corners.resize( 0 ) ;
                facets_ptr.resize( 0 ) ;
                facets_ptr.push_back( 0 ) ;

                std::stack< index_t > S ;
                S.push( i ) ;
                while( !S.empty() ) {
                    index_t f = S.top() ;
                    S.pop() ;
                    visited[ f ] = true ;

                    for( index_t c = mesh.facets.corners_begin( f );
                         c < mesh.facets.corners_end( f );
                         ++c )
                    {
                        corners.push_back( mesh.facet_corners.vertex( c ) ) ;
                        index_t n = mesh.facet_corners.adjacent_facet( c ) ;
                        if( n != NO_ID && !visited[ n ] ) {
                            visited[ n ] = true ;
                            S.push( n ) ;
                        }
                    }
                    facets_ptr.push_back( corners.size() ) ;
                }

                // Create the surface and set its geometry - adjacencies are computed
                set_surface_geometry_bis( create_surface(), corners, facets_ptr ) ;
            }
        }
    }
}

#endif
