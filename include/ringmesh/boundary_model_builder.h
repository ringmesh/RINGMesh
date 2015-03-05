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
 *     Ecole Nationale Supérieure de Géologie - Georessources
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
    /**
     * @brief Base class for the class building a BoundaryModel.
     * @details Derive from this class
     */
    class RINGMESH_API BoundaryModelBuilder {
    public:
        typedef BoundaryModelElement BME ;

        BoundaryModelBuilder( BoundaryModel& model )
              : model_( model ) {}
        virtual ~BoundaryModelBuilder() {}

        // High level functions
        void copy_macro_topology( const BoundaryModel& from ) ;

        void make_vertices_unique() ;

        // Set model attributes
        void set_model_name( const std::string& name )
        {
            model_.name_ = name ;
        }

        void reserve_vertices( index_t nb )
        {
            model_.vertices_.reserve( nb ) ;
        }

        index_t add_vertex( const vec3& vertex )
        {
            model_.vertices_.push_back( vertex ) ;
            return model_.nb_vertices() - 1 ;
        }

        index_t add_vertex( const double* vertex )
        {
            return add_vertex( vec3( vertex[ 0 ], vertex[ 1 ], vertex[ 2 ] ) ) ;
        }

        /*!
         * @brief Generic accessor to a modifiable BoundaryModelElement
         */
        BoundaryModelElement& element(
            BME::TYPE t,
            index_t index )
        {
            return const_cast< BoundaryModelElement& >( model_.element( t, index ) ) ;
        }

        /** @}
         * \name Filling BoundaryModelElement attributes.
         * @{
         */
        void set_model(
            BME::TYPE e_type,
            index_t e_index,
            BoundaryModel* m )
        {
            element( e_type, e_index ).set_model( m ) ;
        }

        void set_element_index(
            BME::TYPE e_type,
            index_t e_index )
        {
            element( e_type, e_index ).set_id( e_index ) ;
        }

        void set_element_name(
            BME::TYPE e_type,
            index_t e_index,
            const std::string& name )
        {
            element( e_type, e_index ).set_name( name ) ;
        }

        void set_element_geol_feature(
            BME::TYPE e_type,
            index_t e_index,
            BME::GEOL_FEATURE geol )
        {
            element( e_type, e_index ).set_geological_feature( geol ) ;
        }

        void add_element_boundary(
            BME::TYPE e_type,
            index_t e_index,
            index_t boundary,
            bool side = false )
        {
            if( e_type == BoundaryModelElement::REGION || e_type ==
                BoundaryModelElement::LAYER )
            {
                element( e_type, e_index ).add_boundary( boundary, side ) ;
            } else { element( e_type, e_index ).add_boundary( boundary ) ;}
        }

        void add_element_in_boundary(
            BME::TYPE e_type,
            index_t e_index,
            index_t in_boundary )
        {
            element( e_type, e_index ).add_in_boundary( in_boundary ) ;
        }

        void set_parent(
            BME::TYPE e_type,
            index_t e_index,
            index_t parent_index )
        {
            element( e_type, e_index ).set_parent( parent_index ) ;
        }

        void add_child(
            BME::TYPE e_type,
            index_t e_index,
            index_t child_index )
        {
            element( e_type, e_index ).add_child( child_index ) ;
        }

        void set_element_vertex(
            BME::TYPE e_type,
            index_t e_index,
            index_t v,
            index_t model_v_id )
        {
            element( e_type, e_index ).set_vertex( v, model_v_id ) ;
        }

        /** @}
         * \name Find and/or create one BoundaryModelElement.
         * @{
         */
        index_t create_element( BME::TYPE e_type ) ;

        void erase_element(
            BME::TYPE type,
            index_t id ) ;

        void resize_elements(
            BME::TYPE e_type,
            index_t nb ) ;

        // Corner
        index_t find_corner( index_t ) const ;
        index_t create_corner( index_t ) ;
        index_t find_or_create_corner( index_t ) ;

        // Line
        index_t find_line( const std::vector< index_t >& vertices ) const ;

        index_t create_line( const std::vector< index_t >& vertices ) ;

        index_t find_or_create_line( const std::vector< index_t >& vertices ) ;

        // Surface
        index_t create_surface() ;

        // Contact
        index_t find_contact( const std::vector< index_t >& interfaces ) const ;

        index_t create_contact( const std::vector< index_t >& interfaces ) ;

        index_t find_or_create_contact( const std::vector< index_t >& interfaces ) ;

        // Interface
        index_t find_interface( const std::string& name ) const ;

        index_t create_interface(
            const std::string& name,
            BME::GEOL_FEATURE type = BME::NO_GEOL ) ;

        // Region
        index_t create_region() ;

        index_t create_region(
            const std::string& name,
            const std::vector< std::pair< index_t, bool > >& boundaries ) ;

        // Layers
        index_t create_layer( const std::string& name ) ;

        // Universe
        void set_universe( const std::vector< std::pair< index_t,
                                                         bool > >& boundaries ) ;

        /** @}
         * \name Set element geometry
         * @{
         */
        void set_corner(
            index_t corner_id,
            index_t vertex_id ) ;

        void set_line(
            index_t id,
            const std::vector< index_t >& vertices ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr,
            const std::vector< index_t >& surface_adjacencies = empty_index_vector ) ;

        void set_surface_geometry_bis(
            index_t surface_id,
            const std::vector< index_t >& corners,
            const std::vector< index_t >& facet_ptr,
            const std::vector< index_t >& corner_adjacent_facets =
                empty_index_vector ) ;

        void set_surface_adjacencies( index_t surface_id ) ;

        /** @}
         * \name Fix model - Check validity and fill missing stuff
         * @{
         */
        bool end_model() ;

        void update_all_ids() ;

        void init_global_model_element_access() ;

        bool complete_element_connectivity() ;

        bool check_basic_element_validity( const BoundaryModelElement& E ) const ;

        bool check_element_connectivity( const BoundaryModelElement& E ) const ;

        void fill_element_geological_feature( BoundaryModelElement& E ) ;

        void fill_elements_boundaries( BME::TYPE type ) ;

        void fill_elements_in_boundaries( BME::TYPE type ) ;

        void fill_elements_parent( BME::TYPE type ) ;

        void fill_elements_children( BME::TYPE type ) ;

        /**
         * @}
         */

    protected:
        BoundaryModel& model_ ;
    } ;

    /*!
     * \brief Build a BoundaryModel from a Gocad Model3D (file_model.ml)
     */
    class RINGMESH_API BoundaryModelBuilderGocad : public BoundaryModelBuilder {
    public:
        BoundaryModelBuilderGocad( BoundaryModel& model )
              : BoundaryModelBuilder( model ) {}
        virtual ~BoundaryModelBuilderGocad() {}

        void load_ml_file( const std::string& ml_file_name ) ;

        index_t determine_line_vertices(
            const Surface& S,
            index_t first_vertex,
            index_t second_vertex,
            std::vector< index_t >& border_vertex_model_ids ) const ;

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

        index_t find_corner( const vec3& ) const ;

        void build_contacts() ;

        void create_surface(
            const std::string& interface_name,
            const std::string& type,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2 ) ;

        /**
         * \brief Check if the surface triangle orientations match the one of the key facet
         */
        bool check_key_facet_orientation( index_t surface ) ;

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
     * \brief Build a BoundaryModel from a file_model.bm
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
     *  - nb_vertices()
     *  - double* vertex_ptr( index_t i )
     *  - index_t nb_corners()
     *  - index_t nb_facets()
     *  - index_t facet_begin( index_t f )
     *  - index_t facet_end( index_t f )
     *  - index_t corner_vertex_index( index_t c )
     *  - signed_index_t corner_adjacent_facet( index_t c )   -1 if no neighbor
     */
    template< class MESH >
    void BoundaryModelBuilderSurface::set_surfaces( const MESH& mesh )
    {
        // / 1. Copy the vertices of the input mesh to the model
        reserve_vertices( mesh.nb_vertices() ) ;
        for( index_t i = 0; i < mesh.nb_vertices(); i++ ) {
            add_vertex( mesh.vertex_ptr( i ) ) ;
        }

        // / 2. Propagate on the input mesh facet to determine its surface connected components
        // Vectors used in the loops
        std::vector< index_t > corners ;
        std::vector< index_t > facets_ptr ;

        corners.reserve( mesh.nb_corners() ) ;
        facets_ptr.reserve( mesh.nb_facets() ) ;

        std::vector< bool > visited( mesh.nb_facets(), false ) ;
        for( index_t i = 0; i < mesh.nb_facets(); i++ ) {
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

                    for( index_t c = mesh.facet_begin( f );
                         c < mesh.facet_end( f );
                         ++c )
                    {
                        corners.push_back( mesh.corner_vertex_index( c ) ) ;
                        index_t n = mesh.corner_adjacent_facet( c ) ;
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
