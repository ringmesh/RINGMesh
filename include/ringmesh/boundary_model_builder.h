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

#ifndef __RINGMESH_BOUNDARY_MODEL_BUILDER__
#define __RINGMESH_BOUNDARY_MODEL_BUILDER__

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
    public:
        BoundaryModelBuilder( BoundaryModel& model )
              : model_( model ) {}
        virtual ~BoundaryModelBuilder() {}

        /*! @}
        * \name Access - Modification of the BoundaryModel
        * @{
        */

        /*!
         *@brief Set the name of the model  
         */
        void set_model_name( const std::string& name )
        {
            model_.name_ = name ;
        }

        /*! 
         *@brief Copy elements and element connectivity of model @param from 
         */
        void copy_macro_topology( const BoundaryModel& from )
        {
            model_.copy_macro_topology( from ) ;
        }

        /*!
         *@brief The model under construction
         */
        const BoundaryModel& model() const 
        {
            return model_ ;
        }

        bool end_model() ;
        bool complete_element_connectivity() ;

        /*! @}
        * \name Creation - Deletion - Access to BoundaryModelElements .
        * @{
        */

        BME::bme_t create_element( BME::TYPE e_type ) ;

        /*!
        * @brief Set the element of the model to the given element.
        * @details It is on purpose that no checking whatsoever is performed.
        *          This way, nil pointers can be set for a following element removal.
        */
        void set_element( const BME::bme_t& id, BoundaryModelElement* E ) const
        {
            if( id.type < BME::NO_TYPE ) {
                model_.modifiable_elements( id.type )[ id.index ] = E ;
            } else {
                ringmesh_assert_not_reached ;
            }
        }      

        /*!
         * @brief Reference to a modifiable element of the model
         */
        BoundaryModelElement& element(
            const BME::bme_t& t ) const
        {
            return *element_ptr(t) ;
        }
        
        /*!
         * @brief Modifiable pointer to an element of the model
         */
        BoundaryModelElement* element_ptr( const BME::bme_t& id ) const
        {
            if( id.type < BME::NO_TYPE ) {
                return model_.elements( id.type )[ id.index ] ;
            } else if( id.type == BME::ALL_TYPES ) {
                return element_ptr( model_.global_to_typed_id( id ) ) ;
            } else {
                ringmesh_assert_not_reached ;
                return &model_.universe_ ;
            }            
        }

        bool get_dependent_elements( std::set< BME::bme_t >& elements ) const ;
        void remove_elements( const std::set< BME::bme_t >& elements ) ;

      
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
            if( t.type == BME::REGION ) {
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


        // Universe
        void set_universe( const std::vector<
                           std::pair< index_t, bool > >& boundaries ) ;
     
        /*! @}
         * \name Set element geometry from geometrical positions   
         * @{
         */
        void set_element_vertex(
            BME::bme_t t,
            index_t v,
            const vec3& point )
        {
            element( t ).set_vertex( v, point, false ) ;
        }

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
            const std::vector< index_t >& surface_facet_ptr ) ;


        /*! @}
        * \name Set element geometry using BoundaryModel vertices
        * @{
        */
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
            const std::vector< index_t >& surface_facet_ptr ) ;

        void set_surface_geometry(
            const BME::bme_t& surface_id,
            const std::vector< index_t >& corners,
            const std::vector< index_t >& facet_ptr ) ;

        void set_surface_adjacencies( const BME::bme_t& surface_id ) ;

        /*!
         * @}
         */

    protected:
        void delete_elements(
            std::vector< std::vector< index_t > >& to_erase ) ;
        void init_global_model_element_access() ;
        void resize_elements( BME::TYPE type, index_t nb ) ;

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
        
        bool load_ml_file( const std::string& ml_file_name ) ;

    protected:
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
        void build_contacts() ;

        void create_surface(
            const std::string& interface_name,
            const std::string& type,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2 ) ;

        /*!
         * @brief Triangle that set the orientation of a TFACE
         *        in a .ml file
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
     */
    class RINGMESH_API BoundaryModelBuilderSurface : public BoundaryModelBuilder {
    public:
        BoundaryModelBuilderSurface( BoundaryModel& model  ) :
            BoundaryModelBuilder( model ) {}
        virtual ~BoundaryModelBuilderSurface() {}

        template< class MESH >
        void set_surfaces( const MESH& mesh ) ;

        bool build_model( bool build_regions = true ) ;
    } ;


    /*!
     * @brief Create the model surfaces from the connected components of the input surfacic mesh
     * @details Add the separately the connected components of the mesh 
     *          as Surface of the model to create 
     *   
     *  The class MESH should implement the following functions
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
        // Propagate on the input mesh facet to determine its surface connected components
        // Vectors used in the loops
        std::vector< index_t > corners ;
        std::vector< index_t > facets_ptr ;
        std::vector< vec3 > vertices ;
        std::vector< index_t > cc_vertex( mesh.vertices.nb(), NO_ID ) ;

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
                vertices.resize( 0 ) ;
                cc_vertex.resize( mesh.vertices.nb(), NO_ID ) ;

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
                        index_t v = mesh.facet_corners.vertex( c ) ;
                        if( cc_vertex[ v ] == NO_ID ) {
                            cc_vertex[ v ] = vertices.size() ;
                            vertices.push_back( mesh.vertices.point( v ) ) ;
                        }
                        corners.push_back( cc_vertex[v] ) ;
                        index_t n = mesh.facet_corners.adjacent_facet( c ) ;
                        if( n != NO_ID && !visited[ n ] ) {
                            visited[ n ] = true ;
                            S.push( n ) ;
                        }
                    }
                    facets_ptr.push_back( corners.size() ) ;
                }

                // Create the surface and set its geometry
                set_surface_geometry( create_element( BME::SURFACE ), vertices, corners, facets_ptr ) ;
            }
        }
    }
}

#endif
