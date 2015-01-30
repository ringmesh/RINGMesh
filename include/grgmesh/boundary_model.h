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
*/

/*! \author Jeanne Pellerin and Arnaud Botella */


#ifndef __GRGMESH_BOUNDARY_MODEL__
#define __GRGMESH_BOUNDARY_MODEL__

#include <grgmesh/common.h>
#include <grgmesh/boundary_model_element.h>
#include <grgmesh/attribute.h>

#include <vector> 
#include <string>

namespace GRGMesh {    
    class BoundaryModelBuilder ;
}

namespace GRGMesh {

    // To move somewhere else
    static std::vector< vec3 > empty_vector ;
    static std::vector< index_t > empty_index_vector ;

    /**
     * \brief The class to describe a volumetric model represented by its boundary surfaces
     *     
     */
    class GRGMESH_API BoundaryModel {       
        friend class BoundaryModelBuilder ;

    public:           
        enum AttributeLocation {
            VERTEX,
            EDGE,
            FACET
        } ;       
        typedef AttributeManager< VERTEX > PointAttributeManager ;
        typedef AttributeManager< FACET > FacetAttributeManager ;
        typedef BoundaryModelElement BME ;
               
        const static index_t NO_ID = index_t( -1 ) ;

        /**
         * \brief Construct an empty BoundaryModel
         */
        BoundaryModel() { } ;
        /**
         * \brief Destroy a BoundaryModel
         */
        virtual ~BoundaryModel(){} ;
        void clear() ;
        const std::string& name() const { return name_ ; }

        /**
         * \name Global access to model vertices and facets
         * @{
         */
        index_t nb_vertices() const { return vertices_.size() ; }        
        const vec3& vertex( index_t p ) const { return vertices_.at(p) ; }
        index_t vertex_index( const vec3& p ) const ;

        index_t nb_facets() const { return nb_facets_in_surfaces_.back() ; }
        void surface_facet( index_t model_facet_id, index_t& surface_id, index_t& surf_facet_id ) const ;     
        index_t model_facet( index_t surface_id, index_t surf_facet_id ) const ;      

        /**
         * \name Accessor to elements
         * @{
         */

        /*!
        * @brief Returns the number of elements of the given type
        * By default returns 0.
        */ 
        inline index_t nb_elements( BME::TYPE type ) const 
        {
             switch( type ){
                case BoundaryModelElement::CORNER    : return corners_.size() ;
                case BoundaryModelElement::LINE      : return lines_.size() ;
                case BoundaryModelElement::SURFACE   : return surfaces_.size() ;
                case BoundaryModelElement::REGION    : return regions_.size() ;
                case BoundaryModelElement::CONTACT   : return contacts_.size() ;
                case BoundaryModelElement::INTERFACE : return interfaces_.size() ;
                case BoundaryModelElement::LAYER     : return layers_.size() ;            
                case BoundaryModelElement::ALL_TYPES : 
                    grgmesh_assert( nb_elements_per_type_.size() > 0 ) ;
                    grgmesh_debug_assert( nb_elements_per_type_.back() == 
                            corners_.size() + lines_.size() + surfaces_.size() + regions_.size() +
                            contacts_.size() + interfaces_.size() + layers_.size() ) ;
                    return nb_elements_per_type_.back() ;                
                default:                
                    return 0 ;
            }
        }
        /*!
         * \brief Returns a const reference the identified BoundaryModelElement   
         *
         * @param[in] t Type of the element 
         * @param[in] index Index of the element
         */
        inline const BoundaryModelElement& 
            element( BME::TYPE type, index_t index ) const 
        {
            grgmesh_assert( index < nb_elements( type ) ) ;
            switch( type ){
                case BoundaryModelElement::CORNER    : return corners_   [ index ] ;
                case BoundaryModelElement::LINE      : return lines_     [ index ] ;
                case BoundaryModelElement::SURFACE   : return surfaces_  [ index ] ;
                case BoundaryModelElement::REGION    : return regions_   [ index ] ;
                case BoundaryModelElement::CONTACT   : return contacts_  [ index ] ;
                case BoundaryModelElement::INTERFACE : return interfaces_[ index ] ;
                case BoundaryModelElement::LAYER     : return layers_    [ index ] ;
                case BoundaryModelElement::ALL_TYPES :
                    {
                        // This must synchro with what is done in the builder
                        index_t t = NO_ID ;
                        for( index_t i = 1; i < nb_elements_per_type_.size(); i++ ) {
                            if( index >= nb_elements_per_type_[i-1] && index < nb_elements_per_type_[i] ) {
                                t = i-1 ; break ;
                            }
                        }
                        grgmesh_assert( t < BME::NO_TYPE ) ;
                        return element( (BME::TYPE) t, index - nb_elements_per_type_[t] ) ;
                    }
                default:
                    grgmesh_assert_not_reached ;
                    return dummy_element ;
            }
        }

        /// Remove these functions ? 
        index_t nb_corners()    const { return nb_elements( BME::CORNER    ) ; }
        index_t nb_lines()      const { return nb_elements( BME::LINE      ) ; }
        index_t nb_surfaces()   const { return nb_elements( BME::SURFACE   ) ; }
        index_t nb_regions()    const { return nb_elements( BME::REGION    ) ; }
        index_t nb_contacts()   const { return nb_elements( BME::CONTACT   ) ; }
        index_t nb_interfaces() const { return nb_elements( BME::INTERFACE ) ; }
        index_t nb_layers()     const { return nb_elements( BME::LAYER     ) ; }

        const Corner& corner( index_t index ) const { return corners_.at(index) ; }
        const Line& line( index_t index ) const { return lines_.at(index) ; }
        const Surface& surface( index_t index ) const { return surfaces_.at(index) ; }

        const BoundaryModelElement& region       ( index_t index ) const { return element( BME::REGION   , index ) ; }
        const BoundaryModelElement& contact      ( index_t index ) const { return element( BME::CONTACT  , index ) ; }
        const BoundaryModelElement& one_interface( index_t index ) const { return element( BME::INTERFACE, index ) ; }
        const BoundaryModelElement& layer        ( index_t index ) const { return element( BME::LAYER    , index ) ; }

        const BoundaryModelElement& universe() const { return universe_ ; }        
       
        PointAttributeManager* vertex_attribute_manager() const
        {
            return const_cast< PointAttributeManager* >( &vertex_attribute_manager_ ) ;
        }
        FacetAttributeManager* facet_attribute_manager() const
        {
            return const_cast< FacetAttributeManager* >( &facet_attribute_manager_ ) ;
        }
         
        index_t find_region( index_t surface_part_id, bool side ) const ;
                
        /// \todo Write a proper IO class for Boundary models
        bool save_gocad_model3d( std::ostream& out ) ;
        void save_as_eobj_file( const std::string& file_name ) ;

        void set_vertex_coordinates( index_t id, const vec3& p ) {
            grgmesh_assert( id < nb_vertices() ) ;
            vertices_[id] = p ;
        }

    private:        
        bool load_gocad_model3d( const std::string& in ) ;

        bool check_model3d_compatibility() ;
        static void save_type( std::ostream& out, BoundaryModelElement::GEOL_FEATURE t ) ;   
        
    private:
        std::string name_ ;

        /** 
         * \brief Coordinates of the vertices of the model elements
         * Storage of vertices is unique for the whole model.
         */
        std::vector< vec3 >                 vertices_ ;

        // Base manifold elements of a model
        std::vector< Corner >               corners_ ;
        std::vector< Line >                 lines_ ;
        std::vector< Surface >              surfaces_ ;
        std::vector< BoundaryModelElement > regions_ ;

        /// The region including all the other regions
        BoundaryModelElement universe_ ;

        /// Sum of the number of facets in all previous surfaces
        /// Must be updated when a Surface is modified !!
        // Size = nb_surface()+1
        std::vector< index_t > nb_facets_in_surfaces_ ;
    
        /** 
         * \brief Contacts between Intefaces
         * Parent of a set of Line
         */
        std::vector< BoundaryModelElement >  contacts_ ;
        /** 
         * \brief Interfaces between layers 
         * Parent of a set of Surface
         */
        std::vector< BoundaryModelElement >  interfaces_ ;

        /** 
         * \brief Rock layers 
         * Parent of a set of Region
         */
        std::vector< BoundaryModelElement >  layers_ ;

        // For a global access to any of the BME
        // MUST be updated if one element is added !!!
        std::vector< index_t > nb_elements_per_type_ ;


        // Attribute managers 
        PointAttributeManager vertex_attribute_manager_ ;
        FacetAttributeManager facet_attribute_manager_ ;
    } ;   


    template< class ATTRIBUTE >
    class BoundaryModelVertexAttribute: public Attribute< BoundaryModel::VERTEX, ATTRIBUTE > {
    public:
        typedef Attribute< BoundaryModel::VERTEX, ATTRIBUTE > superclass ;

        void bind( const BoundaryModel* model, const std::string& name )
        {
            superclass::bind( model->vertex_attribute_manager(), model->nb_vertices(),
                name ) ;
        }

        void bind( const BoundaryModel* model )
        {
            superclass::bind( model->vertex_attribute_manager(),
                model->nb_vertices() ) ;
        }

        BoundaryModelVertexAttribute()
        {
        }

        BoundaryModelVertexAttribute( const BoundaryModel* model )
        {
            bind( model ) ;
        }

        BoundaryModelVertexAttribute( const BoundaryModel* model, const std::string& name )
        {
            bind( model, name ) ;
        }

        static bool is_defined( const BoundaryModel* model, const std::string& name )
        {
            return superclass::is_defined( model->vertex_attribute_manager(), name ) ;
        }
    } ;

    template< class ATTRIBUTE >
    class BoundaryModelFacetAttribute: public Attribute< BoundaryModel::FACET, ATTRIBUTE > {
    public:
        typedef Attribute< BoundaryModel::FACET, ATTRIBUTE > superclass ;

        void bind( const BoundaryModel* model, const std::string& name )
        {
            superclass::bind( model->facet_attribute_manager(), model->nb_facets(),
                name ) ;
        }

        void bind( const BoundaryModel* model )
        {
            superclass::bind( model->facet_attribute_manager(),
                model->nb_facets() ) ;
        }

        BoundaryModelFacetAttribute()
        {
        }

        BoundaryModelFacetAttribute( const BoundaryModel* model )
        {
            bind( model ) ;
        }

        BoundaryModelFacetAttribute( const BoundaryModel* model, const std::string& name )
        {
            bind( model, name ) ;
        }

        static bool is_defined( const BoundaryModel* model, const std::string& name )
        {
            return superclass::is_defined( model->facet_attribute_manager(), name ) ;
        }
    } ;


    /**
    * \brief Build a BoundaryModel
    */ 
    class GRGMESH_API BoundaryModelBuilder {        
    public:
        typedef BoundaryModelElement BME ;

        const static index_t NO_ID = index_t( -1 ) ;

        BoundaryModelBuilder( BoundaryModel& model )
            : model_( model ){}
        virtual ~BoundaryModelBuilder(){} ;

        // High level functions
        bool rebuild() ;
        void copy_macro_topology( const BoundaryModel& from ) ;        
        void update_all_ids() ;
        void make_vertices_unique() ;        

        // Set model attributes
        void set_model_name( const std::string& name ) {
            model_.name_ = name ;
        }
        void reserve_vertices( index_t nb ) {
            model_.vertices_.reserve( nb ) ; 
        }
        index_t add_vertex( const vec3& vertex ) {            
            model_.vertices_.push_back( vertex ) ;
            return model_.nb_vertices()-1 ;
        }
        index_t add_vertex( double* vertex ) {
            return add_vertex( vec3( vertex[0], vertex[1], vertex[2] ) ) ;
        }  

        BoundaryModelElement& element( BME::TYPE t, index_t index ) {
            return const_cast< BoundaryModelElement& >( model_.element( t, index ) ) ;
        }

        /**
        * \name Filling BoundaryModelElement attributes.
        * @{
        */
        void set_model( BME::TYPE e_type, index_t e_index, BoundaryModel* m ) {
            element( e_type, e_index ).set_model( m ) ;
        }
        void set_element_name( BME::TYPE e_type, index_t e_index, const std::string& name ) {
            element( e_type, e_index ).set_name( name ) ;
        }      
        void set_element_geol_feature( BME::TYPE e_type, index_t e_index, BME::GEOL_FEATURE geol ) {
            element( e_type, e_index ).set_geological_feature( geol ) ; 
        }
        void add_element_boundary( BME::TYPE e_type, index_t e_index, index_t boundary, bool side = false ){
            if( e_type == BoundaryModelElement::REGION || e_type == BoundaryModelElement::LAYER ) 
                element( e_type, e_index ).add_boundary( boundary, side ) ;
            else element( e_type, e_index ).add_boundary( boundary ) ;
        }
        void add_element_in_boundary( BME::TYPE e_type, index_t e_index, index_t in_boundary ) {
            element( e_type, e_index ).add_in_boundary( in_boundary ) ;
        }
        void set_parent( BME::TYPE e_type, index_t e_index, index_t parent_index ) {
            element( e_type, e_index ).set_parent( parent_index ) ;
        }
        void add_child( BME::TYPE e_type, index_t e_index, index_t child_index ) {
            element( e_type, e_index ).add_child( child_index ) ;
        }

        /**
        * \name Find and/or create one BoundaryModelElement.
        * @{
        */                     
        index_t create_element( BME::TYPE e_type ) ;

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
        index_t create_interface( const std::string& name, BME::GEOL_FEATURE type = BME::NO_GEOL ) ;

        // Region           
        index_t create_region() ;        
        index_t create_region( 
            const std::string& name,
            const std::vector< std::pair< index_t, bool > >& boundaries ) ;           

        // Layers
        index_t create_layer( const std::string& name ) ;

        // Universe
        void set_universe( const std::vector< std::pair< index_t, bool > >& boundaries ) ;        
        void remove_universe_from_regions( index_t id ) ;

        /**
        * \name Set element geometry 
        * @{
        */
        void set_corner( index_t corner_id, index_t vertex_id ) ;
        void set_line( index_t id, const std::vector< index_t >& vertices ) ;
        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr,
            const std::vector< index_t >& surface_adjacencies = empty_index_vector ) ;

        void set_surface_adjacencies( index_t surface_id ) ;

        void set_surface_first_triangle_as_key( index_t id ) {
            model_.surfaces_[id].set_first_triangle_as_key() ;
        }
        void set_surface_key_facet( index_t id, const Surface::KeyFacet& key ) {
            model_.surfaces_[id].set_key_facet( key ) ;
        } 

        /**
        * \name Fix model - Check validity und fill missing stuff
        * @{
        */
        bool end_model() ;

        bool complete_element_connectivity() ;
        bool check_basic_element_validity( const BoundaryModelElement& E ) const ;        
        void fill_elements_boundaries   ( BME::TYPE type ) ;
        void fill_elements_in_boundaries( BME::TYPE type ) ;
        void fill_elements_parent       ( BME::TYPE type ) ;
        void fill_elements_children     ( BME::TYPE type ) ;      

    protected:
        BoundaryModel& model_ ;
    };



    
    /*!
     * \brief Build a BoundaryModel from a Gocad Model3D (file_model.ml)
     */ 
    class GRGMESH_API BoundaryModelBuilderGocad : public BoundaryModelBuilder {
    public :
        /**
         * \brief Structure used to build contacts when loading a BoundaryModel from .ml file 
         */
        struct Border {
            Border( index_t part, index_t corner, index_t p0, index_t p1):
            part_id_(part), corner_id_(corner), p0_(p0), p1_(p1) {};

            // Id of the Surface owning this Border
            index_t part_id_ ;
            // Id of p0 in the BoundaryModel corner vector
            index_t corner_id_ ;

            // Ids of the starting corner and second vertex on the border in the Surface
            // to which this Border belong
            index_t p0_ ;
            index_t p1_ ;
        } ;


        BoundaryModelBuilderGocad( BoundaryModel& model )
            : BoundaryModelBuilder( model ){}

        void load_ml_file( std::istream& in ) ;   
        
        index_t determine_line_vertices( 
            const Surface& S, 
            index_t first_vertex, 
            index_t second_vertex,            
            std::vector< index_t >& border_vertex_model_ids ) const ;

    private:
        // Geometrical research of points
        index_t find_corner( const vec3& ) const ;

        void build_lines( const std::vector< Border >& borders ) ;
        

        void create_surface(
            const std::string& interface_name = "",
            const std::string& type = "",
            const Surface::KeyFacet& key = Surface::KeyFacet() ) ;
        

        void end_surfaces( const std::vector< index_t >& change_orientation ) ;

        void build_contacts() ;
        /**
         * \brief Check if the surface triangle orientations match the one of the key facet 
         */
        bool check_key_facet_orientation( index_t surface ) ;

        index_t find_key_facet( index_t surface_id, const vec3& p0, const vec3& p1, const vec3& p2, 
            bool& same_orientation ) const ;  
        

    } ;
}

#endif
