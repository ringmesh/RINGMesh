/*[
* NOT ONLY ASGA 
*/
/*! \author Jeanne Pellerin and Arnaud Botella */


#ifndef __GRGMESH_BOUNDARY_MODEL__
#define __GRGMESH_BOUNDARY_MODEL__

#include <grgmesh/common.h>
#include <grgmesh/boundary_model_element.h>
#include <grgmesh/attribute.h>

#include <vector> 
#include <string>

#include <vector> 
#include <map>
#include <string>
#include <algorithm> 


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
        friend class LineMutator ;
        friend class SurfaceMutator ;

    public:
         enum AttributeLocation {
            POINT,
            EDGE,
            FACET
        } ;       
        typedef AttributeManager< POINT > PointAttributeManager ;
        typedef AttributeManager< EDGE >  EdgeAttributeManager ; // Edges of what ?? lines
        typedef AttributeManager< FACET > FacetAttributeManager ;
        
        /**
         * \brief Construct an empty BoundaryModel
         */
        BoundaryModel() ;
        /**
         * \brief Destroy a BoundaryModel
         */
        virtual ~BoundaryModel() ;

        // Accessors to model points - edges - facets
        index_t nb_points() const { return points_.size() ; }
        const vec3& point( index_t p ) const { return points_.at(p) ; }
        index_t point_index( const vec3& p ) const ;

        inline index_t nb_facets() const ;
        void surface_facet( index_t model_facet_id, index_t& surface_id, index_t& surf_facet_id ) const ;     
        index_t model_facet( index_t surface_id, index_t surf_facet_id ) const ;

        // Accessors to model elements
        const std::string& name() const { return name_ ; }

        index_t nb_corners() const { return corners_.size() ; }
        index_t nb_lines() const { return lines_.size() ; }
        index_t nb_surfaces() const { return surfaces_.size() ; }
        index_t nb_regions() const { return regions_.size() ; }
        
        index_t nb_contacts() const { return contacts_.size() ; }        
        index_t nb_interfaces() const { return interfaces_.size() ; }
        index_t nb_layers() const { return layers_.size() ; }       
        
        const Corner& corner( index_t index ) const { return corners_.at(index) ; }
        const Line& line( index_t index ) const { return lines_.at(index) ; }
        const Surface& surface( index_t index ) const { return surfaces_.at(index) ; }
        const BoundaryModelElement& region( index_t index ) const { return regions_.at(index) ; }
        const BoundaryModelElement& universe() const { return universe_ ; }

        const BoundaryModelElement& element( index_t dim, index_t index ) const
        {
            switch( dim ){
                case 0: return corner( index ) ;
                case 1: return line( index ) ;
                case 2: return surface( index ) ;
                case 3: return region( index ) ;
                default:
                    grgmesh_assert_not_reached ;
                    return dummy_element ;
            }
         }

        const BoundaryModelElement& contact( index_t index ) const { return contacts_.at(index) ; }
        const BoundaryModelElement& one_interface( index_t index ) const { return interfaces_.at(index) ; }
        const BoundaryModelElement& layer( index_t index ) const { return layers_.at(index) ; }
        
        /// \todo Write a proper IO class for Boundary models
        bool save_gocad_model3d( std::ostream& out ) ;        
        void save_as_eobj_file( const std::string& file_name ) ;

        // Accessors to attribute managers
        PointAttributeManager* point_attribute_manager() const
        {
            return const_cast< PointAttributeManager* >( &point_attribute_manager_ ) ;
        }
        EdgeAttributeManager* edge_attribute_manager() const
        {
            return const_cast< EdgeAttributeManager* >( &edge_attribute_manager_ ) ;
        }
        FacetAttributeManager* facet_attribute_manager() const
        {
            return const_cast< FacetAttributeManager* >( &facet_attribute_manager_ ) ;
        }

    private:
        bool load_gocad_model3d( const std::string& in ) ;

        bool check_model3d_compatibility() ;
        static void save_type( std::ostream& out, GEOL_FEATURE t ) ;

    private:
        std::string name_ ;

        /** 
         * \brief Coordinates of the vertices of the model elements
         * Storage of points is unique for the whole model.
         */
        std::vector< vec3 >                 points_ ;

        // Base manifold elements of a model
        std::vector< Corner >               corners_ ;
        std::vector< Line >                 lines_ ;
        std::vector< Surface >              surfaces_ ;
        std::vector< BoundaryModelElement > regions_ ;

        /// The biggest volumetric region, defined by Box surfaces
        BoundaryModelElement universe_ ;

        // Number of facets in all previous surfaces
        // This has to be updated when a surface is modified !!
        // How can this be guaranteed ? Which class can modify a surface ?
        // Size = nb_surface()+1
        std::vector< index_t > nb_facets_ ;
    
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

        // Attribute managers 
        PointAttributeManager point_attribute_manager_ ;
        EdgeAttributeManager  edge_attribute_manager_ ;
        FacetAttributeManager facet_attribute_manager_ ;

    } ;   

    template< class ATTRIBUTE >
    class GRGMESH_API BoundaryModelPointAttribute: public Attribute< BoundaryModel::POINT, ATTRIBUTE > {
    public:
        typedef Attribute< BoundaryModel::POINT, ATTRIBUTE > superclass ;

        void bind( BoundaryModel* model, const std::string& name )
        {
            superclass::bind( model->point_attribute_manager(), model->nb_points(),
                name ) ;
        }

        void bind( BoundaryModel* model )
        {
            superclass::bind( model->point_attribute_manager(),
                model->nb_points() ) ;
        }

        BoundaryModelPointAttribute()
        {
        }

        BoundaryModelPointAttribute( BoundaryModel* model )
        {
            bind( model ) ;
        }

        BoundaryModelPointAttribute( BoundaryModel* model, const std::string& name )
        {
            bind( model, name ) ;
        }

        static bool is_defined( BoundaryModel* model, const std::string& name )
        {
            return superclass::is_defined( model->point_attribute_manager(), name ) ;
        }
    } ;

//    template< class ATTRIBUTE >
//    class GRGMESH_API BoundaryModelEdgeAttribute: public Attribute< BoundaryModel::EDGE, ATTRIBUTE > {
//    public:
//        typedef Attribute< BoundaryModel::EDGE, ATTRIBUTE > superclass ;
//
//        void bind( BoundaryModel* model, const std::string& name )
//        {
//            superclass::bind( model->edge_attribute_manager(), model->nb_edges(),
//                name ) ;
//        }
//
//        void bind( BoundaryModel* model )
//        {
//            superclass::bind( model->edge_attribute_manager(),
//                model->nb_edges() ) ;
//        }
//
//        BoundaryModelEdgeAttribute()
//        {
//        }
//
//        BoundaryModelEdgeAttribute( BoundaryModel* model )
//        {
//            bind( model ) ;
//        }
//
//        BoundaryModelEdgeAttribute( BoundaryModel* model, const std::string& name )
//        {
//            bind( model, name ) ;
//        }
//
//        static bool is_defined( BoundaryModel* model, const std::string& name )
//        {
//            return superclass::is_defined( model->edge_attribute_manager(), name ) ;
//        }
//    } ;

    
    template< class ATTRIBUTE >
    class GRGMESH_API BoundaryModelFacetAttribute: public Attribute< BoundaryModel::FACET, ATTRIBUTE > {
    public:
        typedef Attribute< BoundaryModel::FACET, ATTRIBUTE > superclass ;

        void bind( BoundaryModel* model, const std::string& name )
        {
            superclass::bind( model->facet_attribute_manager(), model->nb_facets(),
                name ) ;
        }

        void bind( BoundaryModel* model )
        {
            superclass::bind( model->facet_attribute_manager(),
                model->nb_facets() ) ;
        }

        BoundaryModelFacetAttribute()
        {
        }

        BoundaryModelFacetAttribute( BoundaryModel* model )
        {
            bind( model ) ;
        }

        BoundaryModelFacetAttribute( BoundaryModel* model, const std::string& name )
        {
            bind( model, name ) ;
        }

        static bool is_defined( BoundaryModel* model, const std::string& name )
        {
            return superclass::is_defined( model->facet_attribute_manager(), name ) ;
        }
    } ;



    /**
     * \brief Structure used to build contact when loading a BoundaryModel from .ml file 
     */
    struct Border {
        Border( index_t part, index_t corner, index_t p0, index_t p1):
        part_id_(part), corner_id_(corner), p0_(p0), p1_(p1) {};

        // Id of the Surface owning this Border
        index_t part_id_ ;
        // Id of p0 in the BoundaryModel corner vector
        index_t corner_id_ ;

        // Ids of the starting corner and second point on the border in the Surface
        // to which this Border belong
        index_t p0_ ;
        index_t p1_ ;
    } ;

    /**
     * \brief Build a BoundaryModel
     */ 
    class GRGMESH_API BoundaryModelBuilder {
    public:
        BoundaryModelBuilder( BoundaryModel& model )
            : model_( model ){}
        virtual ~BoundaryModelBuilder(){} ;
        void load_file( const std::string& in ) ;
        bool rebuild() ;
        void copy_macro_topology( const BoundaryModel* from ) ;

        index_t create_interface(
            const std::string& name,
            signed_index_t id = -1,
            GEOL_FEATURE type = default_type ) ;

        void add_interface_child( index_t id, index_t child ) {
            model_.interfaces_[id].add_child( child ) ;
        }

        /** A VOIR SI ion garde tout public l� dedans
         * Pas sur que ce soit passionnant et hyper int�ressant 
         */
    public:
        void reserve_nb_points( index_t size ) {
            model_.points_.reserve( size ) ; 
        }
        void reserve_nb_corners( index_t size ) {
            model_.corners_.reserve( size ) ;
        }
        void reserve_nb_lines( index_t size ) {
            model_.lines_.reserve( size ) ;
        }
        void reserve_nb_surfaces( index_t size ) {
            model_.surfaces_.reserve( size ) ;
        }
        void reserve_nb_interfaces( index_t size ) {
            model_.interfaces_.reserve( size ) ;
        }
        void reserve_nb_contacts( index_t size ) {
            model_.contacts_.reserve( size ) ;
        }
        void reserve_nb_regions( index_t size ) {
            model_.regions_.reserve( size ) ;
        }

        signed_index_t interface_id( const std::string& name ) const ;

        index_t find_or_create_corner( index_t index ) ;
        index_t find_or_create_line( index_t corner0, index_t corner1, std::vector< index_t >& points ) ;
        index_t find_or_create_contact( std::vector< index_t >& interfaces, GEOL_FEATURE type ) ;
        
        signed_index_t find_corner( const vec3& ) const ;
        signed_index_t find_corner( index_t ) const ;
        signed_index_t find_contact( const std::vector< index_t >& interfaces ) const ;
        signed_index_t find_line( index_t corner0, index_t corner1, const std::vector< index_t >& points ) const ;

        signed_index_t find_key_facet( index_t surface_id, const vec3& p0, const vec3& p1, const vec3& p2, 
            bool& same_orientation ) const ;  
         
        /**
         * \brief Check if the surface triangle orientations match the one of the key facet 
         */
        bool check_key_facet_orientation( index_t surface ) const ;
     
        index_t add_point( const vec3& point ) {            
            model_.points_.push_back( point ) ;
            return model_.nb_points()-1 ;
        }
        index_t add_point( double* point) {
            return add_point( vec3( point[0], point[1], point[2] ) ) ;
        }

        void add_corner_boundary( index_t id, index_t b ) {
            model_.corners_[id].add_boundary( b ) ;
        }
        void add_corner_in_boundary( index_t id, index_t b ) {
            model_.corners_[id].add_in_boundary( b ) ;
        }
        void add_line_boundary( index_t id, index_t b ) {
            model_.lines_[id].add_boundary( b ) ;
        }
        void add_line_in_boundary( index_t id, index_t b ) {
            model_.lines_[id].add_in_boundary( b ) ;
        }        
        void add_surface_boundary( index_t id, index_t b ) {
            model_.surfaces_[id].add_boundary( b ) ;
        }
        void add_surface_in_boundary( index_t id, index_t b ) {
            model_.surfaces_[id].add_in_boundary( b ) ;
        }
        
        void add_interface_boundary( index_t id, index_t b ) {
            model_.interfaces_[id].add_boundary( b ) ;
        }
         void add_region_in_boundary( index_t id, index_t b ) {
            model_.regions_[id].add_in_boundary( b ) ;
        }
        void add_region_oriented_boundary( index_t id, index_t b, bool side ) {
            model_.regions_[id].add_boundary( b, side ) ;
        }

        void add_contact_child( index_t id, index_t child ) {
            model_.contacts_[id].add_child( child ) ;
        }
        void add_contact_boundary( index_t id, index_t b ) {
            model_.contacts_[id].add_boundary( b ) ;
        }
        void add_contact_in_boundary( index_t id, index_t b ) {
            model_.contacts_[id].add_in_boundary( b ) ;
        }                     
        
        void add_layer_child( index_t id, index_t child ) {
            model_.layers_[id].add_child( child ) ;
        }

        void set_surface_first_triangle_as_key( index_t id ) {
            model_.surfaces_[id].set_first_triangle_as_key() ;
        }
        
        void set_name( BoundaryModelElement& e, const std::string& name ) {
            e.set_name( name ) ;
        }
        void set_name( const std::string& name ) {
            model_.name_ = name ;
        }

        void set_parent( BoundaryModelElement& e, index_t id ) {
            e.set_parent( id ) ;
        }
            
        index_t create_line( 
            signed_index_t id = -1, 
            const std::vector< index_t >& points = empty_index_vector ) ;
        
        index_t create_surface(
            signed_index_t id = -1,
            signed_index_t parent = -1,
            GEOL_FEATURE type = default_type ) ;

        void create_surface(
            const std::string& interface_name,
            const std::string& type,
            const KeyFacet& key ) ;
        
        index_t create_region( signed_index_t id = -1 ) ;

        index_t create_region(
            const std::string& name,
            const std::vector< std::pair< index_t, bool > >& boundaries,
            signed_index_t id = -1 ) ;   

        index_t create_layer( const std::string& name, signed_index_t id = -1 ) ;
      
        void set_corner( index_t  corner_id, index_t point_id ) ;
        void set_corner( index_t  corner_id, const vec3& point ) ;
        void set_line( index_t id, const std::vector< index_t >& vertices ) ;
        void set_line( index_t id, const std::vector< vec3 >& vertices ) ;
        
        index_t determine_line_vertices( 
            const Surface& S, 
            index_t first_point, 
            index_t second_point,            
            std::vector< index_t >& border_point_model_ids ) const ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& surface_points,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr,
            const std::vector< index_t >& surface_adjacencies = empty_index_vector ) ;
        void set_surface_adjacencies( index_t surface_id ) ;

        void cut_surface_by_line( index_t surface_id, index_t line_id ) ;
             
        void set_universe( const std::vector< std::pair< index_t, bool > >& boundaries ) ;        
        void remove_universe_from_regions( index_t id ) ;

        // This is a big mess  !!!!
        /// \todo Trade the end_something functions in the BoundaryModelBuilder
        /// functions for a smart end_model function
        /// that checks model validity and complete all missing parts
        void make_points_unique() ;
        void build_lines( const std::vector< Border >& borders ) ;
        void build_contacts() ;
        void end_contacts() ;
        void end_lines() ;
        void end_interfaces() ;
        void end_surfaces() ;
        void end_surfaces( const std::vector< index_t >& change_orientation ) ;
        void end_corners() ;
        void end_layers() ; 
        void end_model() ;


        void update_all_ids() ;

        static GEOL_FEATURE determine_geological_type( const std::string& in ) ;
        static GEOL_FEATURE determine_type( const std::vector< GEOL_FEATURE >& types ) ;
       
    protected:
        BoundaryModel& model_ ;
    } ;
}

#endif
