/*[
* NOT ONLY ASGA 
*/
/*! \author Jeanne Pellerin and Arnaud Botella */


#ifndef __GRGMESH_BOUNDARY_MODEL__
#define __GRGMESH_BOUNDARY_MODEL__

#include <grgmesh/boundary_model_element.h>
#include <grgmesh/common.h>
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
    static std::vector< unsigned int > empty_uint_vector ;
    static vec3 dummy = vec3( 0, 0, 0 ) ;

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
        unsigned int nb_points() const { return points_.size() ; }
        const vec3& point( unsigned int p ) const { return points_.at(p) ; }
        unsigned int point_index( const vec3& p ) const ;

        inline uint32 nb_facets() const ;
        void surface_facet( uint32 model_facet_id, uint32& surface_id, uint32& surf_facet_id ) const ;     
        uint32 model_facet( uint32 surface_id, uint32 surf_facet_id ) const ;

        // Accessors to model elements
        const std::string& name() const { return name_ ; }

        unsigned int nb_corners() const { return corners_.size() ; }
        unsigned int nb_lines() const { return lines_.size() ; }
        unsigned int nb_surfaces() const { return surfaces_.size() ; }
        unsigned int nb_regions() const { return regions_.size() ; }
        
        unsigned int nb_contacts() const { return contacts_.size() ; }        
        unsigned int nb_interfaces() const { return interfaces_.size() ; }
        unsigned int nb_layers() const { return layers_.size() ; }       
        
        const Corner& corner( int index ) const { return corners_.at(index) ; }
        const Line& line( int index ) const { return lines_.at(index) ; }
        const Surface& surface( int index ) const { return surfaces_.at(index) ; }
        const BoundaryModelElement& region( int index ) const { return regions_.at(index) ; }

        const BoundaryModelElement& element( int dim, int index ) const 
        {
            switch( dim ){
                case 0: return corner( index ) ;
                case 1: return line( index ) ;
                case 2: return surface( index ) ;
                case 3: return region( index ) ;
                default: grgmesh_assert_not_reached ;
            }
         }

        const BoundaryModelElement& contact( int index ) const { return contacts_.at(index) ; }
        const BoundaryModelElement& one_interface( int index ) const { return interfaces_.at(index) ; }
        const BoundaryModelElement& layer( int index ) const { return layers_.at(index) ; }                      
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
        bool load_gocad_model3d( std::istream& in ) ;

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
        std::vector< uint32 > nb_facets_ ;
    
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
        Border( int part, int corner, int p0, int p1):
        part_id_(part), corner_id_(corner), p0_(p0), p1_(p1) {};

        // Id of the Surface owning this Border
        int part_id_ ; 
        // Id of p0 in the BoundaryModel corner vector
        int corner_id_ ;

        // Ids of the starting corner and second point on the border in the Surface
        // to which this Border belong
        int p0_ ; 
        int p1_ ; 
    } ;

    /**
     * \brief Build a BoundaryModel
     */ 
    class GRGMESH_API BoundaryModelBuilder {
    public:
        BoundaryModelBuilder( BoundaryModel& model )
            : model_( model ){}
        virtual ~BoundaryModelBuilder(){} ;
        void load_file( std::istream& in ) ;
        bool rebuild() ;
        void copy_macro_topology( const BoundaryModel* from ) ;

        unsigned int create_interface(
            const std::string& name,
            int id = -1,
            GEOL_FEATURE type = default_type ) ;

        void add_interface_child( unsigned int id, unsigned int child ) {
            model_.interfaces_[id].add_child( child ) ;
        }

        /** A VOIR SI ion garde tout public l� dedans
         * Pas sur que ce soit passionnant et hyper int�ressant 
         */
    public:
        void reserve_nb_points( unsigned int size ) {
            model_.points_.reserve( size ) ; 
        }
        void reserve_nb_corners( unsigned int size ) {
            model_.corners_.reserve( size ) ;
        }
        void reserve_nb_lines( unsigned int size ) {
            model_.lines_.reserve( size ) ;
        }
        void reserve_nb_surfaces( unsigned int size ) {
            model_.surfaces_.reserve( size ) ;
        }
        void reserve_nb_interfaces( unsigned int size ) {
            model_.interfaces_.reserve( size ) ;
        }
        void reserve_nb_contacts( unsigned int size ) {
            model_.contacts_.reserve( size ) ;
        }
        void reserve_nb_regions( unsigned int size ) {
            model_.regions_.reserve( size ) ;
        }

        int interface_id( const std::string& name ) const ;

        int find_or_create_corner( unsigned int index ) ;
        int find_or_create_line( int corner0, int corner1, std::vector< uint32 >& points ) ;
        int find_or_create_contact( std::vector< int >& interfaces, GEOL_FEATURE type ) ;
        
        int find_corner( const vec3& ) const ;
        int find_corner( uint32 ) const ;
        int find_contact( const std::vector< int >& interfaces ) const ;
        int find_line( int corner0, int corner1, const std::vector< uint32 >& points ) const ;        

        int find_key_facet( uint32 surface_id, const vec3& p0, const vec3& p1, const vec3& p2, 
            bool& same_orientation ) const ;  
         
        /**
         * \brief Check if the surface triangle orientations match the one of the key facet 
         */
        bool check_key_facet_orientation( uint32 surface ) const ;
     
        unsigned int add_point( const vec3& point ) {            
            model_.points_.push_back( point ) ;
            return model_.nb_points()-1 ;
        }       
        void add_corner_boundary( unsigned int id, unsigned int b ) {
            model_.corners_[id].add_boundary( b ) ;
        }
        void add_corner_in_boundary( unsigned int id, unsigned int b ) {
            model_.corners_[id].add_in_boundary( b ) ;
        }
        void add_line_boundary( unsigned int id, unsigned int b ) {
            model_.lines_[id].add_boundary( b ) ;
        }
        void add_line_in_boundary( unsigned int id, unsigned int b ) {
            model_.lines_[id].add_in_boundary( b ) ;
        }        
        void add_surface_boundary( unsigned int id, unsigned int b ) {
            model_.surfaces_[id].add_boundary( b ) ;
        }
        void add_surface_in_boundary( unsigned int id, unsigned int b ) {
            model_.surfaces_[id].add_in_boundary( b ) ;
        }
        
        void add_interface_boundary( unsigned int id, unsigned int b ) {
            model_.interfaces_[id].add_boundary( b ) ;
        }
         void add_region_in_boundary( unsigned int id, unsigned int b ) {
            model_.regions_[id].add_in_boundary( b ) ;
        }
        void add_region_oriented_boundary( unsigned int id, unsigned int b, bool side ) {
            model_.regions_[id].add_boundary( b, side ) ;
        }

        void add_contact_child( unsigned int id, unsigned int child ) {
            model_.contacts_[id].add_child( child ) ;
        }
        void add_contact_boundary( unsigned int id, unsigned int b ) {
            model_.contacts_[id].add_boundary( b ) ;
        }
        void add_contact_in_boundary( unsigned int id, unsigned int b ) {
            model_.contacts_[id].add_in_boundary( b ) ;
        }                     
        
        void add_layer_child( unsigned int id, unsigned int child ) {
            model_.layers_[id].add_child( child ) ;
        }

        void set_surface_first_triangle_as_key( uint32 id ) {
            model_.surfaces_[id].set_first_triangle_as_key() ;
        }
        
        void set_name( BoundaryModelElement& e, const std::string& name ) {
            e.set_name( name ) ;
        }
        void set_name( const std::string& name ) {
            model_.name_ = name ;
        }

        void set_parent( BoundaryModelElement& e, uint32 id ) {
            e.set_parent( id ) ;
        }
            
        unsigned int create_line( 
            int id = -1, 
            const std::vector< unsigned int >& points = empty_uint_vector ) ;
        
        unsigned int create_surface(
            int id = -1,
            int parent = -1,
            GEOL_FEATURE type = default_type ) ;

        void create_surface(
            const std::string& interface_name,
            const std::string& type,
            const KeyFacet& key ) ;
        
        unsigned int create_region( int id = -1 ) ;

        unsigned int create_region(
            const std::string& name,
            const std::vector< std::pair< int, bool > >& boundaries,
            int id = -1 ) ;   

        unsigned int create_layer( const std::string& name, int id = -1 ) ;
      
        void set_corner( unsigned int  corner_id, unsigned int point_id ) ;
        void set_line( unsigned int id, const std::vector< uint32 >& vertices ) ;
        void set_line_geometry( unsigned int id, const std::vector< unsigned int >& line_points ) ;
        
        unsigned int determine_line_vertices( 
            const Surface& S, 
            unsigned int first_point, 
            unsigned int second_point,            
            std::vector< unsigned int >& border_point_model_ids ) const ;

        void set_surface_geometry(
            unsigned int surface_id,
            const std::vector< unsigned int >& surface_points,
            const std::vector< unsigned int >& surface_facets,
            const std::vector< unsigned int >& surface_facet_ptr ) ;
        void set_surface_adjacencies( unsigned int surface_id ) ;

        void cut_surface_by_line( uint32 surface_id, uint32 line_id ) ;
             
        void set_universe( const std::vector< std::pair< int, bool > >& boundaries ) ;        
        void remove_universe_from_regions( unsigned int id ) ;

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
        void end_surfaces( const std::vector< int >& change_orientation ) ;
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



#ifdef TOTOTO
    /**
     * \brief Measure things on a BoundaryModel
     * Output topological measures - and geometrical measures 
     */
    class BoundaryModelMeasure {
    public :
        BoundaryModelMeasure( const BoundaryModel& model ):model_(model){} ;
        ~BoundaryModelMeasure() ;

        /*! Output global information on the model and its elements */
        void print_topology( std::ofstream& out ) const ;
        void print_element_info( std::ofstream& out ) const ;

        static void print_type( std::ostream& out, GEOL_FEATURE t, int dim ) ;

    private:        
        unsigned int nb_surface_inside() const ;
        unsigned int nb_line_inside() const ;
        unsigned int nb_real_corners_inside() const ;
    
        int nb_real_corners() const ;
        int nb_non_manifold_lines() const ;
        int nb_surface_with_free_boundary() const ;    

    private:
        const BoundaryModel& model_ ;
    } ;

#endif     

#endif
