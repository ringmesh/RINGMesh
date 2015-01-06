/*! File part of phD work of Jeanne Pellerin ? copyright ASGA ? 
 *  
 * 
 *  Developed with Arnaud Botella
 *
 *  Will GRGMesh go open source ?
 */ 

#ifndef __GRGMESH_BOUNDARY_MODEL_ELEMENT__
#define __GRGMESH_BOUNDARY_MODEL_ELEMENT__

#include <grgmesh/common.h>
#include <grgmesh/grgmesh_assert.h>

#include <vector> 
#include <map>
#include <string>
#include <algorithm> 
#include <cmath>
#include <cassert>

namespace GRGMesh {
    class BoundaryModel ;
}


namespace GRGMesh {


    /********************************************************************************************/
    /******             General stuff               *********************************************/
    /********************************************************************************************/   
    
   /**
     * \brief A triangle that set the orientation of one Surface's facets
     * 
     */
    struct KeyFacet {
        KeyFacet( const vec3& p0, const vec3& p1, const vec3& p2 ):
            p0_(p0), p1_(p1), p2_(p2){} ;

        KeyFacet():p0_(), p1_(), p2_() {} ;
    
        bool is_default() const {
            if( p0_.x == 1 && p0_.y == -1 && p0_.z == -1 &&
                p1_.x == 1 && p1_.y == -1 && p1_.z == -1  && 
                p2_.x == 1 && p2_.y == -1 && p2_.z == -1  ) return true ;
            else {
                // la compilation ne passe pas ça m'énerve
                /// \todo put this assert back
                //assert( !(p0_==p1_) && !(p0_==p2_) && !(p1_==p2_) ) ;
                return false ;
            }
        }
    public:
        vec3 p0_ ;
        vec3 p1_ ;
        vec3 p2_ ;
    } ;

    /**
     * Types for BoundaryModelElement 
     * \todo Read all types, this is not sufficient
     */
    enum GEOL_FEATURE {
        ALL,
        STRATI,
        FAULT,
        VOI,
        STRATI_FAULT,
        STRATI_VOI,
        FAULT_VOI
    } ;
    /// Default type is all
    static GEOL_FEATURE default_type = ALL ;

    /********************************************************************************************/
    /***********       Components of a Model declaration                 ************************/
    /********************************************************************************************/

    /**
     * \brief Generic class describing the elements of a BoundaryModel
    */
    class GRGMESH_API BoundaryModelElement {
        friend class BoundaryModelBuilder ;

    public:
        BoundaryModelElement(
            BoundaryModel* model = nil, int dim = -1, int id = -1, int parent = -1,
            GEOL_FEATURE type = default_type )
            : model_( model ),  name_( "" ), id_( id ), dim_( dim ), type_( type ),
            parent_( parent )
        {
        }
        
        virtual ~BoundaryModelElement() { }

        const BoundaryModel& model() const { return *model_ ; }
        
        const std::string& name() const { return name_ ; }
        int id() const { return id_ ; }
        uint8 dim() const { return dim_ ; }  
        GEOL_FEATURE type() const { return type_ ; }

        bool is_on_voi() const ;
        bool side( int i ) const { return sides_[i] ; }

        bool has_parent() const { return parent_ != -1 ; }
        const BoundaryModelElement& parent() const ;
        int parent_id() const { return parent_ ; }
        
        unsigned int nb_boundaries() const { return boundaries_.size() ; }
        unsigned int boundary_id( unsigned int x ) const { return boundaries_[x] ; }
        const BoundaryModelElement& boundary( unsigned int x ) const ;
        
        unsigned int nb_in_boundary() const { return in_boundary_.size() ; }
        unsigned int in_boundary_id( unsigned int x ) const { return in_boundary_[x] ; }
        const BoundaryModelElement& in_boundary( unsigned int x ) const ;
        
        unsigned int nb_children() const { return children_.size() ; }
        unsigned int child_id( unsigned int x ) const { return children_[x] ; }
        const BoundaryModelElement& child( unsigned int x ) const ;
        

        virtual unsigned int nb_cells() const ;
        virtual unsigned int nb_points() const {
            grgmesh_assert_not_reached ;  return 0 ;
        }
        virtual unsigned int model_point_id( unsigned int p = 0 ) const {
            grgmesh_assert_not_reached ; return 0 ;
        }
  
        virtual const vec3& point( unsigned int p = 0 ) const {
            grgmesh_assert_not_reached ; return vec3( 0, 0, 0 ) ;
        }
  
        
        
        /*virtual bool is_triangulated() const {
            grgmesh_assert( dim_ == 3 ) ;
            for( unsigned int s = 0; s < nb_boundaries(); s++ ) {
                if( !boundary( s ).is_triangulated() ) return false ;
            }
            return true ;
        }*/
        
    protected:
        void copy_macro_topology(
            const BoundaryModelElement& rhs, BoundaryModel& model ) ;
        
        void set_parent( int p ){ parent_ = p ; }
        void set_name( const std::string& name ) { name_ = name ; }
        void set_type( GEOL_FEATURE type ) { type_ = type ; } 
        void set_dim( int dim ) { dim_ = dim ; }
        void set_id( int id ) { id_ = id ; }
        
        void add_boundary( unsigned int b ) { boundaries_.push_back( b ) ; }
        void add_boundary( unsigned int b, bool side ) { boundaries_.push_back(b) ; sides_.push_back(side) ; }
        virtual void add_in_boundary( unsigned int e ) { in_boundary_.push_back(e) ; }
        void add_child( unsigned int e ){ children_.push_back( e ) ; }
       

    protected :
        /// \todo Can we have something else than a POINTER ?? default constructor is needed...
        BoundaryModel* model_ ;

        /// Name of the element, empty string if none nothing
        std::string name_ ;

        /// Id of this element in the appropriate vector of the BoundaryModel owning it
        int id_ ;
        
        /// Dimension of the element 0 corner; 1 line; 2 surface; 3 region
        uint8 dim_ ;
        
        /// Geological type for this object, default is ALL 
        GEOL_FEATURE type_ ;

        /// Elements on the boundary of this element - their dimension is dim_-1
        std::vector< unsigned int > boundaries_ ;

        /// Flag on which side of the boundary is this element 
        /// Filled for volumetric regions only + (true) or - (false)
        std::vector< bool > sides_ ; 
        
        /// Elements in which boundary this element is - their dimension is dim_+1
        std::vector< unsigned int > in_boundary_ ;

        /// Index of the parent (group of elements to which belong this).
        /// Default value is -1 - no parent.
        int parent_ ;

        /// The group elements making up this one, empty for basic elements
        std::vector< unsigned int > children_ ;
    } ;

    /*-----------------------------------------------------------------------------------------*/
    /*! Corners of the model
    *  Dimension 0
    */
    class GRGMESH_API Corner : public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
    public:
        Corner(
            BoundaryModel* model,
            int id = -1,
            unsigned int p = 0 )
            : BoundaryModelElement( model, 0, id ), p_( p )
        {
        }
        virtual ~Corner()
        {
        }
        bool is_real() const { return in_boundary_.size() > 1 ; }

        virtual const vec3& point( unsigned int p = 0 ) const ;
        
        virtual unsigned int nb_cells() const { return 1 ; }
        virtual unsigned int nb_points() const { return 1 ; }
        virtual unsigned int model_point_id( unsigned int p = 0 ) const { return p_ ; } 
        
    private:
        //void copy_macro_topology( const Corner& rhs, BoundaryModel& model ) ;

        void set_point( unsigned int p ) { p_ = p ; }
    private:
        unsigned int p_ ;
    };
    /*-----------------------------------------------------------------------------------------*/        
    /*! A part of the contact between 2 interfaces
    * 
    *  boundaries -> 1 or 2 corners
    *  in_boundary -> several surface parts
    */
    class GRGMESH_API Line: public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
        friend class LineMutator ;
    public:
        Line( BoundaryModel* model, int id = -1 ) ;
        Line(
            BoundaryModel* model,
            int id,
            const std::vector< unsigned int >& points ) ;
        Line(
            BoundaryModel* model,
            int id,
            unsigned int corner0,
            unsigned int corner1,
            const std::vector< unsigned int >& points ) ;
        virtual ~Line(){} ;

   
        virtual unsigned int nb_cells() const {            
            return points_.size()-1 ; 
        }
        // If the line is closed the last point is equal to the first one
        virtual unsigned int nb_points() const {             
            return points_.size() ; 
        }
        virtual unsigned int model_point_id( unsigned int p ) const {
            return points_.at(p) ;
        }

        //bool contains( const vec3& p ) const ;
        //int find( const vec3& p ) const ;
        //int facet_point_id( const vec3& p ) const ;
        bool is_closed () const { 
            return (boundaries_[0]!= nil ) && (boundaries_[0] == boundaries_[1]) ; 
        }  

        // Returns true if this line is twice in the boundary of the surface e
        bool is_inside_border( const BoundaryModelElement& e ) const ;
            
        virtual const vec3& point( unsigned int line_point_id ) const ;
        inline vec3 segment_barycenter( uint32 s ) const ;

        inline double segment_length( uint32 s ) const ;
            

    private:
        /*void copy_macro_topology(
            const Line& rhs,
            BoundaryModel& model ) ;*/
        void set_vertices( const std::vector< unsigned int >& vertices ) {
            points_ = vertices ;
        }
        //virtual void add_in_boundary( unsigned int e ) ;
        //void set_is_inside_border( bool x ) { is_inside_border_.push_back(x) ; }
    private:
        /// In case of a closed line, the last vertex (equal to the first) is not stored.
        /// ATTENTION ÇA A CHANGÉ SINON C'EST TROP CHIANT, et en fait c'était pas bon
        std::vector< unsigned int > points_ ;
                 
        //std::vector< bool > is_inside_border_ ;
    } ;

    class GRGMESH_API LineMutator {
    public:
        LineMutator( Line& M )
            : M_( M )
        {
        }
        LineMutator( const Line& M )
            : M_( const_cast< Line& >( M ) )
        {
        }
        void set_point( unsigned int id, const vec3& p ) ;
        vec3& point( unsigned int p ) const ;
        std::vector< unsigned int >& points() const { return M_.points_ ; }

        void clear()
        {
            M_.points_.clear() ;
        }

    private:
        Line& M_ ;
    };
    /*-----------------------------------------------------------------------------------------*/
       
    /**
    * \brief A polygonal manifold surface part of a BoundaryModel
    * 
    * This is a BoundaryModelElement of dimension 2.
    * Its boundaries_ is made of Line - It is in_boundary_ of 1 or 2 Region
    *
    */
    class GRGMESH_API Surface : public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
        friend class SurfaceMutator ;
    public:
        Surface(
            BoundaryModel* model,
            int id = -1,
            int parent = -1,
            const GEOL_FEATURE& type = default_type )
            :
            BoundaryModelElement( model, 2, id, parent, type ),
            is_triangulated_( false )
        {
        }
        virtual ~Surface(){} ;
       
        virtual uint32 nb_cells() const { return facets_.empty() ? 0 : facet_ptr_.size() - 1 ; }
        virtual uint32 nb_points() const { return points_.size() ; }               
        virtual uint32 model_point_id( unsigned int p ) const {
            return points_[p] ;
        }
        
        //bool key_facet_orientation() ;
        const KeyFacet& key_facet() const { return key_facet_ ; }       

        bool is_triangulated() const { return is_triangulated_ ; }
        
        unsigned int facet_begin( uint32 f ) const { return facet_ptr_.at(f) ; }
        unsigned int facet_end( uint32 f ) const { return facet_ptr_.at(f+1) ; }
        unsigned int nb_points_in_facet( uint32 f ) const { return facet_end( f ) - facet_begin( f ) ; }
        unsigned int next_in_facet( uint32 f, uint32 v ) const { 
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            if( v != nb_points_in_facet(f)-1 ) return v+1 ;
            else return 0 ;
        }
        unsigned int prev_in_facet( uint32 f, uint32 v ) const {
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            if( v > 0 ) return v-1 ;
            else return nb_points_in_facet(f)-1 ;
        }            
        
        inline const vec3& point( uint32 f, uint32 v ) const ;

        /** 
         * Access to the ids of the points in the model
         * The same point can appear several times - boundaries inside 
         * Necessary for export
         * Should be inline - but linking fails because implementation is in cpp
         * and cannot be there because including boundary_model.h results in circular includes
         */        
        const vec3& point( uint32 surf_point_id ) const ;
        
        /** Returns the id of point \param v in facet \param f in this surface */
        unsigned int surf_point_id( uint32 f, uint32 v ) const {
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            return facets_[facet_begin(f)+v] ; 
        }
        /** Returns the id of point  \param v, in facet \param f in the parent BoundaryModel */ 
        unsigned int model_point_id( uint32 f, uint32 v ) const { 
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            return points_[ facets_[facet_begin(f)+v] ] ; 
        }     

        int surf_point_id( uint32 model_point_id ) const {
            for( uint32 i = 0; i < points_.size() ; ++i ){
                if ( points_[i] == model_point_id ) return i ;
            }
            return -1 ;
        }
        
        /** Returns the id of the adjacent facet of \param f in this surface along the edge starting at \param v */
        int adjacent( uint32 f, uint32 v ) const { 
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            return adjacent_[facet_begin(f)+v] ; 
        }
        //int adjacent_in_neighbor( unsigned int f, unsigned int e ) const ;
        
        bool is_on_border( uint32 f, uint32 v ) const {
            grgmesh_debug_assert( v < nb_points_in_facet(f) ) ;
            return adjacent( f, v ) == -1 ; 
        }
        bool is_on_border( uint32 f ) const {
            for( unsigned int adj = 0; adj < nb_points_in_facet(f); adj++ ) {
                if( is_on_border( f, adj ) ) return true ;
            }
            return false ;
        }
        
        /**
         * For a bidirectional border traversal 
         * From the input facet f, get the facet that share point v and 
         * get the ids of point v and of the following point in this next facet
         * next facet may be the same, from is required to avoid going back 
         */
        void next_on_border(
            uint32 f, uint32 from, uint32 v,
            int32& next_f, int32& v_in_next, int32& to  ) const ;
        
        /**
         * One directional border traversal
         */
        void next_on_border( uint32 f, uint32 e, int32& next_f, int32& next_e ) const ;

        //FacetEdge next_on_border( const FacetEdge& te ) const ;
       
        bool is_triangle( uint32 f ) const { return nb_points_in_facet( f ) == 3 ; }

//        int find_facet( int id0, int id1 ) const ;
        int facet_from_model_point_ids( uint32 i0, uint32 i1 ) const ;
        void edge_from_model_point_ids( uint32 i0, uint32 i1, int32& f, int32& e ) const ;
        void oriented_edge_from_model_point_ids( uint32 i0, uint32 i1, int32& facet, int32& edge ) const ;

        
        //int find_facet( const vec3& p0, const vec3& p1 ) const ; // TO DO A VIRER !!
        //int find_edge( int id0, int id1 ) const ;

        /** Returns the id of the point with surface id p0 in the given facet */
        int facet_point_id( uint32 t, uint32 surf_point_id ) const ;
        /** Returns the id of the points with the given coordinates in the given facet */
     //   int facet_point_id( int t, const vec3& p ) const ;
        //int edge_id( int t, int p0, int p1 ) const ;

        //Box3d bbox() const ;
        //bool contains( const vec3& p ) const ;
        //int find( const vec3& p ) const ;

        int facets_around_point( 
            uint32 surf_point_id, 
            std::vector< uint32 >& result, 
            bool border_only ) const ;
        
        inline vec3 facet_barycenter( int f ) const ;  
        inline double facet_area( int f ) const ;
        inline vec3 facet_normal( int f ) const ;
        void point_normals( std::vector< vec3 >& normals ) const ;
        
        int facets_around_point(
            uint32 surf_point_id,
            std::vector< uint32 >& result,
            bool border_only,
            uint32 first_facet ) const ;

        unsigned int closest_point_in_facet( uint32 f, const vec3& point ) const ;

    private:
        void set_key_facet( const KeyFacet& key ) { key_facet_ = key ; }
        void set_first_triangle_as_key() ;
        
        // On en a besoin à la construction
        // On ne peut pas utiliser model_point ids - car on ne différencie pas
        // alors les arrêtes sur des bords internes = 1 bord - mais points identiques sur
        // triangles en face dans la même surfaace
        int facet_from_surface_point_ids( uint32 i0, uint32 i1 ) const ;


        void set_adjacent( int f, int e, int adjacent ) { adjacent_[facet_begin(f)+e] = adjacent ; }

        void compute_is_triangulated() {
            for( unsigned int f = 0; f < nb_cells(); f++ ) {
                if( !is_triangle( f ) ) {
                    is_triangulated_ = false ;
                    return ;
                }
            }
            is_triangulated_ = true ;
        }       
        //int has_edge( int f, int v0, int v1 ) const ;
        //int has_oriented_edge( int f, int v0, int v1 ) const ;

        /*bool same_point( int i, int j ) const {
            if( i == j ) return true ;
            else if( points_[i] == points_[j] ) return true ;
            else return false ;
        }*/


        void set_geometry(
            const std::vector< uint32 >& points,
            const std::vector< uint32 >& facets,
            const std::vector< uint32 >& facet_ptr )
        {
            // Are these copies parallelized ?
            points_ = points ;
            facets_ = facets ;
            facet_ptr_ = facet_ptr ;

            compute_is_triangulated() ;
        }

        void set_adjacent( const std::vector< int >& adjacent ){
            grgmesh_assert( adjacent.size() == facets_.size() ) ;
            adjacent_ = adjacent ;
        }

       
    private:
        KeyFacet key_facet_ ;

        std::vector< uint32 > points_ ;
        std::vector< uint32 > facets_ ;
        std::vector< uint32 > facet_ptr_ ;

    
        // The adjacent facet is given for each vertex of each facet for the edge
        // starting at this vertex.
        // -1 if the edge is along a Line (Surface boundary)    
        std::vector< int > adjacent_ ;

        bool is_triangulated_ ;
        
    };  


    class GRGMESH_API SurfaceMutator {
    public:
        SurfaceMutator( Surface& M )
            : M_( M )
        {
        }
        SurfaceMutator( const Surface& M )
            : M_( const_cast< Surface& >( M ) )
        {
        }
        void set_point( unsigned int id, const vec3& p ) ;
        vec3& point( unsigned int p ) const ;
        std::vector< unsigned int >& points() const { return M_.points_ ; }
        std::vector< unsigned int >& facets() const { return M_.facets_ ; }
        std::vector< unsigned int >& facet_ptr() const { return M_.facet_ptr_ ; }
        std::vector< int >& adjacents() const { return M_.adjacent_ ; }

        void clear()
        {
            M_.points_.clear() ;
            M_.facets_.clear() ;
            M_.facet_ptr_.clear() ;
            M_.adjacent_.clear() ;
        }

    private:
        Surface& M_ ;
    };



    /**
     * \brief Class to answer geometrical request on BoundaryModelElement
     */
    class GRGMESH_API BoundaryModelElementMeasure {
    public:
        static double size( const BoundaryModelElement* E ) ;
        static double cell_size( const BoundaryModelElement* E, uint32 cell ) ;
        static double distance( const BoundaryModelElement* from,  const vec3& p ) ;
        static double distance( const BoundaryModelElement* from, const BoundaryModelElement* to ) ;
        static vec3 barycenter ( const BoundaryModelElement* E ) ;
        static vec3 barycenter ( const BoundaryModelElement* E, const std::vector< uint32 >& cells ) ;
    } ;



# ifdef TOTOTO
    class BoundaryModelElementMeasure {
    public:
        BoundaryModelElementMeasure( const BoundaryModelElement& in ) :
          element_ ( in ) {};

        // Connectivity measures
        std::vector< const BoundaryModelElement* > compute_neighbors(
            GEOL_FEATURE through, bool exclude_boundaries = true ) const ;
        
        

        // Geometrical measures
        virtual double size () const ;
        virtual double cell_size( int /*i*/ ) const { return 0. ; }

        virtual double distance( const vec3& p ) const ;
        virtual double distance( const BoundaryModelElement& e ) const ;
        virtual double distance( int /*simplex*/, const BoundaryModelElement& /*to*/ ) const { return -1 ; }
        
        virtual double min_angle( const BoundaryModelElement& e ) const ;
        virtual void angles(
            const BoundaryModelElement& /*with*/,
            std::vector< std::pair< double, double > >& /*values*/,
            bool /*same_side*/ ) const {} ;
        
        virtual vec3 average_orientation() const ;      
        double average_angle_to_y() const ;
        double average_dip() const ;

        static void print_categories( std::ostream& out ) ;
        void print( std::ostream& out ) const ;      

        static void print_complexity_categories( std::ostream& out ) ;
        virtual void print_complexity( std::ostream& out ) const ;

        void compute_distances( std::vector< std::pair< double, double > >& values ) const ;
        void compute_angles( std::vector< std::pair< double, double > >& values ) const ;
        int nb_incident_elements( bool exclude_voi = true ) const ;
        int nb_boundary_elements( bool exclude_voi = true ) const ;

    protected:
        const BoundaryModelElement& element_ ;


    } ;

    // In corners 

    virtual void print_complexity( std::ostream& out ) const ;
        virtual double size() const { return 0. ; }
        virtual double distance( const vec3& p ) const { return ::Geomesh::distance( p, point() ) ; }
        virtual double distance( BoundaryModelElement* e ) const { return e->distance( point() ) ; }

        // In Lines
           virtual double size() const ;
        virtual double cell_size( int i ) const ;
        virtual double distance( const vec3& p ) const ;
        virtual double distance( BoundaryModelElement* e ) const ;
        virtual double distance( int simplex, BoundaryModelElement* to ) const ;        
        virtual double min_angle( BoundaryModelElement* e ) const ;
        virtual vec3 average_orientation() const ;      
        
        virtual void angles(
            BoundaryModelElement* with,
            std::vector< std::pair< double, double > >& values,
            bool same_side = false ) const ;

        virtual void print_complexity( std::ostream& out ) const ;

        // On facets
           vec3 closest_normal( const vec3& p, vec3& b ) const ;
            virtual double size() const ;
        virtual double cell_size( int i ) const ;
        virtual double distance( const vec3& p ) const ;
        virtual double distance( BoundaryModelElement* e ) const ;
        virtual double distance( int facet, BoundaryModelElement* to ) const ;
        virtual double min_angle( BoundaryModelElement* e ) const ;
        virtual void angles( 
            BoundaryModelElement* with,
            std::vector< std::pair< double, double > >& values,
            bool same_side = false ) const ;
        virtual vec3 average_orientation() const ;
          virtual void print_complexity( std::ostream& out ) const ;
        void print_mesh( const std::string& filename ) const ;

                vec3 facet_normal( int t ) const ;
        void point_normal( std::vector< vec3 >& normals ) const ;
#endif // Code à porter



} // namespace

#endif

