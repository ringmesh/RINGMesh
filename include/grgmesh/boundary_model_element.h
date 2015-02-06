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


/*! \author Jeanne Pellerin and Arnaud Botella */

#ifndef __GRGMESH_BOUNDARY_MODEL_ELEMENT__
#define __GRGMESH_BOUNDARY_MODEL_ELEMENT__

#include <grgmesh/common.h>
#include <grgmesh/attribute.h>
#include <grgmesh/utils.h>

#include <vector> 
#include <string>

namespace GRGMesh {
    class BoundaryModel ;
}


namespace GRGMesh {     

    /*!
     * \brief Generic class describing one element of a BoundaryModel
     */
    class GRGMESH_API BoundaryModelElement {

    public:
        enum AttributeLocation {
            VERTEX, FACET
        } ;
        typedef AttributeManager< VERTEX > VertexAttributeManager ;
        typedef AttributeManager< FACET > FacetAttributeManager ;       

        /*!
        * @brief Geological feature types for BoundaryModelElement 
        * \todo Read all types, this is not sufficient
        */
        enum GEOL_FEATURE {
            ALL_GEOL,
            NO_GEOL,
            STRATI,
            FAULT,
            VOI,
            STRATI_FAULT,
            STRATI_VOI,
            FAULT_VOI
        } ;
         
        /*!
         * @brief Each BoundaryModelElement has a type
         * @details When no type is defined NO_TYPE should be used
         * There is two main categories of elements
         *   - low-level elements (CORNER, LINE, SURFACE, REGION) have a geometry and connectivity relationships
         *   - high-level elements (CONTACT, INTERFACE, LAYER) that are constuted of low-level elements
         *
         * DO NOT MODIFY THIS ENUM 
         */
        enum TYPE {
            CORNER = 0,
            LINE,
            SURFACE,
            REGION,         
            CONTACT,        
            INTERFACE,      
            LAYER,          
            NO_TYPE,        
            ALL_TYPES       
        } ;

        const static index_t NO_ID = index_t( -1 ) ;

        static GEOL_FEATURE determine_geological_type( const std::string& in ) ;
        static GEOL_FEATURE determine_type( const std::vector< GEOL_FEATURE >& types ) ;
        static std::string geol_name( GEOL_FEATURE ) ;
        static std::string type_name( TYPE t ) ;
       
        // Key functions - They determine which element of which type
        // can fill the different class attributes
        static TYPE parent_type      ( TYPE t ) ;
        static TYPE child_type       ( TYPE t ) ;
        static TYPE boundary_type    ( TYPE t ) ;
        static TYPE in_boundary_type ( TYPE t ) ;      
        static index_t dimension     ( TYPE t ) ;       

        static bool parent_allowed      ( TYPE t ) { return parent_type     (t)!= NO_TYPE ; }
        static bool child_allowed       ( TYPE t ) { return child_type      (t)!= NO_TYPE ; }
        static bool boundary_allowed    ( TYPE t ) { return boundary_type   (t)!= NO_TYPE ; }
        static bool in_boundary_allowed ( TYPE t ) { return in_boundary_type(t)!= NO_TYPE ; } 
        /*!
         * @brief Constructs a BoundaryModelElement
         * 
         * @param[in] model Pointer to the parent model.
         * @param[in] element_type Type of the element to create
         * @param[in] id Index of the element in the corresponding vector in the model
         */
        BoundaryModelElement(
            BoundaryModel* model = NULL, 
            TYPE element_type = NO_TYPE, 
            index_t id = NO_ID 
            )
            : model_( model ), type_( element_type ), id_( id ),
              name_( "" ), geol_feature_( NO_GEOL ), parent_( NO_ID )
        {
        }
        
        virtual ~BoundaryModelElement() { }

        bool operator==( const BoundaryModelElement& rhs ) const ;

        /**
         * \name Accessors to basic information
         * @{
         */
        bool has_model() const { return model_ != NULL ; }
        const BoundaryModel& model() const { return *model_ ; } 
        bool has_name() const { return name_ != "" ; }
        const std::string& name() const { return name_ ; }
        bool has_id() const { return id_ != NO_ID ; } 
        index_t id() const { return id_ ; }
        bool has_type() const { return type_ != NO_TYPE ; }
        TYPE element_type() const { return type_ ; } 
        bool has_geological_feature() const { return geol_feature_ != NO_GEOL ; }
        GEOL_FEATURE geological_feature() const { return geol_feature_ ; }
        bool is_on_voi() const ;

         /**@}
         * \name Connectivity - boundary and in_boundary
         * @{
         */
        index_t nb_boundaries() const { return boundaries_.size() ; }
        index_t boundary_id( index_t x ) const { return boundaries_[x] ; }
        const BoundaryModelElement& boundary( index_t x ) const ;     
        bool side( index_t i ) const { return sides_[i] ; }

        index_t nb_in_boundary() const { return in_boundary_.size() ; }
        index_t in_boundary_id( index_t x ) const { return in_boundary_[x] ; }
        const BoundaryModelElement& in_boundary( index_t x ) const ;
        
        /**@}
         * \name Parent - children relationships
         * @{
         */
        bool has_parent() const { return parent_ != NO_ID ; }
        const BoundaryModelElement& parent() const ;
        index_t parent_id() const { return parent_ ; }

        index_t nb_children() const { return children_.size() ; }
        index_t child_id( index_t x ) const { return children_[x] ; }
        const BoundaryModelElement& child( index_t x ) const ;
        
        
        /**@}
         * \name Accessors to geometry - Reimplemented in Corner, Line, Surface classes
         * @{
         */
        virtual index_t nb_cells() const {
            return 0 ;
        }
        virtual index_t nb_vertices() const {
            return 0 ;
        }      
        virtual index_t model_vertex_id( index_t p = 0 ) const {
            return NO_ID ;
        }  
        virtual const vec3& vertex( index_t p = 0 ) const {
            grgmesh_assert_not_reached ; return dummy_vec3 ;
        }
        virtual void set_vertex( index_t index, index_t model_vertex_id ) {
            return ;
        }
                     
         /**@}
         * \name Accessors to attribute managers
         * @{
         */
        VertexAttributeManager* vertex_attribute_manager() const
        {
            return const_cast< VertexAttributeManager* >( &vertex_attribute_manager_ ) ;
        }
        FacetAttributeManager* facet_attribute_manager() const
        {
            return const_cast< FacetAttributeManager* >( &facet_attribute_manager_ ) ;
        }

        /**@}
         * \name Modification of the element
         * @{
         */
        void copy_macro_topology( const BoundaryModelElement& rhs, BoundaryModel& model ) ;
        
        void set_model( BoundaryModel* m ){ model_ = m  ; } 
        void set_element_type( TYPE t ) { type_ = t ; }
        void set_id( index_t id ) { id_ = id ; }
        void set_name( const std::string& name ) { name_ = name ; }
        void set_geological_feature( GEOL_FEATURE type ) { geol_feature_ = type ; } 
      
        void add_boundary( index_t b ) { 
            grgmesh_assert( boundary_allowed( type_ ) ) ;
            boundaries_.push_back( b ) ; 
        }
        void set_boundary( index_t id, index_t b ) { 
            grgmesh_assert( id < nb_boundaries() ) ;
            boundaries_[id] = b ; 
        }
        void add_boundary( index_t b, bool side ) {
            grgmesh_assert( boundary_allowed( type_ ) ) ;
            boundaries_.push_back(b) ;
            sides_.push_back(side) ; 
        }
        void set_boundary( index_t id, index_t b, bool side ) {
            grgmesh_assert( id < nb_boundaries() && id < sides_.size() ) ;
            boundaries_[id] = b ; 
            sides_[id] = side ; 
        }        
        void add_in_boundary( index_t e ) { 
            grgmesh_assert( in_boundary_allowed( type_ ) ) ;
            in_boundary_.push_back(e) ; 
        }
        void set_in_boundary( index_t id, index_t in_b ) { 
            grgmesh_assert( id < nb_in_boundary() ) ;
            in_boundary_[id] = in_b ; 
        }
        void set_parent( index_t p ){
            grgmesh_assert( parent_allowed( type_ ) ) ;
            parent_ = p ; 
        }       
        void add_child( index_t e ){ 
            grgmesh_assert( child_allowed( type_ ) ) ;
            children_.push_back( e ) ; 
        }
        void set_child( index_t id, index_t c ) {
            grgmesh_assert( id < nb_children() ) ;
            children_[id] = c ;
        }

        /**
         * @}
         */

    protected :
        /// Pointer to the BounadyModel owning this element
        BoundaryModel* model_ ;

        /// Type of the element
        TYPE type_ ;
       
        /// Elements are uniquely identified in the BoundaryModel
        /// the pair TYPE + index
        index_t id_ ;        
       
        /// Name of the element - by default it is an empty string
        std::string name_ ;

        /// Geological feature of this object - default is NO_GEOL 
        GEOL_FEATURE geol_feature_ ;

        /// Elements on the boundary of this element - see boundary_type( TYPE )
        std::vector< index_t > boundaries_ ;

        /// Additional information for oriented boundaries
        /// Side: + (true) or - (false)
        std::vector< bool > sides_ ; 
        
        /// Elements in which boundary this element is - see in_boundary_type( TYPE )
        std::vector< index_t > in_boundary_ ;

        /// Index of the parent - see parent_type( TYPE ) - default value is NO_ID.
        index_t parent_ ;

        /// Elements constituting this one - see child_type( TYPE )
        std::vector< index_t > children_ ;

        // Attribute managers
        VertexAttributeManager vertex_attribute_manager_ ;
        FacetAttributeManager facet_attribute_manager_ ;
    } ;

    /// Element to return when a method failed - to avoid compilation warnings
    const static BoundaryModelElement dummy_element = BoundaryModelElement( nil, BoundaryModelElement::NO_TYPE ) ;

    /*! 
    * @brief A Corner
    * 
    * Element of type CORNER. Its geometry is determined by one vertex.
    * Most corners are at the intersections of at least two Line, but some
    * are in the boundary of a closed Line. 
    */
    class GRGMESH_API Corner : public BoundaryModelElement {       
    public:
        Corner(
            BoundaryModel* model,
            index_t id = NO_ID,
            index_t model_vertex_id = NO_ID )
            : BoundaryModelElement( model, CORNER, id ), vertex_( model_vertex_id )
        {
        }
        virtual ~Corner() { } ;

        virtual index_t nb_cells() const { return 0 ; }
        virtual index_t nb_vertices() const { return 1 ; }
        virtual index_t model_vertex_id( index_t id = 0 ) const { return vertex_ ; } 
        virtual const vec3& vertex( index_t p = 0 ) const ;        
        virtual void set_vertex( index_t toto, index_t model_vertex_id ) { vertex_ = model_vertex_id ; }
        
        void set_vertex( index_t model_vertex_id ) { vertex_ = model_vertex_id ; }

    private:
        index_t vertex_ ;
    };


    /*! 
     * @brief A boundary Line of a Surface
     * 
     * It has 1 or 2 Corners on its boundary and is in the boundary of a least one Surface.
     * It is defined by a set of vertices. Its segments are implicitely defined between vertices
     * vertex(n) and vertex(n+1) for n between 0 and nb_cells()
     *
     * There is no LineMutator since hardly nothing can be performed on a Line without modifying the model
     */
    class GRGMESH_API Line: public BoundaryModelElement {
    public:
        Line( BoundaryModel* model, index_t id = NO_ID ) ;
        Line(
            BoundaryModel* model,
            index_t id,
            const std::vector< index_t >& vertices ) ;
        Line(
            BoundaryModel* model,
            index_t id,
            index_t corner0,
            index_t corner1,
            const std::vector< index_t >& vertices ) ;

        virtual ~Line(){} ;
   
        /*! @brief Returns the number of segments of the Line */       
        virtual index_t nb_cells() const { return vertices_.size()-1 ; }
        virtual index_t nb_vertices() const { return vertices_.size() ; }
        virtual index_t model_vertex_id( index_t p ) const { return vertices_.at(p) ; }
        virtual const vec3& vertex( index_t line_vertex_id ) const ;
        virtual void set_vertex( index_t index, index_t model_vertex_id ) {
            grgmesh_assert( index < nb_vertices() ) ;
            vertices_[index] = model_vertex_id ;
        }

        /*! @brief A Line is closed if its two extremities are identitcal */
        bool is_closed () const {
            grgmesh_assert( nb_boundaries() == 2 ) ;
            return ( boundaries_[0] != NO_ID ) && ( boundaries_[0] == boundaries_[1] ) ; 
        }  
        bool is_inside_border( const BoundaryModelElement& e ) const ;
        bool equal( const std::vector< index_t >& rhs_vertices ) const ;
                
        void set_vertices( const std::vector< index_t >& model_vertex_ids ) { 
            vertices_.resize(0) ;
            vertices_.insert( vertices_.begin(), model_vertex_ids.begin(), model_vertex_ids.end() ) ;
        }
        

        vec3 segment_barycenter( index_t s ) const ;
        double segment_length( index_t s ) const ;
        double total_length() const ;            

    private:
        /// Indices of the model vertices in the line
        /// If the line is closed, the last vertex is equal to the first.
        std::vector< index_t > vertices_ ;                
    } ;


    /*! 
     * @brief Attribute on the vertices of a Line
     */
    template< class ATTRIBUTE >
    class LineVertexAttribute: public Attribute< BoundaryModelElement::VERTEX, ATTRIBUTE > {
    public:
        typedef Attribute< BoundaryModelElement::VERTEX, ATTRIBUTE > superclass ;

        void bind( const Line* line, const std::string& name )
        {
            superclass::bind( line->vertex_attribute_manager(), line->nb_vertices(),
                name ) ;
        }

        void bind( const Line* line )
        {
            superclass::bind( line->vertex_attribute_manager(),
                line->nb_vertices() ) ;
        }

        LineVertexAttribute()
        {
        }

        LineVertexAttribute( const Line* line )
        {
            bind( line ) ;
        }

        LineVertexAttribute( const Line* line, const std::string& name )
        {
            bind( line, name ) ;
        }

        static bool is_defined( const Line* line, const std::string& name )
        {
            return superclass::is_defined( line->vertex_attribute_manager(), name ) ;
        }
    } ;

    /*! 
     * @brief Attribute on the segments of a Line
     */
    template< class ATTRIBUTE >
    class LineFacetAttribute: public Attribute< BoundaryModelElement::FACET, ATTRIBUTE > {
    public:
        typedef Attribute< BoundaryModelElement::FACET, ATTRIBUTE > superclass ;

        void bind( const Line* line, const std::string& name )
        {
            superclass::bind( line->facet_attribute_manager(), line->nb_vertices(),
                name ) ;
        }

        void bind( const Line* line )
        {
            superclass::bind( line->facet_attribute_manager(),
                line->nb_vertices() ) ;
        }

        LineFacetAttribute()
        {
        }

        LineFacetAttribute( const Line* line )
        {
            bind( line ) ;
        }

        LineFacetAttribute( const Line* line, const std::string& name )
        {
            bind( line, name ) ;
        }

        static bool is_defined( const Line* line, const std::string& name )
        {
            return superclass::is_defined( line->facet_attribute_manager(), name ) ;
        }
    } ;

    
    /*!
    * @brief A polygonal manifold surface
    * 
    * This is a BoundaryModelElement of type SURFACE.
    * It is defined by a set of vertices and a set of polygonal facets.
    * Its boundaries are several Lines and it is on the boundary of 1 or 2 Region
    */
    class GRGMESH_API Surface : public BoundaryModelElement {
        friend class SurfaceMutator ;
    public:
        const static index_t NO_ADJACENT = index_t( -1 ) ;

        Surface( BoundaryModel* model, index_t id = NO_ID )
            : BoundaryModelElement( model, SURFACE, id ),
              is_triangulated_( false )
        {
        }
        virtual ~Surface(){} ;
  
        bool is_triangulated() const { return is_triangulated_ ; }
               
        /*!
         * @brief Returns the number of facets 
         */         
        virtual index_t nb_cells() const { return facets_.empty() ? 0 : facet_ptr_.size() - 1 ; }
        /*! 
         * @brief Returns the number of vertices 
         */
        virtual index_t nb_vertices() const { return vertices_.size() ; }
        /*!
         * @brief Get the vertex in the model from a vertex index in the Surface 
         */
        virtual index_t model_vertex_id( index_t p ) const { return vertices_[p] ; }
        /*!
         * @brief Returns the coordinates of the point at the given index in the surface 
         */
        virtual const vec3& vertex( index_t surf_vertex_id ) const ;
        /*! 
         * @brief Returns the coordinates of point \param v in facet \param f 
         */
        const vec3& vertex( index_t f, index_t v ) const ;
        
        virtual void set_vertex( index_t index, index_t new_model_index ) {
            grgmesh_assert( index < nb_vertices() ) ;
            vertices_[index] = new_model_index ;
        }

        /**
         * \name Accessors to facet and vertices
         * @{
         */
        index_t facet_begin( index_t f ) const { return facet_ptr_.at(f) ; }
        index_t facet_end( index_t f ) const { return facet_ptr_.at(f+1) ; }
        index_t nb_vertices_in_facet( index_t f ) const { return facet_end( f ) - facet_begin( f ) ; }
        bool is_triangle( index_t f ) const { return nb_vertices_in_facet( f ) == 3 ; }

        index_t next_in_facet( index_t f, index_t v ) const { 
            grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
            if( v != nb_vertices_in_facet(f)-1 ) return v+1 ;
            else return 0 ;
        }
        index_t prev_in_facet( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
            if( v > 0 ) return v-1 ;
            else return nb_vertices_in_facet(f)-1 ;
        }
        index_t nb_corners() const { return facets_.size() ; }
        index_t model_vertex_id_at_corner( index_t corner ) const { return vertices_[ facets_[corner] ]; }
        
        /*!
         * @brief Convert the facet index in the surface to a facet index in the BoundaryModel 
         */
        index_t model_facet_id( index_t f ) const ;        
        /*! 
         *@brief Returns the surface index of vertex \param v in facet \param f
         */
        index_t surf_vertex_id( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
            return facets_[facet_begin(f)+v] ; 
        }
        /*! 
         * @brief Returns the index of vertex \param v in facet \param f in the parent BoundaryModel
         */ 
        index_t model_vertex_id( index_t f, index_t v ) const { 
            grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
            return vertices_[ surf_vertex_id( f, v ) ] ;
        }     
        /*! 
         * @brief Returns a vertex surface index from its model index 
         * @details Returns the first one only or NO_ID if no point is found 
         */
        index_t surf_vertex_id( index_t model_vertex_id ) const {
            for( index_t i = 0; i < vertices_.size() ; ++i ){
                if ( vertices_[i] == model_vertex_id ) return i ;
            }
            return NO_ID ;
        }
        
        index_t facet_vertex_id( index_t t, index_t surf_vertex_id ) const ;  
        index_t facet_id_from_model( index_t f, index_t model_vertex_id ) const ;
        index_t facet_from_surface_vertex_ids( index_t in0, index_t in1 ) const ;
        index_t facet_from_model_vertex_ids( index_t i0, index_t i1 ) const ;    
        void edge_from_model_vertex_ids(
            index_t i0, index_t i1,
            index_t& f, index_t& e ) const ;
        void oriented_edge_from_model_vertex_ids(
            index_t i0, index_t i1,
            index_t& facet, index_t& edge ) const ;
       
        index_t facets_around_vertex(
            index_t surf_vertex_id, 
            std::vector< index_t >& result, 
            bool border_only ) const ;             
        
        index_t facets_around_vertex(
            index_t surf_vertex_id,
            std::vector< index_t >& result,
            bool border_only,
            index_t first_facet ) const ;      
        /** @}
         * \name Geometrical request on facets
         * @{
         */
        vec3 facet_barycenter( index_t f ) const ;
        double facet_area( index_t f ) const ;
        vec3 facet_normal( index_t f ) const ;
        void vertex_normals( std::vector< vec3 >& normals ) const ;
        index_t closest_vertex_in_facet( index_t f, const vec3& vertex ) const ;
        vec3 edge_barycenter( index_t c ) const ;

        /** @}
         * \name Adjacencies request
         * @{
         */
        /*! @brief Returns the index of the adjacent facet of \param f in this surface 
         *  along the edge starting at \param v */
        index_t adjacent( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
            return adjacent_[facet_begin(f)+v] ; 
        }
        /*! @brief Retruns the index of the adjacent facet at the given corner
        */
        index_t adjacent( index_t c ) const {
            grgmesh_assert( c < adjacent_.size() ) ;
            return adjacent_[c] ; 
        }

        bool is_on_border( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
            return adjacent( f, v ) == NO_ADJACENT ; 
        }
        bool is_on_border( index_t f ) const {
            for( index_t adj = 0; adj < nb_vertices_in_facet(f); adj++ ) {
                if( is_on_border( f, adj ) ) return true ;
            }
            return false ;
        }
              
        void next_on_border(
            index_t f, index_t from, index_t v,
            index_t& next_f, index_t& v_in_next, index_t& to = dummy_index_t ) const ;
        
        void next_on_border(
            index_t f, index_t e,
            index_t& next_f, index_t& next_e ) const ;       


         /** @}
         * \name Modifiers
         * @{
         */
        void set_adjacent( index_t f, index_t e, index_t adjacent ) {
            adjacent_[facet_begin(f)+e] = adjacent ;
        }
        void compute_is_triangulated() {
            for( index_t f = 0; f < nb_cells(); f++ ) {
                if( !is_triangle( f ) ) {
                    is_triangulated_ = false ;
                    return ;
                }
            }
            is_triangulated_ = true ;
        }               
        void set_geometry(
            const std::vector< index_t >& vertices,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) {
            vertices_ = vertices ;
            facets_ = facets ;
            facet_ptr_ = facet_ptr ;
            compute_is_triangulated() ;
        }
        void set_geometry(
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) ;

        void set_adjacent( const std::vector< index_t >& adjacent ){
            grgmesh_assert( adjacent.size() == facets_.size() ) ;
            adjacent_ = adjacent ;
        }
        
      
        /**
         * @}
         */
       
    private:
        /// Indices (in the BoundaryModel) of the vertices defining the surface
        /// The same index can appear several times when there is an Line boundary inside          
        std::vector< index_t > vertices_ ;
        /// Indices (in the Surface) of each vertex in each facet
        std::vector< index_t > facets_ ;
        /// Beginning (and end) of one facet in the facets_ vector
        std::vector< index_t > facet_ptr_ ;
            
        /// Adjacent facet for each vertex in each facet along 
        /// the edge starting at this vertex.
        /// When the edge is along a Surface boundary it is set at NO_ADJACENT 
        std::vector< index_t > adjacent_ ;

        bool is_triangulated_ ;
    };  

    /*! 
     * @brief Attribute on the vertices of a Surface
     */
    template< class ATTRIBUTE >
    class SurfaceVertexAttribute: public Attribute< BoundaryModelElement::VERTEX, ATTRIBUTE > {
    public:
        typedef Attribute< BoundaryModelElement::VERTEX, ATTRIBUTE > superclass ;

        void bind( const Surface* surface, const std::string& name )
        {
            superclass::bind( surface->vertex_attribute_manager(), surface->nb_vertices(),
                name ) ;
        }

        void bind( const Surface* surface )
        {
            superclass::bind( surface->vertex_attribute_manager(),
                surface->nb_vertices() ) ;
        }

        SurfaceVertexAttribute()
        {
        }

        SurfaceVertexAttribute( const Surface* surface )
        {
            bind( surface ) ;
        }

        SurfaceVertexAttribute( const Surface* surface, const std::string& name )
        {
            bind( surface, name ) ;
        }

        static bool is_defined( const Surface* surface, const std::string& name )
        {
            return superclass::is_defined( surface->vertex_attribute_manager(), name ) ;
        }
    } ;
    
    /*! 
     * @brief Attribute on the facets of a Surface
     */
    template< class ATTRIBUTE >
    class SurfaceFacetAttribute: public Attribute< BoundaryModelElement::FACET, ATTRIBUTE > {
    public:
        typedef Attribute< BoundaryModelElement::FACET, ATTRIBUTE > superclass ;

        void bind( const Surface* surface, const std::string& name )
        {
            superclass::bind( surface->facet_attribute_manager(), surface->nb_cells(),
                name ) ;
        }

        void bind( const Surface* surface )
        {
            superclass::bind( surface->facet_attribute_manager(),
                surface->nb_cells() ) ;
        }

        SurfaceFacetAttribute()
        {
        }

        SurfaceFacetAttribute( const Surface* surface )
        {
            bind( surface ) ;
        }

        SurfaceFacetAttribute( const Surface* surface, const std::string& name )
        {
            bind( surface, name ) ;
        }

        static bool is_defined( const Surface* surface, const std::string& name )
        {
            return superclass::is_defined( surface->facet_attribute_manager(), name ) ;
        }
    } ;

    /*! 
     * @brief Class to perform modifications of a Surface
     */
    class GRGMESH_API SurfaceMutator {
    public:
        SurfaceMutator( Surface& S )
            : S_( S )
        {
        }
        SurfaceMutator( const Surface& S )
            : S_( const_cast< Surface& >( S ) )
        {
        }

        std::vector< index_t >& vertices() const  { return S_.vertices_  ; }
        std::vector< index_t >& facets() const    { return S_.facets_    ; }
        std::vector< index_t >& facet_ptr() const { return S_.facet_ptr_ ; }
        std::vector< index_t >& adjacents() const { return S_.adjacent_  ; }

        void clear()
        {
            S_.vertices_.clear() ;
            S_.facets_.clear() ;
            S_.facet_ptr_.clear() ;
            S_.adjacent_.clear() ;
        }
        
        void cut_by_line( const Line& L ) ;

    private:
        Surface& S_ ;
    };


    /*!
     * @brief Class to answer geometrical requests on BoundaryModelElement
     */
    class GRGMESH_API BoundaryModelElementMeasure {
    public:
        static double size( const BoundaryModelElement* E ) ;
        static double cell_size( const BoundaryModelElement* E, index_t cell ) ;
        static double distance( const BoundaryModelElement* from,  const vec3& p ) ;
        static double distance( const BoundaryModelElement* from, const BoundaryModelElement* to ) ;
        static vec3 barycenter ( const BoundaryModelElement* E ) ;
        static vec3 barycenter ( const BoundaryModelElement* E, const std::vector< index_t >& cells ) ;
    } ;

} // namespace

#endif

