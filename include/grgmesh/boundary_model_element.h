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
     * @brief Types for BoundaryModelElement 
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

    /*!
     * \brief Generic class describing one element of a BoundaryModel
     */
    class GRGMESH_API BoundaryModelElement {
        friend class BoundaryModelBuilder ;

    public:
        enum AttributeLocation {
            VERTEX, FACET
        } ;
        typedef AttributeManager< VERTEX > VertexAttributeManager ;
        typedef AttributeManager< FACET > FacetAttributeManager ;
        
        const static index_t NO_ID = index_t( -1 ) ;

        /*!
         * @brief Constructs a BoundaryModelElement
         * 
         * @param[in] model Pointer to the parent model
         * @param[in] dim Dimension of the element 0 for Corners, 1 for Line and Contact,
         *            2 for Surface and Interface, 3 for Region and Layer
         * @param[in] id Index of the element in the corresponding vector in the model
         * @param[in] parent Index of the parent (a Contact, Interface, or Layer) of the element
         * @param[in] type Geological type of the element
         */
        BoundaryModelElement(
            BoundaryModel* model = NULL, 
            index_t dim = NO_ID, 
            index_t id = NO_ID, 
            index_t parent = NO_ID,
            GEOL_FEATURE type = default_type )
            : model_( model ),  name_( "" ), id_( id ),
            dim_( dim ), type_( type ), parent_( parent )
        {
        }
        
        virtual ~BoundaryModelElement() { }

        const BoundaryModel& model() const { return *model_ ; }
        
        // Access to fundamental information
        const std::string& name() const { return name_ ; }
        index_t id() const { return id_ ; }
        index_t dim() const { return dim_ ; }
        GEOL_FEATURE type() const { return type_ ; }
        bool is_on_voi() const ;

        // Parent - only for Line, Surface, and Region
        bool has_parent() const { return parent_ != NO_ID ; }
        const BoundaryModelElement& parent() const ;
        index_t parent_id() const { return parent_ ; }
        
        // Element on the boundary - invalid for Corner
        index_t nb_boundaries() const { return boundaries_.size() ; }
        index_t boundary_id( index_t x ) const { return boundaries_[x] ; }
        const BoundaryModelElement& boundary( index_t x ) const ;
        /*! On which sie of the boundary are we ? Filled for region */
        bool side( index_t i ) const { return sides_[i] ; }

        // Elements containing this element
        index_t nb_in_boundary() const { return in_boundary_.size() ; }
        index_t in_boundary_id( index_t x ) const { return in_boundary_[x] ; }
        const BoundaryModelElement& in_boundary( index_t x ) const ;
        
        // Children - only for Contact, Interface, and Layer
        index_t nb_children() const { return children_.size() ; }
        index_t child_id( index_t x ) const { return children_[x] ; }
        const BoundaryModelElement& child( index_t x ) const ;
        
        // Accessors to the vertices (indices and coordinates) and cells
        // Only valid for Corner, Line, Surface
        virtual index_t nb_cells() const {
            grgmesh_assert_not_reached ;  return 0 ;
        }
        virtual index_t nb_vertices() const {
            grgmesh_assert_not_reached ;  return 0 ;
        }      
        virtual index_t model_vertex_id( index_t p = 0 ) const {
            grgmesh_assert_not_reached ; return 0 ;
        }  
        virtual const vec3& vertex( index_t p = 0 ) const {
            grgmesh_assert_not_reached ; return dummy_vec3 ;
        }
                     
        // Accessors to attribute managers
        VertexAttributeManager* vertex_attribute_manager() const
        {
            return const_cast< VertexAttributeManager* >( &vertex_attribute_manager_ ) ;
        }
        FacetAttributeManager* facet_attribute_manager() const
        {
            return const_cast< FacetAttributeManager* >( &facet_attribute_manager_ ) ;
        }

    protected:        
        void copy_macro_topology(
            const BoundaryModelElement& rhs, BoundaryModel& model ) ;
        
        void set_parent( index_t p ){ parent_ = p ; }
        void set_name( const std::string& name ) { name_ = name ; }
        void set_type( GEOL_FEATURE type ) { type_ = type ; } 
        void set_dim( index_t dim ) { dim_ = dim ; }
        void set_id( index_t id ) { id_ = id ; }
        
        void add_boundary( index_t b ) { boundaries_.push_back( b ) ; }
        void add_boundary( index_t b, bool side ) { boundaries_.push_back(b) ; sides_.push_back(side) ; }
        virtual void add_in_boundary( index_t e ) { in_boundary_.push_back(e) ; }
        void add_child( index_t e ){ children_.push_back( e ) ; }
       

    protected :
        /// \todo Can we have something else than a POINTER ?? default constructor is needed...
        BoundaryModel* model_ ;

        /// Name of the element, empty string if none nothing
        std::string name_ ;

        /// Id of this element in the appropriate vector of the BoundaryModel owning it
        index_t id_ ;
        
        /// Dimension of the element 0 corner; 1 line; 2 surface; 3 region
        index_t dim_ ;
        
        /// Geological type for this object, default is ALL 
        GEOL_FEATURE type_ ;

        /// Elements on the boundary of this element - their dimension is dim_-1
        std::vector< index_t > boundaries_ ;

        /// Flag on which side of the boundary is this element 
        /// Filled for volumetric regions only + (true) or - (false)
        std::vector< bool > sides_ ; 
        
        /// Elements in which boundary this element is - their dimension is dim_+1
        std::vector< index_t > in_boundary_ ;

        /// Index of the parent (group of elements to which belong this).
        /// Default value is NO_ID.
        index_t parent_ ;

        /// The group elements making up this one, empty for basic elements
        std::vector< index_t > children_ ;

        // Attribute managers
        VertexAttributeManager vertex_attribute_manager_ ;
        FacetAttributeManager facet_attribute_manager_ ;
    } ;

    const static BoundaryModelElement dummy_element = BoundaryModelElement( nil, 0 ) ;

    /*! 
    * @brief A Corner - point at the intersection of at least 2 Lines 
    * 
    * Element of dimension 0. Its geometry is determined by one vertex.
    */
    class GRGMESH_API Corner : public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
    public:
        Corner(
            BoundaryModel* model,
            index_t id = NO_ID,
            index_t p_id = 0 )
            : BoundaryModelElement( model, 0, id ), vertex_( p_id )
        {
        }
        virtual ~Corner() { } ;
        /*!
         * @brief A corner in only one (so closed) Line is not real
         */
        bool is_real() const { return in_boundary_.size() > 1 ; }
        
        virtual index_t nb_cells() const { return 1 ; }
        virtual index_t nb_vertices() const { return 1 ; }
        virtual index_t model_vertex_id( index_t id = 0 ) const { return vertex_ ; } 
        virtual const vec3& vertex( index_t p = 0 ) const ;
        
    private:
        void set_vertex( index_t vertex ) { vertex_ = vertex ; }

    private:
        index_t vertex_ ;
    };


    /*! 
     * @brief A line of intersection between two Surfaces
     * 
     * It has 1 or 2 Corners ont its boundary
     * It is in the boundary of several Surface
     */
    class GRGMESH_API Line: public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
        friend class LineMutator ;
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
   
        virtual index_t nb_cells() const {            
            return vertices_.size()-1 ; 
        }
        virtual index_t nb_vertices() const {             
            return vertices_.size() ; 
        }
        virtual index_t model_vertex_id( index_t p ) const {
            return vertices_.at(p) ;
        }
        virtual const vec3& vertex( index_t line_vertex_id ) const ;

        bool is_closed () const { 
            return (boundaries_[0]!= nil ) && (boundaries_[0] == boundaries_[1]) ; 
        }  
        bool is_inside_border( const BoundaryModelElement& e ) const ;
            
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
    class LineVertexAttribute: public Attribute< Line::VERTEX, ATTRIBUTE > {
    public:
        typedef Attribute< Line::VERTEX, ATTRIBUTE > superclass ;

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
    class LineFacetAttribute: public Attribute< Line::FACET, ATTRIBUTE > {
    public:
        typedef Attribute< Line::FACET, ATTRIBUTE > superclass ;

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
     * @brief Class to perform modifications of a Line
     */
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
        void set_vertex( index_t id, const vec3& p ) ;
        vec3& vertex( index_t p ) const ;
        std::vector< index_t >& vertices() const { return M_.vertices_ ; }

        void clear()
        {
            M_.vertices_.clear() ;
        }

    private:
        Line& M_ ;
    };


    /*!
    * @brief A polygonal manifold surface
    * 
    * This is a BoundaryModelElement of dimension 2.
    * Its boundaries are several Lines and it is on the boundary of 1 or 2 Region
    */
    class GRGMESH_API Surface : public BoundaryModelElement {
        friend class BoundaryModelBuilder ;
        friend class SurfaceMutator ;
    public:
        const static index_t NO_ADJACENT = index_t( -1 ) ;

        /*!
         * @brief The triangle that set the orientation of a Surface's facets
         */
        struct KeyFacet {
            KeyFacet( const vec3& p0, const vec3& p1, const vec3& p2 ):
                p0_(p0), p1_(p1), p2_(p2){ } ;

            KeyFacet() {  
                vec3 d( -1, -1, -1) ;
                p0_ = d ; p1_ = d; p2_ = d ; 
            } 
    
            bool is_default() const {
                vec3 d( -1, -1, -1) ;
                if( p0_== d && p1_== d && p2_== d  ) return true ;
                else {
                    grgmesh_assert( p0_!=p1_ && p0_!=p2_ && p1_!=p2_ ) ;
                    return false ;
                }
            }
        public:
            vec3 p0_ ;
            vec3 p1_ ;
            vec3 p2_ ;
        } ;


        Surface(
            BoundaryModel* model,
            index_t id = NO_ID,
            index_t parent = NO_ID,
            const GEOL_FEATURE& type = default_type )
            :
            BoundaryModelElement( model, 2, id, parent, type ),
            is_triangulated_( false )
        {
        }
        virtual ~Surface(){} ;
        
        const KeyFacet& key_facet() const { return key_facet_ ; }       
        bool is_triangulated() const { return is_triangulated_ ; }
               
        /*! @brief Returns the number of facets */         
        virtual index_t nb_cells() const { return facets_.empty() ? 0 : facet_ptr_.size() - 1 ; }
        /*! @brief Returns the number of vertices */
        virtual index_t nb_vertices() const { return vertices_.size() ; }               
        /*! @brief Get the vertex in the model from a vertex index in the Surface */
        virtual index_t model_vertex_id( index_t p ) const { return vertices_[p] ; }
        /*! @brief Returns the coordinates of the point at the given index in the surface */
        virtual const vec3& vertex( index_t surf_vertex_id ) const ;
         /*! @brief Returns the coordinates of point \param v in facet \param f */
        const vec3& vertex( index_t f, index_t v ) const ;
        
        // Facet information
        index_t facet_begin( index_t f ) const { return facet_ptr_.at(f) ; }
        index_t facet_end( index_t f ) const { return facet_ptr_.at(f+1) ; }
        index_t nb_vertices_in_facet( index_t f ) const { return facet_end( f ) - facet_begin( f ) ; }
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
        bool is_triangle( index_t f ) const { return nb_vertices_in_facet( f ) == 3 ; }
        
        /*! @brief Convert the facet index in the surface to a facet index in the BoundaryModel */
        index_t model_facet_id( index_t f ) const ;
        
        /*! @brief Returns the surface index of vertex \param v in facet \param f */
        index_t surf_vertex_id( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
            return facets_[facet_begin(f)+v] ; 
        }
        /*! @brief Returns the index of vertex \param v in facet \param f in the parent BoundaryModel */ 
        index_t model_vertex_id( index_t f, index_t v ) const { 
            grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
            return vertices_[ surf_vertex_id( f, v ) ] ;
        }     
        /*! @brief Returns a vertex surface index from its model index 
         *  @details Returns only the first. Returns NO_ID if no corrsponding point
         *  is found on the surface */
        index_t surf_vertex_id( index_t model_vertex_id ) const {
            for( index_t i = 0; i < vertices_.size() ; ++i ){
                if ( vertices_[i] == model_vertex_id ) return i ;
            }
            return NO_ID ;
        }
        
        index_t facet_vertex_id( index_t t, index_t surf_vertex_id ) const ;  
        index_t facet_from_model_vertex_ids( index_t i0, index_t i1 ) const ;        
        void edge_from_model_vertex_ids(
            index_t i0, index_t i1,
            index_t& f, index_t& e ) const ;
        void oriented_edge_from_model_vertex_ids(
            index_t i0, index_t i1,
            index_t& facet, index_t& edge ) const ;
       
        vec3 facet_barycenter( index_t f ) const ;
        double facet_area( index_t f ) const ;
        vec3 facet_normal( index_t f ) const ;
        void vertex_normals( std::vector< vec3 >& normals ) const ;
        index_t closest_vertex_in_facet( index_t f, const vec3& vertex ) const ;

        
        /*! @brief Returns the index of the adjacent facet of \param f in this surface 
         *  along the edge starting at \param v */
        index_t adjacent( index_t f, index_t v ) const {
            grgmesh_debug_assert( v < nb_vertices_in_facet(f) ) ;
            return adjacent_[facet_begin(f)+v] ; 
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

        index_t facets_around_vertex(
            index_t surf_vertex_id, 
            std::vector< index_t >& result, 
            bool border_only ) const ;             
        
        index_t facets_around_vertex(
            index_t surf_vertex_id,
            std::vector< index_t >& result,
            bool border_only,
            index_t first_facet ) const ;      


    private:
        void set_key_facet( const KeyFacet& key ) { key_facet_ = key ; }
        void set_first_triangle_as_key() ;
        
        // Needed at model building
        index_t facet_from_surface_vertex_ids( index_t i0, index_t i1 ) const ;

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
            const std::vector< index_t >& facet_ptr )
        {
            // Are these copies parallelized ?
            vertices_ = vertices ;
            facets_ = facets ;
            facet_ptr_ = facet_ptr ;

            compute_is_triangulated() ;
        }

        void set_adjacent( const std::vector< index_t >& adjacent ){
            grgmesh_assert( adjacent.size() == facets_.size() ) ;
            adjacent_ = adjacent ;
        }
       
    private:
        KeyFacet key_facet_ ;

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
    class SurfaceVertexAttribute: public Attribute< Surface::VERTEX, ATTRIBUTE > {
    public:
        typedef Attribute< Surface::VERTEX, ATTRIBUTE > superclass ;

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
    class SurfaceFacetAttribute: public Attribute< Surface::FACET, ATTRIBUTE > {
    public:
        typedef Attribute< Surface::FACET, ATTRIBUTE > superclass ;

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
     * @brief Class to perform modification of a Surface
     */
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
        void set_vertex( index_t id, const vec3& p ) ;
        vec3& vertex( index_t p ) const ;
        std::vector< index_t >& vertices() const { return M_.vertices_ ; }
        std::vector< index_t >& facets() const { return M_.facets_ ; }
        std::vector< index_t >& facet_ptr() const { return M_.facet_ptr_ ; }
        std::vector< index_t >& adjacents() const { return M_.adjacent_ ; }

        void clear()
        {
            M_.vertices_.clear() ;
            M_.facets_.clear() ;
            M_.facet_ptr_.clear() ;
            M_.adjacent_.clear() ;
        }

    private:
        Surface& M_ ;
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

