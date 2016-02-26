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

/*! \author Jeanne Pellerin and Arnaud Botella */

#ifndef __RINGMESH_GEO_MODEL_ELEMENT__
#define __RINGMESH_GEO_MODEL_ELEMENT__

#include <ringmesh/common.h>

#include <geogram/mesh/mesh.h>

#include <vector>
#include <string>

namespace GEO {
    class MeshFacetsAABB ;
    class MeshCellsAABB ;
}

namespace RINGMesh {
    class GeoModel ;
    class ColocaterANN ;
}

namespace RINGMesh {
   
    /*!
     * @brief Generic class describing one element of a GeoModel
     */
    class RINGMESH_API GeoModelElement {
        ringmesh_disable_copy( GeoModelElement ) ;
        friend class GeoModelEditor ;
    public:
        /*!
         * @brief Geological feature types for GeoModelElement
         * @todo Read all possible geological features used in RESQML.
         */
        enum GEOL_FEATURE {
            /// All geological features 
            ALL_GEOL,
            /// Default value - No geological feature defined
            NO_GEOL,
            /// Stratigraphic surface - an horizon
            STRATI,
            /*!
             * Unconformity
             * \todo: distinguish between as erosive, baselap or intrusive
             * unconformities --GC
             */
            UNCONFORMITY,
            /// A normal fault
            NORMAL_FAULT,
            /// A reverse fault 
            REVERSE_FAULT,
            /// An unspecified fault 
            FAULT,
            /// Volume Of Interest
            VOI
        } ;

        /*!
         * @brief Type of the GeoModelElement
         * @details When no type is defined NO_TYPE should be used
         * There are two categories of elements
         *   - low-level elements (CORNER, LINE, SURFACE, REGION) which have a 
         *     geometry and connectivity relationships
         *   - high-level elements (CONTACT, INTERFACE, LAYER) 
         *     that are constituted of low-level elements
         * TYPE is used extensively to manage elements, iterate on them, etc.
         *  
         * @warning DO NOT MODIFY THIS ENUM.
         * OK but why ? --GC
         * 
         * @todo Add fault blocks.
         * @todo Encapsulate in functions for testing groups or range of
         * features (See previous comment)--GC
         */
        enum TYPE {
            /// Points at LINE extremities
            CORNER = 0,
            /// One connected component of the intersection of at least 2 SURFACE 
            LINE,
            /// One 2-manifold connected component
            SURFACE,
            /// One volumetric region defined by its boundary SURFACE and
            /// is optionally meshed
            REGION,
            /// A group of LINE, intersection of at least 2 INTERFACE
            CONTACT,
            /// A group of SURFACE, delimit maximum 2 LAYER 
            INTERFACE,
            /// A group of REGION 
            LAYER,
            /// Default TYPE
            NO_TYPE,
            /// Any type - to access generically elements
            ALL_TYPES
        } ;

        /*! 
         * @brief Unique identification of a GeoModelElement in a GeoModel
         * @details Stores the TYPE of the element and its index in the GeoModel.
         *          Default values are NO_TYPE and NO_ID
         */
        struct gme_t {
            gme_t()
                : type( NO_TYPE ), index( NO_ID )
            {
            }
            gme_t( TYPE t, index_t id )
                : type( t ), index( id )
            {
            }
            bool operator!=( const gme_t& rhs ) const
            {
                return type != rhs.type || index != rhs.index ;
            }
            bool operator==( const gme_t& rhs ) const
            {
                return type == rhs.type && index == rhs.index ;
            }
            /*!
             * @brief Sort BME identifiers
             * @details Compare first types, then compare indices, 
             *          beginning with NO_ID indices. 
             * @note In a sorted vector v of gme_t one can find the first surface with
             *       std::lower_bound( v.begin(), v.end(), gme_t( SURFACE, NO_ID ) ) ;
             */
            bool operator<( const gme_t& rhs ) const
            {
                if( type != rhs.type ) {
                    return type < rhs.type ;
                } else {
                    if( index == NO_ID ) return true ;
                    if( rhs.index == NO_ID ) return false ;
                    return index < rhs.index ;
                }
            }
            friend std::ostream& operator<<( std::ostream& os, const gme_t& in )
            {
                os << GeoModelElement::type_name( in.type ) << " " << in.index ;
                return os ;
            }

            bool is_defined() const
            {
                return type != NO_TYPE && type != ALL_TYPES && index != NO_ID ;
            }
            /*!
             * TYPE of the GeoModelElement
             * \todo Should be type_ to be consistent with style guidelines --GC
             */
            TYPE type ;
            /*!
             * Index of the element in the GeoModel
             * \todo Should be index_ to be consistent with style guidelines --GC
             */
            index_t index ;
        } ;

        static const index_t NO_ID = index_t( -1 ) ;

        static GEOL_FEATURE determine_geological_type( const std::string& in ) ;
        static GEOL_FEATURE determine_type(
            const std::vector< GEOL_FEATURE >& types ) ;
        static std::string geol_name( GEOL_FEATURE ) ;
        static bool is_fault( GEOL_FEATURE T )
        {
            return T == FAULT || T == REVERSE_FAULT || T == NORMAL_FAULT ;
        }
        static bool is_stratigraphic_limit( GEOL_FEATURE T )
        {
            return T == STRATI || T == UNCONFORMITY ;
        }

        /*!
         * \name Key functions to access relationships between TYPE s 
         * @{
         */
        static std::string type_name( TYPE t ) ;

        static TYPE parent_type( TYPE t ) ;
        static TYPE child_type( TYPE t ) ;
        static TYPE boundary_type( TYPE t ) ;
        static TYPE in_boundary_type( TYPE t ) ;
        static index_t dimension( TYPE t ) ;
        static bool has_mesh( TYPE t ) ;

        static bool parent_allowed( TYPE t )
        {
            return parent_type( t ) != NO_TYPE ;
        }
        static bool child_allowed( TYPE t )
        {
            return child_type( t ) != NO_TYPE ;
        }
        static bool boundary_allowed( TYPE t )
        {
            return boundary_type( t ) != NO_TYPE ;
        }
        static bool in_boundary_allowed( TYPE t )
        {
            return in_boundary_type( t ) != NO_TYPE ;
        }

        /*!@}
         */

        virtual ~GeoModelElement()
        {
        }


        /*!@}
         * \name Validity checks
         * @{
         */

        /*!
         * @brief Global validity check of the BME
         */
        virtual bool is_valid() const
        {
            return is_connectivity_valid() ;
        }

        /*!
         * @brief Basic checks on the minimum required information 
         * @details Required connectivity information depends on the TYPE.   
         */
        bool is_connectivity_valid() const ;

        /*!
         * \name Accessors to basic information
         * @{
         */

        const GeoModel& model() const
        {
            return model_ ;
        }
        bool has_name() const
        {
            return name() != "" ;
        }
        const std::string& name() const
        {
            return name_ ;
        }
        const gme_t& gme_id() const
        {
            return id_ ;
        }
        index_t index() const
        {
            return gme_id().index ;
        }
        TYPE type() const
        {
            return gme_id().type ;
        }
        bool has_geological_feature() const
        {
            return geol_feature_ != NO_GEOL ;
        }
        GEOL_FEATURE geological_feature() const
        {
            return geol_feature_ ;
        }
        bool is_on_voi() const ;

        /*!@}
         * \name Connectivity - boundary and in_boundary
         * @{
         */
        index_t nb_boundaries() const
        {
            return boundaries_.size() ;
        }
        const gme_t& boundary_gme( index_t x ) const
        {
            return boundaries_[x] ;
        }
        const GeoModelElement& boundary( index_t x ) const ;

        index_t nb_in_boundary() const
        {
            return in_boundary_.size() ;
        }
        const gme_t& in_boundary_gme( index_t x ) const
        {
            return in_boundary_[x] ;
        }
        const GeoModelElement& in_boundary( index_t x ) const ;

        bool is_inside_border( const GeoModelElement& e ) const ;
        bool has_inside_border() const ;

        /*!@}
         * \name Parent - children relationships
         * @{
         */
        bool has_parent() const
        {
            return parent_id().is_defined() ;
        }
        const gme_t& parent_id() const
        {
            return parent_ ;
        }
        const GeoModelElement& parent() const ;

        index_t nb_children() const
        {
            return children_.size() ;
        }
        const gme_t& child_id( index_t x ) const
        {
            return children_[x] ;
        }
        const GeoModelElement& child( index_t x ) const ;

    protected:
        /*!
         * @brief Constructs a GeoModelElement
         * Client code should only create GeoModelElements through
         * GeoModelEditor derived classes.
         *
         * @param[in] model Constant reference to the parent model of this element.
         * @param[in] element_type Type of the element to create
         * @param[in] id Index of the element in the corresponding vector in the model
         * @param[in] name Name of the element, empty by default.
         * @param[in] geological_feature Feature of the element, none by default.
         */
        GeoModelElement(
            const GeoModel& model,
            TYPE element_type,
            index_t id,
            const std::string& name = "",
            GEOL_FEATURE geological_feature = NO_GEOL
            )
            :
                model_( model ),
                id_( element_type, id ),
                name_( name ),
                geol_feature_( geological_feature )
        {
        }


    protected:
        /// Reference to the GeoModel owning this element
        const GeoModel& model_ ;

        /// Unique identifier of the GeoModelElement in the model
        gme_t id_ ;

        /// Name of the element - default is an empty string
        std::string name_ ;

        /// Geological feature of this object - default is NO_GEOL
        GEOL_FEATURE geol_feature_ ;

        /// Elements on the boundary of this element - see boundary_type( TYPE )
        std::vector< gme_t > boundaries_ ;

        /// Elements in which boundary this element is - see in_boundary_type( TYPE )
        std::vector< gme_t > in_boundary_ ;

        /// Parent identification - see parent_type( TYPE )
        gme_t parent_ ;

        /// Elements constituting this one - see child_type( TYPE )
        std::vector< gme_t > children_ ;
    } ;

    typedef GeoModelElement GME ;

    /*!
    * @brief Vertex in a GeoModelElement
    */
    struct GMEVertex {
        GMEVertex( GME::gme_t t, index_t vertex_id_in )
            : gme_id( t ), v_id( vertex_id_in )
        {}
        GMEVertex()
            : gme_id(), v_id( NO_ID )
        {}
        bool operator<( const GMEVertex& rhs ) const
        {
            if( gme_id != rhs.gme_id ) {
                return gme_id < rhs.gme_id ;
            } else {
                return v_id < rhs.v_id ;
            }
        }
        bool operator==( const GMEVertex& rhs ) const
        {
            return gme_id == rhs.gme_id && v_id == rhs.v_id ;
        }
        bool is_defined() const
        {
            return gme_id.is_defined() && v_id != NO_ID ;
        }
        /// GeoModelElement id in its GeoModel
        GME::gme_t gme_id ;
        /// Index of the vertex in the GeoModelElement
        index_t v_id ;
    } ;


    /*!
     * @brief Abstract base class for GeoModelElement 
     *        which have a geometrical representation
     * @details This representation is stored as a GEO::Mesh
     */
    class RINGMESH_API GeoModelMeshElement: public GeoModelElement {
    ringmesh_disable_copy( GeoModelMeshElement ) ;
    friend class GeoModelEditor ;
    friend class GeoModelBuilder ;
    public:

        /*!
         * @brief Name of the attribute storing the index of a vertex in the model
         */
        static const std::string model_vertex_id_att_name() ;

        GeoModelMeshElement(
            const GeoModel& model,
            TYPE element_type,
            index_t id,
            const std::string& name = "",
            GEOL_FEATURE geological_feature = NO_GEOL )
            : GeoModelElement( model, element_type, id, name, geological_feature )
        {
            model_vertex_id_.bind( mesh_.vertices.attributes(),
                model_vertex_id_att_name() ) ;
        }

        virtual ~GeoModelMeshElement() ;

        /*!
         * @brief Global validity of the element
         */
        virtual bool is_valid() const
        {
            return is_connectivity_valid() && is_mesh_valid() ;
            /// \todo Test and add the model vertex validity test            
            /// (no time right now JP)
            // are_model_vertex_indices_valid() ;
        }

        /*!
         * @brief Number of mesh elements
         * vertices for Corner, edges for Line,
         * facets for Surface, cells for Region.
         */
        virtual index_t nb_cells() const = 0 ;

        /*! 
         * @brief Correspondance of vertex index \param lv in a mesh element \param me
         * and vertex index in the GMME
         * @todo changer le nom
         */
        virtual index_t gmme_vertex_index( index_t me, index_t lv ) const = 0 ;
        
        /*!
         * @brief Number of vertices of the mesh
         */
        index_t nb_vertices() const
        {
            return mesh_.vertices.nb() ;
        }

        /*!
        * @brief Global index in GeoModelMesh from the
        * the local one
        * @param[in] v Vertex index in the GeoModelMeshElement
        */
        index_t model_vertex_id( index_t v = 0 ) const
        {
            ringmesh_debug_assert( v < nb_vertices() ) ;            
            return model_vertex_id_[ v ] ;
        }

        index_t model_vertex_id( index_t me, index_t lv ) const
        {            
            return model_vertex_id( gmme_vertex_index( me, lv) ) ;
        }

        /*!
        * @brief Coordinates of a GeoModelMeshElement
        * @param[in] v Index of the vertex in the GeoModelMeshElement
        */
        const vec3& vertex( index_t v = 0 ) const
        {
            ringmesh_debug_assert( v < nb_vertices() ) ;
            return mesh_.vertices.point( v ) ;
        }
      
        const vec3& vertex( index_t me, index_t lv ) const
        {
            return vertex( gmme_vertex_index( me, lv ) ) ;
        }

        /*!
        * @brief Index of the first vertex corresponding to the input model index
        * @details Returns NO_ID if no matching point is found.
        *
        * @param model_vertex_id Index of a vertex in GeoModelMeshVertices
        * @todo changer le nom
        */
        index_t gmme_vertex_index_from_model( index_t model_vertex_id ) const ; 

        std::vector<index_t> gme_vertex_indices( index_t model_vertex_id ) const ;

        /*!
         * @}
         * \name Attribute management
         * @{
         */
        GEO::AttributesManager& vertex_attribute_manager() const
        {
            return mesh_.vertices.attributes() ;
        }
        GEO::AttributesManager& facet_attribute_manager() const
        {
            return mesh_.facets.attributes() ;
        }
        GEO::AttributesManager& cell_attribute_manager() const
        {
            return mesh_.cells.attributes() ;
        }
        void bind_attributes() ;
        void unbind_attributes() ;

        /*! @}
         */

        /*!  
         * @brief Open-bar access to the mesh of the element
         * @details ONLY for internal use.
         * @warning DO NOT directly call this function to modify the facets, edges,
         * or vertices of the element.
         */
        GEO::Mesh& mesh() const
        {
            return const_cast< GEO::Mesh& >( mesh_ ) ;
        }

    protected:
        /*!
         * @brief Check if the mesh stored is valid.
         */
        virtual bool is_mesh_valid() const = 0 ;

        /*!
         * @brief Check that model vertex ids are consistent 
         * with info stored by the model vertices
         */
        bool are_model_vertex_indices_valid() const ;

    protected:
        /// Mesh of the element
        GEO::Mesh mesh_ ;
        /*! Attribute on the Mesh vertices storing the index of
         *  the vertex in the the GeoModel owning this element 
         */
        GEO::Attribute< index_t > model_vertex_id_ ;
    } ;

    /*!
     * @brief A GeoModelElement of type CORNER
     *
     * @details It is a unique point.
     */
    class RINGMESH_API Corner: public GeoModelMeshElement {
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
    public:
        /*! \brief Creates a Corner. A point is added to its Mesh.
         */
        Corner( const GeoModel& model, index_t id )
            : GeoModelMeshElement( model, CORNER, id )
        {
            mesh_.vertices.create_vertex() ;
        }

        ~Corner()
        {
        }

        virtual index_t nb_cells() const
        {
            // A Corner (dim 0) has points
            return mesh_.vertices.nb() ;
        }

        // Only one possible local index 
        virtual index_t gmme_vertex_index( index_t /*me*/, index_t /*lv*/ ) const 
        {
            return 0 ;
        }

    protected:
        virtual bool is_mesh_valid() const ;

    } ;

    /*!
     * @brief A GeoModelElement of type LINE
     *
     * @details One connected component of a 1-manifold.
     *         
     */
    class RINGMESH_API Line: public GeoModelMeshElement {
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
    public:
        Line( const GeoModel& model, index_t id ) ;

        ~Line()
        {
        }

        /*!
         * @brief Number of edges
         */
        virtual index_t nb_cells() const
        {
            return mesh_.edges.nb() ;
        }

        /*!
        * @brief Correspondance of vertex index \param lv in a mesh element \param me
        * and vertex index in the GMME
        */
        virtual index_t gmme_vertex_index( index_t me, index_t lv ) const
        {
            ringmesh_debug_assert( me < nb_cells() ) ;
            ringmesh_debug_assert( lv < 2 ) ;
            return mesh_.edges.vertex( me, lv ) ;
        }

        /*!
         * @brief A Line is closed if its two extremities are identitical 
         */
        bool is_closed() const
        {
            ringmesh_debug_assert( nb_boundaries() == 2 ) ;
            return ( boundaries_[0].is_defined() )
                && ( boundaries_[0] == boundaries_[1] ) ;
        }

    private:
        virtual bool is_mesh_valid() const ;

    } ;
    

    // Forward declaration
    class Surface ;

    /*
    * @todo Comment
    */
    class RINGMESH_API SurfaceTools {
    public:
        SurfaceTools( const Surface& surface ) ;
        ~SurfaceTools() ;

        const GEO::MeshFacetsAABB& aabb() const ;
        const ColocaterANN& ann() const ;

    private:
        const Surface& surface_ ;

        mutable GEO::MeshFacetsAABB* aabb_ ;
        mutable ColocaterANN* ann_ ;
    } ;


    /*!
     * @brief A GeoModelElement of type SURFACE
     *
     * @details One 2-manifold connected component .
     */
    class RINGMESH_API Surface: public GeoModelMeshElement {
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
    public:
        static const index_t NO_ADJACENT = index_t( -1 ) ;

        Surface( const GeoModel& model, index_t id )
            : GeoModelMeshElement( model, SURFACE, id ), tools( *this )
        {}

        ~Surface(){}

        /*!
         * @brief Number of facets
         */
        virtual index_t nb_cells() const
        {
            return mesh_.facets.nb() ;
        }

        /*!
        * @brief Correspondance of vertex index \param lv in a mesh element \param me
        * and vertex index in the GMME
        */
        virtual index_t gmme_vertex_index( index_t me, index_t lv ) const
        {
            ringmesh_debug_assert( me < nb_cells() ) ;
            ringmesh_debug_assert( lv < mesh_.facets.nb_vertices( me ) ) ;
            return mesh_.facets.vertex( me, lv ) ;
        }
         
        bool is_simplicial() const
        {
            return mesh_.facets.are_simplices() ;
        }
      
        /*!
         * \name Accessors to facet and vertices
         * @{
         */
        index_t facet_begin( index_t f ) const
        {
            return mesh_.facets.corners_begin( f ) ;
        }
        index_t facet_end( index_t f ) const
        {
            return mesh_.facets.corners_end( f ) ;
        }
        index_t nb_vertices_in_facet( index_t f ) const
        {
            return mesh_.facets.nb_vertices( f ) ;
        }

        bool is_triangle( index_t f ) const
        {
            return nb_vertices_in_facet( f ) == 3 ;
        }

        index_t next_in_facet( index_t f, index_t v ) const
        {
            ringmesh_debug_assert( v < nb_vertices_in_facet( f ) ) ;
            if( v != nb_vertices_in_facet( f ) - 1 ) {
                return v + 1 ;
            } else {
                return 0 ;
            }
        }

        index_t prev_in_facet( index_t f, index_t v ) const
        {
            ringmesh_debug_assert( v < nb_vertices_in_facet( f ) ) ;
            if( v > 0 ) {
                return v - 1 ;
            } else {
                return nb_vertices_in_facet( f ) - 1 ;
            }
        }

        index_t nb_facet_corners() const
        {
            return mesh_.facet_corners.nb() ;
        }

        index_t model_vertex_id_at_facet_corner( index_t corner ) const
        {
            return GeoModelMeshElement::model_vertex_id(
                mesh_.facet_corners.vertex( corner ) ) ;
        }

        /*!
         * @brief Returns the surface index of vertex \param v in facet \param f
         */
        index_t surf_vertex_id( index_t f, index_t v ) const
        {
            return gmme_vertex_index( f, v ) ;
        }
        
        /*!
         * @brief Returns a vertex surface index from its model index \param model_vertex_id
         * @details If there are two points, returns the first one.
         *          Returns NO_ID if no point is found
         */
        index_t surf_vertex_id( index_t model_vertex_id ) const
        {
            return gmme_vertex_index_from_model( model_vertex_id ) ;
        }

        index_t facet_vertex_id( index_t t, index_t surf_vertex_id ) const ;

        index_t facet_id_from_model( index_t f, index_t model_vertex_id ) const ;

        index_t facet_from_surface_vertex_ids( index_t in0, index_t in1 ) const ;

        index_t facet_from_model_vertex_ids( index_t i0, index_t i1 ) const ;

        void edge_from_model_vertex_ids(
            index_t i0,
            index_t i1,
            index_t& f,
            index_t& e ) const ;

        void oriented_edge_from_model_vertex_ids(
            index_t i0,
            index_t i1,
            index_t& facet,
            index_t& edge ) const ;

        index_t facets_around_vertex(
            index_t surf_vertex_id,
            std::vector< index_t >& result,
            bool border_only ) const ;

        index_t facets_around_vertex(
            index_t surf_vertex_id,
            std::vector< index_t >& result,
            bool border_only,
            index_t first_facet ) const ;

        /*! @}
         * \name Geometrical request on facets
         * @{
         */
        vec3 facet_normal( index_t facet_index ) const ;
        vec3 facet_barycenter( index_t facet_index ) const ;
        double facet_area( index_t facet_index ) const ;
        index_t closest_vertex_in_facet( index_t facet_index, const vec3& to_point ) const ;

        /*! @}
         * \name Adjacencies request
         * @{
         */
        /*! @brief Returns the index of the adjacent facet of \param f in this surface
         *  along the edge starting at \param v */
        index_t adjacent( index_t f, index_t v ) const
        {
            ringmesh_debug_assert( v < nb_vertices_in_facet( f ) ) ;
            return mesh_.facets.adjacent( f, v ) ;
        }

        /*! @brief Returns the index of the adjacent facet at the given corner
         */
        index_t adjacent( index_t c ) const
        {
            ringmesh_assert( c < mesh_.facet_corners.nb() ) ;
            return mesh_.facet_corners.adjacent_facet( c ) ;
        }

        bool is_on_border( index_t f, index_t v ) const
        {
            ringmesh_debug_assert( v < nb_vertices_in_facet( f ) ) ;
            return adjacent( f, v ) == GEO::NO_CELL ;
        }

        bool is_on_border( index_t f ) const
        {
            for( index_t adj = 0; adj < nb_vertices_in_facet( f ); adj++ ) {
                if( is_on_border( f, adj ) ) {
                    return true ;
                }
            }
            return false ;
        }

        void next_on_border(
            index_t f,
            index_t from,
            index_t v,
            index_t& next_f,
            index_t& v_in_next,
            index_t& to ) const ;

        void next_on_border(
            index_t f,
            index_t e,
            index_t& next_f,
            index_t& next_e ) const ;

    public:
        SurfaceTools tools ;

    private:
        virtual bool is_mesh_valid() const ;
        	    
    } ;



    // Defined just after this class
    class Region ;
    /*
     * @todo Comment 
     */
    class RINGMESH_API RegionTools {
    public:
        RegionTools( const Region& region ) ;
        ~RegionTools() ;

        const GEO::MeshCellsAABB& aabb() const ;
        const ColocaterANN& ann() const ;

    private:
        const Region& region_ ;

        mutable GEO::MeshCellsAABB* aabb_ ;
        mutable ColocaterANN* ann_ ;
    } ;

    /*!
     * @brief A GeoModelElement of type REGION
     *
     * @details The Region can be defined only defined by its boundary
     * Surfaces. Its volumetric mesh is optional.
     */
    class RINGMESH_API Region: public GeoModelMeshElement {
        friend class GeoModelEditor ;
        friend class GeoModelBuilder ;
    public:
        Region( const GeoModel& model, index_t id )
            : GeoModelMeshElement( model, REGION, id ), tools( *this )
        {
        }

        Region( 
            const GeoModel& model, 
            index_t id,
            const std::string& name,
            GEOL_FEATURE geological_feature
            )
            :
            GeoModelMeshElement( model, REGION, id, name, geological_feature ),
            tools( *this )
        {
        }

        ~Region()
        {
        }

        /*!
         * Get the number of cells
         */
        virtual index_t nb_cells() const
        {
            return mesh_.cells.nb() ;
        }

        index_t nb_vertices_in_cell( index_t c ) const
        {
            return mesh_.cells.nb_vertices( c ) ;
        }
        virtual index_t gmme_vertex_index( index_t me, index_t lv ) const
        {
            ringmesh_debug_assert( me < nb_cells() ) ;
            ringmesh_debug_assert( lv < mesh_.cells.nb() ) ;
            return mesh_.cells.vertex( me, lv ) ;
        }

        bool is_on_border( index_t cell, index_t facet ) const
        {
            return adjacent_cell(cell, facet) == GEO::NO_CELL;
        }

        index_t adjacent_cell( index_t cell, index_t facet ) const
        {
            return mesh_.cells.adjacent( cell, facet );
        }

        bool is_meshed() const
        {
            return mesh().cells.nb() > 0 ;
        }

        bool side( index_t i ) const
        {
            return sides_[i] ;
        }

        bool is_simplicial() const
        {
            return mesh().cells.are_simplices() ;
        }

        /*!
         * \todo Is connectivity valid should be virtual and we should 
         * reimplement here to check 
         *     if( nb_boundaries() != sides_.size() ) {
         *       GEO::Logger::err( "GeoModelElement" )
         *       << gme_id() << " boundary sides are invalid "
         *       << std::endl ;
         *       return false ;
         *     }

         */

    private:
        virtual bool is_mesh_valid() const ;

    private:
        /*! Additional information to store oriented boundary Surfaces
         * Side: + (true) or - (false)
         * The size of this vector must be the same than boundary_
         */
        std::vector< bool > sides_ ;

    public:
        RegionTools tools ;
    } ;



} // namespace

#endif
