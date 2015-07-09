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

/*! \author Jeanne Pellerin and Arnaud Botella */

#ifndef __RINGMESH_BOUNDARY_MODEL__
#define __RINGMESH_BOUNDARY_MODEL__

#include <ringmesh/common.h>
#include <ringmesh/boundary_model_element.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/file_system.h>
#include <geogram/points/kd_tree.h>

#include <vector>
#include <string>
#include <iterator>
#include <algorithm>

namespace RINGMesh {

    /*!
     * @brief Unique storage of the vertices of a BoundaryModel
     * @details Each instance is unique, unlike vertices in 
     *          the model Corner, Line, and Surface meshes.
     *          Attributes may be defined on the vertices.
     */
    class RINGMESH_API BoundaryModelVertices {
    ringmesh_disable_copy( BoundaryModelVertices ) ;
    public:

        /*!
         * @brief Identification of a vertex in a BoundaryModelElement
         */
        struct VertexInBME {
            VertexInBME( BME::bme_t t, index_t vertex_id_in )
                : bme_id( t ), v_id( vertex_id_in )
            {
            }
            VertexInBME()
                : bme_id(), v_id( NO_ID )
            {
            }
            bool operator<( const VertexInBME& rhs ) const
            {
                if( bme_id != rhs.bme_id ) {
                    return bme_id < rhs.bme_id ;
                } else {
                    return v_id < rhs.v_id ;
                }
            }
            bool operator==( const VertexInBME& rhs ) const
            {
                return bme_id == rhs.bme_id && v_id == rhs.v_id ;
            }
            bool is_defined() const
            {
                return bme_id.is_defined() && v_id != NO_ID ;
            }
            /// Unique identifier of the associated BoundaryModelElement
            BME::bme_t bme_id ;
            /// Index of the vertex in the BME
            index_t v_id ;
        } ;

        /*!
         * @brief Constructor from an existing BoundaryModel
         */
        BoundaryModelVertices( const BoundaryModel& bm )
            : bm_( bm ), kdtree_( nil ), lock_( 0 ), kdtree_to_update_( true )
        {
        }

        ~BoundaryModelVertices()
        {
        }

        /// Maybe this is not the best test, maybe we should store a bool
        /// to store this
        bool is_initialized() const
        {
            return mesh_.vertices.nb() > 0 ;
        }

        /*!
         * @brief Number of vertices stored. 
         * @details Calls initialize() if vertices are not filled yet
         * @todo remove
         */
        index_t nb_unique_vertices() const ;

        /*!
         * @brief Number of vertices stored.
         * @details Calls initialize() if vertices are not filled yet
         */
        index_t nb() const ;

        /*!
         * @brief Coordinates of a vertex of the BoundaryModel
         * @pre v < nb()
         */
        const vec3& unique_vertex( index_t v ) const ;

        /*!
         * @brief Returns the index of the given vertex in the model
         * @param[in] p input point coordinates
         * @return index of the vertex in the model if found 
         * (distance < epsilon), otherwise NO_ID
         */
        index_t vertex_index( const vec3& p ) const ;

        /*!
         * @brief Get the vertices in BME corresponding to the given unique vertex
         */
        const std::vector< VertexInBME >& bme_vertices( index_t v ) const ;

        /*!
         * @brief To use when building the model by first adding its vertices
         * @warning The client is responsible for setting the mapping between the points
         * of the BME and the unique vertex 
         */
        index_t add_unique_vertex( const vec3& point ) ;

        /*!
         * @brief Add a vertex in a BoundaryModelElement 
         *        corresponding to an existing vertex of the model
         */
        void add_unique_to_bme( index_t v, const VertexInBME& v_bme ) ;

        /*! 
         * @brief Change one of the BME vertex associated to a vertex
         * @param v Index of the vertex
         * @param i Index of the BME vertex
         * @param v_bme Id of BME and of the vertex in that BME
         */
        void set_bme( index_t v, index_t i, const VertexInBME& v_bme ) ;

        /*!
         * @brief Set the point coordinates of all the vertices that are 
         *        share this unique vertex, including the unique vertex itself.
         * @param[in] v Index of the vertex
         * @param[in] point New coordinates
         */
        void update_point( index_t v, const vec3& point ) ;

        /*!
         * @brief Clear the vertices - unbind bme_vertices_ - 
         *        set attribute to NO_ID in BME
         * @warning Not stable - crashes because of issues in 
         * Mesh attributes clearing
         * @todo Unbind all attributes !!!! 
         */
        void clear() ;

        /*!
         * @brief Returns the Geogram attribute manager on these vertices
         */
        GEO::AttributesManager& attribute_manager() const
        {
            return mesh_.vertices.attributes() ;
        }

        /*!
         * @brief Remove colocated vertices
         */
        void remove_colocated() ;

        /*! 
         * @brief Remove all invalid VertexInBME and delete the vertices 
         * that are not anymore in any BoundaryModelElement
         */
        void erase_invalid_vertices() ;

        /*! 
         * @brief Delete vertices for which to_delete[i] != i ;
         * @param[in,out] to_delete can be NO_ID or give the index of a 
         *  kept vertex with wich information should be merged.
         *  It is recyled to give the mapping between old and new vertex indices        
         * @pre to_delete[ v ] is either NO_ID, or is equal or inferior to v
         */
        void erase_vertices( std::vector< index_t >& to_delete ) ;

    private:
        /*!
         * @brief Initialize the vertices from the vertices 
         *        of the BoundaryModel Corner s, Line s, and Surface s
         *        Duplicates are not removed yet.
         * @details Fills unique_vertices_, set model vertex id in BME and
         *          initialize the reverse information.
         */
        void initialize() ;

        /*!
         * @brief Set a vertex as invalid
         * @details Put all the ids of its VertexInBME to NO_ID
         */
        void set_invalid_vertex( index_t v ) ;

        /*!
         * @brief True if the vertex is not valid
         */
        bool is_invalid_vertex( index_t v ) const ;

        /*!
         * @brief Delete the KdTree and set the pointer to nil.         
         */
        void set_to_update() ;

        /*!
         * @brief Build the KdTree of the vertices 
         * @pre In debug mode, assert that ann_ pointer is nil,
         *      and that there is no colocated vertices.
         *
         * @note The function is const to be called when accessing a point index
         *  from coordinated without an ugly const-cast.
         */
        void initialize_kdtree() const ;

        /*!
         * @brief Update the model_vertex_id in the corners, lines, surfaces
         *        from new to old index
         */
        void update_bme_model_ids( const GEO::vector< index_t >& old2new ) const ;

    private:
        /// Attached BoundaryModel that owns the vertices
        const BoundaryModel& bm_ ;

        /*! 
         * @brief Mesh storing the coordinates of the vertices that are not colocated
         * @details Each vertex is unique. 
         * With a GEO::Mesh we have attributes on the points without any effort
         */
        GEO::Mesh mesh_ ;

        /*! 
         * Vertices in BoundaryModelElements corresponding to each vertex
         * @todo Change this extremely expensive storage !!!
         */
        std::vector< std::vector< VertexInBME > > bme_vertices_ ;

        /// Kd-tree of the model vertices
        mutable GEO::NearestNeighborSearch_var kdtree_ ;
        mutable bool kdtree_to_update_ ;

        /// Lock to protect from multi-threading during clear()
        GEO::Process::spinlock lock_ ;
    } ;

    /*!
     * @brief The class to describe a volumetric model represented by its boundary surfaces
     */
    class RINGMESH_API BoundaryModel {
    ringmesh_disable_copy( BoundaryModel ) ;
        friend class BoundaryModelBuilder ;

    public:
        const static index_t NO_ID = index_t( -1 ) ;

        /*!
         * @brief Constructs an empty BoundaryModel
         */
        BoundaryModel()
            :
                vertices( *this ),
                debug_directory_( GEO::FileSystem::get_current_working_directory() )
        {
        }

        /*!
         * @brief Delete all BoundaryModelElements stored and owned by the BoundaryModel
         */
        virtual ~BoundaryModel() ;

        void copy( const BoundaryModel& from ) ;

        /*!
         * @brief Name of the model
         */
        const std::string& name() const
        {
            return name_ ;
        }

        /*!
         * @brief Get the directory for debug information
         */
        const std::string& debug_directory() const
        {
            return debug_directory_ ;
        }

        /*!
         * @brief Set the directory where debugging information shall be stored
         * @details Test that this directory exists, if not
         *          keep the previous value.
         *          Default directory is executable directory .
         */
        void set_debug_directory( const std::string& directory ) ;

        /*!
         * @brief Number of unique vertices
         */
        index_t nb_vertices() const
        {
            return vertices.nb() ;
        }

        index_t nb_facets() const ;

        /*!
         * \name Generic BoundaryModelElement accessors
         * @{
         */

        /*!
         * @brief Returns the number of elements of the given type
         * @details By default returns 0.
         */
        inline index_t nb_elements( BME::TYPE type ) const
        {
            if( type < BME::NO_TYPE ) {
                return elements( type ).size() ;
            } else if( type == BME::ALL_TYPES ) {
                ringmesh_assert( !nb_elements_per_type_.empty() ) ;
                return nb_elements_per_type_.back() ;
            } else {
                return 0 ;
            }
        }

        /*!
         * @brief Returns a const reference the identified BoundaryModelElement
         *
         * @param[in] id Id of the element
         *
         */
        inline const BoundaryModelElement& element( BME::bme_t id ) const
        {
            ringmesh_assert( id.index < nb_elements( id.type ) ) ;
            if( id.type < BME::NO_TYPE ) {
                return *elements( id.type )[id.index] ;
            } else if( id.type == BME::ALL_TYPES ) {
                return element( global_to_typed_id( id ) ) ;
            } else {
                return universe_ ;
            }
        }

        /*!
         * @brief Returns a const reference the identified BoundaryModelElement
         *
         * @param[in] type Type of the element
         * @param[in] index Index of the element
         *
         */
        inline const BoundaryModelElement& element(
            BME::TYPE bme_type,
            index_t index ) const
        {
            ringmesh_debug_assert( bme_type < BME::NO_TYPE ) ;
            ringmesh_debug_assert( index < elements( bme_type ).size() ) ;
            return *( elements( bme_type )[index] ) ;
        }

        /*! @}
         * \name Specicalized accessors.
         * @{
         */
        index_t nb_corners() const
        {
            return nb_elements( BME::CORNER ) ;
        }
        index_t nb_lines() const
        {
            return nb_elements( BME::LINE ) ;
        }
        index_t nb_surfaces() const
        {
            return nb_elements( BME::SURFACE ) ;
        }
        index_t nb_regions() const
        {
            return nb_elements( BME::REGION ) ;
        }
        index_t nb_contacts() const
        {
            return nb_elements( BME::CONTACT ) ;
        }
        index_t nb_interfaces() const
        {
            return nb_elements( BME::INTERFACE ) ;
        }
        index_t nb_layers() const
        {
            return nb_elements( BME::LAYER ) ;
        }

        const Corner& corner( index_t index ) const
        {
            return dynamic_cast< const Corner& >( *corners_.at( index ) ) ;
        }
        const Line& line( index_t index ) const
        {
            return dynamic_cast< const Line& >( *lines_.at( index ) ) ;
        }
        const Surface& surface( index_t index ) const
        {
            return dynamic_cast< const Surface& >( *surfaces_.at( index ) ) ;
        }

        const BoundaryModelElement& region( index_t index ) const
        {
            return element( BME::bme_t( BME::REGION, index ) ) ;
        }

        const BoundaryModelElement& contact( index_t index ) const
        {
            return element( BME::bme_t( BME::CONTACT, index ) ) ;
        }

        const BoundaryModelElement& one_interface( index_t index ) const
        {
            return element( BME::bme_t( BME::INTERFACE, index ) ) ;
        }

        const BoundaryModelElement& layer( index_t index ) const
        {
            return element( BME::bme_t( BME::LAYER, index ) ) ;
        }

        const BoundaryModelElement& universe() const
        {
            return universe_ ;
        }

        GEO::AttributesManager& vertex_attribute_manager()
        {
            return vertices.attribute_manager() ;
        }

        index_t find_region( index_t surf_id, bool side ) const ;

        /*! @}
         * \name To save the BoundaryModel.
         * @{
         */
        bool save_gocad_model3d( std::ostream& out ) ;
        void save_as_eobj_file( const std::string& file ) const ;
        void save_surface_as_obj_file( index_t s, const std::string& file ) const ;
        void save_bm_file( const std::string& file_name ) const ;
        void save_smesh_file( const std::string& file_name ) const ;

        /*!
         * @}
         */

        index_t find_element( BME::TYPE type, const std::string& name ) const ;

        void translate( const vec3& ) ;
        void rotate(
            const vec3& origin,
            const vec3& axis,
            float64 angle,
            bool degrees = false ) ;

    private:
        bool check_model_validity() const ;
        bool check_elements_validity() const ;
        bool check_geology_validity() const ;
        bool check_gocad_validity() const ;

        void copy_macro_topology( const BoundaryModel& from ) ;
        void copy_meshes( const BoundaryModel& from ) ;

        /*! 
         * @brief Convert a global BME index into a typed index
         * @details Relies on the nb_elements_per_type_ vector that 
         *          must be updodate 
         *          See the BoundaryModelBuilder::end_model() function
         * @param[in] global A BME id of TYPE - ALL_TYPES
         * @return A BME id of an element of the model, or a invalid one if nothing found
         */
        inline BME::bme_t global_to_typed_id( const BME::bme_t& global ) const
        {
            index_t t = NO_ID ;
            for( index_t i = 1; i < nb_elements_per_type_.size(); i++ ) {
                if( global.index >= nb_elements_per_type_[i - 1]
                    && global.index < nb_elements_per_type_[i] ) {
                    t = i - 1 ;
                    break ;
                }
            }
            if( static_cast< BME::TYPE >( t ) < BME::NO_TYPE ) {
                BME::TYPE T = static_cast< BME::TYPE >( t ) ;
                index_t i = global.index - nb_elements_per_type_[t] ;
                return BME::bme_t( T, i ) ;
            } else {
                return BME::bme_t() ;
            }
        }

        /*!
         * @brief Generic accessor to the storage of elements of the given type
         * @pre The type must be valid NO_TYPE or ALL_TYPES will throw an assertion
         */
        std::vector< BME* >& modifiable_elements( BME::TYPE type )
        {
            return const_cast< std::vector< BME* >& >( elements( type ) ) ;
        }

        const std::vector< BME* >& elements( BME::TYPE type ) const
        {
            switch( type ) {
                case BME::CORNER:
                    return corners_ ;
                case BME::LINE:
                    return lines_ ;
                case BME::SURFACE:
                    return surfaces_ ;
                case BME::REGION:
                    return regions_ ;
                case BME::CONTACT:
                    return contacts_ ;
                case BME::INTERFACE:
                    return interfaces_ ;
                case BME::LAYER:
                    return layers_ ;
                default:
                    ringmesh_assert_not_reached;
                return surfaces_ ; ;
            }
        }

    public:
        BoundaryModelVertices vertices ;

    private:
        // Name of the model
        std::string name_ ;

        // Base manifold elements of a model
        std::vector< BoundaryModelElement* > corners_ ;
        std::vector< BoundaryModelElement* > lines_ ;
        std::vector< BoundaryModelElement* > surfaces_ ;
        std::vector< BoundaryModelElement* > regions_ ;

        /// The region including all the other regions
        BoundaryModelElement universe_ ;

        /*!
         * @brief Contacts between Intefaces
         * Parent of a set of Line
         */
        std::vector< BoundaryModelElement* > contacts_ ;
        /*!
         * @brief Interfaces between layers
         * Parent of a set of Surface
         */
        std::vector< BoundaryModelElement* > interfaces_ ;

        /*!
         * @brief Rock layers
         * Parent of a set of Region
         */
        std::vector< BoundaryModelElement* > layers_ ;

        /// Allow global access to BME. It MUST be updated if one element is added.
        std::vector< index_t > nb_elements_per_type_ ;

        /// Name of the debug directory in which to save stuff 
        /// @note Maybe move this in another class
        std::string debug_directory_ ;
    } ;

}

#endif
