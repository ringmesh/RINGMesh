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

#include <vector>
#include <string>
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
            VertexInBME(
                BME::bme_t t,
                index_t vertex_id_in )
                : bme_id( t ), v_id( vertex_id_in )
            {
            }
            bool operator<( const VertexInBME& rhs ) const
            {
                if( bme_id != rhs.bme_id ) {
                    return bme_id < rhs.bme_id ;
                }
                else {
                    return v_id < rhs.v_id ;
                }
            }
            /// Type of the BME and index
            BME::bme_t bme_id ;
            /// Index of the vertex in the BME
            index_t v_id ;
        } ;

        /*!
         * @brief Vertices are defined for a BoundaryModel
         */
        BoundaryModelVertices( const BoundaryModel& bm )
            : bm_( bm ), ann_( nil ), lock_( 0 )
        {
        }

        ~BoundaryModelVertices() ; 

        /*!
         * @brief Number of vertices stored. 
         * @details Calls initialize_unique_vertices() if they are not filled yet
         */
        index_t nb_unique_vertices() const ;

        /*!
         * @brief Get the index of the BM vertex corresponding to vertex v
         *        in the BME of type T and index id.
         * @details Calls initialize_unique_vertices(), if no vertices yet
         *          The unique_id is stored as an attribute on the vertices of the BME
         * 
         * @todo Remove ? This function is not related to this class
         * This is managed by the BME itself.
         */
        index_t unique_vertex_id(
            BoundaryModelElement::bme_t T, index_t v ) const ;

        /*!
         * @brief Coordinates of a vertex of the BoundaryModel
         * @pre unique_id < nb_unique_vertices()
         */
        const vec3& unique_vertex( index_t unique_id ) const ;        
        
        /*!
         * @brief Returns the index of the given vertex in the model
         * @param[in] p input point coordinates
         * @return index of the vertex in the model if found, otherwise NO_ID
         */
        index_t vertex_index( const vec3& p ) const ;

        /*!
         * @brief Get the vertices in BME that correspond to the given unique vertex
         */
        const std::vector< VertexInBME >& bme_vertices( index_t unique_id ) const ;
        
        /*!
         * @brief To use when building the model by first adding its vertices
         * @warning The client is responsible for setting the mapping between the points
         * of the BME and the unique vertex 
         */
        index_t add_unique_vertex( const vec3& point ) ;

        /*!
         * @brief Add a vertex in a BoundaryModelElement 
         *        corresponding to an existing unique_vertex
         */
        void add_unique_to_bme( 
            index_t unique_id, 
            BME::bme_t bme_id,
            index_t v_id ) ;        

        /*!
         * @brief Set the point coordinates of all the vertices that are 
         *        share this unique vertex, including the unique vertex itself.
         * @param[in] unique_id Index of the unique vertex in the BoundaryModel
         * @param[in] point New coordinates of the vertex 
         */
        void update_point( index_t unique_id, const vec3& point ) ;

        /*!
         * @brief Clear the vertices - unbind unique2bme_ - 
         *        set attribute to NO_ID in BME
         * @warning Not stable - crashes because of issues in 
         * Mesh attributes clearing
         */  
        void clear() ;

        /*!
         * @brief Returns the Geogram attribute manager on these vertices
         */
        GEO::AttributesManager& attribute_manager() const {
            return const_cast<GEO::AttributesManager&> 
                ( unique_vertices_.vertices.attributes() );
        }
        
    private:
        /*!
         * @brief Determine the unique vertices from the vertices 
         *        of the BoundaryModel Corner s, Line s, and Surface s
         * @details Fills unique_vertices_ and set the attributes the global index on
         *          the BoundaryModel Corner, Line and Surface. 
         */
        void initialize_unique_vertices() ;

        /*!
         * @brief Fills the unique2bme vector
         * @details Call initialize_unique_vertices() if unique_vertices_ is empty
         */
        void initialize_reverse() ;

        /*!
        * @brief Get the index of the BM vertex corresponding to 
        *        a given vertex in a BME.
        */
        index_t unique_vertex_id( const VertexInBME& v ) const ;       

        /*!
         * @brief Delete the KdTree and set the pointer to nil.         
         */
        void set_ann_to_update() ;
        
        /*!
         * @brief Build the KdTree. 
         * @pre In debug mode, assert that ann_ pointer is nil.
         * @note Function is const to be called in accessors to point index.
         *  without ugly const-cast.
         */
        void initialize_ann() const ;
       
    private:
        /// Attached BoundaryModel to which belong the vertices
        const BoundaryModel& bm_ ;
        
        /*! 
         * @brief Mesh storing the coordinates of the vertices that are not colocated
         * @details Each point instance is unique. 
         * With a GEO::Mesh we have attributes on the points without any effort
         */
        GEO::Mesh unique_vertices_ ;
               
        /*! 
         * Mapping of a unique vertex to the vertices in the 
         * BoundaryModelElements that have the same coordinates
         */
        GEO::Attribute< std::vector< VertexInBME > > unique2bme_ ;

        /// Kd-tree of the model vertices
        mutable ColocaterANN* ann_ ;

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
        const static index_t NO_ID = index_t( - 1 ) ;

        /*!
         * @brief Constructs an empty BoundaryModel
         */
        BoundaryModel() :
            vertices( *this ),
            debug_directory_( 
                GEO::FileSystem::get_current_working_directory() )
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
        const std::string& name() const { return name_ ; }

        /*!
         * @brief Get the directory for debug information
         */
        const std::string& debug_directory() const { return debug_directory_ ; }

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
        index_t nb_vertices() const { return vertices.nb_unique_vertices() ; }

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
            switch( type ) {
                 case BoundaryModelElement::CORNER    :   return corners_.size() ;
                 case BoundaryModelElement::LINE      :   return lines_.size() ;
                 case BoundaryModelElement::SURFACE   :   return surfaces_.size() ;
                 case BoundaryModelElement::REGION    :   return regions_.size() ;
                 case BoundaryModelElement::CONTACT   :   return contacts_.size() ;
                 case BoundaryModelElement::INTERFACE :   return interfaces_.size() ;
                 case BoundaryModelElement::LAYER     :   return layers_.size() ;
                 case BoundaryModelElement::ALL_TYPES :
                     {
                    ringmesh_assert( nb_elements_per_type_.size() > 0 ) ;
                    ringmesh_debug_assert(
                        nb_elements_per_type_.back()
                            == corners_.size() + lines_.size() + surfaces_.size()
                                + regions_.size() + contacts_.size()
                                + interfaces_.size() + layers_.size() ) ;
                    return nb_elements_per_type_.back() ;
                    }
                 default :  
                     return 0 ;
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
            BME::bme_t id ) const
        {
            ringmesh_assert( id.index < nb_elements( id.type ) ) ;
            switch( id.type ) {
                 case BME::CORNER    :  return *corners_[ id.index ] ;
                 case BME::LINE      :  return *lines_[ id.index ] ;
                 case BME::SURFACE   :  return *surfaces_[ id.index ] ;
                 case BME::REGION    :  return *regions_[ id.index ] ;
                 case BME::CONTACT   :  return *contacts_[ id.index ] ;
                 case BME::INTERFACE :  return *interfaces_[ id.index ] ;
                 case BME::LAYER     :  return *layers_[ id.index ] ;
                 case BME::ALL_TYPES :  return element( global_to_typed_id( id ) ) ;                     
                 default :
                     ringmesh_assert_not_reached ;
                     // By default, return the universe
                     return universe_ ;
            }
        }

        /*! @}
         * \name Specicalized accessors.
         * @{
         */
        index_t nb_corners()    const { return nb_elements( BME::CORNER )    ; }
        index_t nb_lines()      const { return nb_elements( BME::LINE )      ; }
        index_t nb_surfaces()   const { return nb_elements( BME::SURFACE )   ; }
        index_t nb_regions()    const { return nb_elements( BME::REGION )    ; }
        index_t nb_contacts()   const { return nb_elements( BME::CONTACT )   ; }
        index_t nb_interfaces() const { return nb_elements( BME::INTERFACE ) ; }
        index_t nb_layers()     const { return nb_elements( BME::LAYER )     ; }

        const Corner& corner( index_t index ) const { return *corners_.at( index ) ; }
        const Line& line( index_t index ) const { return *lines_.at( index ) ; }
        const Surface& surface( index_t index ) const { return *surfaces_.at( index ) ;}

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

        const BoundaryModelElement& universe() const { return universe_ ; }

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

        signed_index_t find_interface( const std::string& name) const ;
        signed_index_t find_region( const std::string& name) const ;
        

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
        inline BME::bme_t global_to_typed_id(
            const BME::bme_t& global ) const
        {     
            index_t t = NO_ID ;
            for( index_t i = 1; i < nb_elements_per_type_.size(); i++ ) {
                if( global.index >= nb_elements_per_type_[ i - 1 ]
                    && global.index < nb_elements_per_type_[ i ] 
                  ) {
                    t = i - 1 ;
                    break ;
                 }                
            }
            if( ( BME::TYPE ) t < BME::NO_TYPE ) {
                BME::TYPE T = ( BME::TYPE ) ( t ) ;
                index_t i = global.index - nb_elements_per_type_[ t ] ;
                return BME::bme_t( T, i ) ;
            }
            else {
                return BME::bme_t() ;                
            }            
        }

    public:
        BoundaryModelVertices vertices ;

    private:
        // Name of the model
        std::string name_ ;

        // Base manifold elements of a model
        std::vector< Corner* > corners_ ;
        std::vector< Line* > lines_ ;
        std::vector< Surface* > surfaces_ ;
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
