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
     * \todo Implement a BoundaryModelMutator
     */
    class GRGMESH_API BoundaryModel {       
        friend class BoundaryModelBuilder ;

    public:           
        enum AttributeLocation {
            VERTEX
        } ;       
        typedef AttributeManagerImpl< VERTEX > VertexAttributeManager ;
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

        index_t nb_vertices() const { return vertices_.size() ; }        
        index_t vertex_index( const vec3& p ) const ;
        const vec3& vertex( index_t p ) const { return vertices_.at(p) ; }
        void set_vertex_coordinates( index_t id, const vec3& p ) {
            grgmesh_assert( id < nb_vertices() ) ;
            vertices_[id] = p ;
        }
        index_t nb_facets() const ;             

        /** 
        * \name Generic BoundaryModelElement accessors
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

        /** @}
        * \name Specicalized accessors.
        * @{
        */
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

        VertexAttributeManager* vertex_attribute_manager() const
        {
            return const_cast< VertexAttributeManager* >( &vertex_attribute_manager_ ) ;
        }

        index_t find_region( index_t surface_part_id, bool side ) const ;

        /** @}
        * \name To save the BoundaryModel.
        * @{
        */
        bool save_gocad_model3d( std::ostream& out ) ;
        void save_as_eobj_file( const std::string& file_name ) ;
        void save_surface_as_obj_file( index_t s, const std::string& file_name ) const ;
        void save_bm_file( const std::string& file_name ) ;

        /**
        * @}
        */
    private:        
        bool check_model3d_compatibility() ;
        
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

        /// Allow global access to BME. It MUST be updated if one element is added.
        std::vector< index_t > nb_elements_per_type_ ;


        // Attribute manager 
        VertexAttributeManager vertex_attribute_manager_ ;     
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
       
}

#endif
