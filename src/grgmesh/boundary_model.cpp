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


#include <grgmesh/boundary_model.h>
#include <grgmesh/utils.h>

#include <geogram/basic/logger.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <set>


namespace GRGMesh {

    /*!
     * @brief Clear memory (useful to reduce memory impact)
     */
    void BoundaryModel::clear()
    {
        vertices_.clear() ;

        corners_.clear() ;
        lines_.clear() ;
        surfaces_.clear() ;
        regions_.clear() ;

        nb_facets_in_surfaces_.clear() ;
        contacts_.clear() ;
        interfaces_.clear() ;
        layers_.clear() ;
    }

    /*!
     * @brief Write in the stream the geological feature name
     * @param[in,out] out stream in which the name is written
     * @param[in] t geological feature
     *
     *  \todo Implement a real map between these types and the enum - Add more geological types
     */
    void BoundaryModel::save_type(
        std::ostream& out,
        BoundaryModelElement::GEOL_FEATURE t )
    {
        switch( t ) {
            case BoundaryModelElement::STRATI:
                out << "top" ;
                break ;
            case BoundaryModelElement::FAULT:
                out << "fault" ;
                break ;
            case BoundaryModelElement::VOI:
                out << "boundary" ;
                break ;
            case BoundaryModelElement::ALL:
                out << "none" ;
                break ;
            default:
                out << "none" ;
                break ;
        }
    }

    /*!
     * \brief Returns the index of the given vertex in the model
     * \todo Implement the function - Add a KdTree for geometrical request on model vertices
     *
     * @param[in] p input point coordinates
     * @return NO_ID
     */
    index_t BoundaryModel::vertex_index( const vec3& p ) const {
       
        grgmesh_assert_not_reached ;
        return NO_ID ;
    }

     /*!
     * \brief Returns a const reference the identified BoundaryModelElement   
     *
     * @param[in] t Type of the element 
     * @param[in] index Index of the element
     * @return A const reference to the BoundaryModelElement.      
     *
     */
    const BoundaryModelElement& BoundaryModel::element(
        BoundaryModelElement::BM_TYPE t,
        index_t index ) const
    {
        grgmesh_assert( index < nb_elements( t ) ) ;
        switch( t ){
            case BoundaryModelElement::BM_CORNER    : return corners_   [ index ] ;
            case BoundaryModelElement::BM_LINE      : return lines_     [ index ] ;
            case BoundaryModelElement::BM_SURFACE   : return surfaces_  [ index ] ;
            case BoundaryModelElement::BM_REGION    : return regions_   [ index ] ;
            case BoundaryModelElement::BM_CONTACT   : return contacts_  [ index ] ;
            case BoundaryModelElement::BM_INTERFACE : return interfaces_[ index ] ;
            case BoundaryModelElement::BM_LAYER     : return layers_    [ index ] ;
            case BoundaryModelElement::BM_ALL_TYPES :
                {
                    // This must synchro with what is done in the builder
                    index_t t = NO_ID ;
                    for( index_t i = 1; i < nb_elements_per_type_.size(); i++ ) {
                        if( index >= nb_elements_per_type_[i-1] && index < nb_elements_per_type_[i] ) {
                            t = i-1 ; break ;
                        }
                    }
                    grgmesh_assert( t < BME::BM_NO_TYPE ) ;
                    return element( (BME::BM_TYPE) t, index - nb_elements_per_type_[t] ) ;
                }
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
        }
    }

     /*!
     * \brief Private generic accessor to a BoundaryModelElement.
     *        For internal use only.
     *
     * @param[in] element_type Type of the element 
     * @param[in] index Index of the element
     * @return A reference to the BoundaryModelElement .
     *
     *
    BoundaryModelElement& BoundaryModel::element( BoundaryModelElement::BM_TYPE element_type, index_t index ) 
    { 
        grgmesh_assert( index < nb_elements( element_type ) ) ;
        switch( element_type ){
            case BoundaryModelElement::BM_CORNER    : return corners_   [ index ] ;
            case BoundaryModelElement::BM_LINE      : return lines_     [ index ] ;
            case BoundaryModelElement::BM_SURFACE   : return surfaces_  [ index ] ;
            case BoundaryModelElement::BM_REGION    : return regions_   [ index ] ;
            case BoundaryModelElement::BM_CONTACT   : return contacts_  [ index ] ;
            case BoundaryModelElement::BM_INTERFACE : return interfaces_[ index ] ;
            case BoundaryModelElement::BM_LAYER     : return layers_    [ index ] ;
            default:
                grgmesh_assert_not_reached ;
                return corners_.at( 0 ) ;
        }
    }*/

    /*!
     * @brief Returns the number of elements of the given type
     * By default returns 0.
     */ 
    index_t BoundaryModel::nb_elements( BoundaryModelElement::BM_TYPE type ) const 
    {
         switch( type ){
            case BoundaryModelElement::BM_CORNER    : return corners_.size() ;
            case BoundaryModelElement::BM_LINE      : return lines_.size() ;
            case BoundaryModelElement::BM_SURFACE   : return surfaces_.size() ;
            case BoundaryModelElement::BM_REGION    : return regions_.size() ;
            case BoundaryModelElement::BM_CONTACT   : return contacts_.size() ;
            case BoundaryModelElement::BM_INTERFACE : return interfaces_.size() ;
            case BoundaryModelElement::BM_LAYER     : return layers_.size() ;            
            case BoundaryModelElement::BM_ALL_TYPES : 
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
     * \brief Returns the index of the region neighboring the surface.
     * @param[in] surface_part_id Index of the Surface
     * @param[in] side Side of the Surface
     * @return The region index or NO_ID if none found.
     */
    index_t BoundaryModel::find_region( index_t surface_part_id, bool side ) const
    {
        grgmesh_debug_assert( surface_part_id < nb_surfaces() ) ;
        for( index_t r = 0; r < nb_regions(); r++ ) {
            const BoundaryModelElement& cur_region = region( r ) ;
            for( index_t s = 0; s < cur_region.nb_boundaries(); s++ ) {
                if( cur_region.side( s ) == side
                    && cur_region.boundary_id( s ) == surface_part_id ) {
                    return r ;
                }
            }
        }
        return BoundaryModelElement::NO_ID ;
    }

    /*!
     * @brief Convert a global facet index to an index in a Surface
     *
     * @param[in] model_facet_id Facet index in the BoundaryModel
     * @param[out] surface_id Surface containing the facet
     * @param[out] surf_facet_id Index of the facet in this Surface
     */
    void BoundaryModel::surface_facet(
        index_t model_facet_id, index_t& surface_id, index_t& surf_facet_id 
    ) const {
        index_t s = Surface::NO_ID ;
        for( index_t i = 1; i < nb_facets_in_surfaces_.size(); i++ ) {
            if( model_facet_id >= nb_facets_in_surfaces_[i-1] && 
                model_facet_id < nb_facets_in_surfaces_[i] )
            {
                s = i-1 ;
                break ;
            }
        }       
        grgmesh_debug_assert( s != Surface::NO_ID ) ;
            
        surface_id = s ;
        surf_facet_id = model_facet_id - nb_facets_in_surfaces_[s] ;
    }

    /*!
     * @brief Convert a facet index in a surface to its index in the BoundaryModel
     * @param[in] surf_id Surface index
     * @param[in] facet_id Facet index in surf_id
     * @return The global facet index
     */
    index_t BoundaryModel::model_facet( index_t surf_id, index_t facet_id ) const 
    {
        grgmesh_debug_assert( surf_id < nb_surfaces() && 
            facet_id < surface(surf_id).nb_cells() )  ;

        return nb_facets_in_surfaces_[surf_id] + facet_id ;
    }

    /*! 
     *  No support for properties right now
     *  no advance check of coordinate system
     *
     *  Build all the elements of the structural model and their relationships
     */
    /*!
     * @brief Load a .ml file, Gocad Model3d (b-rep file) and fill the model      
     * @details The BoundaryModelBuilder does the job
     * 
     * @param[in] in Name of the file
     * @return false if file opening failed
     */
    bool BoundaryModel::load_gocad_model3d( const std::string& in )
    {
        std::ifstream input( in.c_str() ) ;
        if( !input ) {
            std::cout << "cannot open file:" << in << std::endl ;
            return false ;
        }

        BoundaryModelBuilderGocad builder( *this ) ;
        builder.load_ml_file( input ) ;
        return true ;
    }

    /*!
     * @brief Modify the model so that it is compatible with a Gocad Model3d
     *   and can be saved in .ml format
     * 
     * @return True if this was a success, False if modifications could not be done.
     */
    bool BoundaryModel::check_model3d_compatibility()
    {
         BoundaryModelBuilder builder( *this ) ;

        /// 1. Check that the Interfaces exist
        if( nb_interfaces() == 0 && nb_surfaces() > 0 ) {
          
            /// If not create one Interface per Surface
            for( index_t i = 0; i < surfaces_.size(); ++i ) {
                // Set name, type, links with other elements
                std::ostringstream name ;
                name << "surface_" << i ;
                index_t id = builder.create_interface( name.str() ) ;
                builder.add_child( BoundaryModelElement::BM_INTERFACE, id, i ) ;
            }

            // Set links from surfaces_ toward interfaces_
            for( index_t i = 0; i < interfaces_.size(); ++i ) {
                builder.set_parent( BoundaryModelElement::BM_SURFACE, interfaces_[i].child( 0 ).id(), i ) ;
            }
            // Is it really useful to have contacts, let's hope not... I am not doing it
        }
        /// 2. Set KeyFacet for Surfaces 
        for( index_t i = 0; i < surfaces_.size(); ++i ) {
            Surface& sp = surfaces_[i] ;
            if( sp.nb_vertices() == 0 ) continue ;
            if( sp.key_facet().is_default() ) {
                builder.set_surface_first_triangle_as_key( sp.id() ) ;
            }
        }

        /// 3. Check that the Universe region exists 
        /// \todo Write some code to create the universe (cf. line 805 to 834 de s2_b_model.cpp)
        if( universe_.name() != "Universe" ) {
            GEO::Logger::err( "" )
                << "The region universe is not defined for the model. IMPLEMENTATION TO DO"
                << std::endl ;
            return false ;
        }

        /// 4. Check that each region has a name and valid surfaces
        for( index_t i = 0; i < regions_.size(); ++i ) {
            BoundaryModelElement& region = regions_[i] ;

            if( region.name() == "" ) {
                std::ostringstream name ;
                name << "region_" << i ;
                builder.set_element_name( BoundaryModelElement::BM_REGION, i, name.str() ) ;
            }
            if( region.nb_boundaries() == 0 ) {
                GEO::Logger::err("") << "The region " << region.name()
                    << " has no Surfaces on its boundary" << std::endl ;
                return false ;
            }
        }

        /// 5. Check that all the surfaces_ of the model are triangulated
        /// \todo Implement a triangulation function in SurfaceMutator         
        for( index_t s = 0; s < nb_surfaces(); s++ ) {
            if( !surfaces_[s].is_triangulated() ) {
                GEO::Logger::err("") << "Surface "<< s << " is not triangulated" << std::endl ;
                return false ;
            }
        }               
        return true ;
    }

    /*!
     * @brief Write a region information in a stream
     * @details Used by function to save the Model in a .ml file
     * 
     * @param[in] count Region index in the file
     * @param[in] region The region to save
     * @param[in,out] out The file output stream
     */
    void save_region(
        index_t count,
        const BoundaryModelElement& region,
        std::ostream& out )
    {
        out << "REGION " << count << "  " << region.name() << " " << std::endl ;
        index_t it = 0 ;

        for( index_t i = 0; i < region.nb_boundaries(); ++i ) {
            out << "  " ;
            if( region.side( i ) ) {
                out << "+" ;
            } else {
                out << "-" ;
            }
            out << region.boundary( i ).id() + 1 ;
            it++ ;
            if( it == 5 ) {
                out << std::endl ;
                it = 0 ;
            }
        }
        out << "  0" << std::endl ;
    }

    /*!
     * @brief Write information for on layer in a stream
     * @details Used by function to save the Model in a .ml file
     * 
     * @param[in] count Index of the layer in the file
     * @param[in] offset Offset of region indices in the file
     * @param[in] layer The layer to write
     * @param[in,out] out The output file stream
     */
    void save_layer(
        index_t count,
        index_t offset,
        const BoundaryModelElement& layer,
        std::ostream& out )
    {
        out << "LAYER " << layer.name() << " " << std::endl ;
        index_t it = 0 ;

        for( index_t i = 0; i < layer.nb_children(); ++i ) {
            out << "  " << layer.child_id( i ) + offset + 1 ;
            it++ ;
            if( it == 5 ) {
                out << std::endl ;
                it = 0 ;
            }
        }
        out << "  0" << std::endl ;
    }

    /*!
     * @brief Write basic header for Gocad coordinate system.
     * @param[in,out] out Output .ml file stream
     */
    void save_coordinate_system( std::ostream& out )
    {
        out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl << "NAME Default"
            << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
            << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl << "ZPOSITIVE Elevation"
            << std::endl << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl ;
    }

    /*!
     * @brief Save the model in a .ml file if it compatible
     *
     * @param[in,out] out Output file stream
     * @return false if the model is not compatible with a Gocad model
     */
    bool BoundaryModel::save_gocad_model3d( std::ostream& out )
    {
        out.precision( 16 ) ;
        if( !check_model3d_compatibility() ) {
            GEO::Logger::err("") << "The BoundaryModel " << name_
                << " cannot be saved in .ml format " << std::endl ;
            return false ;
        }

        // Print Gocad Model3d headers
        out << "GOCAD Model3d 1" << std::endl << "HEADER {" << std::endl << "name:"
            << "model_from_graphite" << std::endl << "}" << std::endl ;

        save_coordinate_system( out ) ;

        // Print the TSurf = Interface information
        for( index_t i = 0; i < interfaces_.size(); ++i ) {
            out << "TSURF " << interfaces_[i].name() << std::endl ;
        }

        index_t count = 1 ;
        // Print the TFace = Surface information
        for( index_t i = 0; i < surfaces_.size(); ++i ) {
            const Surface& s = surfaces_[i] ;
            out << "TFACE " << count << "  " ;
            save_type( out, s.geological_feature() ) ;
            out << " " << s.parent().name() << std::endl ;

            const Surface::KeyFacet kf = s.key_facet() ;
            out << "  " << kf.p0_ << std::endl ;
            out << "  " << kf.p1_ << std::endl ;
            out << "  " << kf.p2_ << std::endl ;

            ++count ;
        }

        index_t offset_layer = count ;
        // Print universe, region, and layer information 
        save_region( count, universe_, out ) ;
        ++count ;

        for( index_t i = 0; i < regions_.size(); ++i ) {
            save_region( count, regions_[i], out ) ;
            ++count ;
        }

        for( index_t i = 0; i < layers_.size(); ++i ) {
            save_layer( count, offset_layer, layers_[i], out ) ;
            ++count ;
        }
        out << "END" << std::endl ;

        // Save the geometry of the Surfaces (TFace), Interface (TSurf) by Interface
        for( index_t i = 0; i < interfaces_.size(); ++i ) {
            const BoundaryModelElement& tsurf = interfaces_[i] ;

            // Header
            out << "GOCAD TSurf 1" << std::endl << "HEADER {" << std::endl << "name:"
                << tsurf.name() << std::endl << "name_in_model_list:" << tsurf.name()
                << std::endl << "}" << std::endl ;
            save_coordinate_system( out ) ;

            out << "GEOLOGICAL_FEATURE " << tsurf.name() << std::endl
                << "GEOLOGICAL_TYPE " ;
            save_type( out, tsurf.geological_feature() ) ;
            out << std::endl ;

            out << "PROPERTY_CLASS_HEADER Z {" << std::endl << "is_z:on" << std::endl
                << "}" << std::endl ;

            // Save surfaces_ geometry
            index_t vertex_count = 1 ;
            index_t offset = vertex_count ;

            std::vector< index_t > bstones ;
            std::vector< index_t > next_vertex ;
            std::set< index_t > set_end_corners ;
            
            bstones.reserve( tsurf.nb_boundaries() ) ;
            next_vertex.reserve( tsurf.nb_boundaries() ) ;

            for( index_t j = 0; j < tsurf.nb_children(); ++j ) {
                offset = vertex_count ;

                const Surface& sp = dynamic_cast< const Surface& >( tsurf.child( j ) ) ;

                out << "TFACE" << std::endl ;
                for( index_t k = 0; k < sp.nb_vertices(); ++k ) {
                    out << "VRTX " << vertex_count << " " << sp.vertex( k ) << std::endl ;
                    vertex_count++ ;
                }

                for( index_t k = 0; k < sp.nb_cells(); ++k ) {
                    out << "TRGL " << sp.surf_vertex_id( k, 0 ) + offset << " "
                        << sp.surf_vertex_id( k, 1 ) + offset << " "
                        << sp.surf_vertex_id( k, 2 ) + offset << std::endl ;
                }

                // Gather information on Corners (BStones) and Lines (getting the next point on the line)
                for( index_t k = 0; k < sp.nb_boundaries(); ++k ) {
                    const Line& cp = dynamic_cast< const Line& >( sp.boundary( k ) ) ;                  
      
                    index_t c = cp.model_vertex_id( 0 ) ;
                    index_t next = cp.model_vertex_id( 1 ) ;

                    // To be sure that we have all corners we need to ensure
                    // that all corners at the end of lines are saved too
                    set_end_corners.insert( 
                        sp.surf_vertex_id( cp.model_vertex_id( cp.nb_vertices()-1 ) ) + offset ) ;                    

                    index_t t = sp.facet_from_model_vertex_ids( c, next ) ;
                    grgmesh_assert( t != NO_ID ) ;

                    index_t i0 = sp.surf_vertex_id( t, 0 ) ;
                    index_t i1 = sp.surf_vertex_id( t, 1 ) ;
                    index_t i2 = sp.surf_vertex_id( t, 2 ) ;

                    index_t p0 = sp.model_vertex_id( i0 ) ;
                    index_t p1 = sp.model_vertex_id( i1 ) ;
                    index_t p2 = sp.model_vertex_id( i2 ) ;

                    index_t c_id = NO_ID ;
                    index_t next_id = NO_ID ;

                    if( p0 == c ) c_id = i0 ;
                    else if( p0 == next ) next_id = i0 ;

                    if( p1 == c ) c_id = i1 ;
                    else if( p1 == next ) next_id = i1 ;

                    if( p2 == c ) c_id = i2 ;
                    else if( p2 == next ) next_id = i2 ;

                    grgmesh_assert( c_id != NO_ID && next_id != NO_ID ) ;

                    bstones.push_back( c_id + offset ) ;
                    next_vertex.push_back( next_id + offset ) ;
                }
            }

            // Print Corners and Lines
            std::vector< index_t > end_corners( 
                set_end_corners.begin(), set_end_corners.end() ) ;
            std::vector< bool > end_corner_to_print( end_corners.size(), true ) ;
            
            for( index_t j = 0; j < bstones.size(); ++j ) {
                out << "BSTONE " << bstones[j] << std::endl ;
                // Determine the corners at the end of the lines that are not saved                
                for( index_t k = 0; k < end_corners.size(); k++ ) {
                    if( bstones[j] == end_corners[k] ) {
                        end_corner_to_print[k] = false ;
                        break ;
                    }
                }
            }
            // Print the corners that were at the beginning of none of the contacts
            // in this Interface
            for( index_t j= 0; j < end_corners.size(); j++ ) {
                if( end_corner_to_print[j] ) {
                    out << "BSTONE " << end_corners[j] << std::endl ;
                }
            }
            // Print the the information to build the lines :
            // index of the vertex at the corner and index of the second vertex on the line
            for( index_t j = 0; j < bstones.size(); ++j ) {
                out << "BORDER " << vertex_count << " " << bstones[j] << " "
                    << next_vertex[j] << std::endl ;
                vertex_count++ ;
            }
            out << "END" << std::endl ;
        }
        return true ;
    }    

    /*!
     * @brief Save the surfaces of the model with their facet attributes into an .eobj file.         
     * @details  WARNING It is supposed that these attributes are integer
     * 
     * @param[in] file_name Name of the file
     *
     * \todo Write generic I/O for attributes on a BoundaryModel
     */
    void BoundaryModel::save_as_eobj_file( const std::string& file_name ) {
        std::ofstream out ;
        out.open( file_name.c_str() );
        if( out.bad() ){
            std::cout << "Error when opening the file: " << file_name.c_str() <<std::endl ;
            return ;
        }
        out.precision( 16 ) ;
        std::vector< index_t > offset( nb_surfaces(), 0 ) ;
        index_t cur_offset = 0 ; 

        // Write vertices once for each surface
        for( index_t s = 0; s < nb_surfaces() ; s++ ) {
            const Surface& S = surface(s) ;
            offset[s] = cur_offset ;
            for( index_t p = 0; p < S.nb_vertices(); p++ ){
                const vec3& V = S.vertex(p) ;
                out << "v" 
                    << " " << V.x 
                    << " " << V.y 
                    << " " << V.z 
                    << std::endl ;
            }
            cur_offset += S.nb_vertices() ;
        }

        // Write the facets for a each surface 
        for( index_t s = 0; s < nb_surfaces() ; s++ ) {
            const Surface& S = surface(s) ;
            for( index_t f = 0; f < S.nb_cells(); f++ ){
                out << "f" << " " ;
                for( index_t v = 0; v < S.nb_vertices_in_facet(f); v++ ){                    
                    out << offset[s] + S.surf_vertex_id( f, v )+1 << " " ;            
                }
                out << std::endl ;
            }         
        }

        // Write facet attributes -- WARNING It is supposed that these attributes are integer       
        {
            std::vector< std::string > names ;        
            facet_attribute_manager_.list_named_attributes(names) ;
            std::vector< AttributeStore* > attribute_stores( names.size(), nil ) ;
                    
            for( index_t i=0; i < names.size(); i++ ) 
            {
                // Crash if we cannot get the
                attribute_stores[i] = facet_attribute_manager_.resolve_named_attribute_store( names[i] ) ;
              
                // Output global information on the attribute
                out << "# attribute "<< names[i] << " facet "
                    // In Graphite - which will read this file - there is a stupid map 
                    // betwen the type name and the name in the file 
                    << "integer" 
                    //<< attribute_stores[i]->attribute_type_id().name() // real attribute type
                    << std::endl ;  
            }
            if( names.size() > 0 ) {
                index_t count = 0 ;
                for( index_t s = 0; s < nb_surfaces() ; s++ ) {
                    const Surface& S = surface(s) ;
                    for( index_t f = 0; f < S.nb_cells(); f++ ){            
                        out << "# attrs f " << count+1 ;
                        for( index_t j = 0; j < names.size(); j++ ) {
                            out << " " << *reinterpret_cast<index_t*>( attribute_stores[j]->data(count) ) ;
                        }
                        out << std::endl ;
                        count++ ;
                    }
                }
            }
        }

    }


     /*!
     * @brief Rebuild a model ???
     * \todo Comment rebuild()... 
     * \todo Valgrind finds errors !!!!!!!
     * @return 
     */
    bool BoundaryModelBuilder::rebuild()
    {
   /*     std::vector< index_t > sp_to_remove ;
        std::vector< index_t > new_sp_id( model_.nb_surfaces() ) ;
        index_t offset = 0 ;
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            if( model_.surfaces_[sp].nb_cells() == 0 ) {
                offset++ ;
                sp_to_remove.push_back( sp ) ;
            } else {
                model_.surfaces_[sp - offset] = model_.surfaces_[sp] ;
            }
            new_sp_id[sp] = sp - offset ;
        }
        if( offset == 0 ) return false ;
        model_.surfaces_.erase( model_.surfaces_.end() - offset,
            model_.surfaces_.end() ) ;

        offset = 0 ;
        std::vector< index_t > cp_to_remove ;
        std::vector< index_t > new_cp_id( model_.nb_lines() ) ;
        for( index_t cp = 0; cp < model_.nb_lines(); cp++ ) {
            BoundaryModelElement& line = model_.lines_[cp] ;
            index_t nb_sp_removed = 0 ;
            for( index_t sp = 0; sp < line.nb_in_boundary(); sp++ ) {
                if( Utils::contains( sp_to_remove,
                    line.in_boundary_id( sp ) ) ) {
                    nb_sp_removed++ ;
                } else {
                    line.in_boundary_[sp - nb_sp_removed] =
                        new_sp_id[line.in_boundary_[sp]] ;
                }
            }
            if( nb_sp_removed > 0 ) {
                if( nb_sp_removed == line.nb_in_boundary() ) {
                    offset++ ;
                    cp_to_remove.push_back( cp ) ;
                    continue ;
                } else {
                    line.in_boundary_.erase(
                        line.in_boundary_.end() - nb_sp_removed,
                        line.in_boundary_.end() ) ;
                }
            }
            if( offset > 0 ) {
                model_.lines_[cp - offset] = model_.lines_[cp] ;
            }
            new_cp_id[cp] = cp - offset ;
        }
        if( offset > 0 ) {
            model_.lines_.erase( model_.lines_.end() - offset,
                model_.lines_.end() ) ;
        }

        // Build the contacts
        // Update surfaces
        offset = 0 ;
        std::vector< index_t > s_to_remove ;
        std::vector< index_t > new_s_id( model_.nb_interfaces() ) ;
        for( index_t s = 0; s < model_.nb_interfaces(); s++ ) {
            BoundaryModelElement& surface = model_.interfaces_[s] ;
            surface.boundaries_.clear() ;
            index_t nb_sp_removed = 0 ;
            for( index_t sp = 0; sp < surface.nb_children(); sp++ ) {
                if( Utils::contains( sp_to_remove, surface.child_id( sp ) ) ) {
                    nb_sp_removed++ ;
                } else {
                    surface.children_[sp - nb_sp_removed] =
                        new_sp_id[surface.children_[sp]] ;
                }
            }
            if( nb_sp_removed > 0 ) {
                if( nb_sp_removed == surface.nb_children() ) {
                    offset++ ;
                    s_to_remove.push_back( s ) ;
                    continue ;
                } else {
                    surface.children_.erase( surface.children_.end() - nb_sp_removed,
                        surface.children_.end() ) ;
                }
            }
            if( offset > 0 ) {
                model_.interfaces_[s - offset] = model_.interfaces_[s] ;
            }
            new_s_id[s] = s - offset ;
        }
        if( offset > 0 ) {
            model_.interfaces_.erase( model_.interfaces_.end() - offset,
                model_.interfaces_.end() ) ;
        }
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            BoundaryModelElement& surface = model_.surfaces_[sp] ;
            surface.parent_ = new_s_id[surface.parent_] ;
            for( index_t cp = 0; cp < surface.nb_boundaries(); cp++ ) {
                surface.boundaries_[cp] =
                    new_cp_id[surface.boundaries_[cp]] ;
            }
        }
        model_.contacts_.clear() ;
        build_contacts() ;

        // Then finish the job (the order matters)
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            model_.surfaces_[sp].boundaries_.clear() ;
        }
        offset = 0 ;
        for( index_t r = 0; r < model_.nb_regions(); r++ ) {
            BoundaryModelElement& region = model_.regions_[r] ;
            index_t nb_sp_removed = 0 ;
            for( index_t sp = 0; sp < region.nb_boundaries(); sp++ ) {
                if( Utils::contains( sp_to_remove, region.boundary_id( sp ) ) ) {
                    nb_sp_removed++ ;
                } else {
                    region.boundaries_[sp - nb_sp_removed] =
                        new_sp_id[region.boundaries_[sp]] ;
                    region.sides_[sp - nb_sp_removed] = region.sides_[sp] ;
                }
            }
            if( nb_sp_removed > 0 ) {
                if( nb_sp_removed == region.nb_boundaries() ) {
                    offset++ ;
                    continue ;
                } else {
                    region.sides_.erase( region.sides_.end() - nb_sp_removed,
                        region.sides_.end() ) ;
                    region.boundaries_.erase(
                        region.boundaries_.end() - nb_sp_removed,
                        region.boundaries_.end() ) ;
                }
            }
            if( offset > 0 ) {
                model_.regions_[r - offset] = model_.regions_[r] ;
            }
        }
        if( offset > 0 ) {
            model_.regions_.erase( model_.regions_.end() - offset,
                model_.regions_.end() ) ;
        }
        end_surfaces() ;

        for( index_t i = 0; i < model_.nb_contacts(); ++i ) {
            for( index_t j = 0; j < model_.contacts_[i].nb_in_boundary();
                ++j ) {
                index_t b = model_.contacts_[i].in_boundary_id( j ) ;
                add_interface_boundary( b, i ) ;
            }
        }

        end_lines() ;
        end_contacts() ;

        offset = 0 ;
        std::vector< index_t > co_to_remove ;
        std::vector< index_t > new_co_id( model_.nb_corners() ) ;
        for( index_t co = 0; co < model_.nb_corners(); co++ ) {
            BoundaryModelElement& corner = model_.corners_[co] ;
            index_t nb_cp_removed = 0 ;
            for( index_t cp = 0; cp < corner.nb_in_boundary(); cp++ ) {
                if( Utils::contains( cp_to_remove, corner.in_boundary_id( cp ) ) ) {
                    nb_cp_removed++ ;
                } else {
                    corner.in_boundary_[cp - nb_cp_removed] =
                        new_cp_id[corner.in_boundary_[cp]] ;
                }
            }
            if( nb_cp_removed > 0 ) {
                if( nb_cp_removed == corner.nb_in_boundary() ) {
                    offset++ ;
                    co_to_remove.push_back( co ) ;
                    continue ;
                } else {
                    corner.in_boundary_.erase(
                        corner.in_boundary_.end() - nb_cp_removed,
                        corner.in_boundary_.end() ) ;
                }
            }
            if( offset > 0 ) {
                model_.corners_[co - offset] = model_.corners_[co] ;
            }
            new_co_id[co] = co - offset ;
        }
        if( offset > 0 ) {
            model_.corners_.erase( model_.corners_.end() - offset,
                model_.corners_.end() ) ;
        }

        for( index_t cp = 0; cp < model_.nb_lines(); cp++ ) {
            BoundaryModelElement& line = model_.lines_[cp] ;
            line.boundaries_[0] = new_co_id[line.boundaries_[0]] ;
            line.boundaries_[1] = new_co_id[line.boundaries_[1]] ;
        }

        for( index_t c = 0; c < model_.nb_contacts(); c++ ) {
            BoundaryModelElement& contact = model_.contacts_[c] ;
            for( index_t co = 0; co < contact.nb_boundaries(); co++ ) {
                contact.boundaries_[co] = new_co_id[contact.boundaries_[co]] ;
            }
        }

        for( index_t l = 0; l < model_.nb_layers(); l++ ) {
            BoundaryModelElement& layer = model_.layers_[l] ;
            index_t nb_sp_removed = 0 ;
            for( index_t sp = 0; sp < layer.nb_boundaries(); sp++ ) {
                if( Utils::contains( sp_to_remove, layer.boundary_id( sp ) ) ) {
                    nb_sp_removed++ ;
                } else {
                    layer.boundaries_[sp - nb_sp_removed] =
                        new_sp_id[layer.boundaries_[sp]] ;
                }
            }
            if( nb_sp_removed > 0 ) {
                layer.boundaries_.erase( layer.boundaries_.end() - nb_sp_removed,
                    layer.boundaries_.end() ) ;
            }
        }

        update_all_ids() ;
        return true ;*/

        grgmesh_assert_not_reached ;
        /// \todo Implement the rebuild function for the BoundaryModel
        return false ;
    }

    /*!
     * @brief Copy macro information from a model
     * @details Copy the all the model elements and their relationship ignoring their geometry
     * 
     * @param[in] from Model to copy the information from
     */
    void BoundaryModelBuilder::copy_macro_topology( const BoundaryModel& from )
    {
        model_.name_ = from.name_ ;
        model_.corners_.resize( from.nb_corners(), Corner( &model_ ) ) ;
        model_.lines_.resize( from.nb_lines(), Line( &model_ ) ) ;
        model_.surfaces_.resize( from.nb_surfaces(), Surface( &model_ ) ) ;
        model_.regions_.resize( from.nb_regions(), BoundaryModelElement( &model_, BoundaryModelElement::BM_REGION ) ) ;
        model_.layers_.resize( from.nb_layers(), BoundaryModelElement( &model_, BoundaryModelElement::BM_LAYER ) ) ;
        model_.contacts_.resize( from.nb_contacts(), BoundaryModelElement( &model_, BoundaryModelElement::BM_CONTACT ) ) ;
        model_.interfaces_.resize( from.nb_interfaces(), BoundaryModelElement( &model_, BoundaryModelElement::BM_INTERFACE ) ) ;
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_corners(); i++ ) {
            model_.corners_[i].copy_macro_topology( from.corner( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_lines(); i++ ) {
            model_.lines_[i].copy_macro_topology( from.line( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_surfaces(); i++ ) {
            model_.surfaces_[i].copy_macro_topology( from.surface( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_layers(); i++ ) {
            model_.layers_[i].copy_macro_topology( from.layer( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_regions(); i++ ) {
            model_.regions_[i].copy_macro_topology( from.region( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_contacts(); i++ ) {
            model_.contacts_[i].copy_macro_topology( from.contact( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_interfaces(); i++ ) {
            model_.interfaces_[i].copy_macro_topology( from.one_interface( i ), model_ ) ;
        }
        model_.universe_.copy_macro_topology( from.universe_, model_ ) ;
    }

    /*!
     * @brief Update the indices stored by each element of the model \
     * according to its actual position in the corresponding vector in the model
     */
    void BoundaryModelBuilder::update_all_ids()
    {
        for( index_t co = 0; co < model_.nb_corners(); co++ ) {
            model_.corners_[co].set_id( co ) ;
        }
        for( index_t cp = 0; cp < model_.nb_lines(); cp++ ) {
            model_.lines_[cp].set_id( cp ) ;
        }
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            model_.surfaces_[sp].set_id( sp ) ;
        }
        for( index_t c = 0; c < model_.nb_contacts(); c++ ) {
            model_.contacts_[c].set_id( c ) ;
        }
        for( index_t s = 0; s < model_.nb_interfaces(); s++ ) {
            model_.interfaces_[s].set_id( s ) ;
        }
        for( index_t r = 0; r < model_.nb_regions(); r++ ) {
            model_.regions_[r].set_id( r );
        }
        for( index_t l = 0; l < model_.nb_layers(); l++ ) {
            model_.layers_[l].set_id( l ) ;
        }
    }

    /*!
     * @brief Remove duplicates in the model vertices 
     * @details When reading the file the vertices are duplicated between the different Surfaces,
     * and new vertices are added for Corners. 
     * Compute the duplicates inside the vertices_ vector - Update the vertex vector - 
     * and update the reference to vertices for all model Corners and Surfaces
     * 
     * \todo Check the compute_unique_kdtree function - It is not correct the number of neighbors 
     * is predefined 
     */
    void BoundaryModelBuilder::make_vertices_unique()
    {
        MakeUnique unique( model_.vertices_ ) ;
        unique.unique( 5 ) ;
        model_.vertices_.resize(0) ;
        unique.unique_points( model_.vertices_ ) ;
        const std::vector< index_t >& old2new = unique.indices() ;

        for( index_t s = 0; s < model_.nb_surfaces(); s++ ) {
            Surface& surface = model_.surfaces_[s] ;
            for( index_t p = 0; p < surface.nb_vertices(); p++ ) {
                surface.set_vertex( p, old2new[surface.model_vertex_id(p)] ) ;
            }
        }
        for( index_t l = 0; l < model_.nb_lines(); l++ ) {
            Line& line = model_.lines_[l] ;
            for( index_t p = 0; p < line.nb_vertices(); p++ ) {
                line.set_vertex( p, old2new[line.model_vertex_id(p)] ) ;
            }
        }
        for( index_t co = 0; co < model_.nb_corners(); co++ ) {
            Corner& corner = model_.corners_[co] ;
            corner.set_vertex( old2new[corner.model_vertex_id()] ) ;
        }       
    }

    /*!
     * @brief Creates a element of the given type and add it to the correct vector
     * The BoundaryModelElement is created from its type and its index
     *
     * @param[in] type Type of the element to create
     * @return The index of the created element
     */
    index_t BoundaryModelBuilder::create_element( BoundaryModelElement::BM_TYPE type ) {
        index_t id = model_.nb_elements( type ) ;
        grgmesh_assert( id != NO_ID ) ;
        switch( type ) {
        case BoundaryModelElement::BM_CORNER: 
            {
                model_.corners_.push_back( Corner( &model_, id ) ) ;
                break ;
            }
        case BoundaryModelElement::BM_LINE:
            {
                model_.lines_.push_back( Line( &model_, id ) ) ;
                break ;
            }
        case BoundaryModelElement::BM_SURFACE:
            {
                model_.surfaces_.push_back( Surface( &model_, id ) ) ;
                break ;
            }
        case BoundaryModelElement::BM_REGION:
            {
                model_.regions_.push_back( BoundaryModelElement( &model_, BoundaryModelElement::BM_REGION, id ) ) ;
                break ;
            }
        case BoundaryModelElement::BM_CONTACT:
            {
                model_.contacts_.push_back( BoundaryModelElement( &model_,BoundaryModelElement::BM_CONTACT, id ) ) ;
                break ;
            }
        case BoundaryModelElement::BM_INTERFACE:
            {
                model_.interfaces_.push_back( BoundaryModelElement( &model_,BoundaryModelElement::BM_INTERFACE, id ) ) ;
                break ;
            }
        case BoundaryModelElement::BM_LAYER:
            {
                model_.layers_.push_back( BoundaryModelElement( &model_, BoundaryModelElement::BM_LAYER, id ) ) ;
                break ;
            }
        default:   
            grgmesh_assert_not_reached ;
        }
        return id ;
    }

    /*!
     * @brief Get the index of the Corner at a given model point
     * @param[in] p_id Index of the point
     * @return NO_ID or the index of the Corner
     */
    index_t BoundaryModelBuilder::find_corner( index_t p_id ) const 
    {
        for( index_t i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).model_vertex_id() == p_id ) return i ;
        }
        return NO_ID ;
    }
   
     /*!
     * @brief Create a corner at a vertex.
     * 
     * @param[in] index Index of the vertex in the model
     * @return Index of the Corner
     */
    index_t BoundaryModelBuilder::create_corner( index_t index )
    {
       index_t id = create_element( BoundaryModelElement::BM_CORNER ) ; 
       set_corner( id, index ) ;       
       return id ;
    }

     /*!
     * @brief Find or create a corner at a vertex.
     * 
     * @param[in] index Index of the vertex in the model
     * @return Index of the Corner
     */
    index_t BoundaryModelBuilder::find_or_create_corner( index_t index )
    {
        index_t result = find_corner( index ) ;
        if( result != NO_ID ) return result ;
        else return create_corner( index ) ;
    }

    /*!
     * @brief Looks for a line in the model
     * 
     * @param[in] corner0 Starting corner index
     * @param[in] corner1 Ending corner index
     * @param[in] vertices Indices of the vertices on the line
     * @return NO_ID or the index of the Line
     */
    index_t BoundaryModelBuilder::find_line(
        const std::vector< index_t >& vertices ) const
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            if( model_.line(i).equal( vertices ) ) 
                return i ;
        }
        return NO_ID ;
    }

    /*! 
     * @brief Add a Line knowing only the indices of its points and set its boundary corners 
     * The corners MUST already exist.
     * Used in Geomodeling to convert a surface to a model
     */
    index_t BoundaryModelBuilder::create_line( const std::vector< index_t >& points ) {
        index_t id = create_element( BoundaryModelElement::BM_LINE ) ;
        set_line( id, points ) ;
               
        // Find the indices of the corner at both extremities
        index_t c0 = find_corner( points.front() ) ;
        index_t c1 = find_corner( points.back() ) ;
        grgmesh_assert( c0 != NO_ID && c1 != NO_ID ) ; // Mouais on pourrait peut être le virer celui là
        add_element_boundary( BoundaryModelElement::BM_LINE, id, c0 ) ;
        if( c1 != c0 ) add_element_boundary( BoundaryModelElement::BM_LINE, id, c1 ) ;         

        return id ;
    }

     /*!
     * @brief Find or create a line 
     * 
     * @param[in] corner0 Starting corner index
     * @param[in] corner1 Ending corner index
     * @param[in] vertices Indices of the vertices on this Line
     * @return Index of the Line
     */
    index_t BoundaryModelBuilder::find_or_create_line(
        const std::vector< index_t >& vertices )
    {      
        index_t result = find_line( vertices ) ;
        if( result != NO_ID ) return result ;
        else return create_line( vertices ) ;
    }

  
    /*!
     * @brief Create a contact between the given Interfaces
     * The name of the contact is determined from the names of the interfaces.
     * 
     * @return Index of the Surface in the surfaces_ vector 
     */
    index_t BoundaryModelBuilder::create_surface()
    {      
        return create_element( BoundaryModelElement::BM_SURFACE ) ; 
    }

    /*!
     * @brief Find a Contact
     * @param[in] interfaces Indices of the Interfaces determining the contact
     * @return NO_ID or index of the contact
     */
    index_t BoundaryModelBuilder::find_contact( const std::vector< index_t >& interfaces ) const
    {
        std::vector< const BoundaryModelElement* > comp( interfaces.size() ) ;
        for( index_t i = 0; i < interfaces.size(); ++i ) {
            comp[i] = &model_.one_interface(interfaces[i]) ;
        }

        for( index_t i = 0; i < model_.nb_contacts(); ++i ) {
            if( comp.size() == model_.contact(i).nb_in_boundary() ) {
                bool equal = true ;
                for( index_t j = 0; j < model_.contact(i).nb_in_boundary(); j++ ) {
                    if( comp[j] != &model_.contact(i).in_boundary( j ) ) {
                        equal = false ; break ;
                    }
                }
                if( equal ) return i ;
            }
        }
        return NO_ID ;
    }

    /*!
     * @brief Create a contact between the given Interfaces
     * The name of the contact is determined from the names of the interfaces.
     * 
     * @param[in] interfaces Indices of the intersecting interfaces
     * @return Index of the Contact
     */
    index_t BoundaryModelBuilder::create_contact( const std::vector< index_t >& interfaces ) {        
        // Create a name for this contact
        std::string name = "contact_" ;
        for( index_t i = 0; i < interfaces.size(); ++i ) {
            name += model_.interfaces_[interfaces[i]].name() ;
            name += "_" ;
        }
        
        index_t id = create_element( BoundaryModelElement::BM_CONTACT ) ;
        set_element_name( BoundaryModelElement::BM_CONTACT, id, name ) ;
        
        /*for( index_t i = 0; i < interfaces.size(); ++i ) {
            add_element_in_boundary( BoundaryModelElement::BM_CONTACT, id, interfaces[i] ) ;
        }*/        
        return id ;
    }

    
    /*!
     * @brief Find or create a contact between given Interfaces
     * 
     * @param[in] interfaces Indices of the intersecting interfaces
     * @return Index of the Contact
     */
    index_t BoundaryModelBuilder::find_or_create_contact(
        const std::vector< index_t >& interfaces )
    {
        index_t result = find_contact( interfaces ) ;
        if( result != NO_ID ) return result ;
        else return create_contact( interfaces ) ;        
    }

     /*!
     * @brief Get the index of an Interface from its name
     * 
     * @param[in] name Name of the Interface
     * @return Index of the interface in the model, NO_ID if not found.
     */
    index_t BoundaryModelBuilder::find_interface( const std::string& name ) const
    {
        for( index_t i = 0; i < model_.nb_interfaces(); ++i ) {
            if( model_.one_interface(i).name() == name ) {
                return i ;
            }
        }
        return NO_ID ;
    }


    
    /*!
     * @brief Create a new Interface 
     * 
     * @param[in] name Name of the interface     
     * @param[in] type Type of the interface
     * @return The Interface index.
     */
    index_t BoundaryModelBuilder::create_interface(
        const std::string& name,
        BoundaryModelElement::GEOL_FEATURE type )
    {
        index_t id = create_element( BoundaryModelElement::BM_INTERFACE ) ;
        set_element_geol_feature( BoundaryModelElement::BM_INTERFACE, id, type ) ;
        set_element_name( BoundaryModelElement::BM_INTERFACE, id, name ) ;
        return id ;
    }

     /*
    * @brief Adds an empty region to the model
    *
    *  Used in Geomodeling to convert a surface to a model
    */
    index_t BoundaryModelBuilder::create_region() {
        return create_element( BoundaryModelElement::BM_REGION ) ;
    }

    /*!
     * @brief Adds a new region to the model
     *
     * @param[in] name Name of the region
     * @param[in] boundaries Indices of the surfaces on the region boundary, plus indication on which
     *            side of the surface is the region
     * @return Index of the created region
     */
    index_t BoundaryModelBuilder::create_region(
        const std::string& name,
        const std::vector< std::pair< index_t, bool > >& boundaries )
    {
        index_t id = create_element( BoundaryModelElement::BM_REGION ) ;
        set_element_name( BoundaryModelElement::BM_REGION, id, name ) ;
        for( index_t i = 0; i < boundaries.size(); ++i ) {            
            add_element_boundary( BoundaryModelElement::BM_REGION, id, boundaries[i].first, boundaries[i].second ) ;
        }
        return id ;
    }

    
    /*!
     * @brief Creates a new empty Layer with the given name 
     *
     * @param[in] name Name of the layer
     * @return The layer index
     */
    index_t BoundaryModelBuilder::create_layer(
        const std::string& name )
    {
        index_t id = create_element( BoundaryModelElement::BM_LAYER ) ;
        set_element_name( BoundaryModelElement::BM_LAYER, id, name ) ;
        return id ;
    }
    
    
    /*!
     * @brief Fill the model universe_
     *
     * @param[in] boundaries Indices of the surfaces on the model boundary
     * plus indication on which side of the surface is universe_ ( outside of the model )
     */
    void BoundaryModelBuilder::set_universe(
        const std::vector< std::pair< index_t, bool > >& boundaries )
    {
        model_.universe_.set_name( "Universe" ) ;
        model_.universe_.set_element_type( BoundaryModelElement::BM_REGION ) ;
        model_.universe_.set_model( &model_ ) ;

        for( index_t i = 0; i < boundaries.size(); ++i ) {
            grgmesh_assert( boundaries[i].first < model_.nb_surfaces() ) ;
            model_.universe_.add_boundary( boundaries[i].first,
                boundaries[i].second ) ;
            // If this surface have no type, set it at BoundaryModelElement::VOI
            model_.surfaces_[boundaries[i].first].set_geological_feature( BoundaryModelElement::VOI ) ;
        }
    }

    /*!
     * @brief Remove the universe from the vector of regions and update their indices
     *
     * @param[in] id Index of the Universe region in the regions_ vector
     */
    void BoundaryModelBuilder::remove_universe_from_regions( index_t id ) 
    {
        for( index_t i = 0; i < model_.nb_regions(); ++i ) {
            index_t cur_id = model_.region(i).id() ;            
            if( i > id ) model_.regions_[i].set_id( cur_id-1 ) ;
        }
        model_.regions_.erase( model_.regions_.begin() + id ) ;
    }

   
    /*!
     * @brief Set the vertex for a Corner
     *
     * @param[in] corner_id Index of the corner
     * @param[in] vertex_id Index of the vertex in the model
     */
    void BoundaryModelBuilder::set_corner( index_t corner_id,  index_t vertex_id ) 
    {
        grgmesh_debug_assert( vertex_id < model_.nb_vertices() ) ;
        model_.corners_.at(corner_id).set_vertex( vertex_id ) ;
    }


    /*!
     * @brief Set one Line vertices
     *
     * @param[in] id Line index 
     * @param[in] vertices Indices of the vertices on the line
     */
    void BoundaryModelBuilder::set_line( 
        index_t id, 
        const std::vector< index_t >& vertices )
    {
        grgmesh_assert( id < model_.nb_lines() ) ;
        model_.lines_[id].set_vertices( vertices );
    }

   
    /*!
     * @brief Set the vertices and facets for a surface
     * @details Facet adjencies can be set too.
     *
     * @param[in] surface_id Index of the surface
     * @param[in] vertices Model indices of the vertices 
     * @param[in] facets Indices in the vertices vector to build facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     * @param[in] surface_adjacencies Adjacent facet (size of facet_ptr)
     */
    void BoundaryModelBuilder::set_surface_geometry(
        index_t surface_id,
        const std::vector< index_t >& vertices,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr,
        const std::vector< index_t >& surface_adjacencies )
    {        
        if( facets.size() == 0 ) return ;

        model_.surfaces_[surface_id].set_geometry( vertices, facets, facet_ptr ) ;

        if( surface_adjacencies.empty() )
            set_surface_adjacencies( surface_id ) ;
        else
            model_.surfaces_[surface_id].set_adjacent( surface_adjacencies ) ;
    }
   
    /*!
     * @brief Compute and set the adjacencies between the facets 
     * @details The adjacent facet is given for each vertex of each facet for the edge
     * starting at this vertex.
     * If there is no neighbor inside the same Surface adjacent is set to NO_ADJACENT
     *
     * @param[in] surface_id Index of the surface
     */
    void BoundaryModelBuilder::set_surface_adjacencies( index_t surface_id ) 
    { 
        Surface& S = model_.surfaces_[surface_id] ;
        grgmesh_assert( S.nb_cells() > 0  ) ;

        std::vector< index_t > adjacent ;
        adjacent.resize( S.facet_end( S.nb_cells()-1 ), Surface::NO_ADJACENT ) ;
              
        index_t nb_facets = S.nb_cells() ;
        index_t nb_vertices = S.nb_vertices() ;        
    
        // Allocate some space to store the ids of facets around each vertex
        std::vector< index_t > toto ;
        toto.reserve( 10 ) ;
        std::vector< std::vector< index_t > > 
            vertex_to_facets( nb_vertices, toto ) ;              

        for( index_t f = 0; f < nb_facets; ++f )
        {
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                vertex_to_facets[S.surf_vertex_id( f, v )].push_back( f ) ;
            }
        }
        for( index_t p = 0; p < nb_vertices; ++p ){
            std::sort( vertex_to_facets[p].begin(), vertex_to_facets[p].end() ) ;
        }

        for( index_t f = 0; f < nb_facets; ++f )
        {
            for( index_t v = 0; v < S.nb_vertices_in_facet(f); v++ )
            {                
                index_t cur = S.surf_vertex_id(f, v) ;
                index_t prev = S.surf_vertex_id(f, S.prev_in_facet(f,v)) ;

                const std::vector< index_t >& f_prev = vertex_to_facets[prev] ;
                const std::vector< index_t >& f_cur = vertex_to_facets[cur] ;

                std::vector< index_t > inter(
                    std::min( f_prev.size(), f_cur.size() ) ) ;
                index_t end = std::set_intersection( f_prev.begin(),
                    f_prev.end(), f_cur.begin(), f_cur.end(), inter.begin() )
                    - inter.begin() ;

                if( end == 2 ) {
                    // There is one neighbor
                    index_t f2 = inter[0] == f ? inter[1] : inter[0] ;
                    adjacent[ S.facet_begin(f) + S.prev_in_facet(f,v) ] = f2 ;
                } else {
                    grgmesh_debug_assert( end == 1 ) ;
                }
            }
        }
        S.set_adjacent( adjacent ) ;
    }

    /*! 
     * @brief Basic check of the validity of a BoundaryModelElement     
     * 
     * \todo Write meaningful message when the test fails ? 
     */
    bool BoundaryModelBuilder::check_basic_element_validity( const BoundaryModelElement& E ) const 
    {
        /// Verify that E points to the right BoundaryModel
        /// that its index and type are the right one.
        if( &E.model() != &model_ ) return false ;
        if( E.element_type() == BoundaryModelElement::BM_NO_TYPE ) return false ;
        if( E.id() == NO_ID ) return false ;
        if( E.id() >= model_.nb_elements( E.element_type() ) ) return false ;
        if( !(model_.element(E.element_type(), E.id()) == E) ) return false ;


        /// Verify that the stored model vertex indices are in a valid range
        for( index_t i = 0; i < E.nb_vertices(); ++i ){            
            if( E.model_vertex_id(i) == NO_ID && 
                E.model_vertex_id(i) >= model_.nb_vertices() ) return false ;
        }
        return true ;
    }

    /*!
     * @brief Complete missing information in BoundaryModelElements
     * boundaries - in_boundary - parent - children 
     *
     * WARNING : Une certaine cohérence est supposée exister pour le remplissage des éléments
     * tous les éléments d'un même type ont été initialisé avec les mêmes informations.
     */
    bool BoundaryModelBuilder::complete_element_connectivity() {             
        // Lines
        if( model_.nb_lines() > 0 ) {
            if( model_.line(0).nb_boundaries() == 0 ){
                fill_elements_boundaries( BoundaryModelElement::BM_LINE ) ;
            }
            if( model_.line(0).nb_in_boundary() == 0 ){
                fill_elements_in_boundaries( BoundaryModelElement::BM_LINE ) ;    
            }
            if( model_.line(0).parent_id() == NO_ID && model_.nb_contacts() > 0 ){
                fill_elements_parent( BoundaryModelElement::BM_LINE ) ;
            }
        }
        // Corners
        if( model_.nb_corners() > 0 && model_.corner(0).nb_in_boundary() == 0 ) {
            // Info from line boundaries is used here and should be available
            fill_elements_in_boundaries( BoundaryModelElement::BM_CORNER ) ;
        }    
        // Surfaces - There MUST be at least one
        if( model_.surface(0).nb_boundaries() == 0 ) {            
            fill_elements_boundaries( BoundaryModelElement::BM_SURFACE ) ;
        }        
        if( model_.surface(0).nb_in_boundary() == 0 ) {            
            fill_elements_in_boundaries( BoundaryModelElement::BM_SURFACE ) ;    
        }
        if( model_.surface(0).parent_id() == NO_ID ) {            
            fill_elements_parent( BoundaryModelElement::BM_SURFACE ) ;    
        }
        // Regions
        if( model_.nb_regions() > 0 ) {
            if( model_.region(0).nb_boundaries() == 0 ) {
                fill_elements_boundaries( BoundaryModelElement::BM_REGION ) ;
            }
            if( model_.region(0).parent_id() == NO_ID && model_.nb_layers() > 0 ){
                fill_elements_parent( BoundaryModelElement::BM_REGION ) ;
            }
        }
        // Contacts
        if( model_.nb_contacts() > 0 &&
            model_.contact(0).nb_children() == 0 ) {
            fill_elements_children( BoundaryModelElement::BM_CONTACT ) ;
        }
        // Interfaces
        if( model_.nb_interfaces() > 0 &&
            model_.one_interface(0).nb_children() == 0 ) {
            fill_elements_children( BoundaryModelElement::BM_INTERFACE ) ;
        }
        // Layers
        if( model_.nb_layers() > 0 &&
            model_.layer(0).nb_children() == 0 ) {
                fill_elements_children( BoundaryModelElement::BM_LAYER ) ;
        }
        return true ;
    }

    /*!
     * @brief
     *
     * We have a problem if this is called for regions - No way yet to know the surface orientation
     * 
     * @param type 
     *
     */
    void BoundaryModelBuilder::fill_elements_boundaries( BoundaryModelElement::BM_TYPE type ) 
    {
        BME::BM_TYPE b_type = BME::boundary_type( type ) ;
        if( b_type != BME::BM_NO_TYPE ) {            
            for( index_t i = 0; i < model_.nb_elements( b_type ); ++i ) {
                const BME& b = model_.element( b_type, i ) ;
                for( index_t j = 0; j < b.nb_in_boundary(); ++j ) {
                    add_element_boundary( type, b.in_boundary_id(j), i ) ;
                }
            }
        }
    }

    void BoundaryModelBuilder::fill_elements_in_boundaries( BoundaryModelElement::BM_TYPE type ) 
    {
        BME::BM_TYPE in_b_type = BME::in_boundary_type( type ) ;
        if( in_b_type != BME::BM_NO_TYPE ) {            
            for( index_t i = 0; i < model_.nb_elements( in_b_type ); ++i ) {
                const BME& in_b = model_.element( in_b_type, i ) ;
                for( index_t j = 0; j < in_b.nb_boundaries(); ++j ) {
                    add_element_in_boundary( type, in_b.boundary_id(j), i ) ;
                }
            }
        }        
    }

    void BoundaryModelBuilder::fill_elements_parent( BoundaryModelElement::BM_TYPE type ) 
    {
        BME::BM_TYPE p_type = BME::parent_type( type ) ;
        if( p_type != BME::BM_NO_TYPE ) {
            for( index_t i = 0; i < model_.nb_elements( p_type ); ++i ) {
                const BME& p = model_.element( p_type, i ) ;
                for( index_t j = 0; j < p.nb_children(); ++j ) {
                    set_parent( type, p.child_id(j), i ) ;
                }
            }
        }              
    }

    void BoundaryModelBuilder::fill_elements_children( BoundaryModelElement::BM_TYPE type ) 
    {
        BME::BM_TYPE c_type = BME::child_type( type ) ;
        if( c_type != BME::BM_NO_TYPE ) {
            for( index_t i = 0; i < model_.nb_elements( c_type ); ++i ) {
                index_t parent = model_.element( c_type, i ).parent_id() ;
                if( parent != NO_ID ) add_child( type, parent, i ) ;
            }
        }
    }                                                                          


     /*!
     * @brief Last function to call when building a model
     *
     * 
     * @details Output information on the model and initialise the
     * nb_facets_in_surfaces_ vector
     * 
     * @return False if the model is not valid and cannot be fixed
     * otherwise returns true.
     *
     *  \todo FULL CHECK AND FIX OF THE MODEL CORRECTNESS !!
        \todo Trade the end_something functions in the BoundaryModelBuilder
        functions for a smart end_model function
        that checks model validity and complete all missing parts
     */
    bool BoundaryModelBuilder::end_model()
    {
        // The name should exist
        if( model_.name() == "" ) set_model_name("model_default_name") ;
        // There must be at least 3 vertices
        if( model_.nb_vertices() == 0 ) return false ;
        // And at least one surface
        if( model_.nb_surfaces() == 0 ) return false ;
        
        // The Universe
        /// \todo Write some code to create the universe (cf. line 805 to 834 de s2_b_model.cpp)



        // Fill the nb_elements_per_type_
        // There is 7 different vectors storing the elements        
        {
            index_t count = 0 ;
            model_.nb_elements_per_type_.push_back( count ) ;
            // UNSAFE - mais tant pis
            for( index_t type = BME::BM_CORNER; type < BME::BM_NO_TYPE; type++ ) {
                count += model_.nb_elements( (BME::BM_TYPE) type ) ;
                model_.nb_elements_per_type_.push_back( count ) ;
            }
        }
        complete_element_connectivity() ;
        

        /// 1. Check that the basics to have a valid element
        ///    and that required connectivity relationships are filled

        /// For all the elements of the BoundaryModel check that they have all 
        /// the required attributes for the types

        /// See the static functions   ***_type( BM_TYPE ) in class 
        
        // Sans doute un peu longuet vu l'implémentation mais en même 
        // temps beaucoup mois de lignes de code
        for( index_t i = 0; i < model_.nb_elements( BME::BM_ALL_TYPES ); ++i ){
            const BME& E = model_.element( BME::BM_ALL_TYPES, i ) ;

            if( !check_basic_element_validity( E ) ){
                return false ;
            }

            BME::BM_TYPE T = E.element_type() ;

            if( BME::boundary_type( T ) != BME::BM_NO_TYPE ) {
                if( T != BME::BM_SURFACE ) {
                    // A closed surface - bubble might have no boundary
                    // The others Line - and Region must have one
                    if( E.nb_boundaries() == 0 ){
                        return false ;
                    }
                }
            }

            // In_boundary
            if( BME::boundary_type( T ) != BME::BM_NO_TYPE ) {
                if( E.nb_boundaries() == 0 ){
                    return false ;
                }
            }

            // Parent - High level elements are not mandatory
            // But if the model has elements of the parent type, the element must have a parent
            if( BME::parent_type( T ) != BME::BM_NO_TYPE ) {
                if( E.parent_id() == NO_ID && 
                    model_.nb_elements( BME::parent_type(T) ) > 0 ){
                        return false ;
                }
            }

            // Children
            if( BME::child_type( T ) != BME::BM_NO_TYPE ) {
                if( E.nb_children() == 0 ){
                    return false ;
                }
            }

        }
    

        /// 2. \todo Check the consistency of connectivity relationships between the elements



        /// 3. \todo Check the geometrical consistency of the topological relationships
        


        /// 4. Finally fill the nb_facets_in_surfaces_ vector
        model_.nb_facets_in_surfaces_.resize( model_.nb_surfaces()+1, 0 ) ;
        index_t count = 0 ;
        for( index_t i = 1; i < model_.nb_facets_in_surfaces_.size(); ++i ) {
            count += model_.surface( i-1 ).nb_cells() ;
            model_.nb_facets_in_surfaces_[i] = count ;
        }

#ifdef GRGMESH_DEBUG
        std::cout << "Model " << model_.name() <<" has " << std::endl 
            << std::setw(10) << std::left << model_.nb_vertices()   << " vertices "   << std::endl 
            << std::setw(10) << std::left << model_.nb_facets()   << " facets "   << std::endl  
            << std::setw(10) << std::left << model_.nb_regions()  << " regions "  << std::endl
            << std::setw(10) << std::left << model_.nb_surfaces() << " surfaces " << std::endl
            << std::setw(10) << std::left << model_.nb_lines()    << " lines "    << std::endl 
            << std::setw(10) << std::left << model_.nb_corners()  << " corners "  << std::endl
            << std::endl ;
#endif
        return true ;
    }


    /*!
     * @brief Load and build a BoundaryModel from a Gocad .ml file
     * 
     *  @details This is pretty tricky because of the annoying not well adapted file format. 
     * The correspondance between Gocad::Model3D elements and BoundaryModel elements is :
     * Gocad TSurf  <-> BoundaryModel Interface
     * Gocad TFace  <-> BoundaryModel Surface
     * Gocad Region <-> BoundaryModel Region
     * Gocad Layer  <-> BoundaryModel Layer
     *
     * @param[in] in Input .ml file stream
     */
    void BoundaryModelBuilderGocad::load_ml_file( std::istream& in )
    {
        // Clear the model_ 
        //  Not sure this is actually useful
        model_.clear() ;

        time_t start_load, end_load ;
        time( &start_load ) ;        

        // Count the number of TSurf - Interface
        index_t nb_tsurf = 0 ;
        // Count the number of TFace - Surface
        index_t nb_tface = 0 ;

        // Counters identifying the currently read TSurf or TFace
        index_t tsurf_count = 0 ;
        index_t tface_count = 0 ;

        index_t current_nb_tfaces = 0 ;
        index_t nb_tface_in_prev_tsurf = 0 ;

        // The file contains 2 parts and is read in 2 steps
        // 1. Read model info (true)
        // 2. Read surfaces geometry (false)
        bool read_model = true ;

        // The orientation of positive Z
        // can change for each TSurf and need to be read
        int z_sign = 1 ;

        // In the .ml file - vertices are indexed TSurf by Tsurf
        // They can be duplicated inside one TSurf and betweeen TSurfs
        
        // Indices of the vertices of the currently built TSurf in the model
        std::vector< index_t > tsurf_vertex_ptr ;        
        // Where the vertices of a TFace start in the vertices of the TSurf (offest)
        std::vector< index_t > tface_vertex_start ;

        // Indices of vertices in facets (triangles) of the currently built TFace
        std::vector< index_t > tface_facets ;
        // Starting and ending indices of each facet triangle in the tface_facets vector
        std::vector< index_t > tface_facets_ptr ;
        tface_facets_ptr.push_back( 0 ) ;
        
        // Intermediate information for contact parts building
        std::vector< Border > borders_to_build ;

        // Surfaces for which the KeyFacet orientation should be changed
        // because it does not match triangle orientations.
        std::vector< index_t > change_key_facet ;

        // Begin reading the file 
        InputStream lis( in ) ;

        while( !lis.eof() ) {
            lis.get_line() ;
            std::string keyword ;
            lis >> keyword ;

            if( read_model ) {
                if( std::string( keyword, 0, 5 ) == "name:" ) {                    
                    set_model_name( std::string( keyword, 5 ) );
                }
                else if( keyword == "TSURF" ) {
                    /// 1. Read the TSurf information and create 
                    /// the corresponding Interface from its name
                    std::string temp_str ;
                    std::stringstream tsurf_name ;
                    lis >> temp_str ;
                    tsurf_name << temp_str ;
                    while( !lis.eol() ) {
                        lis >> temp_str ;
                        tsurf_name << "_" << temp_str ;
                    }
                    create_interface( tsurf_name.str() ) ;
                    nb_tsurf++ ;
                } else if( keyword == "TFACE" ) {
                    /// 2. Read the TFace info and build the 
                    /// corresponding Surface from its parent Interface, its type, and 
                    /// its key facet - from which + and - side are determined
                    index_t id ;
                    lis >> id ;
                    std::string type ;
                    lis >> type ;

                    std::string temp_str ;
                    std::stringstream tsurf_name ;
                    lis >> temp_str ;
                    tsurf_name << temp_str ;
                    while( !lis.eol() ) {
                        lis >> temp_str ;
                        tsurf_name << "_" << temp_str ;
                    }
                    // Get the key facet that give the orientation of the surface part
                    // Triangles in Gocad clockwise
                    vec3 p0, p1, p2 ;
                    lis.get_line() ;
                    lis >> p0 ;
                    lis.get_line() ;
                    lis >> p1 ;
                    lis.get_line() ;
                    lis >> p2 ;

                    create_surface( tsurf_name.str(), type,
                        Surface::KeyFacet( p0, p1, p2 ) ) ;
                    nb_tface++ ;

                } else if( keyword == "REGION" ) {
                    /// 3. Read Region information and create them from their name,
                    /// the surfaces on their boundary                    
                    index_t id ;
                    std::string name ;
                    lis >> id >> name ;

                    std::vector< std::pair< index_t, bool > > region_boundaries ;
                    bool end_region = false ;

                    while( !end_region ) {
                        lis.get_line() ;
                        for( index_t i = 0; i < 5; ++i ) {
                            signed_index_t tface_id ;
                            lis >> tface_id ;
                            if( tface_id == 0 ) {
                                end_region = true ;
                                break ;
                            } else {
                                // Correction because ids begin at 1 in the file
                                tface_id = tface_id > 0 ? tface_id-1 : -tface_id-1 ;                                            
                                region_boundaries.push_back(
                                    std::pair< index_t, bool >( tface_id, false ) ) ;
                            }
                        }
                    }
                    // The Universe is not a regular region, it is the regions
                    // outside the volume of interest 
                    if( name != "Universe" ) {
                        create_region( name, region_boundaries ) ;
                    }
                    else {
                        set_universe( region_boundaries ) ;
                    }
                } else if( keyword == "LAYER" ) {
                    /// 4. Build the volumetric layers from their name and 
                    /// the regions in them
                    std::string name ;
                    lis >> name ;
                    index_t layer_id = create_layer( name ) ;
                    bool end_layer = false ;
                    while( !end_layer ) {
                        lis.get_line() ;
                        for( index_t i = 0; i < 5; ++i ) {
                            index_t region_id ;
                            lis >> region_id ;
                            if( region_id == 0 ) {
                                end_layer = true ;
                                break ;
                            } else {
                                region_id -= nb_tface+1 ; // Remove Universe region
                                // Correction because ids begin at 1 in the file
                                add_child( BoundaryModelElement::BM_LAYER, layer_id, region_id-1 ) ;
                            }
                        }
                    }
                } else if( keyword == "END" ) {
                    // End of the high level information on the model
                    // Switch to reading the geometry of the model surfaces
                    read_model = false ;
                    continue ;
                }
            } else {
                if( keyword == "GOCAD" ) {
                    // This is the beginning of a new TSurf = Interface
                    tsurf_count++ ;
                }
                if( keyword == "ZPOSITIVE" ) {
                    std::string positive ;
                    lis >> positive ;
                    if( positive == "Elevation" ) z_sign = 1 ;
                    else if( positive == "Depth" ) z_sign = -1 ;
                } else if( keyword == "END" ) {
                    // This the END of a TSurf
                    if( tsurf_count > 0 ) {
                        // End the last TFace - Surface of this TSurf
                        set_surface_geometry( 
                            tface_count-1,
                            std::vector< index_t >(
                                tsurf_vertex_ptr.begin() + tface_vertex_start.back(),
                                tsurf_vertex_ptr.end() ),
                            tface_facets,
                            tface_facets_ptr ) ;

                        if( !check_key_facet_orientation( tface_count-1 ) ){
                            change_key_facet.push_back( tface_count-1 ) ;
                        }

                        tface_facets.clear() ;
                        tface_facets_ptr.clear() ;
                        tface_facets_ptr.push_back( 0 ) ;

                        // End this TSurf - Interface
                        nb_tface_in_prev_tsurf += tface_vertex_start.size() ;
                        tsurf_vertex_ptr.clear() ;
                        tface_vertex_start.clear() ;
                    }
                } else if( keyword == "TFACE" ) {
                    // Beginning of a new TFace - Surface
                    if( tface_vertex_start.size() > 0 ) {
                        // End the previous TFace - Surface

                        set_surface_geometry( 
                            tface_count-1,std::vector< index_t >(
                                tsurf_vertex_ptr.begin() + tface_vertex_start.back(),
                                tsurf_vertex_ptr.end() ),
                            tface_facets,
                            tface_facets_ptr ) ;

                        if( !check_key_facet_orientation( tface_count-1 ) ) {
                            change_key_facet.push_back( tface_count-1 ) ;
                        }
                        
                        tface_facets.clear() ;
                        tface_facets_ptr.clear() ;
                        tface_facets_ptr.push_back( 0 ) ;
                    }
                    // Register where begin the new TFace vertices
                    tface_vertex_start.push_back( tsurf_vertex_ptr.size() ) ;

                    tface_count++ ;
                }
                /// 4. Read the surface vertices and facets (only triangles in Gocad Model3d files)
                else if( keyword == "VRTX" || keyword == "PVRTX" ) {
                    index_t id ;
                    vec3 p ;
                    lis >> id >> p ;
                    p.z *= z_sign ;
                    tsurf_vertex_ptr.push_back( add_vertex( p ) ) ;
                } else if( keyword == "PATOM" || keyword == "ATOM" ) {
                    index_t id ;
                    index_t v_id ;
                    lis >> id >> v_id ;
                    tsurf_vertex_ptr.push_back( tsurf_vertex_ptr[v_id - 1] ) ;
                } else if( keyword == "TRGL" ) {
                    // Ids of the vertices of each triangle in the TSurf
                    index_t p1, p2, p3 ;
                    lis >> p1 >> p2 >> p3 ;
                    // Change to ids in the TFace
                    p1 += -tface_vertex_start.back()-1 ;
                    p2 += -tface_vertex_start.back()-1 ;
                    p3 += -tface_vertex_start.back()-1 ;

                    tface_facets.push_back( p1 ) ;
                    tface_facets.push_back( p2 ) ;
                    tface_facets.push_back( p3 ) ;
                    tface_facets_ptr.push_back( tface_facets.size() ) ;
                }
                /// 5. Build the corners from their position and the surface parts
                ///    containing them
                else if( keyword == "BSTONE" ) {
                    index_t v_id ;
                    lis >> v_id ;
                    // correction to start at 0
                    v_id-- ;

                    // Get the TFace
                    index_t part_id = tface_vertex_start.size() - 1 ;
                    for( index_t i = 0; i < tface_vertex_start.size(); ++i ) {
                        if( v_id < tface_vertex_start[i] ) {
                            part_id = i - 1 ;
                            break ;
                        }
                    }
                    part_id += nb_tface_in_prev_tsurf ;

                    // c'est plus bon ça -compare geometry
                    index_t new_c = find_or_create_corner( tsurf_vertex_ptr[v_id] ) ;
                }
                /// 6. Read the Border information and store it
                else if( keyword == "BORDER" ) {
                    index_t id, p1, p2 ;
                    lis >> id >> p1 >> p2 ;
                    p1-- ;
                    p2-- ;

                    // Get the global corner id
                    index_t corner_id = find_corner( model_.vertex( tsurf_vertex_ptr[p1] ) ) ;
                    grgmesh_assert( corner_id != NO_ID ) ;

                    // Get the surface
                    index_t part_id = NO_ID ;
                    for( index_t i = 0; i < tface_vertex_start.size(); ++i ) {
                        if( p1 < tface_vertex_start[i] ) {
                            grgmesh_assert( p2 < tface_vertex_start[i] ) ;

                            // Get vertices ids in the surface
                            p1 += -tface_vertex_start[i - 1] ;
                            p2 += -tface_vertex_start[i - 1] ;

                            // i-1 is the id of the TFace in this TSurf
                            part_id = i - 1 ;
                            break ;
                        }
                    }
                    if( part_id == NO_ID ) {
                        // It is in the last built Tface
                        p1 += -tface_vertex_start[tface_vertex_start.size() - 1] ;
                        p2 += -tface_vertex_start[tface_vertex_start.size() - 1] ;
                        part_id = tface_vertex_start.size() - 1 ;
                    }
                    // The number of tfaces in previous tsurf is also to add
                    part_id += nb_tface_in_prev_tsurf ;

                    borders_to_build.push_back(
                        Border( part_id, corner_id, p1, p2 ) ) ;
                }
            }
        }

        make_vertices_unique() ;

        /// 7. Build the Lines
        build_lines( borders_to_build ) ;

        /// 8. Build the Contacts
        build_contacts() ;

    
        /// Eventual changes of the key facet orientation
        end_surfaces( change_key_facet ) ;
     
        // Finish up the model - CRASH if this failed
        grgmesh_assert( end_model() ) ;
        
        time( &end_load ) ;
#ifdef GRGMESH_DEBUG
        std::cout << "Info" << " Boundary model loading time"
            << difftime( end_load, start_load ) << " sec" << std::endl ;
#endif
    }
    

    /*!
     * @brief Find the facet which first 3 vertices are given 
     * @param[in] surface_id Index of the surface
     * @param[in] p0 First point coordinates
     * @param[in] p1 Second point coordinates
     * @param[in] p2 Third point coordinates
     * @param[out] same_sign Is true if the found facet has the same orientation than triangle p0p1p2
     * @return Index of the found facet, NO_ID if none found
     */
    index_t BoundaryModelBuilderGocad::find_key_facet( 
        index_t surface_id, const vec3& p0, const vec3& p1,
        const vec3& p2, bool& same_sign ) const 
    {
        const Surface& surface = model_.surface( surface_id ) ;
        same_sign = false ;
    
        for( index_t t = 0; t < surface.nb_cells(); ++t ){            
            const vec3& pp0 = model_.vertex( surface.model_vertex_id( t, 0 ) ) ;
            const vec3& pp1 = model_.vertex( surface.model_vertex_id( t, 1 ) ) ;
            const vec3& pp2 = model_.vertex( surface.model_vertex_id( t, 2 ) ) ;
            
            if( p0 == pp0 ) {
                if( p1 == pp1 && p2 == pp2 ) {
                    same_sign = true ;
                    return t ;
                }
                if( p1 == pp2 && p2 == pp1 ) {
                    same_sign = false ;
                    return t ;
                }
            }
            if( p0 == pp1 ) {
                if( p1 == pp0 && p2 == pp2 ) {
                    same_sign = false ;
                    return t ;
                }
                if( p1 == pp2 && p2 == pp0 ) {
                    same_sign = true ;
                    return t ;
                }
            }
            if( p0 == pp2 ) {
                if( p1 == pp0 && p2 == pp1 ) {
                    same_sign = true ;
                    return t ;
                }
                if( p1 == pp1 && p2 == pp0 ) {
                    same_sign = false ;
                    return t ;
                }
            }
        }
        return NO_ID ;
    }


    /*!
     * @brief Verify that a surface key facet has an orientation consistent with the surface facets.
     * @param[in] surface_id Index of the surface
     * @return False if the key_facet orientation is not the same than the surface facets, else true.
     */
    bool BoundaryModelBuilderGocad::check_key_facet_orientation( index_t surface_id ) 
    {
        const Surface& S = model_.surface( surface_id ) ;
        const Surface::KeyFacet& key_facet = S.key_facet() ;

        if( key_facet.is_default() ) {
            set_surface_first_triangle_as_key( surface_id ) ;
            return true ;
        }
        else {
            const vec3& p0 = key_facet.p0_ ;
            const vec3& p1 = key_facet.p1_ ;
            const vec3& p2 = key_facet.p2_ ;
            bool same_sign = false ;

            index_t t = find_key_facet( surface_id, p0, p1, p2, same_sign ) ;
            if( t == NO_ID ) {
                vec3 p00 = -1*p0 ;
                vec3 p10 = -1*p1 ;
                vec3 p20 = -1*p2 ;
                // It is because of the sign of Z that is not the same 
                t = find_key_facet( surface_id, p00, p10, p20, same_sign ) ;
            }
            grgmesh_assert( t != NO_ID ) ;
            return same_sign ;
        }
    }
       

    /*!
     * @brief Get the points of a Line between two corners on a Surface
     *   
     * WE ASSUME THAT THE STORAGE OF THE POINTS IS UNIQUE IN THE MODEL AND THAT 
     * SURFACES DO SHARE POINTS ON THEIR CONTACT LINES
     * make_vertices_unique() must have been called first
     *
     * @param[in] S Index of the surface
     * @param[in] id0 Index of the starting point( a corner ) in S
     * @param[in] id1 Index of the second point on the Line in S
     * @param[out] border_vertex_model_ids Indices of vertices on the Line (resized at 0 at the beginning)
     * @return Index of the Corner at which the Line ends
     */
    index_t BoundaryModelBuilderGocad::determine_line_vertices( 
        const Surface& S, 
        index_t id0, 
        index_t id1,
        std::vector< index_t >& border_vertex_model_ids 
        ) const 
    {
        grgmesh_debug_assert( id0 < S.nb_vertices() && id1 < S.nb_vertices() ) ;

        border_vertex_model_ids.resize( 0 ) ;
                 
        // Starting facet that contains the two given vertices
        index_t f = S.facet_from_surface_vertex_ids( id0, id1 ) ;
        grgmesh_assert( f != Surface::NO_ID ) ;

        // Global ids at the model level 
        index_t p0 = S.model_vertex_id( id0 ) ;
        index_t p1 = S.model_vertex_id( id1 ) ;

        border_vertex_model_ids.push_back( p0 ) ;
        border_vertex_model_ids.push_back( p1 ) ;
            
        index_t p1_corner = BoundaryModelBuilder::find_corner( p1 ) ;
        while( p1_corner == NO_ID ) {

            index_t next_f = NO_ID ;
            index_t id1_in_next = NO_ID ;
            index_t next_id1_in_next = NO_ID ;

            // We want the next triangle that is on the boundary and share p1
            // If there is no such triangle, the third vertex of the current triangle is to add
            S.next_on_border( f, S.facet_vertex_id(f, id0), S.facet_vertex_id(f,id1), 
                next_f, id1_in_next, next_id1_in_next ) ;

            grgmesh_assert(
                next_f != NO_ID && id1_in_next != NO_ID
                    && next_id1_in_next != NO_ID ) ;
            
            index_t next_id1 =  S.surf_vertex_id( next_f, next_id1_in_next ) ;

            // Update
            f = next_f ;
            id0 = id1 ;
            id1 = next_id1 ;

            p1 = S.model_vertex_id( next_id1 ) ;
            border_vertex_model_ids.push_back( p1 ) ;         
            p1_corner = BoundaryModelBuilder::find_corner( p1 ) ;
        }
        return p1_corner ; 
    }
    
    /*!
     * @brief Creates all Lines for the model
     *
     * @param[in] borders Information on Surface boundaries gathered at .ml file reading
     */
    void BoundaryModelBuilderGocad::build_lines( const std::vector< Border >& borders )
    {      
        std::vector< index_t > global_ids ;

        for( index_t i = 0; i < borders.size(); ++i ) {
            const Border& b = borders[i] ;

            /// 1- Build the boundary : construct the vector
            /// of vertices on the border
            const Surface& S = model_.surface( b.part_id_) ;

            index_t end_corner_id = determine_line_vertices( 
                S, b.p0_, b.p1_, global_ids ) ;           

            /// 2 - Check if this border already exists
            index_t line_id = find_or_create_line( global_ids ) ;

            // Add the surface in which this line is
            add_element_in_boundary( BoundaryModelElement::BM_LINE, line_id, b.part_id_ ) ;
        }
    }

    /*!
     * @brief Build the Contacts
     * @details One contact is a group of lines shared by the same Interfaces     
     */
    void BoundaryModelBuilderGocad::build_contacts()
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            // The surface part in whose boundary is the part
            std::set< index_t > interfaces ;
            //std::vector< GEOL_FEATURE > types ;
            for( index_t j = 0; j < model_.line(i).nb_in_boundary(); ++j ) {
                index_t sp_id = model_.line(i).in_boundary_id( j ) ;
                const BoundaryModelElement& p = model_.surface(sp_id).parent() ;
                interfaces.insert( p.id() ) ;
                //types.push_back( p.geological_feature() ) ;
            }
            std::vector< index_t > toto( interfaces.begin(), interfaces.end() ) ;
            index_t contact_id = find_or_create_contact( toto ) ;
            add_child( BoundaryModelElement::BM_CONTACT, contact_id, i ) ;
        }
    }

    /*!
     * @brief Finish up Surface 
     * @details Calls the end_surfaces() function and switch KeyFacet orientation of the required
     * surfaces.
     *
     * @param[in] change_orientation Indices of the surfaces in which KeyFacet orientation is not
     * consistent with the its facet orientation
     */
    void BoundaryModelBuilderGocad::end_surfaces(
        const std::vector< index_t >& change_orientation )
    {
        //end_surfaces() ;

        for( index_t i = 0; i < change_orientation.size(); i++ ) {
            index_t s_i = change_orientation[i] ; 
            
            // Change the key facet            
            set_surface_first_triangle_as_key( s_i ) ; 
            const Surface& S = model_.surface( s_i ) ;

            // Change the sign of this Surface in all the regions containing it
            for( index_t j = 0; j < S.nb_in_boundary(); ++j ) {                
                BoundaryModelElement& R =
                    const_cast< BoundaryModelElement& >( S.in_boundary(j) ) ;              

                for( index_t b = 0; b < R.nb_boundaries(); ++b ){
                    if( R.boundary_id(b) == s_i ){
                        bool old_side = R.side( b ) ;
                        R.set_boundary( b, R.boundary_id(b), !old_side ) ;
                    }
                }
            }
        }
    }
    

    /*!
     * @brief Add a Surface to the model
     *
     * @param[in] interface_name Name of the parent. The parent MUST exist.
     * @param[in] type Type of the Surface
     * @param[in] key KeyFacet for this Surface
     */
    void BoundaryModelBuilderGocad::create_surface(
        const std::string& interface_name,
        const std::string& type,
        const Surface::KeyFacet& key )
    {
        index_t parent = find_interface( interface_name ) ;
        if( interface_name != "" ) grgmesh_assert( parent != NO_ID ) ;

        index_t id = create_element( BoundaryModelElement::BM_SURFACE ) ;
        set_parent( BoundaryModelElement::BM_SURFACE, id, parent ) ;
        //set_element_geol_feature( BoundaryModelElement::BM_SURFACE, id, t ) ;
        set_surface_key_facet( id, key ) ;
    }


    /*!
     * @brief Get the index of the Corner at given coordinates
     * @param[in] p Coordinates of the vertex
     * @return NO_ID or the index of the Corner
     */
    index_t BoundaryModelBuilderGocad::find_corner( const vec3& p ) const
    {
        for( index_t i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).vertex() == p ) return i ;
        }
        return NO_ID ;
    }

} 
