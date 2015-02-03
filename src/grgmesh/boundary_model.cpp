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
#include <grgmesh/boundary_model_builder.h>
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
            case BoundaryModelElement::NO_GEOL:
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
                builder.add_child( BoundaryModelElement::INTERFACE, id, i ) ;
            }

            // Set links from surfaces_ toward interfaces_
            for( index_t i = 0; i < interfaces_.size(); ++i ) {
                builder.set_parent( BoundaryModelElement::SURFACE, interfaces_[i].child( 0 ).id(), i ) ;
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
                builder.set_element_name( BoundaryModelElement::REGION, i, name.str() ) ;
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
     * \todo Make this function const
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

     void BoundaryModel::save_surface_as_obj_file( index_t s, const std::string& file_name ) const {
        std::ofstream out ;
        out.open( file_name.c_str() );
        if( out.bad() ){
            std::cout << "Error when opening the file: " << file_name.c_str() <<std::endl ;
            return ;
        }
        out.precision( 16 ) ;
        const Surface& S = surface(s) ;        
        for( index_t p = 0; p < S.nb_vertices(); p++ ){
            const vec3& V = S.vertex(p) ;
            out << "v" 
                << " " << V.x 
                << " " << V.y 
                << " " << V.z 
                << std::endl ;
        }
        for( index_t f = 0; f < S.nb_cells(); f++ ){
            out << "f" << " " ;
            for( index_t v = 0; v < S.nb_vertices_in_facet(f); v++ ){                    
                out << S.surf_vertex_id( f, v )+1 << " " ;            
            }
            out << std::endl ;
        }                  
     }




} // namespace
