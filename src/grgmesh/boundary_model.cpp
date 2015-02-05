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

        contacts_.clear() ;
        interfaces_.clear() ;
        layers_.clear() ;
    }

    /*!
     * @brief Total number of facets in the model Surface
     */
    index_t BoundaryModel::nb_facets() const 
    {
        index_t result = 0 ; 
        for( index_t i = 0; i < nb_surfaces(); ++i ){
            result += surface(i).nb_cells() ;
        }
        return result ;
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
        
        /// 2. Check that the Universe region exists 
        /// \todo Write some code to create the universe (cf. line 805 to 834 de s2_b_model.cpp)
        if( universe_.name() != "Universe" ) {
            GEO::Logger::err( "" )
                << "The region universe is not defined for the model. IMPLEMENTATION TO DO"
                << std::endl ;
            return false ;
        }

        /// 3. Check that each region has a name and valid surfaces
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

        /// 4. Check that all the surfaces_ of the model are triangulated
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
            out << BME::geol_name( s.geological_feature() ) ;
            out << " " << s.parent().name() << std::endl ;
            
            // Print the key facet points, whuich are simply the first three
            // vertices of the first facet
            out << "  " << s.vertex( 0, 0 ) << std::endl ;
            out << "  " << s.vertex( 0, 1 ) << std::endl ;
            out << "  " << s.vertex( 0, 2 ) << std::endl ;

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
            out << BME::geol_name( tsurf.geological_feature() ) ;
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
     * @brief DEBUG function - Save the surfaces of the model with their facet attributes into an .eobj file.         
     * @details WARNING It is supposed that these attributes are integer and that all Surface
     * have the same attributes
     * 
     * @param[in] file_name Name of the file
     *
     * \todo Write generic I/O for attributes on a BoundaryModel
     * \todo Make this function const
     *
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
            // OK NOT NICE WE SUPPOSE ALL SURFACE have the same attributes 
            surface(0).facet_attribute_manager()->list_named_attributes(names) ;
                    
            for( index_t i=0; i < names.size(); i++ ) { 
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
                    
                    // Get the Stores for this surface
                    std::vector< AttributeStore* > stores( names.size() ) ;
                    for( index_t i = 0 ; i < names.size(); ++i ) {
                        stores[i] = S.facet_attribute_manager()->resolve_named_attribute_store( names[i] ) ;
                    }
                    // Output attributes values
                    for( index_t f = 0; f < S.nb_cells(); f++ ){            
                        out << "# attrs f " << count+1 ;
                        for( index_t j = 0; j < names.size(); j++ ) {
                            out << " " << *reinterpret_cast<index_t*>( stores[j]->data(f) ) ;
                        }
                        out << std::endl ;
                        count++ ;
                    }
                }
            }
        }

    }
    /*!
     * @brief Debug: Save a Surface of the model in the file OBJ format is used
     */
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

     /*! 
      * @brief Write in the out stream things to save for CONTACT, INTERFACE and LAYERS
      */
     void save_high_level_bme( std::ofstream& out, const BoundaryModelElement& E ) {
         /// First line:  TYPE - ID - NAME - GEOL
         out << BoundaryModelElement::type_name( E.element_type() ) << " " 
             << E.id() << " " ;
         if( E.has_name() ) out << E.name() << " ";
         else out << "no_name " ;
         out <<  BoundaryModelElement::geol_name( E.geological_feature() ) 
             << std::endl ;

         /// Second line:  IDS of children
         for( index_t j = 0; j < E.nb_children(); ++j ) {
             out << " " << E.child_id(j) ;
         }
         out << std::endl ;    
     }


     void BoundaryModel::save_bm_file( const std::string& file_name ) const 
     {
        std::ofstream out ;
        out.open( file_name.c_str() );
        if( out.bad() ){
            std::cout << "Error when opening the file: " << file_name.c_str() <<std::endl ;
            return ;
        }
        out.precision( 16 ) ;
      
        out << "GRGMESH BOUNDARY MODEL"<< std::endl ;
        out << "NAME " << name() << std::endl ; 

        // Number of the different type of elements 
        for( index_t i = BME::CORNER; i < BME::NO_TYPE; i++ ) {
             BME::TYPE type = (BME::TYPE) i;             
             out <<  "NB_"<< BME::type_name( type ) << " " << nb_elements( type ) << std::endl ;
        }
        // Write high-level elements
        for( index_t i = BME::CONTACT; i < BME::NO_TYPE; i++ ) {
             BME::TYPE type = (BME::TYPE) i ;
             index_t nb = nb_elements( type ) ;
             for( index_t j = 0; j < nb; ++j ) {
                save_high_level_bme( out, element( type, j ) ) ;
             }
        }        
        // Regions
        for( index_t i =0 ; i < nb_regions(); ++i ){
            const BoundaryModelElement& E = region(i) ;
            // Save ID - NAME - 
            out << BME::type_name( BME::REGION ) << " " << E.id() << " " ;
            if( E.has_name() ) out << E.name() ;
            else out << "no_name" ;
            out << std::endl ;

            // Second line Signed ids of boundary surfaces 
            for( index_t j = 0 ; j < E.nb_boundaries(); ++j ){
                if( E.side(j) ) out << "+" ;
                else out << "-" ;
                out << E.boundary_id(j) << " " ;
            }
            out << std::endl ;
        }
        // Universe
        out << "UNIVERSE "<< std::endl ;
        for( index_t j = 0 ; j < universe().nb_boundaries(); ++j ){
            if( universe().side(j) ) out << "+" ;
            else out << "-" ;
            out << universe().boundary_id(j) << " " ;
        }
        out << std::endl ;

        // Vertices
        out << "MODEL_VERTICES" << " " << nb_vertices() << std::endl ;                   
        for( index_t i = 0; i < nb_vertices(); ++i ) {
            out << vertex(i) 
                << std::endl ;
        }        
        // Corners
        for( index_t i = 0; i < nb_corners(); ++i ) {
            out << BME::type_name( BME::CORNER ) << " " 
                << corner(i).id() << " " << corner(i).model_vertex_id() << std::endl ;
        }
        // Lines
        for( index_t i = 0; i < nb_lines(); ++i ) {
            const Line& L = line(i) ;
            out << BME::type_name( BME::LINE) << " " << L.id() << std::endl ;
            out << "LINE_VERTICES " << L.nb_vertices() << std::endl ;
            int count = 0 ;
            for( index_t j = 0; j < L.nb_vertices(); ++j ) {
                out << L.model_vertex_id( j ) << "  " ;
                count++ ;
                if( count == 20 && j+1 < L.nb_vertices() ){
                    count = 0 ;
                    out << std::endl ;
                }
            }
            out << std::endl ;
            out << "IN_BOUNDARY " ;
            for( index_t j = 0; j < L.nb_in_boundary(); ++j ){
                out << L.in_boundary_id( j ) << " " ;
            }
            out << std::endl ;
        }
        // Surfaces
        for( index_t i = 0; i < nb_surfaces(); ++i ) {
            const Surface& S = surface(i) ;
            out << BME::type_name( BME::SURFACE ) << " " << S.id() << std::endl ;
            out << "SURFACE_CORNERS "<< S.nb_corners() << std::endl ;
            out << "SURFACE_FACETS " << S.nb_cells() << std::endl ;
            for( index_t j =0; j < S.nb_cells(); ++j ) {
                out << S.nb_vertices_in_facet(j) << " " ;
                for( index_t v= 0; v < S.nb_vertices_in_facet(j); ++v ){
                    out << S.model_vertex_id( j, v ) << " ";
                }
                out << std::endl ;
            }
        }
    }



} // namespace
