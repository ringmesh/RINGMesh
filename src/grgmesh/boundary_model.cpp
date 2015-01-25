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
     */
    void BoundaryModel::save_type( std::ostream& out, GEOL_FEATURE t )
    {
        switch( t ) {
            case STRATI:
                out << "top" ;
                break ;
            case FAULT:
                out << "fault" ;
                break ;
            case VOI:
                out << "boundary" ;
                break ;
            case ALL:
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
     * \brief Returns a reference the identified BoundaryModelElement     
     *
     * @param[in] dim Type of the element 
     * @param[in] index Index of the element
     * @return A reference to the BoundaryModelElement .
     *
     */
    const BoundaryModelElement& BoundaryModel::element( BM_TYPE t, index_t index ) const {
        switch( t ){
            case BM_CORNER    : return corner       ( index ) ;
            case BM_LINE      : return line         ( index ) ;
            case BM_SURFACE   : return surface      ( index ) ;
            case BM_REGION    : return region       ( index ) ;
            case BM_CONTACT   : return contact      ( index ) ;
            case BM_INTERFACE : return one_interface( index ) ;
            case BM_LAYER     : return layer        ( index ) ;
            default:
                grgmesh_assert_not_reached ;
                return dummy_element ;
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
            if( model_facet_id > nb_facets_in_surfaces_[i-1] && 
                model_facet_id < nb_facets_in_surfaces_[i] )
            {
                s = i-1 ;
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

        BoundaryModelBuilder builder( *this ) ;
        builder.load_file( input ) ;
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
                builder.add_interface_child( id, i ) ;
            }

            // Set links from surfaces_ toward interfaces_
            for( index_t i = 0; i < interfaces_.size(); ++i ) {
                builder.set_parent( surfaces_[interfaces_[i].child( 0 ).id()],  i ) ;
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
                builder.set_name( region, name.str() ) ;
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
     * @brief Map a read geological type with the Geological feature enum
     *
     * @param[in] in Name of the feature
     * @return The geological feature index
     *
     * \todo Keep all the information ( add new GEOL_FEATURE) instead of simplfying it.
     */
    GEOL_FEATURE BoundaryModelBuilder::determine_geological_type( const std::string& in )
    {
        if( in == "" ) return ALL ;
        if( in == "reverse_fault" ) return FAULT ;
        if( in == "normal_fault" ) return FAULT ;
        if( in == "fault" ) return FAULT ;
        if( in == "top" ) return STRATI ;
        if( in == "none" ) return STRATI ;
        if( in == "unconformity" ) return STRATI ;
        if( in == "boundary" ) return VOI ;

        std::cout<< "ERROR" << "Unexpected type in the model file " << in
            << std::endl ;
        return ALL ;
    }

    /*!
     * @brief Compute an intersection type
     *
     * @param[in] types Type that intersect
     * @return Intersection type
     */
    GEOL_FEATURE BoundaryModelBuilder::determine_type(
        const std::vector< GEOL_FEATURE >& types )
    {
        if( types.size() == 0 ) return ALL ;

        // Sort and remove duplicates form the in types
        std::vector< GEOL_FEATURE > in = types ;
        std::sort( in.begin(), in.end() ) ;
        index_t new_size = std::unique( in.begin(), in.end() ) - in.begin() ;
        in.resize( new_size ) ;

        if( in.size() == 1 ) return in[0] ;

        if( in.size() == 2 ) {
            if( in[0] == ALL ) return ALL ;
            if( in[0] == STRATI ) {
                if( in[1] == FAULT ) return STRATI_FAULT ;
                if( in[1] == VOI ) return STRATI_VOI ;
            } else if( in[0] == FAULT ) {
                if( in[1] == VOI ) return FAULT_VOI ;
            }
            // Other cases ? for corners ? what is the vertex ?
            return ALL ;
        }
        return ALL ;
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
        model_.lines_[id].vertices_ = vertices ;
    }

     /*!
     * @brief Add the vertices to the model and set them as the Line vertices
     *
     * @param[in] id Line index 
     * @param[in] vertices Coordinates of the vertices on the line
     * 
     * \todo Check when it is called ! Verify that it is before make_vertices_unique()
     */
    void BoundaryModelBuilder::set_line(
        index_t id,
        const std::vector< vec3 >& vertices )
    {
        grgmesh_assert( id < model_.nb_lines() ) ;

        for( index_t p = 0; p < vertices.size(); p++ ) {
            model_.lines_[id].vertices_.push_back( add_vertex( vertices[p] ) ) ;
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
     * @brief Creates a new empty Layer with the given name 
     *
     * @param[in] name Name of the layer
     * @return The layer index
     */
    index_t BoundaryModelBuilder::create_layer(
        const std::string& name )
    {
        index_t id = model_.layers_.size() ;
        model_.layers_.push_back( BoundaryModelElement( &model_, BM_LAYER, id ) ) ;
        model_.layers_[id].set_name( name ) ;
        return id ;
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
     * @brief Add a point to the model and set it as the corner vertex
     *
     * @param[in] corner_id Index of the corner
     * @param[in] vertex Coordinates of the vertex
     */
    void BoundaryModelBuilder::set_corner( index_t corner_id, const vec3& vertex )
    {
        grgmesh_debug_assert( corner_id < model_.nb_corners() ) ;
        model_.corners_[corner_id].set_vertex( add_vertex( vertex ) ) ;
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
        GEOL_FEATURE type )
    {
        index_t id = model_.nb_interfaces() ;
        model_.interfaces_.push_back( BoundaryModelElement( &model_, BM_INTERFACE, id, NO_ID, type ) ) ;
        model_.interfaces_[id].set_name( name ) ;
        return id ;
    }

    /*!
     * @brief Copy macro information from a model
     * @details Copy the all the model elements and their relationship ignoring their geometry
     * 
     * @param[in] from Model to copy the information from
     */
    void BoundaryModelBuilder::copy_macro_topology( const BoundaryModel* from )
    {
        model_.name_ = from->name_ ;
        model_.corners_.resize( from->nb_corners(), Corner( &model_ ) ) ;
        model_.lines_.resize( from->nb_lines(), Line( &model_ ) ) ;
        model_.surfaces_.resize( from->nb_surfaces(), Surface( &model_ ) ) ;
        model_.regions_.resize( from->nb_regions(), BoundaryModelElement( &model_, BM_REGION ) ) ;
        model_.layers_.resize( from->nb_layers(), BoundaryModelElement( &model_, BM_LAYER ) ) ;
        model_.contacts_.resize( from->nb_contacts(), BoundaryModelElement( &model_, BM_CONTACT ) ) ;
        model_.interfaces_.resize( from->nb_interfaces(), BoundaryModelElement( &model_, BM_INTERFACE ) ) ;
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_corners(); i++ ) {
            model_.corners_[i].copy_macro_topology( from->corner( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_lines(); i++ ) {
            model_.lines_[i].copy_macro_topology( from->line( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_surfaces(); i++ ) {
            model_.surfaces_[i].copy_macro_topology( from->surface( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_layers(); i++ ) {
            model_.layers_[i].copy_macro_topology( from->layer( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_regions(); i++ ) {
            model_.regions_[i].copy_macro_topology( from->region( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_contacts(); i++ ) {
            model_.contacts_[i].copy_macro_topology( from->contact( i ), model_ ) ;
        }
#pragma omp parallel for
        for( index_t i = 0; i < model_.nb_interfaces(); i++ ) {
            model_.interfaces_[i].copy_macro_topology( from->one_interface( i ), model_ ) ;
        }
        model_.universe_.copy_macro_topology( from->universe_, model_ ) ;
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
    void BoundaryModelBuilder::load_file( std::istream& in )
    {
        // Clear the model_ 
        //  Not sure this is actually useful
        model_.name_.clear() ;
        model_.vertices_.clear() ;
        model_.corners_.clear() ; 
        model_.lines_.clear() ; 
        model_.surfaces_.clear() ;
        model_.regions_.clear() ;
        model_.nb_facets_in_surfaces_.clear() ;
        model_.interfaces_.clear() ;
        model_.contacts_.clear() ;
        model_.layers_.clear() ;

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
                    model_.name_ = std::string( keyword, 5 ) ;
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
                                add_layer_child( layer_id, region_id-1 ) ;
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
                            tface_count-1,
                            std::vector< index_t >(
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

        // Then finish the job (the order matters)
        end_surfaces( change_key_facet ) ;
        end_interfaces() ;
        end_lines() ;
        end_contacts() ;
        end_corners() ;
        end_layers() ;
        end_model() ;
        
        time( &end_load ) ;
#ifdef GRGMESH_DEBUG
        std::cout << "Info" << " Boundary model loading time"
            << difftime( end_load, start_load ) << " sec" << std::endl ;
#endif
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
            model_.surfaces_[surface_id].adjacent_ = surface_adjacencies ;
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
    index_t BoundaryModelBuilder::find_key_facet( 
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
    bool BoundaryModelBuilder::check_key_facet_orientation( index_t surface_id ) const 
    {
        Surface& surface = model_.surfaces_[surface_id] ;
        Surface::KeyFacet& key_facet = surface.key_facet_ ;

        if( key_facet.is_default() ) {
            surface.set_first_triangle_as_key() ;
            return true ;
        }
        
        vec3& p0 = key_facet.p0_ ;
        vec3& p1 = key_facet.p1_ ;
        vec3& p2 = key_facet.p2_ ;
        bool same_sign = false ;

        index_t t = find_key_facet( surface_id, p0,p1,p2, same_sign ) ;
        if( t == NO_ID ) {
            // It is because of the sign of Z that is not the same 
            p0.z *= -1 ;
            p1.z *= -1 ;
            p2.z *= -1 ;
            t = find_key_facet( surface_id, p0, p1, p2, same_sign ) ; 
        }
        grgmesh_assert( t != NO_ID ) ;

        return same_sign ;
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
        std::vector< index_t >& adjacent = S.adjacent_ ;
        adjacent.resize( S.facets_.size(), Surface::NO_ADJACENT ) ;
       
        grgmesh_assert( S.facets_.size() > 0  ) ;

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
    }

    /*!
     * @brief Cut the Surface along the line 
     * @details First modify to NO_ADJACENT the neighbors the edges that are along the line 
     * and then duplicate the points along this new boundary
     * Corners are not duplicated - maybe they should be in some cases but not in general..
     * 
     * @param[in] surface_id Surface index
     * @param[in] line_id Line index
     */
    void BoundaryModelBuilder::cut_surface_by_line( index_t surface_id, index_t line_id ) {
        Surface& S = model_.surfaces_[surface_id] ;
        const Line& L = model_.line( line_id ) ;

        for( index_t i= 0; i+1 < L.nb_vertices(); ++i ) {
            index_t p0 = L.model_vertex_id( i ) ;
            index_t p1 =  (i == L.nb_vertices()-1) ? L.model_vertex_id(0) : L.model_vertex_id( i+1 ) ;

            index_t f = Surface::NO_ID ;
            index_t v = Surface::NO_ID ;
            S.edge_from_model_vertex_ids(p0, p1, f, v) ;
            grgmesh_debug_assert( f != Surface::NO_ID && v != Surface::NO_ID ) ;

            index_t f2 = S.adjacent( f, v ) ;
            index_t v2 = Surface::NO_ID ;
            grgmesh_debug_assert( f2 != Surface::NO_ADJACENT ) ;
            S.edge_from_model_vertex_ids( p0, p1, f2, v2 ) ;
            grgmesh_debug_assert( v2 != Surface::NO_ID ) ;

            // Virtual cut - set adjacencies to NO_ADJACENT
            S.set_adjacent( f, v, Surface::NO_ADJACENT ) ;
            S.set_adjacent( f2, v2, Surface::NO_ADJACENT ) ;
        }
        
        // Now travel on one side of the "faked" boundary and actually duplicate
        // the vertices in the surface      
        // Get started in the surface - find (again) one of the edge that contains 
        // the first two vertices of the line
        index_t f = Surface::NO_ID ;
        index_t v = Surface::NO_ID ;
        S.oriented_edge_from_model_vertex_ids( L.model_vertex_id( 0 ), L.model_vertex_id( 1 ), f, v ) ;
        grgmesh_assert( f != Surface::NO_ID && v != Surface::NO_ID ) ;

        index_t id0 = S.surf_vertex_id( f, v ) ;
        index_t id1 = S.surf_vertex_id( f, S.next_in_facet(f,v) ) ;
        
        // Stopping criterion
        index_t last_vertex = L.model_vertex_id( L.nb_vertices()-1 ) ;
        // Hopefully we have all the vertices on the Line.. 
        /// \todo Check that all vertices on the line are recovered
        while( S.model_vertex_id( id1 ) != last_vertex ) {
            // Get the next vertex on the border 
            // Same algorithm than in determine_line_vertices function
            index_t next_f = Surface::NO_ID ;
            index_t id1_in_next = Surface::NO_ID ;
            index_t next_id1_in_next = Surface::NO_ID ;

            // Get the next facet and next triangle on this boundary
            S.next_on_border( f, S.facet_vertex_id(f, id0), S.facet_vertex_id(f,id1), 
                next_f, id1_in_next, next_id1_in_next ) ;
            grgmesh_assert(
                next_f != Surface::NO_ID && id1_in_next != Surface::NO_ID
                    && next_id1_in_next != Surface::NO_ID ) ;
            
            index_t next_id1 = S.surf_vertex_id( next_f, next_id1_in_next ) ;
            
            // Duplicate the vertex at id1 
            // After having determined the next 1 we can probably get both at the same time 
            // but I am lazy, and we must be careful not to break next_on_border function (Jeanne)
            std::vector< index_t > facets_around_id1 ;
            S.facets_around_vertex( id1, facets_around_id1, false, f ) ;

            S.vertices_.push_back( S.model_vertex_id(id1) ) ;
            grgmesh_debug_assert( S.nb_vertices() > 0 ) ;
            index_t new_id1 = S.nb_vertices()-1 ;
            
            for( index_t i = 0; i < facets_around_id1.size(); ++i ){
                index_t cur_f = facets_around_id1[i] ;
                for( index_t cur_v = 0; cur_v < S.nb_vertices_in_facet( cur_f ) ; cur_v++ )
                {
                    if( S.surf_vertex_id( cur_f, cur_v ) == id1 ) {
                        S.facets_[ S.facet_begin( cur_f ) + cur_v ] = new_id1 ;
                        break ;
                    }
                }
            }
            // Update
            f = next_f ;
            id0 = new_id1 ;
            id1 = next_id1 ;
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
        MakeUnique unique(  model_.vertices_ ) ;
        unique.unique( 5 ) ;
        model_.vertices_.resize(0) ;
        unique.unique_points( model_.vertices_ ) ;
        const std::vector< index_t >& old2new = unique.indices() ;

        for( index_t s = 0; s < model_.nb_surfaces(); s++ ) {
            Surface& surface = model_.surfaces_[s] ;
            for( index_t p = 0; p < surface.nb_vertices(); p++ ) {
                surface.vertices_[p] = old2new[surface.vertices_[p]] ;
            }
        }
        for( index_t l = 0; l < model_.nb_lines(); l++ ) {
            Line& line = model_.lines_[l] ;
            for( index_t p = 0; p < line.nb_vertices(); p++ ) {
                line.vertices_[p] = old2new[line.vertices_[p]] ;
            }
        }
        for( index_t co = 0; co < model_.nb_corners(); co++ ) {
            Corner& corner = model_.corners_[co] ;
            corner.vertex_ = old2new[corner.vertex_] ;
        }       
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
        model_.universe_.set_element_type( BM_REGION ) ;
        model_.universe_.model_ = &model_ ;

        for( index_t i = 0; i < boundaries.size(); ++i ) {
            grgmesh_assert( boundaries[i].first < model_.nb_surfaces() ) ;
            model_.universe_.add_boundary( boundaries[i].first,
                boundaries[i].second ) ;
            // If this surface have no type, set it at VOI
            model_.surfaces_[boundaries[i].first].set_geological_feature( VOI ) ;
        }
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
        index_t id = model_.regions_.size() ;
        model_.regions_.push_back( BoundaryModelElement( &model_, BM_REGION, id ) ) ;
        model_.regions_[id].set_name( name ) ;

        for( index_t i = 0; i < boundaries.size(); ++i ) {
            grgmesh_assert( boundaries[i].first < model_.nb_surfaces() ) ;
            add_region_oriented_boundary( id, boundaries[i].first,
                boundaries[i].second ) ;
        }
        return id ;
    }

    /*
    * @brief Adds an empty region to the model
    *
    *  Used in Geomodeling to convert a surface to a model
    */
    index_t BoundaryModelBuilder::create_region() {
        index_t id = model_.regions_.size() ;
        model_.regions_.push_back( 
            BoundaryModelElement( &model_, BM_REGION, id ) ) ;
        return id ;
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
    index_t BoundaryModelBuilder::determine_line_vertices( 
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
            
        index_t p1_corner = find_corner( p1 ) ;
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
            p1_corner = find_corner( p1 ) ;
        }
        return p1_corner ; 
    }
    
    /*!
     * @brief Creates all Lines for the model
     *
     * @param[in] borders Information on Surface boundaries gathered at .ml file reading
     */
    void BoundaryModelBuilder::build_lines( const std::vector< Border >& borders )
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
            index_t line_id = find_or_create_line( b.corner_id_,
                end_corner_id, global_ids ) ;

            // Add the surface in which this line is
            add_line_in_boundary( line_id, b.part_id_ ) ;
        }
    }

    /*!
     * @brief Build the Contacts
     * @details One contact is a group of lines shared by the same Interfaces     
     */
    void BoundaryModelBuilder::build_contacts()
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            // The surface part in whose boundary is the part
            std::set< index_t > interfaces ;
            std::vector< GEOL_FEATURE > types ;
            for( index_t j = 0; j < model_.lines_[i].nb_in_boundary(); ++j ) {
                index_t sp_id = model_.lines_[i].in_boundary_id( j ) ;
                const BoundaryModelElement& p = model_.surfaces_[sp_id].parent() ;
                interfaces.insert( p.id() ) ;
                types.push_back( p.geological_feature() ) ;
            }
            std::vector< index_t > toto( interfaces.begin(), interfaces.end() ) ;
            index_t contact_id = find_or_create_contact( toto, determine_type( types ) ) ;
            add_contact_child( contact_id, i ) ;
        }
    }

    /*! 
     * @brief Set the Interfaces delimited by each Contact     
     * @details Warning the Interfaces must be finished first.
     */
    void BoundaryModelBuilder::end_contacts()
    {
        for( index_t i = 0; i < model_.nb_contacts(); ++i ) {
            std::set< const BoundaryModelElement* > corners ;
            for( index_t j = 0; j < model_.contacts_[i].nb_children(); ++j ) {
                const BoundaryModelElement& child = model_.contacts_[i].child( j ) ;
                corners.insert( &child.boundary( 0 ) ) ;
                corners.insert( &child.boundary( 1 ) ) ;
            }
            for( std::set< const BoundaryModelElement* >::iterator it( corners.begin() );
                it != corners.end(); ++it ) {
                add_contact_boundary( i, ( *it )->id() ) ;
            }
        }
    }

    /*! 
     * @brief Set the parent Contact for all the Lines
     */
    void BoundaryModelBuilder::end_lines()
    {
        for( index_t i = 0; i < model_.nb_contacts(); ++i ) {
            for( index_t j = 0; j < model_.contacts_[i].nb_children(); ++j ) {
                index_t child = model_.contacts_[i].child_id( j ) ;
                model_.lines_[child].set_parent( i ) ;
                model_.lines_[child].set_geological_feature(
                    model_.contacts_[i].geological_feature() ) ;
            }
        }
    }

    /*!
     * @brief Fills the children and the in_boundary information for Interfaces 
     *
     * \todo End interface function not correct, in_boundary not well set
     */
    void BoundaryModelBuilder::end_interfaces()
    {
        for( index_t i = 0; i < model_.nb_surfaces(); ++i ) {
            index_t parent = model_.surfaces_[i].parent_id() ;
            add_interface_child( parent, i ) ;
        }
        for( index_t i = 0; i < model_.nb_contacts(); ++i ) {
            for( index_t j = 0; j < model_.contacts_[i].nb_in_boundary(); ++j ) {
                index_t b = model_.contacts_[i].in_boundary_id( j ) ;
                add_interface_boundary( b, i ) ;
            }
        }       
    }

    /*!
     * @brief Finishes up layers
     * \todo What do we want to have in the boundaries of layers ? surfaces or interfaces ?
     */
    void BoundaryModelBuilder::end_layers()
    {
        for( index_t i = 0; i < model_.nb_layers(); ++i ) {
            BoundaryModelElement& layer = model_.layers_[i] ;
            std::set< index_t > boundaries ;

            for( index_t r = 0; r < layer.nb_children(); r++ ) {
                const BoundaryModelElement& region = layer.child( r ) ;
                model_.regions_[region.id()].parent_ = i ;
                for( index_t sp = 0; sp < region.nb_boundaries(); sp++ ) {
                    boundaries.insert( region.boundary_id( sp ) ) ;
                }
            }

            for( std::set< index_t >::const_iterator it( boundaries.begin() );
                it != boundaries.end(); it++ ) {
                layer.add_boundary( *it ) ;
            }
        }
    }

    /*!
     * @brief Sets boundary and in_boundary information on surfaces
     */
    void BoundaryModelBuilder::end_surfaces()
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            for( index_t j = 0; j < model_.lines_[i].nb_in_boundary();
                ++j ) {
                index_t s_id = model_.lines_[i].in_boundary_id( j ) ;
                add_surface_boundary( s_id, i ) ;
            }
        }
        for( index_t i = 0; i < model_.nb_regions(); ++i ) {
            for( index_t j = 0; j < model_.regions_[i].nb_boundaries(); ++j ) {
                index_t s_id = model_.regions_[i].boundary_id( j ) ;
                add_surface_in_boundary( s_id, i ) ;
            }
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
    void BoundaryModelBuilder::end_surfaces(
        const std::vector< index_t >& change_orientation )
    {
        end_surfaces() ;
        for( index_t i = 0; i < change_orientation.size(); i++ ) {
            index_t s_i = change_orientation[i] ; 
            
            // Change the key facet
            Surface& surface = model_.surfaces_[s_i] ;
            surface.set_first_triangle_as_key() ;

            // Change the sign of this Surface in all the regions containing it
            for( index_t j = 0; j < surface.nb_in_boundary(); ++j ) {                
                BoundaryModelElement& region = model_.regions_[ surface.in_boundary_id( j ) ] ;              
                for( index_t b = 0; b < region.nb_boundaries(); ++b ){
                    if( region.boundary(b).id() == s_i ){
                        region.sides_[i] = !region.sides_[i] ;
                    }
                }
            }
        }
    }


    /*!
     * @brief Last function to call when building a model
     * @details Output information on the model and initialise the
     * nb_facets_in_surfaces_ vector
     * 
     */
    void BoundaryModelBuilder::end_model()
    {
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
    }                                                                          


    /*!
     * @brief Set on the boundary of which lines are corners
     */
    void BoundaryModelBuilder::end_corners()
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            index_t c0_id = model_.lines_[i].boundary_id( 0 ) ;
            index_t c1_id = model_.lines_[i].boundary_id( 1 ) ;

            add_corner_in_boundary( c0_id, i ) ;
            if( c1_id != c0_id ) add_corner_in_boundary( c1_id, i ) ;
        }
    }

    /*!
     * @brief Find or create a contact between given Interfaces
     * 
     * @param[in] interfaces Indices of the intersecting interfaces
     * @param[in] type Type of the intersection
     * @return Index of the Contact
     */
    index_t BoundaryModelBuilder::find_or_create_contact(
        std::vector< index_t >& interfaces,
        GEOL_FEATURE type )
    {
        index_t result = find_contact( interfaces ) ;
        if( result == NO_ID ) {
            // Create a name for this contact
            std::string name = "contact_" ;
            for( index_t i = 0; i < interfaces.size(); ++i ) {
                name += model_.interfaces_[interfaces[i]].name() ;
                name += "_" ;
            }
            result = model_.nb_contacts() ;
            model_.contacts_.push_back( BoundaryModelElement( &model_, BM_CONTACT, result ) ) ;
            model_.contacts_[result].set_name( name ) ;
            model_.contacts_[result].set_geological_feature( type ) ;

            for( index_t i = 0; i < interfaces.size(); ++i ) {
                add_contact_in_boundary( result, interfaces[i] ) ;
            }
        }
        return result ;
    }

    /*! 
     * @brief Add a Line knowing only the indices of its points
     *
     * Used in Geomodeling to convert a surface to a model
     */
    void BoundaryModelBuilder::create_line( const std::vector< index_t > points ) {
        index_t id = model_.lines_.size() ;
        model_.lines_.push_back( Line( &model_, id, points ) ) ;
    }

    /*!
     * @brief Add a Surface to the model
     *
     * @param[in] interface_name Name of the parent. The parent MUST exist.
     * @param[in] type Type of the Surface
     * @param[in] key KeyFacet for this Surface
     */
    void BoundaryModelBuilder::create_surface(
        const std::string& interface_name,
        const std::string& type,
        const Surface::KeyFacet& key )
    {
        index_t parent = interface_id( interface_name ) ;
        if( interface_name != "" ) grgmesh_assert( parent != NO_ID ) ;

        index_t id = model_.nb_surfaces() ;
        GEOL_FEATURE t = determine_geological_type( type ) ;
        model_.surfaces_.push_back( Surface( &model_, id, parent, t ) ) ;
        model_.surfaces_[id].set_key_facet( key ) ;

        model_.interfaces_[parent].set_geological_feature( t ) ;
    }

    /*!
     * @brief Find or create a corner at a vertex.
     * 
     * @param[in] index Index of the vertex in the model
     * @return Index of the Corner
     */
    index_t BoundaryModelBuilder::find_or_create_corner( index_t index )
    {
        index_t result = find_corner( model_.vertex( index ) ) ;
        if( result == NO_ID ) {
            // Create the corner
            result = model_.nb_corners() ;
            model_.corners_.push_back( Corner( &model_, result, index ) ) ;
        }
        return result ;
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
        index_t corner0,
        index_t corner1,
        std::vector< index_t >& vertices )
    {
        grgmesh_assert( corner0 != NO_ID ) ;
        grgmesh_assert( corner1 != NO_ID ) ; 
        
        index_t result = find_line( corner0, corner1, vertices ) ;

        if( result == NO_ID ) {
            result = model_.nb_lines() ;
            model_.lines_.push_back(
                Line( &model_, result, corner0, corner1, vertices ) ) ;
        }
        return result ;
    }
  
    /*!
     * @brief Get the index of an Interface from its name
     * 
     * @param[in] name Name of the Interface
     * @return Index of the interface in the model, NO_ID if not found.
     */
    index_t BoundaryModelBuilder::interface_id( const std::string& name ) const
    {
        for( index_t i = 0; i < model_.nb_interfaces(); ++i ) {
            if( model_.one_interface(i).name() == name ) {
                return i ;
            }
        }
        return NO_ID ;
    }


    /*!
     * @brief Get the index of the Corner at given coordinates
     * @param[in] p Coordinates of the vertex
     * @return NO_ID or the index of the Corner
     */
    index_t BoundaryModelBuilder::find_corner( const vec3& p ) const
    {
        for( index_t i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).vertex() == p ) return i ;
        }
        return NO_ID ;
    }

    /*!
     * @brief Get the index of the Corner at a given model point
     * @param[in] p_id Index of the point
     * @return NO_ID or the index of the Corner
     */
    index_t BoundaryModelBuilder::find_corner( index_t p_id ) const 
    {
        for( index_t i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).vertex_ == p_id ) return i ;
        }
        return NO_ID ;
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
        index_t corner0,
        index_t corner1,
        const std::vector< index_t >& vertices ) const
    {
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            const Line& cp = model_.line(i) ;

            // Not the same number of vertices  = not the same Line
            if( cp.nb_vertices() != vertices.size() ) continue ; 

            if( corner0 == cp.boundary_id( 0 ) && corner1 == cp.boundary_id( 1 ) ) {
                if( std::equal( vertices.begin(), vertices.end(), cp.vertices_.begin() ) ){
                    return i ;
                }
            }

            if( corner1 == cp.boundary_id( 0 ) && corner0 == cp.boundary_id( 1 ) ) {
                if( std::equal( vertices.begin(), vertices.end(), cp.vertices_.rbegin() ) ){
                    return i ;
                }               
            }
        }
        return NO_ID ;
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
                if( equal ) {
                    return i ;
                }
            }
        }
        return NO_ID ;
    }


    /*!
     * @brief Rebuild a model ???
     * \todo Comment rebuild()... 
     * \todo Valgrind finds errors !!!!!!!
     * @return 
     */
    bool BoundaryModelBuilder::rebuild()
    {
        std::vector< index_t > sp_to_remove ;
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
        return true ;
    }

    /*!
     * @brief Update the indices stored by each element of the model \
     * according to its actual position in the corresponding vector in the model
     */
    void BoundaryModelBuilder::update_all_ids()
    {
        for( index_t co = 0; co < model_.nb_corners(); co++ ) {
            model_.corners_[co].id_ = co ;
        }
        for( index_t cp = 0; cp < model_.nb_lines(); cp++ ) {
            model_.lines_[cp].id_ = cp ;
        }
        for( index_t sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            model_.surfaces_[sp].id_ = sp ;
        }
        for( index_t c = 0; c < model_.nb_contacts(); c++ ) {
            model_.contacts_[c].id_ = c ;
        }
        for( index_t s = 0; s < model_.nb_interfaces(); s++ ) {
            model_.interfaces_[s].id_ = s ;
        }
        for( index_t r = 0; r < model_.nb_regions(); r++ ) {
            model_.regions_[r].id_ = r ;
        }
        for( index_t l = 0; l < model_.nb_layers(); l++ ) {
            model_.layers_[l].id_ = l ;
        }
    }

} 
