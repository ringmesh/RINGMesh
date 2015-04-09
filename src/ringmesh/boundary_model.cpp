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

#include <ringmesh/boundary_model.h>
#include <ringmesh/boundary_model_builder.h>
#include <ringmesh/utils.h>

#include <geogram/basic/logger.h>
#include <geogram/points/colocate.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <set>

namespace RINGMesh {


    /*! 
     * @brief Trigger an assertion if several vertices of a mesh at the same geometric location
     * @note Code modified from geomgram/mesh/mesh_repair.cpp 
     * @param[in] M the mesh
     * @param[in] colocate_epsilon tolerance      
     */
    void assert_no_colocate_vertices( GEO::Mesh& M, double colocate_epsilon) {
        GEO::vector<index_t> old2new;

        index_t nb_new_vertices = 0;
        if(colocate_epsilon == 0.0) {
            nb_new_vertices = GEO::Geom::colocate_by_lexico_sort(
                M.vertices.point_ptr(0), 3, M.vertices.nb(),
                old2new, M.vertices.dimension()
            );
        } else {
            nb_new_vertices = GEO::Geom::colocate(
                M.vertices.point_ptr(0), 3, M.vertices.nb(),
                old2new, colocate_epsilon, M.vertices.dimension()
            );
        }
        if( nb_new_vertices != M.vertices.nb() ) {
            geo_assert_not_reached ;
        }       
    }

    /*! 
     * @brief Merges the vertices of a mesh that are at the same geometric location
     * @note Code modified from geomgram/mesh/mesh_repair.cpp 
     * @param[in] M the mesh
     * @param[in] colocate_epsilon tolerance for merging vertices
     * @param[out] old2new mapping from previous M.vertices to new M.vertices
     */
    void repair_colocate_vertices(
        GEO::Mesh& M, 
        double colocate_epsilon, 
        GEO::vector<index_t>& old2new) 
    {
        old2new.clear() ;

        index_t nb_new_vertices = 0;
        if(colocate_epsilon == 0.0) {
            nb_new_vertices = GEO::Geom::colocate_by_lexico_sort(
                M.vertices.point_ptr(0), 3, M.vertices.nb(),
                old2new, M.vertices.dimension()
            );
        } else {
            nb_new_vertices = GEO::Geom::colocate(
                M.vertices.point_ptr(0), 3, M.vertices.nb(),
                old2new, colocate_epsilon, M.vertices.dimension()
            );
        }        
        if(nb_new_vertices == M.vertices.nb()) {
            return;
        }        
        for(index_t c = 0; c < M.facet_corners.nb(); c++) {
            M.facet_corners.set_vertex(c, old2new[M.facet_corners.vertex(c)]);
        }

        // Some index magic to flag the point to delete and the right 
        // mapping between old and new vertices of the mesh
        GEO::vector< index_t > to_delete( old2new.size() ) ;
        for(index_t i = 0; i < old2new.size(); i++) {
            if(old2new[i] == i) {
                to_delete[i] = 0;
            } else {
                to_delete[i] = 1;
            }
        }
        M.vertices.delete_elements( to_delete, false );

        // The to_delete vector is used for mapping in the delete_elements function
        // We need it to get the correct mapping
        for( index_t i=0; i< old2new.size(); i++ ) {
            if( to_delete[i] != NO_ID ) {
                old2new[i] = to_delete[i] ;
            }
            else {
                old2new[i] = to_delete[ old2new[i] ] ;
            }
        }
    }


    /************************************************************************/


    void BoundaryModelVertices::initialize_unique_vertices()
    {
        index_t nb_corners = bm_.nb_corners() ;
        index_t nb_lines = bm_.nb_lines() ;
        index_t nb_surfaces = bm_.nb_surfaces() ;

        // Total number of vertices in the Corners - Lines and Surfaces of the model
        index_t nb = bm_.nb_corners() ;

        for( index_t l = 0; l < nb_lines; l++ ) {
            nb += bm_.line( l ).nb_vertices() ;
        }
        for( index_t s = 0; s < nb_surfaces; s++ ) {
            nb += bm_.surface( s ).nb_vertices() ;
        }

        // Get out if the BME have no vertex
        if( nb == 0 ) {
            return ;
        }
        
        std::vector< vec3 > all_vertices( nb ) ;
        index_t index = 0 ;
        for( index_t c = 0; c < nb_corners; c++ ) {
            all_vertices[index++] = bm_.corner( c ).vertex() ;
        }
        for( index_t l = 0; l < nb_lines; l++ ) {
            const Line& line = bm_.line( l ) ;
            for( index_t v = 0; v < line.nb_vertices(); v++ ) {
                all_vertices[index++] = line.vertex( v ) ;
            }
        }
        for( index_t s = 0; s < nb_surfaces; s++ ) {
            const Surface& surface = bm_.surface( s ) ;
            for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
                all_vertices[index++] = surface.vertex( v ) ;
            }
        }
        
        unique_vertices_.vertices.create_vertices( all_vertices.size() ) ;
        unique_vertices_.vertices.assign_points( all_vertices[0].data(), 3, all_vertices.size() ) ;

        GEO::vector< index_t > old2new ;
        repair_colocate_vertices( unique_vertices_, epsilon, old2new ) ;   

        // We do the same loop as above
        index = 0 ;
        for( index_t c = 0; c < nb_corners; c++ ) {
             Corner& C = const_cast< Corner& >( bm_.corner( c ) ) ;
             C.set_model_vertex_id( old2new[index++] ) ;
             // I am crazy paranoid - Jeanne
             ringmesh_debug_assert( C.vertex() == unique_vertex( old2new[index-1] ) ) ;
        }
        for( index_t l = 0; l < nb_lines; l++ ) {
            Line& L = const_cast< Line& >( bm_.line( l ) ) ;
            for( index_t v = 0; v < L.nb_vertices(); v++ ) {
                L.set_model_vertex_id( v, old2new[index++] ) ;                
                ringmesh_debug_assert( L.vertex(v) == unique_vertex( old2new[index-1] ) ) ;
            }
        }
        for( index_t s = 0; s < nb_surfaces; s++ ) {
            Surface& S = const_cast< Surface& >( bm_.surface( s ) ) ;
            for( index_t v = 0; v < S.nb_vertices(); v++ ) {
                S.set_model_vertex_id( v, old2new[index++] ) ;
                ringmesh_debug_assert( S.vertex(v) == unique_vertex( old2new[index-1] ) ) ;
            }
        }

        // Paranoia - Jeanne
        assert_no_colocate_vertices( unique_vertices_, epsilon ) ;
    }

    
    void BoundaryModelVertices::initialize_reverse()
    {
        typedef BoundaryModelElement BME ;

        if( unique_vertices_.vertices.nb() == 0 ) {
            initialize_unique_vertices() ;
        }
        if( !unique2bme_.is_bound() ) {
            unique2bme_.bind( attribute_manager(), "unique2bme") ; 
        }
        
        for( index_t c = 0; c < bm_.nb_corners(); c++ ) {            
            unique2bme_[c].push_back( VertexInBME( BME::CORNER, c, 0 ) ) ;
        }
        for( index_t l = 0; l < bm_.nb_lines(); l++ ) {
            for( index_t v = 0; v <  bm_.line( l ).nb_vertices(); v++ ) {
                VertexInBME cur( BME::LINE, l, v ) ;
                unique2bme_[unique_vertex_id( cur )].push_back( cur ) ;
            }
        }
        for( index_t s = 0; s < bm_.nb_surfaces(); s++ ) {
            for( index_t v = 0; v < bm_.surface( s ).nb_vertices(); v++ ) {
                VertexInBME cur( BME::SURFACE, s, v ) ;
                unique2bme_[unique_vertex_id( cur )].push_back( cur ) ;
            }
        }
    }

    void BoundaryModelVertices::update_point( index_t v, const vec3& point )
    {
        const std::vector< VertexInBME >& bme_v = bme_vertices( v ) ;
        for( index_t i = 0; i < bme_v.size(); i++ ) {
            const VertexInBME& info = bme_v[i] ;
            const_cast< BoundaryModelElement& >( bm_.element( info.bme_type,
                info.bme_id ) ).set_vertex( info.v_id, point, false ) ;
        }
    }


    const std::vector< BoundaryModelVertices::VertexInBME >&
    BoundaryModelVertices::bme_vertices( index_t v ) const
    {
        if( !unique2bme_.is_bound() ) {
            const_cast< BoundaryModelVertices* >( this )->initialize_reverse() ;
        }
        return unique2bme_[v] ;
    }

    index_t BoundaryModelVertices::add_unique_vertex( const vec3& point ) 
    {
        return unique_vertices_.vertices.create_vertex( point.data() ) ;
    }

    void BoundaryModelVertices::add_unique_to_bme( 
        index_t unique_id, 
        BoundaryModelElement::TYPE bme_type, 
        index_t bme_id,
        index_t v_id ) 
    {
        if( !unique2bme_.is_bound() ) {
            unique2bme_.bind( attribute_manager(), "unique2bme") ; 
        }
        ringmesh_assert( unique_id < nb_unique_vertices() ) ;
        unique2bme_[unique_id].push_back( VertexInBME( bme_type, bme_id, v_id ) ) ;
    } 


    index_t BoundaryModelVertices::nb_unique_vertices() const
    {
        if( unique_vertices_.vertices.nb() == 0 ) {
            const_cast< BoundaryModelVertices* >( this )->initialize_unique_vertices() ;
        }
        return unique_vertices_.vertices.nb() ;
    }


    index_t BoundaryModelVertices::unique_vertex_id(
        BoundaryModelElement::TYPE type,
        index_t element,
        index_t v ) const
    {
        if( unique_vertices_.vertices.nb() == 0 ) {
            const_cast< BoundaryModelVertices* >( this )->initialize_unique_vertices() ;
        }
        ringmesh_assert( v < bm_.element( type, element ).nb_vertices() ) ;      
        return bm_.element( type, element ).model_vertex_id( v ) ;
    }


    index_t BoundaryModelVertices::unique_vertex_id(
        const VertexInBME& v ) const 
    {
        return unique_vertex_id( v.bme_type, v.bme_id, v.v_id ) ;
    }


    const vec3& BoundaryModelVertices::unique_vertex( index_t v ) const
    {       
        ringmesh_assert( v < nb_unique_vertices() ) ;
        return unique_vertices_.vertices.point(v) ;
    }

    void BoundaryModelVertices::clear() 
    {
        /// \todo Unbind all attributes !!!! otherwise we'll get a crash
        // For the moment 
        if( unique2bme_.is_bound() ) unique2bme_.unbind() ;
        
        unique_vertices_.clear( true, true ) ; 

        // Clear the information for the Corner - Line - Surface
        for( index_t c = 0; c < bm_.nb_corners(); c++ ) {
             Corner& C = const_cast< Corner& >( bm_.corner( c ) ) ;
             C.set_model_vertex_id( NO_ID ) ;
        }
        for( index_t l = 0; l < bm_.nb_lines(); l++ ) {
            Line& L = const_cast< Line& >( bm_.line( l ) ) ;
            for( index_t v = 0; v < L.nb_vertices(); v++ ) {
                L.set_model_vertex_id( v, NO_ID ) ;                
            }
        }
        for( index_t s = 0; s < bm_.nb_surfaces(); s++ ) {
            Surface& S = const_cast< Surface& >( bm_.surface( s ) ) ;
            for( index_t v = 0; v < S.nb_vertices(); v++ ) {
                S.set_model_vertex_id( v, NO_ID ) ;
            }
        }
    }

    /*******************************************************************************/

    
    BoundaryModel::~BoundaryModel()
    {
        for( index_t i = 0; i < corners_.size(); i++ ) {
            if( corners_[i] ) delete corners_[i] ;
        }
        for( index_t i = 0; i < lines_.size(); i++ ) {
            if( lines_[i] ) delete lines_[i] ;
        }
        for( index_t i = 0; i < surfaces_.size(); i++ ) {
            if( surfaces_[i] ) delete surfaces_[i] ;
        }
        for( index_t i = 0; i < regions_.size(); i++ ) {
            if( regions_[i] ) delete regions_[i] ;
        }
        for( index_t i = 0; i < contacts_.size(); i++ ) {
            if( contacts_[i] ) delete contacts_[i] ;
        }
        for( index_t i = 0; i < interfaces_.size(); i++ ) {
            if( interfaces_[i] ) delete interfaces_[i] ;
        }
        for( index_t i = 0; i < layers_.size(); i++ ) {
            if( layers_[i] ) delete layers_[i] ;
        }
    }

    /*!
     * @brief Total number of facets in the model Surface s
     */
    index_t BoundaryModel::nb_facets() const
    {
        index_t result = 0 ;
        for( index_t i = 0; i < nb_surfaces(); ++i ) {
            result += surface( i ).nb_cells() ;
        }
        return result ;
    }


    /*!
     * @brief Returns the index of the given vertex in the model
     * \todo Implement the function - Add a KdTree for geometrical request on model vertices
     *
     * @param[in] p input point coordinates
     * @return NO_ID
     */
    index_t BoundaryModel::vertex_index( const vec3& p ) const
    {
        ringmesh_assert_not_reached ;
        return NO_ID ;
    }


    /*!
     * @brief Returns the index of the region neighboring the surface.
     * @param[in] surface_part_id Index of the Surface
     * @param[in] side Side of the Surface
     * @return The region index or NO_ID if none found.
     */
    index_t BoundaryModel::find_region(
        index_t surface_part_id,
        bool side ) const
    {
        ringmesh_debug_assert( surface_part_id < nb_surfaces() ) ;
        for( index_t r = 0; r < nb_regions(); r++ ) {
            const BoundaryModelElement& cur_region = region( r ) ;
            for( index_t s = 0; s < cur_region.nb_boundaries(); s++ ) {
                if( cur_region.side( s ) == side
                    && cur_region.boundary_id( s ) == surface_part_id )
                {
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
                builder.set_parent( BoundaryModelElement::SURFACE,
                    one_interface( i ).child_id( 0 ), i ) ;
            }

            // Is it really useful to have contacts, let's hope not... I am not doing it
        }

        /// 2. Check that the Universe region exists
        /// \todo Write some code to create the universe (cf. line 805 to 834 de s2_b_model.cpp)
        if( universe_.name() != "Universe" ) {
            GEO::Logger::err( "" )
            <<
            "The region universe is not defined for the model. IMPLEMENTATION TO DO"
            << std::endl ;
            return false ;
        }

        /// 3. Check that each region has a name and valid surfaces
        for( index_t i = 0; i < regions_.size(); ++i ) {
            const BoundaryModelElement& region = this->region( i ) ;

            if( region.name() == "" ) {
                std::ostringstream name ;
                name << "region_" << i ;
                builder.set_element_name( BoundaryModelElement::REGION, i, name.str() ) ;
            }
            if( region.nb_boundaries() == 0 ) {
                GEO::Logger::err( "" ) << "The region " << region.name()
                                       << " has no Surfaces on its boundary" <<
                std::endl ;
                return false ;
            }
        }

        /// 4. Check that all the surfaces_ of the model are triangulated
        /// \todo Implement a triangulation function in SurfaceMutator
        for( index_t s = 0; s < nb_surfaces(); s++ ) {
            if( !surface( s ).is_triangulated() ) {
                GEO::Logger::err( "" ) << "Surface " << s <<
                " is not triangulated" << std::endl ;
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
            GEO::Logger::err( "" ) << "The BoundaryModel " << name_
                                   << " cannot be saved in .ml format " << std::endl ;
            return false ;
        }

        // Print Gocad Model3d headers
        out << "GOCAD Model3d 1" << std::endl << "HEADER {" << std::endl << "name:"
            << name() << std::endl << "}" << std::endl ;

        save_coordinate_system( out ) ;

        // Print the TSurf = Interface information
        for( index_t i = 0; i < nb_interfaces(); ++i ) {
            out << "TSURF " << one_interface( i ).name() << std::endl ;
        }

        index_t count = 1 ;

        // Print the TFace = Surface information
        for( index_t i = 0; i < nb_surfaces(); ++i ) {
            const Surface& s = surface( i ) ;
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

        for( index_t i = 0; i < nb_regions(); ++i ) {
            save_region( count, region( i ), out ) ;
            ++count ;
        }

        for( index_t i = 0; i < nb_layers(); ++i ) {
            save_layer( count, offset_layer, layer( i ), out ) ;
            ++count ;
        }
        out << "END" << std::endl ;

        // Save the geometry of the Surfaces (TFace), Interface (TSurf) by Interface
        for( index_t i = 0; i < nb_interfaces(); ++i ) {
            const BoundaryModelElement& tsurf = one_interface( i ) ;

            // Header
            out << "GOCAD TSurf 1" << std::endl << "HEADER {" << std::endl <<
            "name:"
                << tsurf.name() << std::endl << "name_in_model_list:" << tsurf.name()
                << std::endl << "}" << std::endl ;
            save_coordinate_system( out ) ;

            out << "GEOLOGICAL_FEATURE " << tsurf.name() << std::endl
                << "GEOLOGICAL_TYPE " ;
            out << BME::geol_name( tsurf.geological_feature() ) ;
            out << std::endl ;

            out << "PROPERTY_CLASS_HEADER Z {" << std::endl << "is_z:on" <<
            std::endl
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
                    out << "VRTX " << vertex_count << " " << sp.vertex( k ) <<
                    std::endl ;
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

                    vec3 c = cp.vertex( 0 ) ;
                    vec3 next = cp.vertex( 1 ) ;

                    // To be sure that we have all corners we need to ensure
                    // that all corners at the end of lines are saved too
                    std::vector< index_t > result ;
                    sp.tools.ann().get_colocated( cp.vertex( cp.nb_vertices() - 1 ), result ) ;
                    ringmesh_debug_assert( !result.empty() ) ;
                    set_end_corners.insert( result[0] + offset ) ;

                    result.clear() ;
                    sp.tools.ann().get_colocated( c, result ) ;
                    ringmesh_debug_assert( !result.empty() ) ;
                    index_t c_id = result[0] ;
                    result.clear() ;
                    sp.tools.ann().get_colocated( next, result ) ;
                    ringmesh_debug_assert( !result.empty() ) ;
                    index_t next_id = result[0] ;

                    ringmesh_assert( c_id != NO_ID && next_id != NO_ID ) ;

                    bstones.push_back( c_id + offset ) ;
                    next_vertex.push_back( next_id + offset ) ;
                }
            }

            // Print Corners and Lines
            std::vector< index_t > end_corners(
                set_end_corners.begin(), set_end_corners.end() ) ;
            std::vector< bool > end_corner_to_print( end_corners.size(), true ) ;

            for( index_t j = 0; j < bstones.size(); ++j ) {
                out << "BSTONE " << bstones[ j ] << std::endl ;

                // Determine the corners at the end of the lines that are not saved
                for( index_t k = 0; k < end_corners.size(); k++ ) {
                    if( bstones[ j ] == end_corners[ k ] ) {
                        end_corner_to_print[ k ] = false ;
                        break ;
                    }
                }
            }

            // Print the corners that were at the beginning of none of the contacts
            // in this Interface
            for( index_t j = 0; j < end_corners.size(); j++ ) {
                if( end_corner_to_print[ j ] ) {
                    out << "BSTONE " << end_corners[ j ] << std::endl ;
                }
            }

            // Print the the information to build the lines :
            // index of the vertex at the corner and index of the second vertex on the line
            for( index_t j = 0; j < bstones.size(); ++j ) {
                out << "BORDER " << vertex_count << " " << bstones[ j ] << " "
                    << next_vertex[ j ] << std::endl ;
                vertex_count++ ;
            }
            out << "END" << std::endl ;
        }
        return true ;
    }


    /*! To save the attributes in a Graphite readable file, we need to write the correct
     * keyword for the attribute type - We restrict ourselves to the 3 types
     * int          "integer"
     * double       "real"
     * float        "real"
     * bool         "boolean"
     */
    inline std::string alias_name( const std::string& in )
    {
        if( in == "int" ) {return "integer" ;} else if( in == "index" ) {
            return "integer" ;
        } else if( in == "double" ) {
            return "real" ;
        } else if( in ==
                   "float" )
        {
            return "real" ;
        } else if( in ==
                   "bool" )
        {
            return "boolean" ;
        }
        ringmesh_assert_not_reached ;
        return "" ;
    }


    /*!
     * @brief DEBUG function - Save the surfaces of the model with their facet attributes into an .eobj file.
     * @details WARNING We assume that all Surface have the same attributes - if not this function will most
     *  certainly crash.
     *
     * @param[in] file_name Name of the file
     *
     * \todo Make this function const
     *
     */
    void BoundaryModel::save_as_eobj_file( const std::string& file_name )
    {
        std::ofstream out ;
        out.open( file_name.c_str() ) ;
        if( out.bad() ) {
            std::cout << "Error when opening the file: " << file_name.c_str() <<
            std::endl ;
            return ;
        }
        out.precision( 16 ) ;
        std::vector< index_t > offset( nb_surfaces(), 0 ) ;
        index_t cur_offset = 0 ;

        // Write vertices once for each surface
        for( index_t s = 0; s < nb_surfaces(); s++ ) {
            const Surface& S = surface( s ) ;
            offset[ s ] = cur_offset ;
            for( index_t p = 0; p < S.nb_vertices(); p++ ) {
                const vec3& V = S.vertex( p ) ;
                out << "v"
                    << " " << V.x
                    << " " << V.y
                    << " " << V.z
                    << std::endl ;
            }
            cur_offset += S.nb_vertices() ;
        }

        // Write the facets for a each surface
        for( index_t s = 0; s < nb_surfaces(); s++ ) {
            const Surface& S = surface( s ) ;
            for( index_t f = 0; f < S.nb_cells(); f++ ) {
                out << "f" << " " ;
                for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                    out << offset[ s ] + S.surf_vertex_id( f, v ) + 1 << " " ;
                }
                out << std::endl ;
            }
        }

        // Write facet attributes
        {
            // Get the attributes that can be saved on the first Surface
//            std::vector< SerializedAttribute< BME::FACET > > facet_attribs ;
//            get_serializable_attributes( surface(
//                    0 ).facet_attribute_manager(), facet_attribs ) ;
//
//            for( index_t i = 0; i < facet_attribs.size(); i++ ) {
//                // Output global information on the attribute
//                out << "# attribute " << facet_attribs[ i ].name() << " facet "
//                    << alias_name( facet_attribs[ i ].type_name() )
//                    << std::endl ;
//            }
//            if( facet_attribs.size() > 0 ) {
//                // Global counter for all the facets of all surfaces
//                index_t count = 0 ;
//                for( index_t s = 0; s < nb_surfaces(); s++ ) {
//                    const Surface& S = surface( s ) ;
//
//                    std::vector< SerializedAttribute< BME::FACET > > cur_attribs ;
//                    get_serializable_attributes(
//                        S.facet_attribute_manager(), cur_attribs ) ;
//
//                    ringmesh_assert( cur_attribs.size() == facet_attribs.size() ) ;
//                    for( index_t i = 0; i < facet_attribs.size(); ++i ) {
//                        ringmesh_assert(
//                            facet_attribs[ i ].type_name() ==
//                            cur_attribs[ i ].type_name() &&
//                            facet_attribs[ i ].name() == cur_attribs[ i ].name() ) ;
//                    }
//
//                    // Output attributes values
//                    for( index_t f = 0; f < S.nb_cells(); f++ ) {
//                        out << "# attrs f " << count + 1 << " " ;
//                        serialize_write_attributes( out, f, cur_attribs ) ;
//                        out << std::endl ;
//                        count++ ;
//                    }
//                }
//            }
        }
    }


    /*!
     * @brief Debug: Save a Surface of the model in the file OBJ format is used
     */
    void BoundaryModel::save_surface_as_obj_file(
        index_t s,
        const std::string& file_name ) const
    {
        std::ofstream out ;
        out.open( file_name.c_str() ) ;
        if( out.bad() ) {
            std::cout << "Error when opening the file: " << file_name.c_str() <<
            std::endl ;
            return ;
        }
        out.precision( 16 ) ;
        const Surface& S = surface( s ) ;
        for( index_t p = 0; p < S.nb_vertices(); p++ ) {
            const vec3& V = S.vertex( p ) ;
            out << "v"
                << " " << V.x
                << " " << V.y
                << " " << V.z
                << std::endl ;
        }
        for( index_t f = 0; f < S.nb_cells(); f++ ) {
            out << "f" << " " ;
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                out << S.surf_vertex_id( f, v ) + 1 << " " ;
            }
            out << std::endl ;
        }
    }


    /*!
     * @brief Write in the out stream things to save for CONTACT, INTERFACE and LAYERS
     */
    void save_high_level_bme(
        std::ofstream& out,
        const BoundaryModelElement& E )
    {
        /// First line:  TYPE - ID - NAME - GEOL
        out << BoundaryModelElement::type_name( E.element_type() ) << " "
            << E.id() << " " ;
        if( E.has_name() ) { out << E.name() << " " ;} else { out << "no_name " ;}
        out <<  BoundaryModelElement::geol_name( E.geological_feature() )
            << std::endl ;

        /// Second line:  IDS of children
        for( index_t j = 0; j < E.nb_children(); ++j ) {
            out << " " << E.child_id( j ) ;
        }
        out << std::endl ;
    }


    /*!
     * @brief Save the BoundaryModel into a dedicated format bm
     */
    void BoundaryModel::save_bm_file( const std::string& file_name )
    {
        std::ofstream out ;
        out.open( file_name.c_str() ) ;
        if( out.bad() ) {
            std::cout << "Error when opening the file: " << file_name.c_str() <<
            std::endl ;
            return ;
        }
        out.precision( 16 ) ;

        out << "RINGMESH BOUNDARY MODEL" << std::endl ;
        out << "NAME " << name() << std::endl ;

        // Number of the different type of elements
        for( index_t i = BME::CORNER; i < BME::NO_TYPE; i++ ) {
            BME::TYPE type = (BME::TYPE) i ;
            out <<  "NB_" << BME::type_name( type ) << " " << nb_elements( type ) <<
            std::endl ;
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
        for( index_t i = 0; i < nb_regions(); ++i ) {
            const BoundaryModelElement& E = region( i ) ;

            // Save ID - NAME -
            out << BME::type_name( BME::REGION ) << " " << E.id() << " " ;
            if( E.has_name() ) {out << E.name() ;} else {out << "no_name" ;}
            out << std::endl ;

            // Second line Signed ids of boundary surfaces
            for( index_t j = 0; j < E.nb_boundaries(); ++j ) {
                if( E.side( j ) ) {out << "+" ;} else {out << "-" ;}
                out << E.boundary_id( j ) << " " ;
            }
            out << std::endl ;
        }

        // Universe
        out << "UNIVERSE " << std::endl ;
        for( index_t j = 0; j < universe().nb_boundaries(); ++j ) {
            if( universe().side( j ) ) {out << "+" ;} else {out << "-" ;}
            out << universe().boundary_id( j ) << " " ;
        }
        out << std::endl ;

//        // Vertices and attributes on vertices
//        out << "MODEL_VERTICES" << " " << nb_vertices() << std::endl ;
//        out << "MODEL_VERTEX_ATTRIBUTES " ;
//        std::vector< SerializedAttribute< VERTEX > > vertex_attribs ;
//        get_serializable_attributes( &vertex_attribute_manager_, vertex_attribs, out ) ;
//        for( index_t i = 0; i < nb_vertices(); ++i ) {
//            out << vertex( i )  << "  " ;
//            serialize_write_attributes( out, i, vertex_attribs ) ;
//            out << std::endl ;
//        }

        // Corners
        for( index_t i = 0; i < nb_corners(); ++i ) {
            out << BME::type_name( BME::CORNER ) << " "
                << corner( i ).id() << " " << corner( i ).vertex() <<
            std::endl ;
        }

        // Lines
        for( index_t i = 0; i < nb_lines(); ++i ) {
            const Line& L = line( i ) ;
            out << BME::type_name( BME::LINE ) << " " << L.id() << std::endl ;
            out << "LINE_VERTICES " << L.nb_vertices() << std::endl ;
            for( index_t j = 0; j < L.nb_vertices(); ++j ) {
                out << L.vertex( j ) << std::endl ;
            }
//            out << "LINE_VERTEX_ATTRIBUTES " ;
//            std::vector< SerializedAttribute< BME::VERTEX > > line_vertex_attribs ;
//            get_serializable_attributes(
//                L.vertex_attribute_manager(), line_vertex_attribs, out ) ;
//            for( index_t j = 0; j < L.nb_vertices(); ++j ) {
//                out << j << "  " ;
//                serialize_write_attributes( out, j, line_vertex_attribs ) ;
//                out << std::endl ;
//            }
//            out << "LINE_SEGMENT_ATTRIBUTES " ;
//            std::vector< SerializedAttribute< BME::FACET > > line_segments_attribs ;
//            get_serializable_attributes(
//                L.facet_attribute_manager(), line_segments_attribs, out ) ;
//            if( line_segments_attribs.size() > 0 ) {
//                for( index_t j = 0; j < L.nb_cells(); ++j ) {
//                    out << j << "  " ;
//                    serialize_write_attributes( out, j, line_segments_attribs ) ;
//                    out << std::endl ;
//                }
//            }
            out << "IN_BOUNDARY " ;
            for( index_t j = 0; j < L.nb_in_boundary(); ++j ) {
                out << L.in_boundary_id( j ) << " " ;
            }
            out << std::endl ;
        }

        // Surfaces
        for( index_t i = 0; i < nb_surfaces(); ++i ) {
            const Surface& S = surface( i ) ;
            out << BME::type_name( BME::SURFACE ) << " " << S.id() << std::endl ;
            out << "SURFACE_VERTICES " << S.nb_vertices() << std::endl ;
            for( index_t j = 0; j < S.nb_vertices(); ++j ) {
                out << S.vertex( j ) << std::endl ;
            }
//            out << "SURFACE_VERTEX_ATTRIBUTES " ;
//            std::vector< SerializedAttribute< BME::VERTEX > > surface_vertex_attribs ;
//            get_serializable_attributes(
//                S.vertex_attribute_manager(), surface_vertex_attribs, out ) ;
//            for( index_t j = 0; j < S.nb_vertices(); ++j ) {
//                out << j << "  " ;
//                serialize_write_attributes( out, j, surface_vertex_attribs ) ;
//                out << std::endl ;
//            }

            out << "SURFACE_CORNERS " << S.nb_corners() << std::endl ;
            out << "SURFACE_FACETS " << S.nb_cells() << std::endl ;
//            out << "SURFACE_FACET_ATTRIBUTES " ;
//            std::vector< SerializedAttribute< BME::FACET > > surface_facet_attribs ;
//            get_serializable_attributes(
//                S.facet_attribute_manager(), surface_facet_attribs, out ) ;
//
            for( index_t j = 0; j < S.nb_cells(); ++j ) {
                out << S.nb_vertices_in_facet( j ) << " " ;
                for( index_t v = 0; v < S.nb_vertices_in_facet( j ); ++v ) {
                    out << S.surf_vertex_id( j, v ) << " " ;
                }
//                serialize_write_attributes( out, j, surface_facet_attribs ) ;
                out << std::endl ;
            }
        }
    }

    signed_index_t BoundaryModel::find_interface( const std::string& name) const {
        for(index_t i = 0 ; i < nb_interfaces() ; i++ ) {
            if( one_interface(i).name() == name ) {
                return i ;
            }
        }
        GEO::Logger::err("") << "Surface name did not match with an actual\
                                interface name of the Boundary Model. Abort.. " 
                                << std::endl ;
        ringmesh_assert_not_reached ;
        return -1 ;
    }

    signed_index_t BoundaryModel::find_region( const std::string& name) const {
        for(index_t r = 0 ; r < nb_regions() ; r++ ) {
            if( region(r).name() == name ) {
                return r ;
            }
        }
        GEO::Logger::err("") << "Region name did not match with an actual\
                                region name of the Boundary Model. Abort.. " 
                                << std::endl ;
        ringmesh_assert_not_reached ;
        return -1 ;
    }


} // namespace
