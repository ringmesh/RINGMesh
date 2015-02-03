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


#ifndef __GRGMESH_BOUNDARY_MODEL_FROM_SURFACE__
#define __GRGMESH_BOUNDARY_MODEL_FROM_SURFACE__

#include <grgmesh/common.h>
#include <grgmesh/boundary_model.h>

#include <vector> 
#include <string>

namespace GRGMesh {

    /*!
     * @brief Builder of a BoundaryModel form a surface mesh
     * The class MESH should implement :
     *
     *
     */
    template< MESH > 
    class GRGMESH_API BoundaryModelFromSurface : public BoundaryModelBuilder {
        BoundaryModelFromSurface( const MESH& surface, BoundaryModel& model  ) ;
        virtual ~BoundaryModelFromSurface() {};

        void build_model() ;    
    
    private:      
        void create_surfaces() ;
        void create_lines_and_corners() ;
        
        void determine_connected_components(
            const std::vector< index_t >& colocalized,  // colocalized Corners ? halfedges ?
            std::vector< index_t >& adjacent_parts ) ;
    

    private:
        const MESH& S_;
        
        // Index of the connected component of one Surface facet
        std::vector< index_t > surface_id_ ;  // sp_id_ ;
        // Index of of the line connected componenet for each Surface halfedge - NO_ID else
        // Sans doute mieux de stocker ceux pour qui on a effectivement une valeur
        std::vector< index_t > line_id_ ; // cp_id_ ;
        // Index of the Corner for each Surface vertex  - NO_ID else
        // Sans doute mieux de stocker ceux pour qui on a effectivement une valeur
        std::vector< index_t > corner_id_ ; // c_id_;

        // Utils
        MapVertexAttribute< int > v_id_ ;
        int nb_surfaces_ ;
    } ;

    
    /*!
     * @brief Class to Detect colocated edges from a MESH
     *
     * EN IMPLEMENTER UNE AUTRE ????
     */ 
    template< MESH > class Colocater {
    public: 
        Colocater::~Colocater() { 
        }        

        Colocater::Colocater( 
            index_t nb_surfaces,
            const std::vector< index_t >& in_edges, 
            MapHalfedgeAttribute< unsigned int >& id 
            ):  halfedges_(in),
            id_(id),
            nb_elements_per_query_( 4*nb_sp )
        {
            if( !halfedges_.empty() ) {
                int nb_vertices = in.size() ;
                ann_points_ = annAllocPts( nb_vertices, 3 ) ;
                // Take the barycenters
                for( unsigned int i = 0; i < nb_vertices; i++ ) {
                    vec3 vv = (in[i]->vertex()->point() + in[i]->opposite()->vertex()->point() ) /2. ;
                    std::copy( vv.data(), vv.data() + 3, &ann_points_[i][0] ) ;
                }
                ann_tree_ = new ANNkd_tree( &ann_points_[0], nb_vertices, 3 ) ;
            }
        }

        void Colocater::get_colocated_halfedges(
            Map::Halfedge* h, 
            std::vector< Map::Halfedge* >& result 
            ) const {
                result.clear() ;
                if( halfedges_.empty() ) return ;

                unsigned int h1 = id_[h] ;
                ogf_debug_assert( h1 < halfedges_.size() ) ;
                std::vector< int > neighbors ;
                int nb_neighbors = ogf_min( nb_elements_per_query_, ann_tree_->nPoints() ) ;
                neighbors.resize( nb_neighbors ) ;
                ANNdistArray dist = new ANNdist[nb_neighbors] ;
                vec3 v = (h->vertex()->point() + h->opposite()->vertex()->point() ) /2. ;
                ann_tree_->annkSearch( v.data(), nb_neighbors, &neighbors[0], dist ) ;
                delete[] dist ;

                for( int i = 0; i < neighbors.size(); ++i ){
                    int h2 = neighbors[i] ;
                    if( equal( h1, h2 ) ){
                        result.push_back( halfedges_[h2] ) ;
                    }
                    else break ;
                }        
                ogf_debug_assert( !result.empty() ) ;
        }

        bool Colocater::equal( int i1, int i2 ) const {
            ogf_assert( i1 <halfedges_.size() ) ;
            ogf_assert( i2 <halfedges_.size() ) ;
            Map::Halfedge* h1 = halfedges_[i1] ;
            Map::Halfedge* h2 = halfedges_[i2] ;

            const vec3& p11 = h1->vertex()->point() ;
            const vec3& p12 = h1->opposite()->vertex()->point() ;
            const vec3& p21 = h2->vertex()->point() ;
            const vec3& p22 = h2->opposite()->vertex()->point() ;

            if( p11 == p21 && p12 == p22 ) return true ;
            if( p11 == p22 && p12 == p21 ) return true ;
            return false ;
        }

    private:
        const std::vector< Map::Halfedge* > halfedges_ ;
        MapHalfedgeAttribute< unsigned int >& id_ ;
        ANNpointArray ann_points_ ;
        ANNkd_tree* ann_tree_ ;
        int nb_elements_per_query_ ;

    } ;


    template< MESH > 
    BoundaryModelFromSurface::BoundaryModelFromSurface( const MESH& surface, BoundaryModel& model  )
        :BoundaryModel( model ), S_( surface )
    {            
        surface_id_.resize( S_.nb_facets(), NO_ID ) ;
        line_id_.resize( S_.nb_corners(), NO_ID ) ;
        corner_id_.resize( S_.nb_vertices(), NO_ID ) ;            
    }
  
    template< MESH > 
    void BoundaryModelFromSurface::build_model() {
        /// 1. Attribute to each facet of the MESH a surface connected component
        /// index and create the Surface in the model
        create_surfaces() ;

        /// 2. Compute the contacts and the corners and the regions
        create_lines_and_corners() ;
    }

    template< MESH > 
    void BoundaryModelFromSurface::create_surfaces() {      
        // Info to store for the vertices of the surface
        std::vector< index_t > id_in_part( S_.nb_vertices(), NO_ID ) ;
        std::vector< index_t > part      ( S_.nb_vertices(), NO_ID ) ;

        std::vector< index_t > nb_vertices_per_part ;        
        {   
            index_t cur_part = 0 ;
            for( index_t it = 0 ; it < S_.nb_vertices(); ++it ) {
                if( part[it] == NO_ID ) {
                    index_t cur_id_in_part = 0 ;
        
                    std::stack< index_t > stack ;
                    stack.push(it) ;                    
                    while( !stack.empty() ) {
                        index_t v = stack.top() ;
                        stack.pop() ;
                 
                        if( part[v] != NO_ID ) continue ;

                        part[v] = cur_part ;
                        id_in_part[v] = cur_id_in_part ;
                        cur_id_in_part++ ;

                        // Get the neighboring vertices 
                        for( index_t j = 0; j < S_.nb_neighbors( v ); ++j ){
                            if( part[ S_.neighbors( v, j ) ] == NO_ID )
                                stack.push( S_.neighbors( v, j ) ) ;
                        }
                    }
                    // Add an empty Surface to the BoundaryModel
                    // Parent, key_facet, type are not known
                    builder_.create_surface( ) ;
                    nb_vertices_per_part.push_back( cur_id_in_part ) ;
                    cur_part++ ;
                }
            }
        }

        nb_surfaces_ = nb_vertices_per_part.size() ;
        // Paranoia
        grgmesh_assert( nb_surfaces_ == model_->nb_surfaces() ) ;

        // Get facet connected components 
        // Fill the Surface
        std::vector< index_t > points ;
        points.reserve( *(std::max_element( nb_vertices_per_part.begin(), nb_vertices_per_part.end() )) ) ;
        
        std::vector< index_t > facets ;
        facets.reserve ( S_.nb_corners() ) ;
        std::vector< index_t > facets_ptr ;
        facets_ptr.reserve ( S_->nb_facets() ) ;


        std::vector< index_t > facet_id_in_part( S_.nb_facets(), NO_ID ) ;
        
        for( index_t it = 0; it < S_.nb_facets(), it++ ) {
            if( surface_id_[it] == NO_ID ) {
                index_t cur_id = part[ S_.vertex( it, 0 ) ] ;
               
                // Get the facets that are in the same connected component than the current facet
                points.resize( nb_vertices_per_part[cur_id] ) ;
                std::vector< bool > visited( nb_vertices_per_part[cur_id] ) ;
                index_t cur_facet_id = 0 ;
                facets.resize( 0 ) ;
                facets_ptr.resize( 0 ) ;
                facets_ptr.push_back( 0 ) ;
                 
                std::stack< index_t > stack ;
                stack.push(it) ;
                while( !stack.empty() ) {
                    index_t f = stack.top() ;
                    stack.pop() ;
                    
                    if( surface_id_[f] != NO_ID ) continue ;
                    surface_id_[f] = cur_id ;
                    facet_id_in_part[f] = cur_facet_id++ ;
                     
                    for( index_t j = 0 ; j < S_.nb_facet_vertices( f ); ++j )
                    {
                        index_t v = S_.vertex( f, j ) ;
                        index_t id = id_in_part[v] ;
                        if( !visited[id] ) {
                            points[id] = builder_.add_vertex( v->point().data() ) ;
                            visited[id] = true ;
                        }
                        facets.push_back( id ) ;
                        index_t neighbor = S_.adjacent( f, j ) ;
                        if( neighbor != NO_ID && surface_id_[neighbor] == NO_ID ) stack.push( neighbor ) ;
                    }
                    facets_ptr.push_back( facets.size() ) ;
                }
                grgmesh_assert( std::count( visited.begin(), visited.end(), false ) == 0 ) ;
              
                // Set the points and triangles of the part, default orientation is
                // kept and the adjacent triangles are not computed
                builder_.set_surface_geometry(
                    cur_id, points, facets, facets_ptr ) ;

                builder_.set_surface_adjacencies( cur_id ) ;
        
                cur_id++ ;
            }
        }
    }

    template< MESH > 
    void BoundaryModelFromSurface::determine_parts( 
        const std::vector< index_t >& colocalized, 
        std::vector< index_t >& adjacent_parts 
    ){
        for( index_t j = 0; j < colocalized.size(); ++j ){
            index_t cur = colocalized[j] ;
            grgmesh_assert( line_id_[cur] == NO_ID ) ;        
            adjacent_parts[j] = surface_id_[ S_.facet( cur ) ]; // Access to facet id from corner id
        }
        std::sort( adjacent_parts.begin(), adjacent_parts.end() ) ;
        int nn = std::unique( adjacent_parts.begin(), adjacent_parts.end() )-adjacent_parts.begin() ;
        adjacent_parts.resize(nn) ;
    }
 
    /*! Computes and creates the corners, contacts and regions of the 
     *  using geometrical tests
     */
    template< MESH > 
    void BoundaryModelFromSurface::create_lines_and_corners() {
      
        /// 1. Gather corners( edge associated to a vertex in a facet) on the border of the surface

        /// 1.1 Flag all the halfedges that are on the borders
        index_t nb = 0 ;
        {
            for( index_t j = 0; j < S_.nb_corners(); ++j ){
                if( S_.adjacent( j ) == NO_ID ) {
                    line_id_[j] = 0 ; 
                    ++nb ;
                }                
            }
        }
        /// 1.2 Fill a vector of the halfedges on the border
        std::vector< index_t > on_border( nb, NO_ID ) ;
        
        // Atttribute to store the id of each halfedge in on_border
        std::vector< index_t > id( S_.nb_corners(), NO_ID ) ;
        {
            index_t default_id = nb+1;         
            index_t cur = 0 ;            
            for( index_t j = 0; j < S_.nb_corners(); ++j ){
                if( line_id_[j] == 0 ){
                    on_border[cur] = j ;
                    id[j] = cur ;
                    ++cur ;
                }     
                else id[j] = default_id ;
            }
        }

        /// 1.3 Build a Colocater to find halfedges that have the same vertices     
        Colocater colocalization( nb_surfaces_, on_border, id ) ;
        
        // Il faut aussi porter ça


        /// 2. Build the contact parts 
        // For each contact part built keep 
        // The ids of the corners (2 per contact)
        std::vector< index_t > contact_corners ;
        // The adjacent surfaces
        std::vector< std::vector < index_t > > surfaces_in_contact ;
        // Info to build the regions
        std::vector< ContactSort > regions_info;

        {
            ///////// Vectors used in the loop, allocated here to avoid multiple reallocation
            std::vector< index_t > colocalized ;
            colocalized.reserve( 2*nb_surfaces_ ) ;  // should be enough

            // Vertices on the contact (used to build the GRGMesh::Line)
            std::vector< index_t > points ;
            points.reserve( nb + 1 ) ;

            std::vector< index_t > one_contact ;
            one_contact.reserve( nb ) ;
            index_t cur_contact_id = 1 ;   //// On commence à 1 au lieu de 0 car à  0 signifie qu'on est au bord
            ///////////////////////////////////////

            for( index_t i = 0; i < nb; ++i ){
                index_t h = on_border[i] ;
                if( line_id_[h] > 0 ) continue ; /// PROBLEME A CAUSE DU 0

                one_contact.clear() ;
                points.clear() ;
                                     
                one_contact.push_back( h ) ;

                index_t vc0 = S_.vertex_at_corner(h) ;
                index_t vc1 = S_.next_vertex_from corner(h) ;
                points.push_back( model_->nb_vertices() ) ;
                builder_.add_vertex( S_.point( vc0 ) ) ;
                points.push_back( model_->nb_vertices() ) ;
                builder_.add_vertex( S_.point( vc1 ) ) ;
                
                // There is at least one colocalized halfedge: opposite
                colocalization.get_colocated_halfedges( h, colocalized ) ;            
                // Get ids of the surfaces in contact along the segment           
                std::vector< index_t > adjacent_parts ( colocalized.size() ) ;
                determine_parts( colocalized, adjacent_parts ) ; 
              
                one_contact.insert( one_contact.end(), colocalized.begin(), colocalized.end() ) ;

                // Get info to sort GRGMesh::Surfaces around the contact
                regions_info.push_back( ContactSort( cur_contact_id ) ) ;
                for( index_t j = 0; j < colocalized.size(); ++j ){
                    index_t h = colocalized[j] ;                    
                    index_t f = S_.facet( h ) 
                    regions_info.back().add_surface( 
                        surface_id_[ f ], 
                        S_.point( S_.next_in_facet( h ) )// ; h->opposite()->vertex()->point(),
                        S_.point( h ) ; // h->vertex()->point(),
                        S_.point( S_.next_in_facet( S_.next_inf_facet(h) ) ) // h->next()->vertex()->point() 
                        ) ;
                }                
               
                // Build the contact propating forward on the border
                // While the adjacent surfaces are the same the points and
                // halfedges are added to the GRGMesh::Line built
                bool same_surfaces = true ;            
                index_t next = S_.next_on_border( h ) ;
                do { 
                    if( line_id_[ next ] == 0 ){
                        colocalization.get_colocated_halfedges( next, colocalized ) ;
                        
                        std::vector< index_t > cur_parts ( colocalized.size() ) ;
                        determine_parts( colocalized, cur_parts ) ;

                        if( cur_parts.size() == adjacent_parts.size() &&
                            std::equal( cur_parts.begin(), cur_parts.end(), adjacent_parts.begin() )
                        ){
                            one_contact.push_back( next ) ;
                            points.insert( points.begin(), model_->nb_vertices() ) ;
                            builder_.add_vertex( S_.point( next) ) ;
                            vc0 = S_.vertex( next ) ;
                            one_contact.insert( one_contact.end(), colocalized.begin(), colocalized.end() ) ;
                        }
                        else same_surfaces = false ;                    
                    }
                    else same_surfaces = false ;            
                    next = next_on_border( next ) ;

                } while( same_surfaces && next != h ) ;

                bool is_a_loop = (next == h) ;

                if( !is_a_loop ) {
                    // Propagate backward to reach the other extremity 
                    same_surfaces = true ;            
                    index_t prev = S_.prev_on_border( h ) ;
                    do { 
                        grgmesh_assert( prev != h ) ;
                        if( line_id_[ prev ] == -1 ){
                            colocalization.get_colocated_halfedges( prev, colocalized ) ;
                            
                            std::vector< index_t > cur_parts ( colocalized.size() ) ;
                            determine_parts( colocalized, cur_parts ) ;

                            if( cur_parts.size() == adjacent_parts.size() &&
                                std::equal(cur_parts.begin(), cur_parts.end(),adjacent_parts.begin() )
                            ){
                                one_contact.push_back( prev ) ;
                                points.push_back( model_->nb_vertices() ) ;
                                builder_.add_vertex( S_.point( prev ) ) ;
                                vc1 = S_.vertex( prev ) ;
                                one_contact.insert( one_contact.end(), colocalized.begin(), colocalized.end() ) ;
                            }
                            else same_surfaces = false ;                    
                        }
                        else same_surfaces = false ;            
                        prev = S_.prev_on_border( h ) ;

                    } while( same_surfaces ) ;
                }
               
                // Now we have all the points on the contact 
                // The first and the last are corners if the contact is not a loop
                // Find or create the corners and fill the corner attribute on the surface
                int corner0 = corner_id_[vc0] ;
                if( corner0 < 0 ) {
                    corner0 = builder_.find_or_create_corner( points[0] ) ;
                    corner_id_[vc0] = corner0 ;
                }
                int corner1 = corner_id_[vc1] ;
                if( corner1 < 0 ) {
                    corner1 = builder_.find_or_create_corner( points.back() ) ;
                    corner_id_[vc1] = corner1 ;
                }
           
                // Create the current Line
                builder_.create_line( points ) ; //  cur_contact_id, points ) ;
                
                contact_corners.push_back(corner0) ;
                contact_corners.push_back(corner1) ;
                surfaces_in_contact.push_back(adjacent_parts) ;

                // Flag all the CORNER that were found along this contact
                for( index_t j = 0; j< one_contact.size(); ++j ) {
                    line_id_[one_contact[j]] = cur_contact_id ;
                }
                cur_contact_id++ ;
            } // End For
        }

        // Add links between elements of the model         
        grgmesh_assert( surfaces_in_contact.size() == contact_corners.size()/2 ) ;
        for( index_t i = 0 ; i < surfaces_in_contact.size(); ++i ) {

            // Corner-Contact and contact corner
            builder_.add_element_in_boundary( GRGMesh::BM_CORNER, contact_corners[2*i], i ) ;
            builder_.add_element_in_boundary( GRGMesh::BM_CORNER, contact_corners[2*i+1], i ) ;
            builder_.add_element_boundary( GRGMesh::BM_LINE, i, contact_corners[2*i] ) ;
            builder_.add_element_boundary( GRGMesh::BM_LINE, i, contact_corners[2*i+1] ) ;
                    
            // Interface->Contact and Contact->Interface
            const std::vector< index_t >& surfs = surfaces_in_contact[i] ;
            for( index_t j = 0; j < surfs.size(); ++j ) {
                // A revérifier c'est n'imp ça j'ai l'impression
                builder_.add_element_in_boundary( GRGMesh::BM_SURFACE, surfs[j], i ) ;
                builder_.add_element_boundary( GRGMesh::BM_LINE, i, j ) ;
            }            
        }

        /// 3. Build the regions 

        // Sort surfaces around the contacts
        for( index_t i = 0 ; i < regions_info.size(); ++ i) {
            regions_info[i].sort() ;
        }
        
        
        if( model_->nb_surfaces() == 1 ) {
            builder_.create_region() ;
            // I am not sure this is the correct side (Jeanne)
            builder_.add_element_boundary( GRGMesh::BM_REGION, 0, 0, true ) ;
            builder_.add_element_in_boundary( GRGMesh::BM_SURFACE, 0, 0 ) ;
        }
        else if( model_->nb_surfaces() > 1 ) {
            index_t cur_region_id  = 0 ;
            // For each side of each GRGMesh::Surface store the region in it is
            std::vector< index_t > surf_2_region ( 2*model_->nb_surfaces(), NO_ID ) ;

            // Start with the first Interface on its + side            
            std::stack< std::pair< index_t, bool > > S ;
            S.push( std::pair< index_t, bool > ( 0, true ) ) ;

            while( !S.empty() ){
                std::pair< index_t, bool > cur = S.top() ;
                S.pop() ;

                // Already visited
                if( surf_2_region[ cur.second == true ? 2*cur.first : 2*cur.first+1 ] != NO_ID ) continue ;
                                                
                // Create a new region
                builder_.create_region() ;

                std::stack< std::pair< index_t, bool > > SR ;
                SR.push( cur ) ;

                while( !SR.empty() ) {
                    std::pair< index_t, bool > s = SR.top() ;
                    SR.pop() ;

                    index_t s_id = s.second == true ? 2*s.first : 2*s.first+1 ;
                    // Already in a region
                    if( surf_2_region[ s_id ] != NO_ID ) continue ;
                    
                    // Add the surface to the current region
                    builder_.add_element_boundary( GRGMesh::BM_REGION, cur_region_id, s.first, s.second ) ;
                    surf_2_region[ s_id ] = cur_region_id ;

                    // Check the other side of the surface and push it in S
                    index_t s_id_opp = !s.second == true ? 2*s.first : 2*s.first+1 ;
                    if( surf_2_region[ s_id_opp ] == -1 ) 
                        S.push( std::pair< index_t, bool >(s.first, !s.second ) );                    

                    // For each contact, push the next oriented surface in the region
                    const GRGMesh::BoundaryModelElement& surface_part = model_->surface( s.first ) ;
                    for( index_t i = 0 ; i < surface_part.nb_boundaries(); ++i ){
                        index_t cp_id = surface_part.boundary( i ).id() ;
                        grgmesh_assert( cp_id < regions_info.size() ) ;
                        const std::pair< index_t, bool >& n = regions_info[ cp_id ].next( s ) ;
                        index_t n_id =  n.second == true ? 2*n.first : 2*n.first+1 ;
                        if( surf_2_region[ n_id ] == -1 )
                            SR.push( n ) ;
                    }                
                }
                cur_region_id ++ ;
            }
            
            // Check if all the surfaces were visited
            // If not, this means that there are additionnal regions included in those built
            // TO DO implement the code to take into account included region building
            grgmesh_assert( std::count( surf_2_region.begin(), surf_2_region.end(), -1 ) == 0 ) ;        
        }

        // We need to remove form the regions_ the one corresponding
        // to the universe_ that is the one with the biggest volume
        double max_volume = -1. ;
        index_t universe_id = NO_ID ;
        for( index_t i = 0; i < model_->nb_regions(); ++i ){
            double cur_volume = BoundaryModelElementMeasure::size( &model_->region(i) ) ;
            if( cur_volume > max_volume ) {
                max_volume = cur_volume ;
                universe_id = i ;
            }
        }
        grgmesh_assert( universe_id != NO_ID ) ;
        
        const BoundaryModelElement& cur_region = model_->region(universe_id) ;
        std::vector< std::pair< index_t, bool > > univ_boundaries( cur_region.nb_boundaries() ) ;
        for( index_t i = 0; i < cur_region.nb_boundaries(); ++i ){
            univ_boundaries[i].first = cur_region.boundary( i ).id() ;
            univ_boundaries[i].second = cur_region.side( i ) ;
        }
        builder_.set_universe( univ_boundaries ) ;

        // Decrease by one the ids of the regions that are after the
        // one converted to the universe
        // Remove the region converted to universe from the regions
        builder_.remove_universe_from_regions( universe_id ) ;
               
        // Set the links surface to region
        for( index_t i = 0; i < model_->nb_regions(); ++i ) {
            const BoundaryModelElement& region = model_->region(i) ;
            for( index_t j = 0; j < region.nb_boundaries(); ++j ) {
                builder_.add_element_in_boundary( SURFACE, region.boundary( j ).id(), i ) ;
            }
        }
      
    } 
}

#endif