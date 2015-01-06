/*[
* Association Scientifique pour la Geologie et ses Applications (ASGA)
* Copyright (c) 1993-2013 ASGA. All Rights Reserved.
*
* This program is a Trade Secret of the ASGA and it is not to be:
* - reproduced, published, or disclosed to other,
* - distributed or displayed,
* - used for purposes or on Sites other than described
*   in the GOCAD Advancement Agreement,
* without the prior written authorization of the ASGA. Licencee
* agrees to attach or embed this Notice on all copies of the program,
* including partial copies or modified versions thereof.
]*/ 

/*! \author Jeanne Pellerin */

#ifndef __GRGMESH_BOUNDARY_MODEL__
#define __GRGMESH_BOUNDARY_MODEL__

#include <grgmesh/common.h>
#include <grgmesh/boundary_model_element.h>

#include <vector> 
#include <string>

namespace GRGMesh {
    struct Border ;
    class BoundaryModelBuilder ;
}

namespace GRGMesh {
    
/********************************************************************************************/
/***********             BoundaryModel class declaration              ************************/
/********************************************************************************************/

    class GRGMESH_API BoundaryModel {
        friend class BoundaryModelBuilder ;
        friend class ContactPartMutator ;
        friend class SurfacePartMutator ;
    public:
        BoundaryModel() ;
        ~BoundaryModel() ;

        Box3d bbox() const ;

        /* To print global information on the Model 
         * and for each one of its components in this stream 
         */
        void print_topology( std::ofstream& out ) const ;
        void print_element_info( std::ofstream& out ) const ;

        /* Accessors, counters */
        signed_index_t surface_id( const std::string& name ) const ;
        signed_index_t find_corner( const vec3& ) const ;
        signed_index_t find_contact( const std::vector< signed_index_t >& interfaces ) const ;
        signed_index_t find_contact_part(
            signed_index_t corner0,
            signed_index_t corner1,
            const std::vector< vec3 >& points ) const ;
        signed_index_t find_contact_part( vec3 point0, vec3 point1 ) const ;

        std::vector< BoundaryModelElement* > find_interface_parts(
            const std::vector< const BoundaryModelElement* >& contact_parts,
            bool consider_free_contacts = false ) const ;

        std::vector< BoundaryModelElement* > find_contacts(
            const BoundaryModelElement* corner0,
            const BoundaryModelElement* corner1 ) const ;

        std::vector< BoundaryModelElement* > get_closed_contacts() const ;
        signed_index_t find_region( signed_index_t surface_part_id, bool side ) const ;

        signed_index_t nb_real_corners() const ;
        signed_index_t nb_non_manifold_contact_parts() const ;
        signed_index_t nb_surface_with_free_boundary() const ;

        //Accessors
        index_t nb_points() const { return points_.size() ; }
        index_t nb_corners() const { return corners_.size() ; }
        index_t nb_contact_parts() const { return contact_parts_.size() ; }
        index_t nb_surface_parts() const { return surface_parts_.size() ; }
        index_t nb_contacts() const { return contacts_.size() ; }
        index_t nb_surfaces() const { return surfaces_.size() ; }
        index_t nb_regions() const { return regions_.size() ; }
        index_t nb_layers() const { return layers_.size() ; }
        index_t nb_surface_part_inside() const ;
        index_t nb_contact_part_inside() const ;
        index_t nb_real_corners_inside() const ;
        const std::string& gocad_name() const { return gocad_name_ ; }

        const vec3& point( index_t p ) const { return points_[p] ; }
        const Corner& corner( signed_index_t index ) const { return corners_[index] ; }
        const ContactPart& contact_part( signed_index_t index ) const { return contact_parts_[index] ; }
        const SurfacePart& surface_part( signed_index_t index ) const { return surface_parts_[index] ; }
        const BoundaryModelElement& contact( signed_index_t index ) const { return contacts_[index] ; }
        const BoundaryModelElement& surface( signed_index_t index ) const { return surfaces_[index] ; }
        const BoundaryModelElement& region( signed_index_t index ) const { return regions_[index] ; }
        const BoundaryModelElement& layer( signed_index_t index ) const { return layers_[index] ; }
        bool is_triangulated_model() const {
            for( index_t s = 0; s < nb_surface_parts(); s++ ) {
                if( !surface_parts_[s].is_triangulated() ) return false ;
            }
            return true ;
        }
        bool is_resolution_set() const {
            for( index_t s = 0; s < nb_surface_parts(); s++ ) {
                if( !surface_parts_[s].is_resolution_set() ) return false ;
            }
            return true ;
        }
        bool is_U_set() const {
            for( index_t s = 0; s < nb_surface_parts(); s++ ) {
                if( !surface_parts_[s].is_U_set() ) return false ;
            }
            return true ;
        }
        bool is_V_set() const {
            for( index_t s = 0; s < nb_surface_parts(); s++ ) {
                if( !surface_parts_[s].is_V_set() ) return false ;
            }
            return true ;
        }
        bool is_W_set() const {
            for( index_t s = 0; s < nb_surface_parts(); s++ ) {
                if( !surface_parts_[s].is_W_set() ) return false ;
            }
            return true ;
        }
        bool is_UVW_set() const {
            for( index_t s = 0; s < nb_surface_parts(); s++ ) {
                if( !surface_parts_[s].is_UVW_set() ) return false ;
            }
            return true ;
        }

        static void print_type( std::ostream& out, GEOL_FEATURE t, signed_index_t dim ) ;
        static void save_type( std::ostream& out, GEOL_FEATURE t ) ;
        static GEOL_FEATURE determine_geological_type( const std::string& in ) ;
        static GEOL_FEATURE determine_type( const std::vector< GEOL_FEATURE >& types ) ;

    protected:
        bool load_gocad_model3d( std::istream& in ) ;
        bool save_gocad_model3d( std::ostream& out ) ;

        bool check_model3d_compatibility() ;

    protected:
        std::string gocad_name_ ;

        // Basic elements
        std::vector< vec3 >                 points_ ;
        std::vector< Corner >               corners_ ;     
        std::vector< ContactPart >          contact_parts_ ;
        std::vector< SurfacePart >          surface_parts_ ; 
        std::vector< BoundaryModelElement > regions_ ;
        
        BoundaryModelElement universe_ ;

        /// \todo Borders and In Boundary of Contact and Surfaces are not properly built

        // Geological interfaces and contact between them
        std::vector< BoundaryModelElement >  contacts_ ;
        std::vector< BoundaryModelElement >  surfaces_ ;
        std::vector< BoundaryModelElement >  layers_ ;
    } ;    

    class GRGMESH_API BoundaryModelBuilder {
    public:
        BoundaryModelBuilder( BoundaryModel* model )
            : model_( model )
        {
        }
        void load_file( std::istream& in ) ;
        bool rebuild() ;
        void copy_macro_topology( const BoundaryModel* from ) ;
        void reserve_nb_corners( index_t size ) {
            model_->corners_.reserve( size ) ;
        }
        void reserve_nb_contact_parts( index_t size ) {
            model_->contact_parts_.reserve( size ) ;
        }
        void reserve_nb_surface_parts( index_t size ) {
            model_->surface_parts_.reserve( size ) ;
        }
        void reserve_nb_surfaces( index_t size ) {
            model_->surfaces_.reserve( size ) ;
        }
        void reserve_nb_contacts( index_t size ) {
            model_->contacts_.reserve( size ) ;
        }
        void reserve_nb_regions( index_t size ) {
            model_->regions_.reserve( size ) ;
        }
        void set_universe( const std::vector< std::pair< int, bool > >& boundaries ) ;

        void add_point( float64* point ) {
            add_point( vec3( point ) ) ;
        }
        void add_point( const vec3& point ) {
            model_->points_.push_back( point ) ;
        }

        void set_corner( index_t id, const float64* p ) {
            set_corner( id, vec3( p ) ) ;
        }
        void set_corner( index_t id, const vec3& p ) {
            grgmesh_debug_assert( id < model_->nb_corners() ) ;
            model_->corners_[id].set_corner( model_->nb_points() ) ;
            add_point( p ) ;
        }
        void add_corner_boundary( index_t id, index_t b ) {
            model_->corners_[id].add_boundary( b ) ;
        }
        void add_corner_in_boundary( index_t id, index_t b ) {
            model_->corners_[id].add_in_boundary( b ) ;
        }

        index_t create_contact_part(
            signed_index_t id = -1,
            const std::vector< index_t >& points = empty_uint_vector )
        {
            if( id == -1 ) id = model_->nb_contact_parts() ;
            grgmesh_debug_assert( id == model_->nb_contact_parts() ) ;
            model_->contact_parts_.push_back( ContactPart( model_, id, points ) ) ;
            return id ;
        }
        void set_contact_part( index_t id, const std::vector< vec3 >& vertices ) {
            grgmesh_debug_assert( id < model_->nb_contact_parts() ) ;
            for( index_t p = 0; p < vertices.size(); p++ ) {
                model_->contact_parts_[id].vertices_.push_back( model_->nb_points() ) ;
                add_point( vertices[p] ) ;
            }
        }
        void add_contact_part_boundary( index_t id, index_t b ) {
            model_->contact_parts_[id].add_boundary( b ) ;
        }
        void add_contact_part_in_boundary( index_t id, index_t b ) {
            model_->contact_parts_[id].add_in_boundary( b ) ;
        }

        index_t create_surface_part(
            signed_index_t id = -1,
            signed_index_t parent = -1,
            GEOL_FEATURE type = default_type )
        {
            if( id == -1 ) id = model_->nb_surface_parts() ;
            grgmesh_debug_assert( id == model_->nb_surface_parts() ) ;
            model_->surface_parts_.push_back( SurfacePart( model_, id, parent, type ) ) ;
            return id ;
        }
        void create_surface_part(
            const std::string& interface_name,
            const std::string& type,
            const KeyFacet& key ) ;

        void add_surface_part_boundary( index_t id, index_t b ) {
            model_->surface_parts_[id].add_boundary( b ) ;
        }
        void add_surface_part_in_boundary( index_t id, index_t b ) {
            model_->surface_parts_[id].add_in_boundary( b ) ;
        }
        bool set_surface_part_points_and_facets(
            index_t sp_id,
            const std::vector< index_t >& points,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr,
            bool compute_adjacent = true )
        {
            return model_->surface_parts_[sp_id].set_points_and_facets( points,
                facets, facet_ptr, compute_adjacent ) ;
        }
        void set_surface_part_adjacent_facets(
            index_t sp_id,
            const std::vector< signed_index_t >& in )
        {
            model_->surface_parts_[sp_id].set_adjacent_facets( in ) ;
        }
        void end_surface( index_t id ) {
            SurfacePart& surface = model_->surface_parts_[id] ;
            if( surface.nb_simplices() > 0 )
                surface.set_first_triangle_as_key() ;
            surface.compute_is_triangulated() ;
        }

        index_t create_surface(
            const std::string& name,
            signed_index_t id = -1,
            GEOL_FEATURE type = default_type )
        {
            if( id == -1 ) id = model_->nb_surfaces() ;
            grgmesh_debug_assert( id == model_->nb_surfaces() ) ;
            model_->surfaces_.push_back( BoundaryModelElement( model_, 2, id, -1, type ) ) ;
            model_->surfaces_[id].set_name( name ) ;
            return id ;
        }
        void add_surface_child( index_t id, index_t child ) {
            model_->surfaces_[id].add_child( child ) ;
        }
        void add_surface_boundary( index_t id, index_t b ) {
            model_->surfaces_[id].add_boundary( b ) ;
        }
        void reserve_surface_points( index_t id, index_t nb ) {
            model_->surface_parts_[id].points_.reserve( nb ) ;
        }

        void add_contact_child( index_t id, index_t child ) {
            model_->contacts_[id].add_child( child ) ;
        }
        void add_contact_boundary( index_t id, index_t b ) {
            model_->contacts_[id].add_boundary( b ) ;
        }
        void add_contact_in_boundary( index_t id, index_t b ) {
            model_->contacts_[id].add_in_boundary( b ) ;
        }

        index_t create_region( signed_index_t id = -1 ) {
            if( id == -1 ) id = model_->regions_.size() ;
            grgmesh_debug_assert( id == model_->regions_.size() ) ;
            model_->regions_.push_back( BoundaryModelElement( model_, 3, id ) ) ;
            return id ;
        }
        index_t create_region(
            const std::string& name,
            const std::vector< std::pair< int, bool > >& boundaries,
            signed_index_t id = -1 ) ;
        void add_region_oriented_boundary(
            index_t id,
            index_t b,
            bool side )
        {
            model_->regions_[id].add_oriented_boundary( b, side ) ;
        }
        void add_region_in_boundary( index_t id, index_t b ) {
            model_->regions_[id].add_in_boundary( b ) ;
        }
        void remove_univers_from_regions( index_t id ) {
            for( signed_index_t i = 0; i < model_->nb_regions(); ++i ) {
                signed_index_t cur_id = model_->region(i).id() ;
                grgmesh_assert( i == cur_id ) ;
                if( i > id ) model_->regions_[i].set_id( cur_id-1 ) ;
            }
            model_->regions_.erase( model_->regions_.begin() + id ) ;
        }

        index_t create_layer( const std::string& name, signed_index_t id = -1 ) {
            if( id == -1 ) id = model_->layers_.size() ;
            grgmesh_debug_assert( id == model_->layers_.size() ) ;
            model_->layers_.push_back( BoundaryModelElement( model_, 3, id ) ) ;
            model_->layers_[id].set_name( name ) ;
            return id ;
        }
        void add_layer_child( index_t id, index_t child ) {
            model_->layers_[id].add_child( child ) ;
        }

        /* Functions to help building the BoundaryModel */

        signed_index_t find_or_create_corner( index_t index ) ;
        signed_index_t find_or_create_contact_part(
            signed_index_t corner0,
            signed_index_t corner1,
            std::vector< vec3 >& points ) ;
        signed_index_t find_or_create_contact(
            std::vector< signed_index_t >& interfaces,
            GEOL_FEATURE type ) ;

        void build_contact_parts( const std::vector< Border >& borders) ;
        void build_contacts() ;
        void end_contacts() ;
        void end_contact_parts() ;
        void end_surfaces() ;
        void end_surface_parts() ;
        void end_surface_parts( const std::vector< signed_index_t >& change_orientation ) ;
        void end_corners() ;
        void end_layers() ;
        void end_model() ;
        void update_all_ids() ;
    protected:
        BoundaryModel* model_ ;
    } ;
}
#endif
