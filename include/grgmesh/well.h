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
 

#ifndef __GRGMESH_WELL_
#define __GRGMESH_WELL_

#include <grgmesh/common.h>

namespace GRGMesh {
    class BoundaryModel ;
    class Edge ;
    class Well ;
}

namespace GRGMesh {

    class GRGMESH_API WellCorner {
    public:
        WellCorner( const vec3& point, signed_index_t id = -1 )
            : point_( point ), surface_id_( id ), resolution_( -1 )
        {
        }

        const vec3& point() const { return point_ ; }
        void set_surface_id( signed_index_t id ) { surface_id_ = id ; }
        signed_index_t surface_id() const { return surface_id_ ; }
        double resolution() const { return resolution_ ; }
        void set_resolution( double d ) { resolution_ = d ; }

    private:
        vec3 point_ ;
        signed_index_t surface_id_ ;
        double resolution_ ;
    } ;

    class GRGMESH_API WellPart {
    public:
        WellPart()
            : well_( nil ), id_( -1 )
        {
            corners_.reserve( 2 ) ;
        }
        void add_corner( signed_index_t c ) { corners_.push_back( c ) ; }
        signed_index_t corner( signed_index_t c ) const {
            grgmesh_debug_assert( c < corners_.size() ) ;
            return corners_[c] ;
        }
        void add_points( const std::vector< vec3 >& p ) { points_ = p ; }
        const vec3& point( index_t p ) const {
            grgmesh_debug_assert( p < points_.size() ) ;
            return points_[p] ;
        }
        index_t nb_points() const { return points_.size() ; }
        const std::vector< vec3 >& points() const { return points_ ; }
        double resolution( index_t p ) const {
            grgmesh_debug_assert( p < resolutions_.size() ) ;
            return resolutions_[p] ;
        }
        std::vector< double >& resolutions() { return resolutions_ ; }
        void set_well( Well* well ) { well_ = well ; }
        const Well* well() const { return well_ ; }
        void set_id( signed_index_t id ) { id_ = id ; }
        signed_index_t id() const { return id_ ; }

    private:
        Well* well_ ;
        signed_index_t id_ ;
        std::vector< vec3 > points_ ;
        std::vector< signed_index_t > corners_ ;
        std::vector< double > resolutions_ ;
    } ;

    class GRGMESH_API Well {
    public:
        Well() {}
        void copy_corners_and_informations( Well& well ) const ;
        void set_well_in_parts() {
            for( index_t p = 0; p < nb_parts(); p++ ) {
                parts_[p].set_well( this ) ;
            }
        }
        void add_part_points( index_t part, const std::vector< vec3 >& p ) {
            parts_[part].add_points( p ) ;
        }
        void get_part_edges( index_t p, std::vector< Edge >& edges ) const ;
        void get_region_edges( index_t p, std::vector< Edge >& edges ) const ;

        void add_part( const WellPart& part, index_t r ) {
            parts_.push_back( part ) ;
            part_region_id_.push_back( r ) ;
            parts_.back().set_id( parts_.size()-1 ) ;
        }
        signed_index_t find_or_create_corner( const vec3& p, signed_index_t id = -1 ) ;
        signed_index_t find_corner( const vec3& p ) const ;
        const WellCorner& corner( index_t c ) const {
            grgmesh_debug_assert( c < corners_.size() ) ;
            return corners_[c] ;
        }
        const WellPart& part( index_t part ) const {
            grgmesh_debug_assert( part < parts_.size() ) ;
            return parts_[part] ;
        }
        signed_index_t part_region_id( index_t part ) const {
            grgmesh_debug_assert( part < nb_parts() ) ;
            return part_region_id_[part] ;
        }
        double part_length( index_t part ) const ;
        void all_part_vertices(
            index_t part,
            std::vector< vec3 >& vertices ) const ;

        index_t nb_corners() const { return corners_.size() ; }
        index_t nb_parts() const { return parts_.size() ; }

        void set_name( const std::string& name ) { name_ = name ; }
        const std::string& name() const { return name_ ; }

    private:
        signed_index_t add_corner( const vec3& p, signed_index_t id ) {
            corners_.push_back( WellCorner( p, id ) ) ;
            return corners_.size()-1 ;
        }

    private:
        std::vector< WellCorner > corners_ ;
        std::vector< WellPart > parts_ ;
        std::vector< index_t > part_region_id_ ;
        std::string name_ ;
    } ;

    class GRGMESH_API WellGroup {
    public:
        WellGroup() ;
        virtual ~WellGroup() ;

        void get_region_edges(
            index_t region,
            std::vector< Edge >& edges ) const ;

        BoundaryModel* model() const { return model_ ; }
        void set_model( GRGMesh::BoundaryModel* model ) { model_ = model ; }
        bool is_well_already_added( const std::string& name ) const ;
        void add_well( const Well& w ) { wells_.push_back( w ) ; }

        index_t nb_wells() const { return wells_.size() ; }
        const Well& well( index_t w ) const { return wells_[w] ; }
        void update_wells() {
            for( index_t w = 0; w < nb_wells(); w++ ) {
                wells_[w].set_well_in_parts() ;
            }
        }

    protected:
        std::vector< Well > wells_ ;
        BoundaryModel* model_ ;
    } ;

}
#endif

