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

#include <grgmesh/boundary_model_element.h>
#include <grgmesh/boundary_model.h>
#include <grgmesh/utils.h>
#include <set>
#include <stack>
#include <fstream>

#define _USE_MATH_DEFINES // Otherwise M_PI not defined interfere with other math stuff ? Jeanne
#include <math.h>

namespace GRGMesh {

    // See http://www.geometrictools.com/LibMathematics/Distance/Distance.html
    float64 point_triangle_squared_distance(
           const vec3& point,
           const vec3& V0,
           const vec3& V1,
           const vec3& V2,
           vec3& closest_point = dummy
       ) {
           vec3 diff = V0 - point;
           vec3 edge0 = V1 - V0;
           vec3 edge1 = V2 - V0;
           float64 a00 = length2(edge0) ;
           float64 a01 = dot(edge0, edge1) ;
           float64 a11 = length2(edge1) ;
           float64 b0 = dot(diff, edge0) ;
           float64 b1 = dot(diff, edge1) ;
           float64 c  = length2(diff) ;
           float64 det = ::fabs(a00*a11 - a01*a01);
           float64 s = a01*b1 - a11*b0;
           float64 t = a01*b0 - a00*b1;
           float64 sqrDistance;

           if (s + t <= det) {
           if (s < 0.0) {
               if (t < 0.0)  { // region 4
               if (b0 < 0.0) {
                   t = 0.0;
                   if (-b0 >= a00) {
                   s = 1.0;
                   sqrDistance = a00 + 2.0*b0 + c;
                   } else {
                   s = -b0/a00;
                   sqrDistance = b0*s + c;
                   }
               } else {
                   s = 0.0;
                   if (b1 >= 0.0) {
                   t = 0.0;
                   sqrDistance = c;
                   } else if (-b1 >= a11) {
                   t = 1.0;
                   sqrDistance = a11 + 2.0*b1 + c;
                   } else {
                   t = -b1/a11;
                   sqrDistance = b1*t + c;
                   }
               }
               } else  { // region 3
               s = 0.0;
               if (b1 >= 0.0) {
                   t = 0.0;
                   sqrDistance = c;
               } else if (-b1 >= a11) {
                   t = 1.0;
                   sqrDistance = a11 + 2.0*b1 + c;
               } else {
                   t = -b1/a11;
                   sqrDistance = b1*t + c;
               }
               }
           } else if (t < 0.0) { // region 5
               t = 0.0;
               if (b0 >= 0.0) {
               s = 0.0;
               sqrDistance = c;
               } else if (-b0 >= a00) {
               s = 1.0;
               sqrDistance = a00 + 2.0*b0 + c;
               } else {
               s = -b0/a00;
               sqrDistance = b0*s + c;
               }
           } else  { // region 0
               // minimum at interior point
               float64 invDet = float64(1.0)/det;
               s *= invDet;
               t *= invDet;
               sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
               t*(a01*s + a11*t + 2.0*b1) + c;
           }
           } else {
               float64 tmp0, tmp1, numer, denom;

           if (s < 0.0)  { // region 2
               tmp0 = a01 + b0;
               tmp1 = a11 + b1;
               if (tmp1 > tmp0) {
               numer = tmp1 - tmp0;
               denom = a00 - 2.0*a01 + a11;
               if (numer >= denom) {
                   s = 1.0;
                   t = 0.0;
                   sqrDistance = a00 + 2.0*b0 + c;
               } else {
                   s = numer/denom;
                   t = 1.0 - s;
                   sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
                   t*(a01*s + a11*t + 2.0*b1) + c;
               }
               } else {
               s = 0.0;
               if (tmp1 <= 0.0) {
                   t = 1.0;
                   sqrDistance = a11 + 2.0*b1 + c;
               }
               else if (b1 >= 0.0) {
                   t = 0.0;
                   sqrDistance = c;
               } else {
                   t = -b1/a11;
                   sqrDistance = b1*t + c;
               }
               }
           } else if (t < 0.0) { // region 6
               tmp0 = a01 + b1;
               tmp1 = a00 + b0;
               if (tmp1 > tmp0) {
               numer = tmp1 - tmp0;
               denom = a00 - 2.0*a01 + a11;
               if (numer >= denom) {
                   t = 1.0;
                   s = 0.0;
                   sqrDistance = a11 + 2.0*b1 + c;
               } else {
                   t = numer/denom;
                   s = 1.0 - t;
                   sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
                   t*(a01*s + a11*t + 2.0*b1) + c;
               }
               } else {
               t = 0.0;
               if (tmp1 <= 0.0) {
                   s = 1.0;
                   sqrDistance = a00 + 2.0*b0 + c;
               } else if (b0 >= 0.0) {
                   s = 0.0;
                   sqrDistance = c;
               } else {
                   s = -b0/a00;
                   sqrDistance = b0*s + c;
               }
               }
           } else { // region 1
               numer = a11 + b1 - a01 - b0;
               if (numer <= 0.0) {
               s = 0.0;
               t = 1.0;
               sqrDistance = a11 + 2.0*b1 + c;
               } else {
               denom = a00 - 2.0*a01 + a11;
               if (numer >= denom) {
                   s = 1.0;
                   t = 0.0;
                   sqrDistance = a00 + 2.0*b0 + c;
               } else {
                   s = numer/denom;
                   t = 1.0 - s;
                   sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
                   t*(a01*s + a11*t + 2.0*b1) + c;
               }
               }
           }
           }

           // Account for numerical round-off error.
           if (sqrDistance < 0.0) {
           sqrDistance = 0.0;
           }

           closest_point = V0 + s*edge0 + t*edge1 ;
           return sqrDistance;
       }

/********************************************************************************************/
/******             BoundaryModelElement implementation     ***********************************/
/********************************************************************************************/

    const BoundaryModelElement* BoundaryModelElement::parent() const
    {
        if( !has_parent() ) return nil ;
        switch( dim() ) {
            case 1:
                return &model_->contact( parent_ ) ;
            case 2:
                return &model_->surface( parent_ ) ;
            default:
                grgmesh_assert_not_reached ;
                return nil ;
        }
    }

    const BoundaryModelElement* BoundaryModelElement::boundary( uint32 x ) const
    {
        if( x >= nb_boundaries() ) return nil ;
        uint32 id = boundaries_[x] ;
        switch( dim() ) {
            case 1:
                return &model_->corner( id ) ;
            case 2:
                if( has_parent() ) return &model_->contact_part( id ) ;
                return &model_->contact( id ) ;
            case 3:
                return &model_->surface_part( id ) ;
            default:
                grgmesh_assert_not_reached ;
                return nil ;
        }
    }

    const BoundaryModelElement* BoundaryModelElement::in_boundary( uint32 x ) const
    {
        if( x >= nb_in_boundary() ) return nil ;
        uint32 id = in_boundary_[x] ;
        switch( dim() ) {
            case 0:
                return &model_->contact_part( id ) ;
            case 1:
                if( has_parent() ) return &model_->surface_part( id ) ;
                return &model_->surface( id ) ;
            case 2:
                return &model_->region( id ) ;
            default:
                grgmesh_assert_not_reached ;
                return nil ;
        }
    }

    const BoundaryModelElement* BoundaryModelElement::child( uint32 x ) const
    {
        if( has_parent() || x >= nb_children() ) return nil ;
        uint32 id = children_[x] ;
        switch( dim() ) {
            case 1:
                return &model_->contact_part( id ) ;
            case 2:
                return &model_->surface_part( id ) ;
            case 3:
                return &model_->region( id ) ;
            default:
                grgmesh_assert_not_reached ;
                return nil ;
        }
    }

    void BoundaryModelElement::copy_macro_topology(
        const BoundaryModelElement& rhs,
        BoundaryModel& model )
    {
        model_ = &model ;
        name_ = rhs.name_ ;
        id_ = rhs.id_ ;
        dim_ = rhs.dim_ ;
        type_ = rhs.type_ ;
        parent_ = rhs.parent_ ;
        boundaries_ = rhs.boundaries_ ;
        sides_ = rhs.sides_ ;
        in_boundary_ = rhs.in_boundary_ ;
        children_ = rhs.children_ ;
    }

    /*! Return the neighbors of the element, that is the one that share one of
     *  its boundaries. Only the boudaries of the given type are considered.
     *  If exclude_voi is true, the negighbors belonging to the VOI are ignored
     */
    std::vector< const BoundaryModelElement* > BoundaryModelElement::compute_neighbors(
        GEOL_FEATURE through,
        bool exclude_voi
    ) const {
        std::vector< const BoundaryModelElement* > result ;
     
        for( int32 i = 0; i < boundaries_.size(); ++i ){
            const BoundaryModelElement* b = boundary( i ) ;
            if( through != ALL && b->type() != through ) continue ;

            for( int32 j = 0; j < b->nb_in_boundary(); ++j ){
                const BoundaryModelElement* bb = b->in_boundary( j ) ;
                if( bb == this ) continue ;
                if( exclude_voi && bb->is_on_voi() ) continue ;
                else result.push_back( bb ) ;
            }            
        }
        std::sort( result.begin(), result.end() ) ;
        int32 size = std::unique(result.begin(), result.end()) - result.begin() ;
        result.resize( size ) ;

        return result ;
    }

    float64 BoundaryModelElement::size() const {
        float64 result = 0. ;
        // If this element has children sum up their sizes
        for( uint32 i = 0; i < nb_children(); ++i ){
            result += child(i)->size() ;
        }

        if( result == 0 ){
            if( dim_ == 3 ) {
                // Compute the volume if this is a region
                for( uint32 i = 0; i < nb_boundaries(); i++ ) {
                    const SurfacePart* surface = dynamic_cast< const SurfacePart* >( boundary( i ) ) ;
                    for( uint32 t = 0; t < surface->nb_simplices(); t++ ) {
                        float64 cur_volume = ( dot( surface->point( t, 0 ),
                            cross( surface->point( t, 1 ), surface->point( t, 2 ) ) ) )
                            / static_cast< float64 >( 6 ) ;
                        sides_[i] ? result -= cur_volume : result += cur_volume ;
                    }
                }
            }
        }
        return fabs( result ) ;
    }

    void BoundaryModelElement::change_boundary_side( int32 id ) {
        for( uint32 i = 0; i < boundaries_.size(); ++i ){
            if( boundaries_[i] == id ){
                sides_[i] = !sides_[i] ;
            }
        }
    }

    float64 BoundaryModelElement::distance( const vec3& p ) const {
        if( nb_children() == 0 ) return big_float64 ;
        else {
            float64 result = big_float64 ;
            for( uint32 i = 0; i < nb_children(); ++i ){
                result = std::min( result, child( i )->distance(p) ) ;
            }
            return result ;
        }
    }
    float64 BoundaryModelElement::distance( BoundaryModelElement* e ) const {
        if( nb_children() == 0 ) return big_float64 ;
        else {
            float64 result = big_float64 ;
            for( uint32 i = 0; i < nb_children(); ++i ){
                result = std::min( result, child( i )->distance(e) ) ;
            }
            return result ;
        }
    }

    float64 BoundaryModelElement::min_angle( BoundaryModelElement* e ) const {
        if( nb_children() == 0 ) return 999 ;
        if( dim_ == 1 || dim_ == 2 ) return 999 ;

        float64 result = 180. ;
        for( uint32 i=0; i < nb_children(); ++i ){
            result = std::min( result, child( i )->min_angle( e ) ) ;
        }
        return result ;        
    }

    uint32 BoundaryModelElement::nb_simplices() const {
        uint32 result = 0 ;
         for( uint32 i=0; i < nb_children(); ++i ){
            result += child( i )->nb_simplices() ;
        }
        return result ;
    }
  
    /*! Mean orientation
     *  For surfaces it is the average Normal
     */
    vec3 BoundaryModelElement::average_orientation() const {
        if( nb_children() == 0 ) return vec3(-99999, -99999, -99999) ;
        if( dim_ == 0 || dim_ == 3 ) return vec3(-99999, -99999, -99999) ;
    
        vec3 result (0,0,0) ;
        float64 total_size = 0 ;
        for( uint32 i=0; i < nb_children(); ++i ){
            float64 s = child( i )->size() ;
            result += s * child( i )->average_orientation() ;
            total_size += s ;
        }
        if(total_size > 10e-30 ) return result/total_size ;
        else return vec3(0.,0.,0.) ;
    }

    /*! Probably not correct 
     */
    float64 BoundaryModelElement::average_angle_to_y() const {
        vec3 n = average_orientation() ;
        if( n.x >= 0 ) return std::acos( n.y / sqrt( n.x*n.x + n.y*n.y ) )*180./M_PI ;
        else return M_PI - std::acos( n.y / sqrt( n.x*n.x + n.y*n.y ) )*180./M_PI ;
    }

    /*! Probably not correct 
     */
    float64 BoundaryModelElement::average_dip() const {
        vec3 n = average_orientation() ;
        float64 phi = std::acos( n.z )*180. / M_PI ;
        if( phi <= 90. ) return 90.-phi ;
        else return phi-90. ;
    }


    /*! What information goes out in the CSV file in which
     *  some complexity will be computed
     */
    void BoundaryModelElement::print_categories( std::ostream& out ) {
        out << "Dimension" << SEP 
            << "Id" << SEP
            << "Name" << SEP
            << "Parent" << SEP
            << "Type" << SEP
            << "Nb boundaries"<< SEP
            << "Nb in boundary of" << SEP 
            << "Number of neighbors" << SEP
            << "Nb children" << SEP
            << "Size" << SEP

            << "Boundary ids" << SEP
            << "In boundary of ids" << SEP
            << "Neighbor ids" << SEP
            << "Children ids " << SEP
            
            << std::endl ;
    }
    
    /* To synchronize with the categories = columns */
    void BoundaryModelElement::print( std::ostream& out ) const {    
    
        out << dim_ << SEP
            << id_ << SEP
            << name_ << SEP
            << parent_ << SEP ;
        
        BoundaryModel::print_type( out, type_, dim_ ) ;

        out << boundaries_.size() << SEP
            << in_boundary_.size() << SEP ;

        std::vector< const BoundaryModelElement* > n = compute_neighbors( ALL, false ) ;
        out << n.size() << SEP ;

        if( children_.size() == 0 ) out << "-1" << SEP ;
        else out << children_.size() << SEP ;
    
        out << size() << SEP ;
        
        print_vector_contents( boundaries_, out ) ; 
        print_vector_contents( in_boundary_, out ) ; 
        print_vector_contents( n, out ) ; 
        print_vector_contents( children_, out ) ; 

        out << std::endl ;
    }


    struct CompPair {
         bool operator()( const std::pair< float64, float64 >& r, const std::pair< float64, float64 >& l ) const {
            return r.first < l.first ;
        }
    } ;

    /* To get some statistics on a vector of pairs
     * Not very efficient.
     */
    void print_stats( std::ostream& out, std::vector< std::pair< float64, float64> >& values,
        float64 min1 = -big_float64, float64 min2 = -big_float64, float64 min3 = -big_float64,
        float64 max1 =  big_float64, float64 max2 =  big_float64, float64 max3 =  big_float64
    ) { 
        if( values.size() == 0 ){
            out << "No values " << SEP ;
            out << "" << SEP
                << "" << SEP
                << "" << SEP  ;            
            if( min1 != -big_float64 ) out << "" << SEP ;
            if( min2 != -big_float64 ) out << "" << SEP ;
            if( min3 != -big_float64 ) out << "" << SEP ;
            if( max1 !=  big_float64 ) out << "" << SEP ;
            if( max2 !=  big_float64 ) out << "" << SEP ;
            if( max3 !=  big_float64 ) out << "" << SEP ;
            return ;
        }
        std::sort( values.begin(), values.end() ) ;

        float64 W = 0 ;
        float64 V = 0 ;

        float64 q1 = 0. ;
        float64 q2 = 0. ;
        float64 q3 = 0. ;

        float64 r1 = 0. ;
        float64 r2 = 0. ;
        float64 r3 = 0. ;
        
        for( uint32 i = 0 ; i < values.size();  ++i ){
            V += values[i].second * values[i].first ;
            W += values[i].second ;

            // Some of the tests are useless if the default values are used
            if( values[i].first < min1 ) q1 += values[i].second ;
            if( values[i].first < min2 ) q2 += values[i].second ;
            if( values[i].first < min3 ) q3 += values[i].second ;

            if( values[i].first > max1 ) r1 += values[i].second ;
            if( values[i].first > max2 ) r2 += values[i].second ;
            if( values[i].first > max3 ) r3 += values[i].second ;
        }
        if( W < 10e-30 ) {
            out << "Nil weights" << SEP << std::endl ;
            out << "" << SEP
                << "" << SEP
                << "" << SEP ;             
            if( min1 != -big_float64 ) out << "" << SEP ;
            if( min2 != -big_float64 ) out << "" << SEP ;
            if( min3 != -big_float64 ) out << "" << SEP ;
            if( max1 !=  big_float64 ) out << "" << SEP ;
            if( max2 !=  big_float64 ) out << "" << SEP ;
            if( max3 !=  big_float64 ) out << "" << SEP ;
            return ;
        }

        float64 mean = V/W ;

        float64 SD = 0 ;
        for( uint32 i = 0 ; i < values.size();  ++i ){
            SD += values[i].second * (values[i].first-mean)*(values[i].first-mean) ;
        }
        SD = sqrt( SD / W ) ;           

        out << values[0].first     << SEP
            << values.back().first << SEP 
            << mean                << SEP
            << SD                  << SEP ;
        if( min1 != -big_float64 )
            out << q1              << SEP ;
        if( min2 != -big_float64 )
            out << q2            << SEP ;
        if( min3 != -big_float64 )
            out << q3            << SEP ;       

        if( max1 !=  big_float64 )
            out << r1            << SEP ;
        if( max2 !=  big_float64 )
            out << r2            << SEP ;
        if( max3 !=  big_float64 )
            out << r3            << SEP ;

        out << W                 << SEP ;

    }
   
    void BoundaryModelElement::print_complexity_categories( std::ostream& out ) {
        out << "ELEMENT"           << SEP 
            << "ID"                << SEP
            << "NAME"              << SEP
            << "BORDER ELEMENTS"   << SEP
            << "INCIDENT ELEMENTS" << SEP
            << "NEIGHBORS"         << SEP 
            << "SIZE"              << SEP
            << "PROJECTED SIZE"    << SEP
            << "ASPECT RATIO"      << SEP
            << "AZIMUT"            << SEP
            << "DIP"               << SEP ;
            // Distance statistics
        out << "MIN DIST"          << SEP
            << "MAX DIST"          << SEP 
            << "W. AV. DIST"       << SEP
            << "W. STD DEV DIST"   << SEP
            << " < 1m"            << SEP 
            << " < 10m"           << SEP 
            << " < 100m"          << SEP 
            << "Total size"      << SEP ;
                 
            // Angle statistics
        out << "MIN ANGLE"         << SEP
            << "MAX ANGLE"         << SEP
            << "W. AV. ANGLE"       << SEP
            << "W. STD DEV ANGLE"   << SEP
            << " < 10 deg"        << SEP 
            << " > 170 deg"       << SEP
            << "Total size"      << SEP
            << std::endl ;
    }
    
    /*! Print the complexity for a region 
     *  \todo consider other possible types 
     */
    void BoundaryModelElement::print_complexity( std::ostream& out ) const {

        out << "Region" << SEP 
            << id_ << SEP
            << name_ << SEP
            << nb_boundary_elements() << SEP
            << "" << SEP
            << compute_neighbors(ALL).size() << SEP 
            << "" << SEP
            << "" << SEP
            << "" << SEP
            << "" << SEP
            << "" << SEP ;

        std::vector< std::pair< float64, float64 > > values ;
  
        compute_distances( values ) ;
        print_stats( out, values, 1., 10., 100. ) ;

        compute_angles( values ) ;
        print_stats( out, values, 10., -big_float64, -big_float64, 170. ) ;
           
        out << std::endl ;
    }

    /*! Returns true if this element or one of the element containing it
     *  is on the Volume Of Interest
     *  This info is strored in the type of the element
     */
    bool BoundaryModelElement::is_on_voi() const
    {
        if( type_ == ALL ) {
            for( int32 j = 0; j < nb_in_boundary(); ++j ) {
                GEOL_FEATURE t = in_boundary( j )->type() ;
                if( t == VOI || t == STRATI_VOI || t == FAULT_VOI ) return true ;
            }
        } else if( type_ == VOI || type_ == STRATI_VOI || type_ == FAULT_VOI )
            return true ;

        return false ;
    }

    void get_all_boundary_elements(
        const BoundaryModelElement* e,
        std::vector< const BoundaryModelElement* >& result )
    {
        result.clear() ;
        std::set< const BoundaryModelElement* > b1 ;
        for( int32 i = 0; i < e->nb_boundaries(); ++i ) {
            const BoundaryModelElement* b = e->boundary( i ) ;
            result.push_back( b ) ;
            for( uint32 j = 0; j < b->nb_boundaries(); b++ ) {
                b1.insert( b->boundary( j ) ) ;
            }
        }

        std::set< const BoundaryModelElement* > b2 ;
        for( std::set< const BoundaryModelElement* >::const_iterator it( b1.begin() );
            it != b1.end(); ++it ) {
            for( uint32 i = 0; i < ( *it )->nb_boundaries(); i++ ) {
                b2.insert( ( *it )->boundary( i ) ) ;
            }
        }        
        result.insert( result.end(), b1.begin(), b1.end() ) ;
        result.insert( result.end(), b2.begin(), b2.end() ) ;        
    }


    
    /*! Do the 2 elements share a boundary ? 
     */
    bool are_connected(
        const BoundaryModelElement* e1,
        const BoundaryModelElement* e2,
        bool check_all )
    {

        if( !check_all ) {
            // Check only if the elements share a boundary
            for( int32 i = 0; i < e1->nb_boundaries(); ++i ) {
                const BoundaryModelElement* b1 = e1->boundary( i ) ;
                for( int32 j = 0; j < e2->nb_boundaries(); ++j ) {
                    if( e2->boundary( i ) == b1 ) return true ;
                }
            }
        }
        else {
            std::vector< const BoundaryModelElement* > b1 ;
            std::vector< const BoundaryModelElement* > b2 ;

            get_all_boundary_elements( e1, b1 ) ;
            get_all_boundary_elements( e2, b2 ) ;
           
            for( int32 i = 0; i < b1.size(); ++i ){
                if( std::count( b2.begin(), b2.end(), b1[i] ) > 0 ) return true ;
            }
        }
        return false ;
    }

    /*! Approximated distances between boundaries 
     *  Distances are measured from the barycenters of the triangles/segments/ points 
     *  to the triangles/segments/points on other boundaries
     *     
     *  For volumes distances are only measured from the horizons to the other surfaces (because smapling is correct)
     *  For surfaces distances are measured from contacts that are not on the VOI
     * 
     *  Distance between intersecting surface tends toward 0.
     */
    void BoundaryModelElement::compute_distances( std::vector< std::pair< float64, float64 > >& values ) const {
    
        if( dim_ == 0 || dim_ == 1 ) return ;

        // If the element has no or only one boundary get out
        if( boundaries_.size() < 2 ) {
            values.clear() ;
            return ;
        }

        // Allocate room for the values vector
        int32 values_size = 0 ;
        std::vector< int32 > nb_simplex ;
        std::vector< const BoundaryModelElement* > from ;
        for( int32 i = 0; i < nb_boundaries(); ++i ){
            const BoundaryModelElement* b = boundary( i ) ;

            // Skip surfaces that are not horizons( unconf. )  when 
            // the element is a region
            if( dim_ == 3 && b->type() != STRATI ) continue ;
            if( dim_ == 2 && b->is_on_voi() ) continue ;

            from.push_back( b ) ;            
            nb_simplex.push_back( b->nb_simplices() ) ;
            values_size += nb_simplex.back() ;
        }
     
        values.resize( values_size ) ;

        int32 count = 0 ;

        std::vector< const BoundaryModelElement* > to ;
        for( int32 i = 0; i < from.size(); ++i ){
            const BoundaryModelElement* b = from[i] ;
            
            to.clear() ;
        
            for( int32 j = 0; j < nb_boundaries(); ++j ){
                const BoundaryModelElement* b1 = boundary( j ) ;

                if( b == b1
                    || (b->is_on_voi() && b1->is_on_voi())
                    || ( dim_ == 3 && b->parent() != nil && b->parent() == b1->parent() )
                  ) continue ;
                else to.push_back( b1 ) ;            
            }

            for( int32 j = 0; j < nb_simplex[i]; ++j ) {
                float64 min = big_float64 ;
                for( int32 k = 0; k < to.size(); ++k ) {
                    min = std::min( min, b->distance(j, to[k]) );
                }
                values[count+j].first = min ;
                values[count+j].second = b->simplex_size( j ) ;
            }

            count += nb_simplex[i] ;
        }
    }

    /*! Fills values with angles measured between the boundaries
     *  A weight (length of a segment or fixed at 1 for corner) is associated to the
     *  angle which is in degree between 0 and 180 
     * 
     *  Angles between two surfaces that have the same parent are ignored
     *  
     *  Returns 999 if the computation fails
     */
    void BoundaryModelElement::compute_angles( std::vector< std::pair< float64, float64 > >& values ) const {

        // It is rather difficult to get the size of the values vector 
        // here. Should not be too big, and greatly inferior than when used for distances just before

        // If the element has no or only one boundary get out
        if( boundaries_.size() < 2 ) {
            values.clear() ;
            return ;
        }

        values.clear() ;
        for( int32 i = 0; i < nb_boundaries(); ++i ){
            const BoundaryModelElement* b = boundary( i ) ;
            for( int32 j = i + 1; j < nb_boundaries(); ++j ) {
                const BoundaryModelElement* b1 = boundary( j ) ;
                bool same_side = false ;
                if( dim_ == 3 && sides_[i] == sides_[j] ) same_side = true ;
                if( are_connected( b, b1, false )
                    && !( dim_ == 3 && b->parent() != nil
                        && b->parent() == b1->parent() ) ) {
                    // Add the values for this contact
                    b->angles( b1, values, same_side ) ;
                }
            }
        }                      
    }

    int32 BoundaryModelElement::nb_boundary_elements( bool exclude_voi ) const
    {
        int32 result = 0 ;
        std::set< const BoundaryModelElement* > b2 ;
        for( int32 i = 0; i < nb_boundaries(); ++i ){
            const BoundaryModelElement* b = boundary( i ) ;
            if( exclude_voi && b->is_on_voi() ) continue ;
            if( b->dim() == 0 ){
                // Remove false corners should be removed on closed lines 
                const Corner* c = dynamic_cast< const Corner* >( b ) ;
                if( !c->is_real() ) continue ;
            }
            // Otherwise add an element and check its boundaries
            result++ ;

            for( int32 j = 0; j < b->nb_boundaries(); ++j ) {
                const BoundaryModelElement* bb = b->boundary( j ) ;
                if( exclude_voi && bb->is_on_voi() ) continue ;
                if( bb->dim() == 0 ){
                    const Corner* c = dynamic_cast< const Corner* >( bb ) ;
                    if( !c->is_real() ) continue ;
                }

                b2.insert( bb ) ;
            }
        }
        result += b2.size() ;

        std::set< const BoundaryModelElement* > b3 ;
        for( std::set< const BoundaryModelElement* >::iterator it( b2.begin() );
            it != b2.end(); ++it ) {
            for( int32 j = 0; j < ( *it )->nb_boundaries(); ++j ) {
                const BoundaryModelElement* bb = ( *it )->boundary( j ) ;
                if( exclude_voi && bb->is_on_voi() ) continue ;
                if( bb->dim() == 0 ) {
                    const Corner* c = dynamic_cast< const Corner* >( bb ) ;
                    if( !c->is_real() ) continue ;
                }
           
                else b3.insert( bb ) ;
            }          
        }
        result += b3.size() ;
        return result ;
    }

    int32 BoundaryModelElement::nb_incident_elements( bool exclude_voi ) const
    {
        int32 result = 0 ;
               
        std::set< const BoundaryModelElement* > b2 ;
        for( int32 i = 0; i < nb_in_boundary(); ++i ) {
            const BoundaryModelElement* b = in_boundary( i ) ;
            if( exclude_voi && b->is_on_voi() ) continue ;
            else {
                result++ ;
                for( int32 j = 0; j < b->nb_in_boundary(); ++j ) {
                    const BoundaryModelElement* bb = b->in_boundary( j ) ;
                    if( exclude_voi && bb->is_on_voi() ) continue ;
                    else b2.insert( bb ) ;
                }
            }
        }            
        result += b2.size() ;

        std::set< const BoundaryModelElement* > b3 ;
        for( std::set< const BoundaryModelElement* >::iterator it( b2.begin() );
            it != b2.end(); ++it ) {
            for( int32 j = 0; j < ( *it )->nb_in_boundary(); ++j ) {
                const BoundaryModelElement* bb = ( *it )->in_boundary( j ) ;
                if( exclude_voi && bb->is_on_voi() ) continue ;
                else b3.insert( bb ) ;
            }           
        }
        result += b3.size() ;
        return result ;
    }

/***********************************************************************************************/

    const vec3& Corner::vertex( uint32 p ) const
    {
        return model_->point( p_ ) ;
    }

    void Corner::copy_macro_topology(
        const Corner& rhs,
        BoundaryModel& model )
    {
        BoundaryModelElement::copy_macro_topology( rhs, model ) ;
    }
    
    void Corner::print_complexity( std::ostream& out ) const {
        if( !is_real() ) return ;

        out << "Corner" << SEP 
            <<  id_ << SEP
            << "" << SEP
            << nb_boundary_elements() << SEP
            << nb_incident_elements() << SEP

            // All other fields are empty
            << std::endl ;
    }
                

/********************************************************************************************/
/******             ContactPart implementation            ***********************************/
/********************************************************************************************/

    ContactPart::ContactPart( BoundaryModel* model, int32 id ):
        BoundaryModelElement( model, 1, id )
    { 
        boundaries_.resize( 2, nil) ; 
    }

    ContactPart::ContactPart(
        BoundaryModel* model,
        int32 id,
        const std::vector< uint32 >& points )
        : BoundaryModelElement( model, 1, id ), vertices_( points )
    {
    }

    ContactPart::ContactPart(
        BoundaryModel* model,
        int32 id,
        uint32 corner0,
        uint32 corner1,
        const std::vector< uint32 >& points
    ):  BoundaryModelElement( model, 1, id ),
        vertices_( points )
    {
        boundaries_.push_back( corner0 ) ;
        boundaries_.push_back( corner1 ) ;
    } ;

    const vec3& ContactPart::vertex( uint32 p ) const
    {
        return model_->point( vertices_[p] ) ;
    }

    void ContactPart::copy_macro_topology(
        const ContactPart& rhs,
        BoundaryModel& model )
    {
        BoundaryModelElement::copy_macro_topology( rhs, model ) ;
        is_inside_border_ = rhs.is_inside_border_ ;
    }
    void ContactPart::add_in_boundary( uint32 e )
    {
        for( uint32 i = 0; i < nb_in_boundary(); i++ ) {
            if( in_boundary_[i] == e ) {
                grgmesh_debug_assert( !is_inside_border_[i] ) ;
                is_inside_border_[i] = true ;
                return ;
            }
        }
        BoundaryModelElement::add_in_boundary( e ) ;
        is_inside_border_.push_back( false ) ;
    }

    bool ContactPart::contains( const vec3& p ) const {
        return find( p ) != -1 ;
    }
    int32 ContactPart::find( const vec3& p ) const
    {
        for( uint32 i = 0; i < vertices_.size(); ++i ) {
            if( vertex( i ) == p ) return i ;
        }
        return -1 ;
    }

    float64 ContactPart::size() const {
        float64 result = 0. ;
        for( uint32 i = 1; i < vertices_.size(); ++i ){
            result += length( vertex( i )-vertex( i-1 ) ) ;
        }
        if( is_closed() ) result += length( model_->point( vertices_.back() )-vertex( 0 ) );
        return result ;
    }

    float64 ContactPart::simplex_size( int32 i ) const {
        if( i < vertices_.size()-1 ) return length(vertex( i+1 ) - vertex( i)) ;
        else {
            grgmesh_debug_assert( i < vertices_.size() ) ;
            return length( model_->point( vertices_.back() )-vertex( 0 ) );
        }
    }


    /* Min distance betweeen the given point and the segments of this contact part
     */
    float64 ContactPart::distance( const vec3& p ) const {
        float64 result = big_float64 ;
        for( uint32 i = 1; i < vertices_.size(); ++i ){
            // Distance betweena a point and a segment COPY from smwh else
            const vec3& p0 = vertex( i-1 ) ;
            const vec3& p1 = vertex( i ) ;

            float64 distance_pt_2_segment  = big_float64 ;
            vec3 c = (p1-p0)/2 ;        
            float64 half = length( c - p1 ) ;
            float64 cp_dot_p0p1 = dot( p-c, p1-p0 ) ;

            if( cp_dot_p0p1 < -half ) distance_pt_2_segment =  length( p - p0 ) ;
            else if( cp_dot_p0p1 > half ) distance_pt_2_segment = length( p - p1 ) ;
            else {
                vec3 projection = c + cp_dot_p0p1*(p1-p0) ;
                distance_pt_2_segment = length( p - projection ) ;
            }    
            result = std::min( distance_pt_2_segment, result ) ;
        }
        if( is_closed() ) {
            // COPY BEUUURKKK !!
            const vec3& p0 = model_->point( vertices_.back() ) ;
            const vec3& p1 = vertex( 0 ) ;

            float64 distance_pt_2_segment  = big_float64 ;
            vec3 c = (p1-p0)/2 ;        
            float64 half = length( c - p1 ) ;
            float64 cp_dot_p0p1 = dot( p-c, p1-p0 ) ;

            if( cp_dot_p0p1 < -half ) distance_pt_2_segment =  length( p - p0 ) ;
            else if( cp_dot_p0p1 > half ) distance_pt_2_segment = length( p - p1 ) ;
            else {
                vec3 projection = c + cp_dot_p0p1*(p1-p0) ;
                distance_pt_2_segment = length( p - projection ) ;
            }    
            result = std::min( distance_pt_2_segment, result ) ;
        }
        return result ;
    }
    vec3 ContactPart::average_orientation() const {
        vec3 s(0., 0., 0. ) ;
        for( uint32 i = 1; i < vertices_.size(); ++i ){
            s += vertex( i-1 ) - vertex( i );
        }
        if( is_closed() ) {
            s += vertex( 0 )-model_->point( vertices_.back() ) ;
        }
        return normalize( s ) ;
    }


    /*! Minimal distance betweeen the vertices of this element to
     *  the triangle of e.
     *  Totally inefficient but not the priority right now
     */
    float64 ContactPart::distance( BoundaryModelElement* e ) const {
        float64 result = big_float64 ;
        for( uint32 i = 0; i < vertices_.size(); ++i ) {
            result = std::min( result, e->distance( vertex( i ) ) ) ;
        }
        return result ;
    }

    /*! Min distance fron the i-th simplex barycenter to the given elements
     */
    float64 ContactPart::distance( int32 s, BoundaryModelElement* to ) const {
        vec3 centroid ;
        if( s < vertices_.size()-1 ) {
            centroid = ( vertex( s+1 )+vertex( s )) * 0.5 ;
        }
        else {
            grgmesh_debug_assert( s == vertices_.size()-1 ) ;
            centroid = (vertex( 0 )+model_->point( vertices_.back())) * 0.5 ;
        }        
        return to->distance( centroid ) ;        
    }

    /*! Add the measured angle with a weight at 1. for the corners at
     *  which this contact intersect "with"
     */
    void ContactPart::angles( 
        BoundaryModelElement* in, 
        std::vector< std::pair< float64, float64 > >& values,
        bool same_side
    ) const {
        ContactPart* cp = dynamic_cast< ContactPart* >( in ) ;
        if( cp == nil ) return ;


        std::vector< const BoundaryModelElement* > shared ;
        for( int32 i = 0; i < nb_boundaries(); ++i ){
            for( int32 j = 0; j < cp->nb_boundaries(); ++j ){
                if( cp->boundaries_[j] == boundaries_[i] ) {
                    shared.push_back( cp->boundary( j ) ) ;
                }
            }
        }

        for( int32 i = 0; i < shared.size(); ++i ){
            const Corner* c = dynamic_cast< const Corner* >( shared[i] ) ;
            grgmesh_debug_assert( c != nil ) ;

            const vec3& p = c->vertex() ; 

            vec3 e1 = vertex( 1 ) - vertex( 0 ) ;
            if( p != vertex( 0 ) ){
                grgmesh_debug_assert( p == model_->point( vertices_.back() ) ) ;
                e1 =  vertex( vertices_.size()-2 )- model_->point( vertices_.back() ) ;
            }
            vec3 e2 = cp->vertex( 1 ) - cp->vertex( 0 ) ;
            if( p != cp->vertex( 0 ) ) {
                grgmesh_debug_assert( p == cp->vertex( cp->nb_vertices()-1 ) ) ;
                e2 = cp->vertex( cp->nb_vertices()-2 ) - cp->vertex( cp->nb_vertices()-1 ) ;
            }
            e1 = normalize( e1 ) ;
            e2 = normalize( e2 ) ;

            float64 a = std::acos( dot(e1,e2) ) * 180. / M_PI ;
            values.push_back( std::pair< float64, float64 >( a, 1. ) ) ;
        }

    }

    float64 ContactPart::min_angle( BoundaryModelElement* in ) const {
        float64 result = 999. ;

        ContactPart* cp = dynamic_cast< ContactPart* >( in ) ;
        if( cp != nil ) {

            std::vector< const BoundaryModelElement* > shared ;
            for( int32 i = 0; i < nb_boundaries(); ++i ){
                for( int32 j = 0; j < cp->nb_boundaries(); ++j ){
                    if( cp->boundaries_[j] == boundaries_[i] ) {
                        shared.push_back( cp->boundary( j ) ) ;
                    }
                }
            }
            
            for( int32 i = 0; i < shared.size(); ++i ){
                const Corner* c = dynamic_cast< const Corner* >( shared[i] ) ;
                grgmesh_debug_assert( c != nil ) ;

                const vec3& p = c->vertex() ; 

                vec3 e1 = vertex( 1 ) - vertex( 0 ) ;
                if( p != vertex( 0 ) ){
                    grgmesh_debug_assert( p == model_->point( vertices_.back() ) ) ;
                    e1 =  vertex( vertices_.size()-2 )- model_->point( vertices_.back() ) ;
                }
                vec3 e2 = cp->vertex( 1 ) - cp->vertex( 0 ) ;
                if( p != cp->vertex( 0 ) ) {
                    grgmesh_debug_assert( p == cp->vertex( cp->nb_vertices()-1 ) ) ;
                    e2 = cp->vertex( cp->nb_vertices()-2 ) - cp->vertex( cp->nb_vertices()-1 ) ;
                }
                e1 = normalize( e1 ) ;
                e2 = normalize( e2 ) ;
                
                float64 a = std::acos( dot(e1,e2) ) * 180. / M_PI ;
                //if( dot(e1,e2) < 0 ) a = std::acos( dot(-e1,e2) ) * 180. / Pi ;

                result = std::min( result, a ) ;
            }                
        }         
        return result ;     
    }

    void ContactPart::print_complexity( std::ostream& out ) const {
        out << "Contact"           << SEP 
            << id_                << SEP
            << ""              << SEP
            << nb_boundary_elements() << SEP
            << nb_incident_elements() << SEP
            << compute_neighbors(ALL).size() << SEP 
            << size()              << SEP ;

        float64 d = 0. ;
        if( !is_closed() ){
            d = length(vertex(0) - model_->point( vertices_.back())) ;
        }
        else {
            // Compute the max distance between 2 point divided by two
            for( uint32 i = 0; i < vertices_.size() ; ++i ) {
                for( uint32 j = 0; j < vertices_.size() ; ++j ) {
                    d = std::max( d, length2(vertex(i)- vertex(j) )) ;
                }
            }
            d = sqrt(d)/2. ;
        }
        float64 aspect = size() ;
        if( d > 10e-30 ) aspect /= d ;
        else aspect = big_float64 ;

        out << d     << SEP 
            << aspect << SEP               
            << average_angle_to_y() << SEP
            << average_dip()        << SEP


            << std::endl ;
    }

    void ContactPartMutator::set_point( uint32 id, const vec3& p ) {
        M_.model_->points_[M_.vertices_[id]] = p ;
    }

    vec3& ContactPartMutator::point( uint32 p ) const
    {
        return M_.model_->points_[ M_.vertices_[p] ] ;
    }

/********************************************************************************************/
/******             SurfacePart implementation            ***********************************/
/********************************************************************************************/

    const vec3& SurfacePart::point( int32 f, int32 v ) const
    {
        return model_->point( points_[facets_[facet_begin( f ) + v]] ) ;
    }
    const vec3& SurfacePart::vertex( uint32 v ) const
    {
        return model_->point( points_[v] ) ;
    }

    void SurfacePart::set_first_triangle_as_key()
    {
        key_facet_ = KeyFacet( model_->point( points_[facets_[0]] ),
            model_->point( points_[facets_[1]] ),
            model_->point( points_[facets_[2]] ) ) ;
    }

    int32 SurfacePart::adjcent_in_neighbor( uint32 f, uint32 e ) const
    {
        int32 adj = adjacent( f, e ) ;
        if( adj == -1 ) return -1 ;

        for( uint32 i = 0; i < nb_points_in_facet( adj ); i++ ) {
            if( adjacent( adj, i ) == f ) return i ;
        }
        return -1 ;
    }

    FacetEdge SurfacePart::next_on_border( const FacetEdge& te ) const
    {
        grgmesh_debug_assert( te.facet_ != -1 ) ;
        grgmesh_debug_assert( te.facet_ < nb_simplices() ) ;
        grgmesh_debug_assert( te.edge_ != -1 ) ;
        grgmesh_debug_assert( te.edge_ < 4 ) ;
        grgmesh_debug_assert( is_on_border( te.facet_, te.edge_ ) ) ;
        uint32 ref_index = point_index( te.facet_,
            edge_vertex( te.facet_, te.edge_, 1 ) ) ;
        std::set< int32 > used_triangles ;
        std::stack< int32 > S ;
        S.push( te.facet_ ) ;
        while( !S.empty() ) {
            int32 cur_t = S.top() ;
            S.pop() ;
            if( used_triangles.find( cur_t ) != used_triangles.end() ) { continue ; }
            used_triangles.insert( cur_t ) ;
            for( uint32 e = 0; e < 3 ; e++ ) {
                if( point_index( cur_t, edge_vertex( cur_t, e, 0 ) ) == ref_index ) {
                    if( is_on_border( cur_t, e ) ) {
                        return FacetEdge( cur_t, e ) ;
                    } else {
                        S.push( adjacent( cur_t, e ) ) ;
                        break ;
                    }
                }
            }
        }
        grgmesh_debug_assert( false ) ;
        return FacetEdge() ;
    }

    int32 SurfacePart::has_edge( int32 facet, int32 in0, int32 in1 ) const {
        int32 p0 = point_index( facet, 0 ) ;
        int32 p1 = point_index( facet, 1 ) ;
        int32 p2 = point_index( facet, 2 ) ;

        if( is_triangle( facet ) ) {

            if( in0 == p0 ) {
                if( in1 == p1 ) return 2 ;
                if( in1 == p2 ) return 1 ;
            }
            if( in0 == p1 ) {
                if( in1 == p0 ) return 2 ;
                if( in1 == p2 ) return 0 ;
            }
            if( in0 == p2 ) {
                if( in1 == p0 ) return 1 ;
                if( in1 == p1 ) return 0 ;
            }
        } else {
            int32 p3 = point_index( facet, 3 ) ;

            if( in0 == p0 ) {
                if( in1 == p1 ) return 0 ;
                if( in1 == p3 ) return 3 ;
            }
            if( in0 == p1 ) {
                if( in1 == p0 ) return 0 ;
                if( in1 == p2 ) return 1 ;
            }
            if( in0 == p2 ) {
                if( in1 == p1 ) return 1 ;
                if( in1 == p3 ) return 2 ;
            }
            if( in0 == p3 ) {
                if( in1 == p2 ) return 2 ;
                if( in1 == p0 ) return 3 ;
            }
        }
        return -1 ;               
    }

    int32 SurfacePart::find_triangle ( int32 in0, int32 in1 ) const {
        for( uint32 f = 0; f < nb_simplices(); ++f ) {
            if( has_edge( f, in0, in1 ) != -1 ) {
                return f ;
            }
        }
        return -1 ; 
    }

    bool SurfacePart::contains( const vec3& p ) const {
        return find( p ) != -1 ;
    }
    int32 SurfacePart::find( const vec3& p ) const
    {
        for( uint32 i = 0; i < points_.size(); ++i ) {
            if( model_->point( points_[i] ) == p ) return i ;
        }
        return -1 ;
    }
    /*! Returns the id of a triangle that has these two points 
     *  if any else returns -1 
     *  
     *  WARNING There might TWO such triangles
     */
    int32 SurfacePart::find_triangle( const vec3& p0, const vec3& p1 ) const {
        // There might be several points with the same coordinates
        // Test all possible pairs
        
        std::vector< int32 > i0 ;
        std::vector< int32 > i1 ;

        for( uint32 v = 0; v < points_.size(); ++v ) {
            if( model_->point( points_[v] ) == p0 ) i0.push_back( v ) ;
            if( model_->point( points_[v] ) == p1 ) i1.push_back( v ) ;
        }
        int32 t = -1 ;
        for( int32 i = 0; i < i0.size(); ++i ){
            for( int32 j = 0; j < i1.size(); ++j ){
                t = find_triangle( i0[i], i1[j] ) ; 
                if( t != -1 ) return t ;
            }
        }
        return -1 ;
    }

    /*! Fills the adjacent_ vector. There is one per triangle. 
     */
    void SurfacePart::compute_adjacent_facets() {
        adjacent_.resize( facets_.size(), -1 ) ;
        std::vector< int32 > facets ;
        facets.reserve( 6 ) ;
        std::vector< std::vector< int32 > > facet_points( nb_vertices(), facets ) ;

        for( uint32 f = 0; f < nb_simplices(); ++f ){
            for( uint32 v = 0; v < nb_points_in_facet( f ); v++ ) {
                facet_points[point_index( f, v )].push_back( f ) ;
            }
        }
        for( uint32 p = 0; p < nb_vertices(); ++p ){
            std::sort( facet_points[p].begin(), facet_points[p].end() ) ;
        }

        for( uint32 f = 0; f < nb_simplices(); ++f ){
            uint32 nb_edges = is_triangle( f ) ? 3 : 4 ;
            for( uint32 e = 0; e < nb_edges; ++e ){
                if( !is_on_border( f, e ) ) continue ;

                int32 v0 = point_index( f, edge_vertex( f, e, 0 ) ) ;
                int32 v1 = point_index( f, edge_vertex( f, e, 1 ) ) ;

                const std::vector< int32 >& facets0 = facet_points[v0] ;
                const std::vector< int32 >& facets1 = facet_points[v1] ;

                std::vector< int32 > inter(
                    std::min( facets0.size(), facets1.size() ) ) ;
                std::vector< int32 >::iterator end = std::set_intersection(
                    facets0.begin(), facets0.end(), facets1.begin(), facets1.end(),
                    inter.begin() ) ;
                int32 nb_intersections = end - inter.begin() ;
                if( nb_intersections == 2 ) {
                    int32 f2 = inter[0] == f ? inter[1] : inter[0] ;
                    int32 e2 = has_edge( f2, v0, v1 ) ;
                    adjacent( f, e ) = f2 ;
                    adjacent( f2, e2 ) = f ;
                } else {
                    grgmesh_debug_assert( nb_intersections == 1 ) ;
                }
            }
        }
    }

    /*! Change the triangles if their orientations do not
     * match the one of the key facet
     */
    bool SurfacePart::key_facet_orientation() {
        
        if( key_facet_.is_default() ) {
            set_first_triangle_as_key() ;
            return true ;
        }
        
        vec3& p0 = key_facet_.p0_ ;
        vec3& p1 = key_facet_.p1_ ;
        vec3& p2 = key_facet_.p2_ ;
        int32 t = find_triangle( p0,p1,p2 ) ;
        if( t == -1 ) {
            // It is because of the sign of Z that is not the same 
            p0.z *= -1 ;
            p1.z *= -1 ;
            p2.z *= -1 ;
            t = find_triangle( p0, p1, p2 ) ; 
        }
        grgmesh_debug_assert( t > -1 ) ;

        int32 i0 = point_id( t, p0 ) ;
        int32 i1 = point_id( t, p1 ) ;
        int32 i2 = point_id( t, p2 ) ;

        bool same_sign = true ;
        if( i0 == 0 && i1 == 2 ) same_sign = false ;      
        else if( i0 == 1 && i1 == 0 ) same_sign = false ;
        else if( i0 == 2 && i1 == 1 ) same_sign = false ;

        return same_sign ;
   }

    int32 SurfacePart::point_id( int32 t, int32 p0 ) const {
        for( uint32 v = 0; v < nb_points_in_facet(t); v++ ) {
            if( same_point( point_index( t, v ), p0) ) return v ;
        }
        return -1 ;
    }

    int32 SurfacePart::point_id( int32 t, const vec3& p ) const {
        for( uint32 v = 0; v < nb_points_in_facet(t); v++ ) {
            if( point( t, v ) == p ) return v ;
        }
        return -1 ;
    }

    // this is a copy
    struct comp_vec3bis {
        bool operator()( const vec3& l, const vec3& r ) const {
            if( l.x != r.x ) return l.x < r.x ;
            if( l.y != r.y ) return l.y < r.y ;
            return l.z < r.z ;
        }
    } ;

    int32 SurfacePart::find_triangle( const vec3& p0, const vec3& p1, const vec3& p2 ) const {
        for( uint32 t = 0; t < nb_simplices(); ++t ){
            const vec3& pp0 = point( t, 0 )   ;
            const vec3& pp1 = point( t, 1 )   ;
            const vec3& pp2 = point( t, 2 )   ;
            
            if( p0 == pp0 ) {
                if( p1 == pp1 && p2 == pp2 ) return t ;
                if( p1 == pp2 && p2 == pp1 ) return t ;
            }
            if( p0 == pp1 ) {
                if( p1 == pp0 && p2 == pp2 ) return t ;
                if( p1 == pp2 && p2 == pp0 ) return t ;
            }
            if( p0 == pp2 ) {
                if( p1 == pp0 && p2 == pp1 ) return t ;
                if( p1 == pp1 && p2 == pp0 ) return t ;
            }
        }
        return -1 ;
    }


    int32 SurfacePart::edge_id( int32 t, int32 p0, int32 p1 ) const {
        int32 t_0 = point_id( t, p0 ) ;
        int32 t_1 = point_id( t, p1 ) ;
       
        if( t_0 > t_1 ) { int32 tmp = t_0 ; t_0 = t_1 ; t_1 = tmp ; }
        
        if     ( t_0 == 0 && t_1 == 1 ) return 2 ;
        else if( t_0 == 0 && t_1 == 2 ) return 1 ;
        else if( t_0 == 1 && t_1 == 2 ) return 0 ;
        else return -1 ;
    }

    int32 SurfacePart::triangles_around_point(
        int32 shared_point,
        std::vector< int32 >& result,
        bool border_only ) const
    {
        result.resize(0) ;
        std::stack< int32 > S ;
        for( uint32 t = 0; t < nb_simplices(); ++t ) {
            for( uint32 v = 0; v < nb_points_in_facet(t); v++ ) {
                if( point_index( t, v ) == shared_point ) {
                    return triangles_around_point_with_hint( shared_point, result,
                        border_only, t ) ;
                }
            }
        }
        grgmesh_assert_not_reached ;
        return -1 ;
    }

    int32 SurfacePart::triangles_around_point_with_hint(
        int32 shared_point,
        std::vector< int32 >& result,
        bool border_only,
        int32 triangle_hint ) const
    {
        result.resize( 0 ) ;
        std::stack< int32 > S ;
        S.push( triangle_hint ) ;

        std::vector< int32 > visited ;
        visited.reserve( 20 ) ;
        do {
            int32 t = S.top() ;
            S.pop() ;
            if( vector_contains( visited, t ) ) continue ;
            visited.push_back( t ) ;
            for( uint32 v = 0; v < nb_points_in_facet(t); ++v ) {
                if( point_index( t, v ) == shared_point ) {

                    for( uint32 adj = 0; adj < 3; adj++ ) {
                        if( !is_on_border( t, adj ) ) {
                            S.push( adjacent( t, adj ) ) ;
                        }
                    }
                    if( border_only ) {
                        if( !is_on_border( t ) ) break ;
                        if( is_on_border( t, edge_vertex( t, v, 0 ) )
                            || is_on_border( t, edge_vertex( t, v, 1 ) ) ) {
                            result.push_back( t ) ;
                        }
                    } else {
                        result.push_back( t ) ;
                    }
                    break ;
                }
            }
        } while( !S.empty() ) ;

        return result.size() ;
    }

    vec3 SurfacePart::barycenter( int32 f ) const {
        vec3 barycenter ;
        for( uint32 i = 0; i < nb_points_in_facet( f ); i++ ) {
            barycenter += point( f, i ) ;
        }
        return barycenter / nb_points_in_facet( f ) ;
    }

    uint32 SurfacePart::closest_point_in_facet(
        uint32 f,
        const vec3& v ) const
    {
       uint32 result = 0 ;
       float64 dist = big_float64 ;
       for( uint32 p = 0; p < nb_points_in_facet( f ); p++ ) {
           float64 distance = length2( v - point( f, p ) ) ;
           if( dist > distance ) {
               dist = distance ;
               result = p ;
           }
       }
       return result ;
    }

    float64 SurfacePart::facet_resolution( uint32 f ) const
    {
        float64 result = 0.0 ;
        for( uint32 p = 0; p < nb_points_in_facet( f ); p++ ) {
            result += resolution( f, p )  ;
        }
        return result / (float64)nb_points_in_facet( f ) ;
    }

    float64 SurfacePart::size() const {
        float64 result = 0. ;
        for( uint32 i = 0; i < nb_simplices(); i++ ) {
            result += Utils::triangle_area( point( i, 0 ), point( i, 1 ),
                point( i, 2 ) ) ;
            if( !is_triangle( i ) ) {
                result += Utils::triangle_area( point( i, 0 ), point( i, 2 ),
                    point( i, 3 ) ) ;
            }
        }
        return result ;
    }

    float64 SurfacePart::simplex_size( int32 t ) const
    {
        grgmesh_debug_assert( t < nb_simplices() ) ;
        float64 result = Utils::triangle_area( point( t, 0 ), point( t, 1 ),
            point( t, 2 ) ) ;
        if( !is_triangle( t ) ) {
            result += Utils::triangle_area( point( t, 0 ), point( t, 2 ),
                point( t, 3 ) ) ;
        }
        return result ;
    }

    float64 SurfacePart::distance( const vec3& p ) const {
        float64 result = big_float64 ;

        for( uint32 i = 0; i < nb_simplices(); i++ ) {
            float64 cur_result = point_triangle_squared_distance( p, point( i, 0 ),
                    point( i, 1 ), point( i, 2 ) ) ;
            if( !is_triangle( i ) ) {
                cur_result += point_triangle_squared_distance( p, point( i, 0 ), point( i, 2 ),
                    point( i, 3 ) ) ;
                cur_result /= static_cast< float64 >( 2.0 ) ;
            }
            result = std::min( cur_result, result ) ;
        }
        if( result != big_float64 ) result = sqrt( result ) ;
        return result ;
    }
   
    /*! Minimal distance betweeen the vertices of this element to
     *  the triangle of e.
     *
     *  Totally inefficient but not the priority right now
     */
    float64 SurfacePart::distance( BoundaryModelElement* e ) const {
        float64 result = big_float64 ;

        for( uint32 i = 0; i < points_.size(); ++i ) {
            result = std::min( result, e->distance( model_->point( points_[i] ) ) ) ;
        }
        return result ;
    }

    float64 SurfacePart::distance( int32 t, BoundaryModelElement* to ) const {
        grgmesh_debug_assert( t < nb_simplices() ) ;
        if( is_triangle( t ) ) {
            return to->distance( point( t, 0 ) + point( t, 1 ) + point( t, 2 ) ) / 3. ;
        } else {
            return to->distance(
                point( t, 0 ) + point( t, 1 ) + point( t, 2 ) + point( t, 3 ) ) / 4. ;
        }
    }

    vec3 SurfacePart::facet_normal( int32 t ) const {
        const vec3& p0 = point( t, 0 )  ;
        const vec3& p1 = point( t, 1 )  ;
        const vec3& p2 = point( t, 2 )  ;
        vec3 c0 = cross(p0-p2, p1-p2) ;
        if( !is_triangle(t) ) {
            const vec3& p3 = point( t, 3 )  ;
            c0 += cross(p0-p3, p2-p3) ;
            c0 /= static_cast< float64 >( 2.0 ) ;
        }
        return normalize( c0 ) ;
    }

    vec3 SurfacePart::average_orientation() const {
        float64 total_a = 0 ;
        vec3 result (0,0,0) ;
        for( uint32 t = 0; t < nb_simplices() ; ++t){
            vec3 n = facet_normal( t ) ;
            float64 a = simplex_size( t ) ;

            result += a*n ;
            total_a += a ;
        }
        if( total_a > 10e-30 ) return result/total_a ;
        else return vec3(0.,0.,0.) ;
    }

    void SurfacePart::angles( 
        BoundaryModelElement* in , 
        std::vector< std::pair< float64, float64 > >& values,
        bool same_side
    ) const {
        if( in == this ) return ;
        SurfacePart* sp = dynamic_cast< SurfacePart* >( in ) ;
        if( sp == nil ) return ;

        std::vector< const BoundaryModelElement* > shared ;
        for( int32 i = 0; i < nb_boundaries(); ++i ) {
            for( int32 j = 0; j < sp->nb_boundaries(); ++j ) {
                if( sp->boundaries_[j] == boundaries_[i] ) {
                    shared.push_back( sp->boundary( j ) ) ;
                }
            }
        }

        for( int32 i = 0; i < shared.size(); ++i ){
            const ContactPart* cp = dynamic_cast< const ContactPart* >( shared[i] ) ;
            grgmesh_debug_assert( cp != nil ) ;

            // Find the triangles sharing each segment 
            for( uint32 j = 1; j < cp->nb_vertices(); ++j ){
                const vec3& p0 = cp->vertex( j-1 ) ;
                const vec3& p1 = cp->vertex( j ) ;

                int32 t1 = find_triangle( p0, p1 ) ;
                int32 t2 = sp->find_triangle( p0, p1 ) ;
                if( t1 == -1 || t2 == -1 ) {
                    // Should not happen
                    continue ;
                }
                // Get the angle between these triangles
                float64 d = dot(facet_normal(t1), sp->facet_normal(t2)) ;
           
                float64 a = std::acos( d ) * 180 / M_PI ;
                if( same_side ) a = 180. - a ;
                values.push_back( std::pair< float64, float64 >(a, length(p0-p1)) ) ;
            }
            // If the line is closed check the closing segment (Copy)
            if( cp->is_closed() ) {
                const vec3& p0 = cp->vertex( 0 ) ;
                const vec3& p1 = cp->vertex( cp->nb_vertices()-2 ) ;

                int32 t1 = find_triangle( p0, p1 ) ;
                int32 t2 = sp->find_triangle( p0, p1 ) ;
                if( t1 == -1 || t2 == -1 ) {
                    // Should not happen
                    continue ;
                }

                float64 d = dot(facet_normal(t1), sp->facet_normal(t2)) ;
                 
                float64 a = std::acos( d ) * 180 / M_PI ;
                if( same_side ) a = 180. - a ;
                values.push_back( std::pair< float64, float64 >( a, length(p0-p1) ) ) ;
            }
        }
    }

    /*! If e is a SurfacePart that share a boundary line with this SurfacePart
     *  returns the min angle between the two
     *  else return 999.
     */
    float64 SurfacePart::min_angle( BoundaryModelElement* in ) const {
        float64 result = 999. ;
        if( in == this ) return result ;

        SurfacePart* sp = dynamic_cast< SurfacePart* >( in ) ;
        if( sp != nil ) {
            // Find the contact shared by the two
            std::vector< const BoundaryModelElement* > shared ;
            for( int32 i = 0; i < nb_boundaries(); ++i ) {
                for( int32 j = 0; j < sp->nb_boundaries(); ++j ) {
                    if( sp->boundaries_[j] == boundaries_[i] ) {
                        shared.push_back( sp->boundary( j ) ) ;
                    }
                }
            }
            
            for( int32 i = 0; i < shared.size(); ++i ){
                const ContactPart* cp = dynamic_cast< const ContactPart* >( shared[i] ) ;
                grgmesh_debug_assert( cp != nil ) ;
                
                // Find the triangle sharing each segment 
                 for( uint32 j = 1; j < cp->nb_vertices(); ++j ){
                    const vec3& p0 = cp->vertex( j-1 ) ;
                    const vec3& p1 = cp->vertex( j ) ;

                    int32 t1 = find_triangle( p0, p1 ) ;
                    int32 t2 = sp->find_triangle( p0, p1 ) ;
                    if( t1 == -1 || t2 == -1 ) {
                        // Should not happen
                        continue ;
                    }
                    // Get the angle between these triangles
                    float64 d = dot(facet_normal(t1), sp->facet_normal(t2)) ;
                    //if( d < 0 ) d = dot(-triangle_normal(t1), sp->triangle_normal(t2)) ;                   
                    float64 a = std::acos( d ) * 180 / M_PI ;

                    result = std::min( result, a ) ;
                 }
                 // If the line is closed check the closing segment (Copy)
                 if( cp->is_closed() ) {
                    const vec3& p0 = cp->vertex( 0 ) ;
                    const vec3& p1 = cp->vertex( cp->nb_vertices()-2 ) ;

                    int32 t1 = find_triangle( p0, p1 ) ;
                    int32 t2 = sp->find_triangle( p0, p1 ) ;
                    if( t1 == -1 || t2 == -1 ) {
                        // Should not happen
                        continue ;
                    }

                    float64 d = dot(facet_normal(t1), sp->facet_normal(t2)) ;
                    //if( d < 0 ) d = dot(-triangle_normal(t1), sp->triangle_normal(t2)) ;                   
                    float64 a = std::acos( d ) * 180 / M_PI ;

                    result = std::min( result, a ) ;
                 }
            }                
        }
         
        return result ;
    }

    void SurfacePart::print_complexity( std::ostream& out ) const {
        out << "Surface"           << SEP 
            << id_                << SEP
            << parent()->name()              << SEP
            << nb_boundary_elements() << SEP
            << nb_incident_elements() << SEP
            << compute_neighbors(ALL).size() << SEP 
            << size()              << SEP 
            << ""    << SEP
            << ""      << SEP
            
            << average_angle_to_y() << SEP
            << average_dip()        << SEP ;
    

        std::vector< std::pair< float64, float64 > > values ;
  
        compute_distances( values ) ;
        print_stats( out, values, 1., 10., 100. ) ;

        compute_angles( values ) ;
        print_stats( out, values, 10., -big_float64, -big_float64, 170. ) ;
           
        out << std::endl ;

    }

    void SurfacePart::print_mesh( const std::string& filename ) const
    {
        std::ofstream file( filename.c_str(), std::ios::trunc | std::ios::out ) ;
        file << "Surface : " << id() << std::endl ;
        file << "========== Points =========" << std::endl ;
        for( uint32 p = 0; p < nb_vertices(); p++ ) {
            file << p << " -> " << vertex( p ) << std::endl ;
        }
        file << "========== Facets =========" << std::endl ;
        for( uint32 p = 0; p < nb_simplices(); p++ ) {
            file << p << " ->" ;
            for( uint32 v = 0; v < nb_points_in_facet(p); v++ ) {
                file << " " << point_index( p, v ) ;
            }
            file << std::endl ;
        }
        file << "========== Facet ptr =========" << std::endl ;
        for( uint32 p = 0; p < nb_simplices(); p++ ) {
            file << p << " -> " << facet_begin(p) << " " << facet_end(p) << std::endl ;
        }
        file << "========== Adjacents =========" << std::endl ;
        for( uint32 p = 0; p < nb_simplices(); p++ ) {
            file << p << " ->" ;
            for( uint32 v = 0; v < nb_points_in_facet(p); v++ ) {
                file << " " << adjacent( p, v ) ;
            }
            file << std::endl ;
        }

    }

    void SurfacePart::point_normal( std::vector< vec3 >& normals ) const
    {
        normals.resize( nb_vertices() ) ;
        for( uint32 f = 0; f < nb_simplices(); f++ ) {
            vec3 normal = facet_normal( f ) ;
            for( uint32 p = 0; p < nb_points_in_facet( f ); p++ ) {
                uint32 id = point_index( f, p ) ;
                normals[id] += normal ;
            }
        }
        for( uint32 p = 0; p < nb_vertices(); p++ ) {
            normals[p] = normalize( normals[p] ) ;
        }
    }

    Box3d SurfacePart::bbox() const
    {
        Box3d result ;
        for( uint32 p = 0; p < nb_vertices(); p++ ) {
            result.add_point( vertex( p ) ) ;
        }
        return result ;
    }

    void SurfacePartMutator::set_point( uint32 id, const vec3& p ) {
        M_.model_->points_[M_.points_[id]] = p ;
    }

    vec3& SurfacePartMutator::point( uint32 p ) const
    {
        return M_.model_->points_[ M_.points_[p] ] ;
    }
}

