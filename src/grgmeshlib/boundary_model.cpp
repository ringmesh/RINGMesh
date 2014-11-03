/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advanceme:nt Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

/*! \author Jeanne Pellerin */

#include <grgmeshlib/boundary_model.h>
#include <grgmeshlib/utils.h>

#include <fstream>
#include <sstream>
#include <cmath>
#include <time.h>
#include <set>

namespace GRGMesh {

    /********************************************************************************************/
    /***********             BoundaryModel class implementation           ************************/
    /********************************************************************************************/

    BoundaryModel::BoundaryModel()
        : universe_( this )
    {
    }

    BoundaryModel::~BoundaryModel()
    {
    }

    Box3d BoundaryModel::bbox() const
    {
        Box3d result ;
        for( uint32 i = 0; i < surface_parts_.size(); ++i ) {
            for( uint32 j = 0; j < surface_parts_[i].nb_vertices(); ++j ) {
                result.add_point( surface_parts_[i].vertex( j ) ) ;
            }
        }
        return result ;
    }

    template< class T > void print_elements(
        const std::vector< T >& v,
        std::ofstream& out,
        bool skip_voi = true )
    {
        for( uint32 i = 0; i < v.size(); ++i ) {
            if( skip_voi && v[i].is_on_voi() )
                continue ;
            else
                v[i].print( out ) ;
        }
    }

    int32 BoundaryModel::find_region( int32 surface_part_id, bool side ) const
    {
        grgmesh_debug_assert( surface_part_id < nb_surface_parts() ) ;
        for( uint32 r = 0; r < nb_regions(); r++ ) {
            const BoundaryModelElement& cur_region = region( r ) ;
            for( uint32 s = 0; s < cur_region.nb_boundaries(); s++ ) {
                if( cur_region.side( s ) == side
                    && cur_region.boundary( s )->id() == surface_part_id ) {
                    return r ;
                }
            }
        }
        return -1 ;
    }

    int32 BoundaryModel::nb_real_corners() const
    {
        int32 result = 0 ;
        for( int32 i = 0; i < corners_.size(); ++i ) {
            if( corners_[i].is_real() ) ++result ;
        }
        return result ;
    }

    int32 BoundaryModel::nb_surface_with_free_boundary() const
    {
        int32 result = 0 ;
        for( uint32 i = 0; i < surface_parts_.size(); ++i ) {
            for( uint32 j = 0; j < surface_parts_[i].nb_boundaries(); ++j ) {
                if( surface_parts_[i].boundary( j )->nb_in_boundary() == 1 ) {
                    ++result ;
                    break ;
                }
            }
        }
        return result ;
    }

    int32 BoundaryModel::nb_non_manifold_contact_parts() const
    {
        int32 result = 0 ;
        for( uint32 i = 0; i < contact_parts_.size(); ++i ) {
            const ContactPart& cp = contact_parts_[i] ;
            if( cp.nb_in_boundary() > 1 ) ++result ;
        }
        return result ;
    }

    uint32 BoundaryModel::nb_surface_part_inside() const
    {
        uint32 result = 0 ;
        for( uint32 i = 0; i < surface_parts_.size(); ++i ) {
            if( surface_parts_[i].type() != VOI ) ++result ;
        }
        return result ;
    }

    uint32 BoundaryModel::nb_contact_part_inside() const
    {
        uint32 result = 0 ;
        for( uint32 i = 0; i < contact_parts_.size(); ++i ) {
            if( !contact_parts_[i].is_on_voi() ) ++result ;
        }
        return result ;
    }

    uint32 BoundaryModel::nb_real_corners_inside() const
    {
        uint32 result = 0 ;
        for( int32 i = 0; i < corners_.size(); ++i ) {
            if( corners_[i].is_real() && !corners_[i].is_on_voi() ) ++result ;
        }
        return result ;
    }

    void BoundaryModel::print_topology( std::ofstream& out ) const
    {

        out << "Topological information of model" << SEP << gocad_name_ << SEP
            << std::endl ;

        out << "Numbers of" << SEP << "Regions" << SEP << "Surfaces Parts Inside"
            << SEP << "Surface Parts Box" << SEP << "Blind Surface parts" << SEP
            << "Contact Parts Inside" << SEP << "Contact Parts Box" << SEP
            << "Corners Inside" << SEP << "Corners Box" << SEP << std::endl ;

        int32 s_in = nb_surface_part_inside() ;
        int32 b_in = nb_contact_part_inside() ;
        int32 c_in = nb_real_corners_inside() ;

        out << SEP << regions_.size() << SEP << s_in << SEP
            << surface_parts_.size() - s_in << SEP << nb_surface_with_free_boundary()
            << SEP << b_in << SEP << contact_parts_.size() - b_in << SEP << c_in
            << SEP << nb_real_corners() - c_in << SEP << std::endl ;

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "REGIONS" << std::endl ;
        print_elements( regions_, out ) ;

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "SURFACES" << std::endl ;
        print_elements( surfaces_, out ) ;

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "CONTACTS" << std::endl ;
        print_elements( contacts_, out ) ;

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "SURFACE PARTS" << std::endl ;
        print_elements( surface_parts_, out ) ;

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "CONTACT PARTS" << std::endl ;
        print_elements( contact_parts_, out ) ;

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "CORNERS" << std::endl ;
        print_elements( corners_, out ) ;

    }

    void BoundaryModel::print_element_info( std::ofstream& out ) const
    {
        out << gocad_name_ << std::endl ;

        BoundaryModelElement::print_complexity_categories( out ) ;

        for( uint32 i = 0; i < regions_.size(); ++i ) {
            regions_[i].print_complexity( out ) ;
        }

        for( uint32 i = 0; i < surface_parts_.size(); ++i ) {
            if( surface_parts_[i].is_on_voi() ) continue ;
            surface_parts_[i].print_complexity( out ) ;
        }
        for( uint32 i = 0; i < contact_parts_.size(); ++i ) {
            if( contact_parts_[i].is_on_voi() ) continue ;
            contact_parts_[i].print_complexity( out ) ;
        }
        for( uint32 i = 0; i < corners_.size(); ++i ) {
            if( corners_[i].is_on_voi() ) continue ;
            corners_[i].print_complexity( out ) ;
        }
    }

    int32 BoundaryModel::surface_id( const std::string& name ) const
    {
        for( uint32 i = 0; i < surfaces_.size(); ++i ) {
            if( surfaces_[i].name() == name ) {
                return i ;
            }
        }
        return -1 ;
    }

    int32 BoundaryModel::find_corner( const vec3& p ) const
    {
        for( int32 i = 0; i < corners_.size(); ++i ) {
            if( corners_[i].vertex() == p ) return i ;
        }
        return -1 ;
    }

    /*!
     * Return the contact part sampled by the segment [point0,point1].
     * If this is an interior segment default value return is -1
     * \author Gaetan FUSS
     * \date 2014
     */
    int32 BoundaryModel::find_contact_part( vec3 point0, vec3 point1 ) const
    {
        for( uint32 k = 0; k < contact_parts_.size(); ++k ) {
            // If the points are on the same contact_part k
            if( contact_parts_[k].contains( point0 )
                && contact_parts_[k].contains( point1 ) ) {
                return k ;
            }
        }
        // Else default value is -1 for interior segment
        return -1 ;
    }

    int32 BoundaryModel::find_contact_part(
        int32 corner0,
        int32 corner1,
        const std::vector< vec3 >& points ) const
    {

        for( uint32 i = 0; i < contact_parts_.size(); ++i ) {
            const ContactPart& cp = contact_parts_[i] ;

            if( corner0 == cp.boundary_id( 0 )
                &&    // same corners
                corner1 == cp.boundary_id( 1 )
                && cp.nb_vertices() == points.size() ) {         // same points
                bool equal = true ;
                for( uint32 p = 0; p < cp.nb_vertices(); p++ ) {
                    if( cp.vertex( p ) != points[p] ) {
                        equal = false ; break ;
                    }
                }
                if( equal ) return i ;
            }

            if( corner1 == cp.boundary_id( 0 )
                && corner0 == cp.boundary_id( 1 )
                && cp.nb_vertices() == points.size() ) {         // same points
                bool equal = true ;
                for( uint32 p = 0; p < cp.nb_vertices(); p++ ) {
                    if( cp.vertex( cp.nb_vertices() - p - 1 ) != points[p] ) {
                        equal = false ; break ;
                    }
                }
                if( equal ) return i ;
            }
        }
        return -1 ;
    }

    std::vector< BoundaryModelElement* > BoundaryModel::find_contacts(
        const BoundaryModelElement* corner0,
        const BoundaryModelElement* corner1 ) const
    {
        std::vector< BoundaryModelElement* > result ;
        for( uint32 i = 0; i < contact_parts_.size(); ++i ) {
            if( ( corner0 == contact_parts_[i].boundary( 0 )
                && corner1 == contact_parts_[i].boundary( 1 ) )
                || ( corner0 == contact_parts_[i].boundary( 1 )
                    && corner1 == contact_parts_[i].boundary( 0 ) ) ) {
                const BoundaryModelElement* cc = &contact_parts_[i] ;
                BoundaryModelElement* c = const_cast< BoundaryModelElement* >( cc ) ;
                result.push_back( c ) ;
            }
        }
        return result ;
    }

    std::vector< BoundaryModelElement* > BoundaryModel::get_closed_contacts() const
    {
        std::vector< BoundaryModelElement* > result ;
        for( uint32 i = 0; i < contact_parts_.size(); ++i ) {
            if( contact_parts_[i].is_closed() ) {
                const BoundaryModelElement* cc = &contact_parts_[i] ;
                BoundaryModelElement* c = const_cast< BoundaryModelElement* >( cc ) ;
                result.push_back( c ) ;
            }
        }
        return result ;
    }

    std::vector< BoundaryModelElement* > BoundaryModel::find_interface_parts(
        const std::vector< const BoundaryModelElement* >& contacts,
        bool consider_free_contacts ) const
    {

        std::set< const BoundaryModelElement* > set( contacts.begin(),
            contacts.end() ) ;

        std::vector< BoundaryModelElement* > exact_match_with_set ;
        std::vector< BoundaryModelElement* > contains_set ;

        for( uint32 i = 0; i < surface_parts_.size(); ++i ) {
            std::set< const BoundaryModelElement* > cur ;

            if( consider_free_contacts ) {
                for( uint32 j = 0; j < surface_parts_[i].nb_boundaries(); j++ ) {
                    cur.insert( surface_parts_[i].boundary( j ) ) ;
                }
            } else {
                for( uint32 j = 0; j < surface_parts_[i].nb_boundaries(); ++j ) {
                    const BoundaryModelElement* b = surface_parts_[i].boundary( j ) ;
                    if( b->nb_in_boundary() > 1 ) cur.insert( b ) ;
                }
            }
            if( set.size() == cur.size()
                && std::equal( set.begin(), set.end(), cur.begin() ) ) {
                const BoundaryModelElement* cc = &surface_parts_[i] ;
                BoundaryModelElement* c = const_cast< BoundaryModelElement* >( cc ) ;
                exact_match_with_set.push_back( c ) ;
            } else if( std::includes( cur.begin(), cur.end(), set.begin(),
                set.end() ) ) {
                const BoundaryModelElement* cc = &surface_parts_[i] ;
                BoundaryModelElement* c = const_cast< BoundaryModelElement* >( cc ) ;
                contains_set.push_back( c ) ;
            }
        }
        if( exact_match_with_set.size() > 0 ) {
            return exact_match_with_set ;
        } else
            return contains_set ;
    }
    int32 BoundaryModel::find_contact( const std::vector< int32 >& interfaces ) const
    {
        std::vector< const BoundaryModelElement* > comp( interfaces.size() ) ;
        for( uint32 i = 0; i < interfaces.size(); ++i ) {
            comp[i] = &surfaces_[interfaces[i]] ;
        }

        for( uint32 i = 0; i < contacts_.size(); ++i ) {
            if( comp.size() == contacts_[i].nb_in_boundary() ) {
                bool equal = true ;
                for( uint32 j = 0; j < contacts_[i].nb_in_boundary(); j++ ) {
                    if( comp[j] != contacts_[i].in_boundary( j ) ) {
                        equal = false ; break ;
                    }
                }
                if( equal ) {
                    return i ;
                }
            }
        }
        return -1 ;
    }

    void BoundaryModel::print_type( std::ostream& out, GEOL_FEATURE t, int32 dim )
    {
        switch( dim ) {
            case 0:
                out << "" ;
                break ;

            case 1:
                switch( t ) {
                    case STRATI:
                        //out << "Horizon-horizon contact" ;
                        out << "H-H" ;
                        break ;
                    case FAULT:
                        //out << "Fault-fault contact" ;
                        out << "F-F" ;
                        break ;
                    case VOI:
                        out << "B-B" ;
                        break ;
                    case STRATI_FAULT:
                        out << "F-H" ;
                        break ;
                    case STRATI_VOI:
                        out << "H-B" ;
                        break ;
                    case FAULT_VOI:
                        out << "F-B" ;
                        break ;
                    default:
                        out << "" ;
                        break ;
                }
                break ;

            case 2:
                switch( t ) {
                    case STRATI:
                        out << "Horizon" ;
                        break ;
                    case FAULT:
                        out << "Fault" ;
                        break ;
                    case VOI:
                        out << "Box" ;
                        break ;
                    default:
                        out << "" ;
                        break ;
                }
                break ;

            case 3:
                out << "" ;
                break ;

            default:
                out << "" ;
                break ;
        }
        out << SEP ;
    }

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

    GEOL_FEATURE BoundaryModel::determine_geological_type( const std::string& in )
    {
        if( in == "" ) return ALL ;
        if( in == "reverse_fault" ) return FAULT ;
        if( in == "normal_fault" ) return FAULT ;
        if( in == "fault" ) return FAULT ;
        if( in == "top" ) return STRATI ;
        if( in == "none" ) return STRATI ;
        if( in == "unconformity" ) return STRATI ;
        if( in == "boundary" ) return VOI ;

        std::cout << "Unexpected type in the model file " << in
            << std::endl ;
        return ALL ;
    }

    /*! Determine the type of a geological feature that is at the intersection
     *  of the input features
     */
    GEOL_FEATURE BoundaryModel::determine_type(
        const std::vector< GEOL_FEATURE >& types )
    {
        if( types.size() == 0 ) return ALL ;

        // Sort and remove duplicates form the in types
        std::vector< GEOL_FEATURE > in = types ;
        std::sort( in.begin(), in.end() ) ;
        int32 new_size = std::unique( in.begin(), in.end() ) - in.begin() ;
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
            // Other cases ? for corners ? what is the point ?
            return ALL ;
        }
        return ALL ;
    }

    /*! Load a .ml file, gocad model3d (b-rep file) 
     *  Works correctly only if there is a UNIQUE Model3d in the file
     * 
     *  No support for properties right now
     *  no advance check of coordinate system
     *
     *  Build all the elements of the structural model and their relationships
     */
    bool BoundaryModel::load_gocad_model3d( std::istream& in )
    {
        BoundaryModelBuilder builder( this ) ;
        builder.load_file( in ) ;
        return true ;
    }

    /*! Check that the BoundaryModel can be saved as a .ml 
     */
    bool BoundaryModel::check_model3d_compatibility()
    {
        // Check that the Surfaces exist
        if( nb_surfaces() == 0 && nb_surface_parts() > 0 ) {
            BoundaryModelBuilder builder( this ) ;
            // Creates one Surface per SurfacePart
            for( uint32 i = 0; i < surface_parts_.size(); ++i ) {
                // set name, type, links
                std::ostringstream name ;
                name << "surface_" << i ;
                uint32 id = builder.create_surface( name.str() ) ;
                builder.add_surface_child( id, i ) ;
            }

            // Set links surface parts toward surfaces
            for( uint32 i = 0; i < surfaces_.size(); ++i ) {
                surface_parts_[surfaces_[i].child( 0 )->id()].set_parent( i ) ;
            }

            // Is it really useful to have contacts, let's hope not... I am not doing it
        }

        for( uint32 i = 0; i < surface_parts_.size(); ++i ) {
            SurfacePart& sp = surface_parts_[i] ;
            if( sp.nb_vertices() == 0 ) continue ;
            if( sp.key_facet().is_default() ) {
                sp.set_first_triangle_as_key() ;
            }
        }

        // the universe should exist, perhaps create a function copying the lines 
        // 805 to 834 de s2_b_model.cpp
        if( universe_.name() != "Universe" ) {
            std::cout
                << "The region universe is not defined for the model. IMPLEMENTATION TO DO"
                << std::endl ;
            return false ;
        }

        // Check that each region has a name and valid surfaces
        for( uint32 i = 0; i < regions_.size(); ++i ) {
            BoundaryModelElement& region = regions_[i] ;

            if( region.name() == "" ) {
                std::ostringstream name ;
                name << "region_" << i ;
                region.set_name( name.str() ) ;
            }
            if( region.nb_boundaries() == 0 ) {
                std::cout << " The region " << region.name()
                    << " has no SurfaceParts on its boundary" << std::endl ;
                return false ;
            }
        }

        if( !is_triangulated_model() ) {
            std::cout << "The model contains quads " << std::endl ;
            return false ;
        }

        return true ;
    }

    void save_region(
        int32 count,
        const BoundaryModelElement& region,
        std::ostream& out )
    {
        out << "REGION " << count << "  " << region.name() << " " << std::endl ;
        int32 it = 0 ;

        for( int32 i = 0; i < region.nb_boundaries(); ++i ) {
            out << "  " ;
            if( region.side( i ) ) {
                out << "+" ;
            } else {
                out << "-" ;
            }
            out << region.boundary( i )->id() + 1 ;
            it++ ;
            if( it == 5 ) {
                out << std::endl ;
                it = 0 ;
            }
        }
        out << "  0" << std::endl ;
    }

    void save_layer(
        int32 count,
        int32 offset,
        const BoundaryModelElement& layer,
        std::ostream& out )
    {
        out << "LAYER " << layer.name() << " " << std::endl ;
        int32 it = 0 ;

        for( int32 i = 0; i < layer.nb_children(); ++i ) {
            out << "  " << layer.child_id( i ) + offset + 1 ;
            it++ ;
            if( it == 5 ) {
                out << std::endl ;
                it = 0 ;
            }
        }
        out << "  0" << std::endl ;
    }

    void save_coordinate_system( std::ostream& out )
    {
        out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl << "NAME Default"
            << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
            << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl << "ZPOSITIVE Elevation"
            << std::endl << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl ;
    }

    bool BoundaryModel::save_gocad_model3d( std::ostream& out )
    {
        out.precision( 16 ) ;
        if( !check_model3d_compatibility() ) {
            std::cout << "The BoundaryModel " << gocad_name_
                << " cannot be saved in .ml format " << std::endl ;
            return false ;
        }

        // Print Model3d headers
        out << "GOCAD Model3d 1" << std::endl << "HEADER {" << std::endl << "name:"
            << "model_from_graphite" << std::endl << "}" << std::endl ;

        save_coordinate_system( out ) ;

        // Print the TSurf = Surfaces info
        for( int32 i = 0; i < surfaces_.size(); ++i ) {
            out << "TSURF " << surfaces_[i].name() << std::endl ;
        }

        int32 count = 1 ;
        // Print the TFace info
        for( int32 i = 0; i < surface_parts_.size(); ++i ) {
            const SurfacePart& s = surface_parts_[i] ;
            out << "TFACE " << count << "  " ;
            save_type( out, s.type() ) ;
            out << " " << s.parent()->name() << std::endl ;

            const KeyFacet kf = s.key_facet() ;
            out << "  " << kf.p0_ << std::endl ;
            out << "  " << kf.p1_ << std::endl ;
            out << "  " << kf.p2_ << std::endl ;

            ++count ;
        }

        int32 offset_layer = count ;
        // Region info + universe
        save_region( count, universe_, out ) ;
        ++count ;

        for( int32 i = 0; i < regions_.size(); ++i ) {
            save_region( count, regions_[i], out ) ;
            ++count ;
        }

        for( int32 i = 0; i < layers_.size(); ++i ) {
            save_layer( count, offset_layer, layers_[i], out ) ;
            ++count ;
        }

        out << "END" << std::endl ;

        // Now save each one of the surfaces
        for( int32 i = 0; i < surfaces_.size(); ++i ) {
            const BoundaryModelElement& tsurf = surfaces_[i] ;

            // Header
            out << "GOCAD TSurf 1" << std::endl << "HEADER {" << std::endl << "name:"
                << tsurf.name() << std::endl << "name_in_model_list:" << tsurf.name()
                << std::endl << "}" << std::endl ;
            save_coordinate_system( out ) ;

            out << "GEOLOGICAL_FEATURE " << tsurf.name() << std::endl
                << "GEOLOGICAL_TYPE " ;
            save_type( out, tsurf.type() ) ;
            out << std::endl ;

            out << "PROPERTY_CLASS_HEADER Z {" << std::endl << "is_z:on" << std::endl
                << "}" << std::endl ;

            // On entre les TFACES
            uint32 vertex_count = 1 ;
            uint32 offset = vertex_count ;

            std::vector< uint32 > bstones ;
            std::vector< uint32 > next_point ;

            for( int32 j = 0; j < tsurf.nb_children(); ++j ) {
                offset = vertex_count ;

                const SurfacePart* sp = dynamic_cast< const SurfacePart* >( tsurf.child( j ) ) ;

                out << "TFACE" << std::endl ;
                for( int32 k = 0; k < sp->nb_vertices(); ++k ) {
                    out << "VRTX " << vertex_count << " " << sp->vertex( k ) << std::endl ;
                    vertex_count++ ;
                }

                for( int32 k = 0; k < sp->nb_simplices(); ++k ) {
                    out << "TRGL " << sp->point_index( k, 0 ) + offset << " "
                        << sp->point_index( k, 1 ) + offset << " "
                        << sp->point_index( k, 2 ) + offset << std::endl ;
                }

                // Info on corners and contacts

                for( int32 k = 0; k < sp->nb_boundaries(); ++k ) {
                    const ContactPart* cp = dynamic_cast< const ContactPart* >( sp->boundary( k ) ) ;

                    const vec3& c = cp->vertex( 0 ) ;
                    const vec3& next = cp->vertex( 1 ) ;

                    int32 t = sp->find_triangle( c, next ) ;
                    grgmesh_debug_assert( t != -1 ) ;

                    int32 i0 = sp->point_index( t, 0 ) ;
                    int32 i1 = sp->point_index( t, 1 ) ;
                    int32 i2 = sp->point_index( t, 2 ) ;

                    const vec3& p0 = sp->vertex( i0 ) ;
                    const vec3& p1 = sp->vertex( i1 ) ;
                    const vec3& p2 = sp->vertex( i2 ) ;

                    int32 c_id = -1 ;
                    int32 next_id = -1 ;

                    if( p0 == c )
                        c_id = i0 ;
                    else if( p0 == next ) next_id = i0 ;

                    if( p1 == c )
                        c_id = i1 ;
                    else if( p1 == next ) next_id = i1 ;

                    if( p2 == c )
                        c_id = i2 ;
                    else if( p2 == next ) next_id = i2 ;

                    grgmesh_debug_assert( c_id != -1 && next_id != -1 ) ;

                    bstones.push_back( c_id + offset ) ;
                    next_point.push_back( next_id + offset ) ;
                }
            }

            // Print corners and contact
            for( int32 j = 0; j < bstones.size(); ++j ) {
                out << "BSTONE " << bstones[j] << std::endl ;
            }
            for( int32 j = 0; j < bstones.size(); ++j ) {
                out << "BORDER " << vertex_count << " " << bstones[j] << " "
                    << next_point[j] << std::endl ;
                vertex_count++ ;
            }
            out << "END" << std::endl ;
        }
        return true ;
    }


    void BoundaryModelBuilder::copy_macro_topology( const BoundaryModel* from )
    {
        model_->corners_.resize( from->nb_corners(), Corner( model_ ) ) ;
        model_->contact_parts_.resize( from->nb_contact_parts(), ContactPart( model_ ) ) ;
        model_->surface_parts_.resize( from->nb_surface_parts(), SurfacePart( model_ ) ) ;
        model_->regions_.resize( from->nb_regions(), BoundaryModelElement( model_ ) ) ;
        model_->layers_.resize( from->nb_layers(), BoundaryModelElement( model_ ) ) ;
        model_->contacts_.resize( from->nb_contacts(), BoundaryModelElement( model_ ) ) ;
        model_->surfaces_.resize( from->nb_surfaces(), BoundaryModelElement( model_ ) ) ;
#pragma omp parallel for
        for( uint32 i = 0; i < model_->nb_corners(); i++ ) {
            model_->corners_[i].copy_macro_topology( from->corner( i ), *model_ ) ;
        }
#pragma omp parallel for
        for( uint32 i = 0; i < model_->nb_contact_parts(); i++ ) {
            model_->contact_parts_[i].copy_macro_topology( from->contact_part( i ),
                *model_ ) ;
        }
#pragma omp parallel for
        for( uint32 i = 0; i < model_->nb_surface_parts(); i++ ) {
            model_->surface_parts_[i].copy_macro_topology( from->surface_part( i ),
                *model_ ) ;
        }
#pragma omp parallel for
        for( uint32 i = 0; i < model_->nb_layers(); i++ ) {
            model_->layers_[i].copy_macro_topology( from->layer( i ), *model_ ) ;
        }
#pragma omp parallel for
        for( uint32 i = 0; i < model_->nb_regions(); i++ ) {
            model_->regions_[i].copy_macro_topology( from->region( i ), *model_ ) ;
        }
#pragma omp parallel for
        for( uint32 i = 0; i < model_->nb_contacts(); i++ ) {
            model_->contacts_[i].copy_macro_topology( from->contact( i ), *model_ ) ;
        }
#pragma omp parallel for
        for( uint32 i = 0; i < model_->nb_surfaces(); i++ ) {
            model_->surfaces_[i].copy_macro_topology( from->surface( i ), *model_ ) ;
        }
        model_->universe_.copy_macro_topology( from->universe_, *model_ ) ;
    }

    void BoundaryModelBuilder::load_file( std::istream& in )
    {
        time_t start_load, end_load ;
        time( &start_load ) ;
        InputStream lis( in ) ;

        // Intermediate for contact parts building
        std::vector< Border > borders_to_build ;

        /******* Counters, parameters and vectors used in the loop reading the file *********/
        int32 nb_tsurf = 0 ;
        int32 nb_tface = 0 ;

        int32 tsurf_count = 0 ;
        int32 tface_count = 0 ;

        int32 current_nb_tfaces = 0 ;
        int32 nb_tface_in_prev_tsurf = 0 ;

        // The reading of the file is in 2 steps
        // 1. Read model info (true)
        // 2. Read surfaces geometry (false)
        bool read_model = true ;

        // Depending on sign of positive Z
        // read for each TSurf of the file
        int32 z_sign = 1 ;

        // Vertices of the currently built TSurf
        std::vector< uint32 > tface_vertex_ptr ;
        // Triangles (3 ids per trgl) of the currently built TFace (part of a TSurf)
        std::vector< uint32 > tface_triangles ;
        // Starting and ending indices of each triangle in tface_triangles
        std::vector< uint32 > tface_ptr ;
        tface_ptr.push_back( 0 ) ;
        // For each TFace of the current TSurf, offest for the first vertex
        std::vector< int32 > v_offset ;

        // Ids of the surface parts on which the KeyFacet should be changed
        std::vector< int32 > change_key_facet ;

        /****** File reading **********************************/
        while( !lis.eof() ) {

            lis.get_line() ;
            std::string keyword ;
            lis >> keyword ;

            if( read_model ) {
                if( keyword == "TSURF" ) {
                    /// 1. Build all the interfaces
                    std::string temp_str ;
                    std::stringstream tsurf_name ;
                    lis >> temp_str ;
                    tsurf_name << temp_str ;
                    while( !lis.eol() ) {
                        lis >> temp_str ;
                        tsurf_name << "_" << temp_str ;
                    }
                    create_surface( tsurf_name.str() ) ;
                    nb_tsurf++ ;
                } else if( keyword == "TFACE" ) {
                    /// 2. Build the surface parts from their parent interface and type
                    int32 id ;
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

                    // Get the key tface that give the orientation of the surface part
                    // Triangles in Gocad clockwise
                    vec3 p0, p1, p2 ;
                    lis.get_line() ;
                    lis >> p0 ;
                    lis.get_line() ;
                    lis >> p1 ;
                    lis.get_line() ;
                    lis >> p2 ;

                    create_surface_part( tsurf_name.str(), type,
                        KeyFacet( p0, p1, p2 ) ) ;
                    nb_tface++ ;
                } else if( keyword == "REGION" ) {
                    /// 3. Build the volumetric regions from their name and
                    /// boundary surface parts
                    int32 id ;
                    std::string name ;
                    lis >> id >> name ;

                    std::vector< std::pair< int, bool > > region_boundaries ;
                    bool end_region = false ;

                    while( !end_region ) {
                        lis.get_line() ;
                        for( uint32 i = 0; i < 5; ++i ) {
                            int32 tface_id ;
                            lis >> tface_id ;
                            if( tface_id == 0 ) {
                                end_region = true ;
                                break ;
                            } else {
                                // Correction because ids begin at 1 in the file
                                if( tface_id > 0 )
                                    region_boundaries.push_back(
                                        std::pair< int, bool >( tface_id - 1,
                                            true ) ) ;
                                else
                                    region_boundaries.push_back(
                                        std::pair< int, bool >( -tface_id - 1,
                                            false ) ) ;
                            }
                        }
                    }
                    /// The Universe region is not added to the Surface Model
                    if( name != "Universe" )
                        create_region( name, region_boundaries ) ;
                    else
                        set_universe( region_boundaries ) ;
                } else if( keyword == "LAYER" ) {
                    /// 4. Build the volumetric layers
                    std::string name ;
                    lis >> name ;
                    uint32 layer_id = create_layer( name ) ;
                    bool end_layer = false ;
                    while( !end_layer ) {
                        lis.get_line() ;
                        for( uint32 i = 0; i < 5; ++i ) {
                            int32 region_id ;
                            lis >> region_id ;
                            if( region_id == 0 ) {
                                end_layer = true ;
                                break ;
                            } else {
                                region_id -= nb_tface + 1 ; // Remove Univers region
                                // Correction because ids begin at 1 in the file
                                add_layer_child( layer_id, region_id - 1 ) ;
                            }
                        }
                    }
                } else if( keyword == "END" ) {
                    read_model = false ;
                    continue ;
                }
            } else {
                if( keyword == "GOCAD" ) {
                    // This is the beginning of a new TSurf
                    tsurf_count++ ;
                }
                if( keyword == "ZPOSITIVE" ) {
                    std::string positive ;
                    lis >> positive ;
                    if( positive == "Elevation" )
                        z_sign = 1 ;
                    else if( positive == "Depth" ) z_sign = -1 ;
                } else if( keyword == "END" ) {
                    if( tsurf_count > 0 ) {
                        // End last surface part
                        if( !model_->surface_parts_[tface_count - 1].set_points_and_facets(
                            std::vector< uint32 >(
                                tface_vertex_ptr.begin() + v_offset.back(),
                                tface_vertex_ptr.end() ), tface_triangles,
                            tface_ptr ) ) {
                            change_key_facet.push_back( tface_count - 1 ) ;
                        }

                        tface_triangles.clear() ;
                        tface_ptr.clear() ;
                        tface_ptr.push_back( 0 ) ;

                        // End the previous surface
                        nb_tface_in_prev_tsurf += v_offset.size() ;
                        tface_vertex_ptr.clear() ;
                        v_offset.clear() ;
                    }
                } else if( keyword == "TFACE" ) {
                    // Beginning of a new surface part
                    if( v_offset.size() > 0 ) {
                        // End last surface part
                        if( !model_->surface_parts_[tface_count - 1].set_points_and_facets(
                            std::vector< uint32 >(
                                tface_vertex_ptr.begin() + v_offset.back(),
                                tface_vertex_ptr.end() ), tface_triangles,
                            tface_ptr ) ) {
                            change_key_facet.push_back( tface_count - 1 ) ;
                        }
                        tface_triangles.clear() ;
                        tface_ptr.clear() ;
                        tface_ptr.push_back( 0 ) ;
                    }
                    // Register where begin the new TFace vertices
                    v_offset.push_back( tface_vertex_ptr.size() ) ;

                    tface_count++ ;
                }
                /// 4. Read the surface parts and set their geometry
                else if( keyword == "VRTX" || keyword == "PVRTX" ) {
                    int32 id ;
                    vec3 p ;
                    lis >> id >> p ;
                    p.z *= z_sign ;
                    tface_vertex_ptr.push_back( model_->nb_points() ) ;
                    add_point( p ) ;
                } else if( keyword == "PATOM" || keyword == "ATOM" ) {
                    int32 id ;
                    int32 v_id ;
                    lis >> id >> v_id ;
                    tface_vertex_ptr.push_back( tface_vertex_ptr[v_id - 1] ) ;
                } else if( keyword == "TRGL" ) {
                    int32 p1, p2, p3 ;
                    lis >> p1 >> p2 >> p3 ;
                    p1 += -v_offset.back() - 1 ;
                    p2 += -v_offset.back() - 1 ;
                    p3 += -v_offset.back() - 1 ;

                    tface_triangles.push_back( p1 ) ;
                    tface_triangles.push_back( p2 ) ;
                    tface_triangles.push_back( p3 ) ;
                    tface_ptr.push_back( tface_triangles.size() ) ;
                }
                /// 5. Build the corners from their position and the surface parts
                ///    containing them
                else if( keyword == "BSTONE" ) {
                    int32 v_id ;
                    lis >> v_id ;
                    // correction to start at 0
                    v_id-- ;

                    // Get the TFace
                    int32 part_id = v_offset.size() - 1 ;
                    for( int32 i = 0; i < v_offset.size(); ++i ) {
                        if( v_id < v_offset[i] ) {
                            part_id = i - 1 ;
                            break ;
                        }
                    }
                    part_id += nb_tface_in_prev_tsurf ;

                    int32 new_c = find_or_create_corner( tface_vertex_ptr[v_id] ) ;
                }
                /// 6. Read the Border information and store it
                else if( keyword == "BORDER" ) {
                    int32 id, p1, p2 ;
                    lis >> id >> p1 >> p2 ;
                    p1-- ;
                    p2-- ;

                    // Get the global corner id
                    int32 corner_id = model_->find_corner( model_->point( tface_vertex_ptr[p1] ) ) ;
                    grgmesh_debug_assert( corner_id > -1 ) ;

                    // Get the surface part
                    int32 part_id = -1 ;
                    for( int32 i = 0; i < v_offset.size(); ++i ) {
                        if( p1 < v_offset[i] ) {
                            grgmesh_debug_assert( p2 < v_offset[i] ) ;

                            // Get points ids in the interface part
                            p1 += -v_offset[i - 1] ;
                            p2 += -v_offset[i - 1] ;

                            // i-1 is the id of the TFace in this TSurf
                            part_id = i - 1 ;

                            break ;
                        }
                    }
                    if( part_id == -1 ) {
                        // It is in the last built Tface
                        p1 += -v_offset[v_offset.size() - 1] ;
                        p2 += -v_offset[v_offset.size() - 1] ;
                        part_id = v_offset.size() - 1 ;
                    }
                    // the number of tfaces in previous tsurf is also to add
                    part_id += nb_tface_in_prev_tsurf ;

                    borders_to_build.push_back(
                        Border( part_id, corner_id, p1, p2 ) ) ;
                }
            }
        }

        /// 7. Build the contact parts
        build_contact_parts( borders_to_build ) ;

        /// 8. Build the contacts
        build_contacts() ;

        // Then finish the job (the order matters)
        end_surface_parts( change_key_facet ) ;
        end_surfaces() ;
        end_contact_parts() ;
        end_contacts() ;
        end_corners() ;
        end_layers() ;
        end_model() ;

        time( &end_load ) ;
        std::cout << " Boundary model loading time"
            << difftime( end_load, start_load ) << " sec" << std::endl ;
    }

    void BoundaryModelBuilder::end_model()
    {
        MakeUnique mu( model_->points_ ) ;
        mu.unique( 2 ) ;
        model_->points_.clear() ;
        mu.unique_points( model_->points_ ) ;

        const std::vector< int32 >& indices = mu.indices() ;
        for( uint32 s = 0; s < model_->nb_surface_parts(); s++ ) {
            SurfacePart& surface = model_->surface_parts_[s] ;
            for( uint32 p = 0; p < surface.nb_vertices(); p++ ) {
                surface.points_[p] = indices[surface.points_[p]] ;
            }
        }
        for( uint32 c = 0; c < model_->nb_contact_parts(); c++ ) {
            ContactPart& contact = model_->contact_parts_[c] ;
            for( uint32 p = 0; p < contact.nb_vertices(); p++ ) {
                contact.vertices_[p] = indices[contact.vertices_[p]] ;
            }
        }
        for( uint32 co = 0; co < model_->nb_corners(); co++ ) {
            Corner& corner = model_->corners_[co] ;
            corner.p_ = indices[corner.p_] ;
        }
    }

    void BoundaryModelBuilder::set_universe(
        const std::vector< std::pair< int, bool > >& boundaries )
    {
        model_->universe_.set_name( "Universe" ) ;
        model_->universe_.set_dim( 3 ) ;

        for( uint32 i = 0; i < boundaries.size(); ++i ) {
            grgmesh_debug_assert( boundaries[i].first < model_->nb_surface_parts() ) ;
            model_->universe_.add_oriented_boundary( boundaries[i].first,
                boundaries[i].second ) ;

            // If this surface_part have no type, set it at VOI
            model_->surface_parts_[boundaries[i].first].set_type( VOI ) ;
        }
    }

    uint32 BoundaryModelBuilder::create_region(
        const std::string& name,
        const std::vector< std::pair< int, bool > >& boundaries,
        int32 id )
    {
        id = create_region( id ) ;
        model_->regions_[id].set_name( name ) ;

        for( uint32 i = 0; i < boundaries.size(); ++i ) {
            grgmesh_debug_assert( boundaries[i].first < model_->nb_surface_parts() ) ;
            add_region_oriented_boundary( id, boundaries[i].first,
                boundaries[i].second ) ;
        }
        return id ;
    }


    /*! Build the ContactParts
     *  from the information collected in the Border structures
     *
     *  fill their in_boundary_ vector
     */
    void BoundaryModelBuilder::build_contact_parts( const std::vector< Border >& borders )
    {

        std::vector< vec3 > border_vertices ;
        std::vector< int32 > triangles_sharing_p1 ;

        for( uint32 i = 0; i < borders.size(); ++i ) {
            const Border& b = borders[i] ;

            /// 1- Build the boundary : construct the vector
            /// of vertices on the border
            const SurfacePart& part = model_->surface_part( b.part_id_) ;

            // Stuff used in while loop
            int32 id0 = b.p0_ ;
            int32 id1 = b.p1_ ;
            grgmesh_debug_assert( id0 < part.nb_vertices() && id1 < part.nb_vertices() ) ;

            vec3 p0 = part.vertex( id0 ) ;
            vec3 p1 = part.vertex( id1 ) ;

            int32 t = part.find_triangle( id0, id1 ) ;
            grgmesh_debug_assert( t != -1 ) ;

            border_vertices.resize( 0 ) ;
            border_vertices.push_back( p0 ) ;
            border_vertices.push_back( p1 ) ;

            // We want the next triangle that is on the boundary and share p1
            // If there is no such triangle, the third point of the current triangle is to add

            int32 p1_corner = model_->find_corner( p1 ) ;
            while( p1_corner == -1 ) {

                int32 nb_t = part.triangles_around_point( id1, triangles_sharing_p1,
                    true ) ;
                grgmesh_debug_assert( nb_t < 3 && nb_t > 0 ) ;

                int32 other_t = triangles_sharing_p1[0] ;
                int32 new_id1 = -1 ;

                if( nb_t == 2 ) {
                    if( other_t == t ) other_t = triangles_sharing_p1[1] ;
                    grgmesh_debug_assert( other_t != t ) ;

                    // Now get the other point that is on the boundary opposite to p1
                    int32 p1_id_other = part.point_id( other_t, id1 ) ;
                    grgmesh_debug_assert( p1_id_other != -1 ) ;

                    // Edges containing p1 in the other_t
                    int32 e0 = part.edge_vertex( other_t, p1_id_other, 0 ) ;
                    int32 e1 = part.edge_vertex( other_t, p1_id_other, 1 ) ;

                    bool be0 = part.is_on_border( other_t, e0 ) ;
                    bool be1 = part.is_on_border( other_t, e1 ) ;

                    // Only one must be on the boundary otherwise there is a corner missing
                    grgmesh_debug_assert( be0 != be1 ) ;

                    int32 border_e = be0 ? e0 : e1 ;
                    int32 new_p1_inother =
                        part.edge_vertex( other_t, border_e, 0 ) == p1_id_other ?
                            part.edge_vertex( other_t, border_e, 1 ) :
                            part.edge_vertex( other_t, border_e, 0 ) ;

                    new_id1 = part.point_index( other_t, new_p1_inother ) ;
                } else if( nb_t == 1 ) {
                    // p1 must be in two border edges in t triangle
                    grgmesh_debug_assert( other_t == t ) ;

                    int32 v1 = part.point_id( t, id1 ) ;
                    int32 v0 = part.point_id( t, id0 ) ;

                    new_id1 =
                        part.edge_vertex( t, v1, 0 ) == v0 ?
                            part.point_index( t, part.edge_vertex( t, v1, 1 ) ) :
                            part.point_index( t, part.edge_vertex( t, v1, 0 ) ) ;
                }

                grgmesh_debug_assert( new_id1 != -1 ) ;
                // Update p1, id0, and id1  and current triangle
                t = other_t ;
                id0 = id1 ;
                id1 = new_id1 ;
                p1 = part.vertex( new_id1 ) ;
                p1_corner = model_->find_corner( p1 ) ;
                border_vertices.push_back( p1 ) ;
            }

            // 2 - Check if this border already exists
            int32 contact_part_id = find_or_create_contact_part( b.corner_id_,
                p1_corner, border_vertices ) ;

            // Add that the surface part in which this line is
            add_contact_part_in_boundary( contact_part_id, b.part_id_ ) ;
        }
    }

    /*! Build the contact form the contact_parts_
     *
     */
    void BoundaryModelBuilder::build_contacts()
    {
        for( uint32 i = 0; i < model_->nb_contact_parts(); ++i ) {
            // The surface part in whose boundary is the part
            std::set< int32 > interfaces ;
            std::vector< GEOL_FEATURE > types ;
            for( uint32 j = 0; j < model_->contact_parts_[i].nb_in_boundary(); ++j ) {
                uint32 sp_id = model_->contact_parts_[i].in_boundary_id( j ) ;
                const BoundaryModelElement* p = model_->surface_parts_[sp_id].parent() ;
                grgmesh_debug_assert( p != nil ) ;
                interfaces.insert( p->id() ) ;
                types.push_back( p->type() ) ;
            }
            std::vector< int32 > toto( interfaces.begin(), interfaces.end() ) ;
            int32 contact_id = find_or_create_contact( toto, model_->determine_type( types ) ) ;
            add_contact_child( contact_id, i ) ;
        }
    }

    /*! Set the surfaces delimited by each contact
     *  Warning the surfaces must be finished first
     */
    void BoundaryModelBuilder::end_contacts()
    {
        for( uint32 i = 0; i < model_->nb_contacts(); ++i ) {
            std::set< const BoundaryModelElement* > corners ;
            for( uint32 j = 0; j < model_->contacts_[i].nb_children(); ++j ) {
                const BoundaryModelElement* child = model_->contacts_[i].child( j ) ;
                corners.insert( child->boundary( 0 ) ) ;
                corners.insert( child->boundary( 1 ) ) ;
            }
            for( std::set< const BoundaryModelElement* >::iterator it( corners.begin() );
                it != corners.end(); ++it ) {
                add_contact_boundary( i, ( *it )->id() ) ;
            }
        }
    }

    /*! Set the parent contact for all the contact parts
     */
    void BoundaryModelBuilder::end_contact_parts()
    {
        for( uint32 i = 0; i < model_->nb_contacts(); ++i ) {
            for( uint32 j = 0; j < model_->contacts_[i].nb_children(); ++j ) {
                uint32 child = model_->contacts_[i].child_id( j ) ;
                model_->contact_parts_[child].set_parent( i ) ;
                model_->contact_parts_[child].set_type( model_->contacts_[i].type() ) ;
            }
        }
    }

    void BoundaryModelBuilder::end_surfaces()
    {
        for( uint32 i = 0; i < model_->nb_surface_parts(); ++i ) {
            int32 parent = model_->surface_parts_[i].parent_id() ;
            add_surface_child( parent, i ) ;
        }

        for( uint32 i = 0; i < model_->nb_contacts(); ++i ) {
            for( uint32 j = 0; j < model_->contacts_[i].nb_in_boundary(); ++j ) {
                uint32 b = model_->contacts_[i].in_boundary_id( j ) ;
                add_surface_boundary( b, i ) ;
            }
        }
    }

    void BoundaryModelBuilder::end_layers()
    {
        for( uint32 i = 0; i < model_->nb_layers(); ++i ) {
            BoundaryModelElement& layer = model_->layers_[i] ;
            std::set< uint32 > boundaries ;
            for( uint32 r = 0; r < layer.nb_children(); r++ ) {
                const BoundaryModelElement* region = layer.child( r ) ;
                model_->regions_[region->id()].parent_ = i ;
                for( uint32 sp = 0; sp < region->nb_boundaries(); sp++ ) {
                    boundaries.insert( region->boundary_id( sp ) ) ;
                }
            }

            for( std::set< uint32 >::const_iterator it( boundaries.begin() );
                it != boundaries.end(); it++ ) {
                layer.add_boundary( *it ) ;
            }
        }
    }

    /*! Set the parent of surface_parts
     */
    void BoundaryModelBuilder::end_surface_parts()
    {
        for( uint32 i = 0; i < model_->nb_contact_parts(); ++i ) {
            for( uint32 j = 0; j < model_->contact_parts_[i].nb_in_boundary();
                ++j ) {
                uint32 s_id = model_->contact_parts_[i].in_boundary_id( j ) ;
                add_surface_part_boundary( s_id, i ) ;
            }
        }
        for( uint32 i = 0; i < model_->nb_regions(); ++i ) {
            for( uint32 j = 0; j < model_->regions_[i].nb_boundaries(); ++j ) {
                uint32 s_id = model_->regions_[i].boundary_id( j ) ;
                add_surface_part_in_boundary( s_id, i ) ;
            }
        }
    }

    void BoundaryModelBuilder::end_surface_parts(
        const std::vector< int32 >& change_orientation )
    {
        end_surface_parts() ;
        for( int32 i = 0; i < change_orientation.size(); i++ ) {
            // Change the key facet
            SurfacePart& sp = model_->surface_parts_[change_orientation[i]] ;
            sp.set_first_triangle_as_key() ;
            grgmesh_debug_assert( change_orientation[i] == sp.id() ) ;

            // Change the sign of this SP in all the regions containing it
            for( int32 j = 0; j < sp.nb_in_boundary(); ++j ) {
                uint32 region = sp.in_boundary_id( j ) ;
                model_->regions_[region].change_boundary_side(
                    change_orientation[i] ) ;
            }
        }
    }

    /*! Set in what surface parts are the corners
     */
    void BoundaryModelBuilder::end_corners()
    {
        for( uint32 i = 0; i < model_->nb_contact_parts(); ++i ) {
            uint32 c0_id = model_->contact_parts_[i].boundary_id( 0 ) ;
            uint32 c1_id = model_->contact_parts_[i].boundary_id( 1 ) ;

            add_corner_in_boundary( c0_id, i ) ;
            if( c1_id != c0_id ) add_corner_in_boundary( c1_id, i ) ;
        }
    }

    int32 BoundaryModelBuilder::find_or_create_contact(
        std::vector< int32 >& interfaces,
        GEOL_FEATURE type )
    {
        int32 result = model_->find_contact( interfaces ) ;
        if( result == -1 ) {
            // Create a name for this contact
            std::string name = "contact_" ;
            for( uint32 i = 0; i < interfaces.size(); ++i ) {
                name += model_->surfaces_[interfaces[i]].name() ;
                name += "_" ;
            }
            result = model_->nb_contacts() ;
            model_->contacts_.push_back( BoundaryModelElement( model_, 1, result ) ) ;
            model_->contacts_[result].set_name( name ) ;
            model_->contacts_[result].set_type( type ) ;

            for( uint32 i = 0; i < interfaces.size(); ++i ) {
                add_contact_in_boundary( result, interfaces[i] ) ;
            }
        }
        return result ;
    }

    /*! Add a surface part to the model
     *  The parent with the given name (interface name) MUST exist
     */
    void BoundaryModelBuilder::create_surface_part(
        const std::string& interface_name,
        const std::string& type,
        const KeyFacet& key )
    {
        int32 parent = model_->surface_id( interface_name ) ;
        grgmesh_debug_assert( parent != -1 ) ;

        int32 id = model_->nb_surface_parts() ;
        GEOL_FEATURE t = model_->determine_geological_type( type ) ;
        model_->surface_parts_.push_back( SurfacePart( model_, id, parent, t ) ) ;
        model_->surface_parts_[id].set_key_facet( key ) ;

        model_->surfaces_[parent].set_type( t ) ;
    }

    int32 BoundaryModelBuilder::find_or_create_corner( uint32 index )
    {
        int32 result = model_->find_corner( model_->point( index ) ) ;
        if( result == -1 ) {
            // Create the corner
            result = model_->nb_corners() ;
            model_->corners_.push_back( Corner( model_, result, index ) ) ;
        }
        return result ;
    }

    int32 BoundaryModelBuilder::find_or_create_contact_part(
        int32 corner0,
        int32 corner1,
        std::vector< vec3 >& points )
    {
        int32 result = model_->find_contact_part( corner0, corner1, points ) ;

        if( result == -1 ) {
            result = model_->nb_contact_parts() ;
            if( corner1 == -1 ) {
                // Closed contact part
                corner1 = corner0 ;
                // Remove the last point
                grgmesh_debug_assert( points[0] == points.back() ) ;
                points.pop_back() ;
            }

            std::vector< uint32 > indices( points.size() ) ;
            for( uint32 p = 0; p < points.size(); p++ ) {
                indices[p] = model_->nb_points() ;
                add_point( points[p] ) ;
            }

            model_->contact_parts_.push_back(
                ContactPart( model_, result, corner0, corner1, indices ) ) ;
        }
        return result ;
    }

    bool BoundaryModelBuilder::rebuild()
    {
        /*! \todo Valgrind finds errors !!!!!!! */
        std::vector< uint32 > sp_to_remove ;
        std::vector< uint32 > new_sp_id( model_->nb_surface_parts() ) ;
        uint32 offset = 0 ;
        for( uint32 sp = 0; sp < model_->nb_surface_parts(); sp++ ) {
            if( model_->surface_parts_[sp].nb_simplices() == 0 ) {
                offset++ ;
                sp_to_remove.push_back( sp ) ;
            } else {
                model_->surface_parts_[sp - offset] = model_->surface_parts_[sp] ;
            }
            new_sp_id[sp] = sp - offset ;
        }
        if( offset == 0 ) return false ;
        model_->surface_parts_.erase( model_->surface_parts_.end() - offset,
            model_->surface_parts_.end() ) ;

        offset = 0 ;
        std::vector< uint32 > cp_to_remove ;
        std::vector< uint32 > new_cp_id( model_->nb_contact_parts() ) ;
        for( uint32 cp = 0; cp < model_->nb_contact_parts(); cp++ ) {
            BoundaryModelElement& contact_part = model_->contact_parts_[cp] ;
            uint32 nb_sp_removed = 0 ;
            for( uint32 sp = 0; sp < contact_part.nb_in_boundary(); sp++ ) {
                if( vector_contains( sp_to_remove,
                    contact_part.in_boundary_id( sp ) ) ) {
                    nb_sp_removed++ ;
                } else {
                    contact_part.in_boundary_[sp - nb_sp_removed] =
                        new_sp_id[contact_part.in_boundary_[sp]] ;
                }
            }
            if( nb_sp_removed > 0 ) {
                if( nb_sp_removed == contact_part.nb_in_boundary() ) {
                    offset++ ;
                    cp_to_remove.push_back( cp ) ;
                    continue ;
                } else {
                    contact_part.in_boundary_.erase(
                        contact_part.in_boundary_.end() - nb_sp_removed,
                        contact_part.in_boundary_.end() ) ;
                }
            }
            if( offset > 0 ) {
                model_->contact_parts_[cp - offset] = model_->contact_parts_[cp] ;
            }
            new_cp_id[cp] = cp - offset ;
        }
        if( offset > 0 ) {
            model_->contact_parts_.erase( model_->contact_parts_.end() - offset,
                model_->contact_parts_.end() ) ;
        }

        // Build the contacts
        // Update surfaces
        offset = 0 ;
        std::vector< uint32 > s_to_remove ;
        std::vector< uint32 > new_s_id( model_->nb_surfaces() ) ;
        for( uint32 s = 0; s < model_->nb_surfaces(); s++ ) {
            BoundaryModelElement& surface = model_->surfaces_[s] ;
            surface.boundaries_.clear() ;
            uint32 nb_sp_removed = 0 ;
            for( uint32 sp = 0; sp < surface.nb_children(); sp++ ) {
                if( vector_contains( sp_to_remove, surface.child_id( sp ) ) ) {
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
                model_->surfaces_[s - offset] = model_->surfaces_[s] ;
            }
            new_s_id[s] = s - offset ;
        }
        if( offset > 0 ) {
            model_->surfaces_.erase( model_->surfaces_.end() - offset,
                model_->surfaces_.end() ) ;
        }
        for( uint32 sp = 0; sp < model_->nb_surface_parts(); sp++ ) {
            BoundaryModelElement& surface_part = model_->surface_parts_[sp] ;
            surface_part.parent_ = new_s_id[surface_part.parent_] ;
            for( uint32 cp = 0; cp < surface_part.nb_boundaries(); cp++ ) {
                surface_part.boundaries_[cp] =
                    new_cp_id[surface_part.boundaries_[cp]] ;
            }
        }
        model_->contacts_.clear() ;
        build_contacts() ;

        // Then finish the job (the order matters)
        for( uint32 sp = 0; sp < model_->nb_surface_parts(); sp++ ) {
            model_->surface_parts_[sp].boundaries_.clear() ;
        }
        offset = 0 ;
        for( uint32 r = 0; r < model_->nb_regions(); r++ ) {
            BoundaryModelElement& region = model_->regions_[r] ;
            uint32 nb_sp_removed = 0 ;
            for( uint32 sp = 0; sp < region.nb_boundaries(); sp++ ) {
                if( vector_contains( sp_to_remove, region.boundary_id( sp ) ) ) {
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
                model_->regions_[r - offset] = model_->regions_[r] ;
            }
        }
        if( offset > 0 ) {
            model_->regions_.erase( model_->regions_.end() - offset,
                model_->regions_.end() ) ;
        }
        end_surface_parts() ;

        for( uint32 i = 0; i < model_->nb_contacts(); ++i ) {
            for( uint32 j = 0; j < model_->contacts_[i].nb_in_boundary();
                ++j ) {
                uint32 b = model_->contacts_[i].in_boundary_id( j ) ;
                add_surface_boundary( b, i ) ;
            }
        }

        end_contact_parts() ;
        end_contacts() ;

        offset = 0 ;
        std::vector< uint32 > co_to_remove ;
        std::vector< uint32 > new_co_id( model_->nb_corners() ) ;
        for( uint32 co = 0; co < model_->nb_corners(); co++ ) {
            BoundaryModelElement& corner = model_->corners_[co] ;
            uint32 nb_cp_removed = 0 ;
            for( uint32 cp = 0; cp < corner.nb_in_boundary(); cp++ ) {
                if( vector_contains( cp_to_remove, corner.in_boundary_id( cp ) ) ) {
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
                model_->corners_[co - offset] = model_->corners_[co] ;
            }
            new_co_id[co] = co - offset ;
        }
        if( offset > 0 ) {
            model_->corners_.erase( model_->corners_.end() - offset,
                model_->corners_.end() ) ;
        }

        for( uint32 cp = 0; cp < model_->nb_contact_parts(); cp++ ) {
            BoundaryModelElement& contact_part = model_->contact_parts_[cp] ;
            contact_part.boundaries_[0] = new_co_id[contact_part.boundaries_[0]] ;
            contact_part.boundaries_[1] = new_co_id[contact_part.boundaries_[1]] ;
        }

        for( uint32 c = 0; c < model_->nb_contacts(); c++ ) {
            BoundaryModelElement& contact = model_->contacts_[c] ;
            for( uint32 co = 0; co < contact.nb_boundaries(); co++ ) {
                contact.boundaries_[co] = new_co_id[contact.boundaries_[co]] ;
            }
        }

        for( uint32 l = 0; l < model_->nb_layers(); l++ ) {
            BoundaryModelElement& layer = model_->layers_[l] ;
            uint32 nb_sp_removed = 0 ;
            for( uint32 sp = 0; sp < layer.nb_boundaries(); sp++ ) {
                if( vector_contains( sp_to_remove, layer.boundary_id( sp ) ) ) {
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

    void BoundaryModelBuilder::update_all_ids()
    {
        for( uint32 co = 0; co < model_->nb_corners(); co++ ) {
            model_->corners_[co].id_ = co ;
        }
        for( uint32 cp = 0; cp < model_->nb_contact_parts(); cp++ ) {
            model_->contact_parts_[cp].id_ = cp ;
        }
        for( uint32 sp = 0; sp < model_->nb_surface_parts(); sp++ ) {
            model_->surface_parts_[sp].id_ = sp ;
        }
        for( uint32 c = 0; c < model_->nb_contacts(); c++ ) {
            model_->contacts_[c].id_ = c ;
        }
        for( uint32 s = 0; s < model_->nb_surfaces(); s++ ) {
            model_->surfaces_[s].id_ = s ;
        }
        for( uint32 r = 0; r < model_->nb_regions(); r++ ) {
            model_->regions_[r].id_ = r ;
        }
        for( uint32 l = 0; l < model_->nb_layers(); l++ ) {
            model_->layers_[l].id_ = l ;
        }
    }

    /*
    void BoundaryModelBuilder::set_surface_part_map( uint32 id, Map* map )
    {
        grgmesh_debug_assert( id < model_->nb_surface_parts() ) ;
        SurfacePart& surface = model_->surface_parts_[id] ;
        MapVertexAttribute< int32 > vertex_id( map ) ;
        surface.points_.reserve( map->size_of_vertices() ) ;
        int32 start = 0 ;
        FOR_EACH_VERTEX( Map, map, v ){
        surface.points_.push_back( model_->nb_points() ) ;
        model_->points_.push_back( v->point() ) ;
        vertex_id[v] = start++ ;
    }

        MapFacetAttribute< int32 > facet_id( map ) ;
        surface.facets_.reserve( 3 * map->size_of_facets() ) ;
        surface.facet_ptr_.reserve( 3 * map->size_of_facets() ) ;
        start = 0 ;
        surface.facet_ptr_.push_back( 0 ) ;
        FOR_EACH_FACET( Map, map, f ){
        Map::Halfedge* h = f->halfedge() ;
        do {
            surface.facets_.push_back( vertex_id[h->vertex()] ) ;
            h = h->next() ;
        }while( h != f->halfedge() ) ;
        facet_id[f] = start++ ;
        surface.facet_ptr_.push_back( surface.facets_.size() ) ;
    }

        surface.adjacent_.resize( surface.facets_.size(), -1 ) ;
        FOR_EACH_FACET( Map, map, f ){
        int32 t = facet_id[f] ;

        Map::Halfedge* h1 = f->halfedge() ;
        Map::Halfedge* h2 = h1->next() ;
        Map::Halfedge* h0 = h1->prev() ;

        // Adjacent facets
        Map::Facet* f0 = h0->opposite()->facet() ;
        Map::Facet* f1 = h1->opposite()->facet() ;
        Map::Facet* f2 = h2->opposite()->facet() ;

        if( f0 ) surface.adjacent_[3*t ] = facet_id[f0] ;
        if( f1 ) surface.adjacent_[3*t+1] = facet_id[f1] ;
        if( f2 ) surface.adjacent_[3*t+2] = facet_id[f2] ;
    }

        if( surface.nb_simplices() > 0 ) surface.set_first_triangle_as_key() ;

        surface.compute_is_triangulated() ;
    }
 */
}

