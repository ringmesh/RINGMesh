/*! This file is originally part of the Geomodeling plugin of Graphite */
/*! \author Jeanne Pellerin */


#include <grgmesh/boundary_model.h>
#include <grgmesh/boundary_model_element.h>
#include <grgmesh/utils.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <time.h>
#include <set>

namespace GRGMesh {

    BoundaryModel::BoundaryModel()
        : point_attribute_manager_(), facet_attribute_manager_()
    { 
        universe_ = BoundaryModelElement( this ) ;       
    }

    BoundaryModel::~BoundaryModel()
    {
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

    /** 
     * Returns the index of the given point in the model
     * \todo Add a KdTree for geometrical request on model points
     */ 
    unsigned int BoundaryModel::point_index( const vec3& p ) const {
       
        grgmesh_assert_not_reached ;
        return -1 ;
    }

    uint32 BoundaryModel::nb_facets() const {
        return nb_facets_.back() ;    
    }

    void BoundaryModel::surface_facet(
        uint32 model_facet_id, uint32& surface_id, uint32& surf_facet_id 
    ) const {
        int s = -1 ;
        for( uint32 i = 1; i < nb_facets_.size(); i++ ) {
            if( model_facet_id > nb_facets_[i-1] && 
                model_facet_id < nb_facets_[i] )
            {
                s = i-1 ;
            }
        }       
        grgmesh_debug_assert( s != -1 ) ;
            
        surface_id = s ;
        surf_facet_id = model_facet_id - nb_facets_[s] ;
    }

    uint32 BoundaryModel::model_facet( uint32 surf_id, uint32 facet_id ) const 
    {
        grgmesh_debug_assert( surf_id < nb_surfaces() && 
            facet_id < surface(surf_id).nb_cells() )  ;

        return nb_facets_[surf_id] + facet_id ;
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
        BoundaryModelBuilder builder( *this ) ;
        builder.load_file( in ) ;
        return true ;
    }

    /*! Check that the BoundaryModel can be saved as a .ml 
     */
    bool BoundaryModel::check_model3d_compatibility()
    {
         BoundaryModelBuilder builder( *this ) ;

        // Check that the Surfaces exist
        if( nb_interfaces() == 0 && nb_surfaces() > 0 ) {
          
            // Creates one Surface per Surface
            for( unsigned int i = 0; i < surfaces_.size(); ++i ) {
                // set name, type, links
                std::ostringstream name ;
                name << "surface_" << i ;
                unsigned int id = builder.create_interface( name.str() ) ;
                builder.add_interface_child( id, i ) ;
            }

            // Set links surface parts toward surfaces
            for( unsigned int i = 0; i < interfaces_.size(); ++i ) {
                builder.set_parent( surfaces_[interfaces_[i].child( 0 ).id()],  i ) ;
            }

            // Is it really useful to have contacts, let's hope not... I am not doing it
        }

        for( unsigned int i = 0; i < surfaces_.size(); ++i ) {
            Surface& sp = surfaces_[i] ;
            if( sp.nb_points() == 0 ) continue ;
            if( sp.key_facet().is_default() ) {
                builder.set_surface_first_triangle_as_key( sp.id() ) ;
            }
        }

        // the universe should exist, perhaps create a function copying the lines 
        // 805 to 834 de s2_b_model.cpp
        if( universe_.name() != "Universe" ) {
            std::cout <<  "Error"
                << "The region universe is not defined for the model. IMPLEMENTATION TO DO"
                << std::endl ;
            return false ;
        }

        // Check that each region has a name and valid surfaces
        for( unsigned int i = 0; i < regions_.size(); ++i ) {
            BoundaryModelElement& region = regions_[i] ;

            if( region.name() == "" ) {
                std::ostringstream name ;
                name << "region_" << i ;
                builder.set_name( region, name.str() ) ;
            }
            if( region.nb_boundaries() == 0 ) {
                std::cout << "Error"  << " The region " << region.name()
                    << " has no Surfaces on its boundary" << std::endl ;
                return false ;
            }
        }

        // Check that all surfaces_ of the model are trinagulated
        /// \todo Implement a triangulation function in SurfaceMutator         
        for( unsigned int s = 0; s < nb_surfaces(); s++ ) {
            if( !surfaces_[s].is_triangulated() ) {
                std::cout<< "Error" << "Surface "<< s << " is not triangulated" << std::endl ;
                return false ;
            }
        }       
        
        return true ;
    }

    void save_region(
        int count,
        const BoundaryModelElement& region,
        std::ostream& out )
    {
        out << "REGION " << count << "  " << region.name() << " " << std::endl ;
        int it = 0 ;

        for( uint32 i = 0; i < region.nb_boundaries(); ++i ) {
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

    void save_layer(
        int count,
        int offset,
        const BoundaryModelElement& layer,
        std::ostream& out )
    {
        out << "LAYER " << layer.name() << " " << std::endl ;
        int it = 0 ;

        for( int i = 0; i < layer.nb_children(); ++i ) {
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
            std::cout << "Error"  << "The BoundaryModel " << name_
                << " cannot be saved in .ml format " << std::endl ;
            return false ;
        }

        // Print Model3d headers
        out << "GOCAD Model3d 1" << std::endl << "HEADER {" << std::endl << "name:"
            << "model_from_graphite" << std::endl << "}" << std::endl ;

        save_coordinate_system( out ) ;

        // Print the TSurf = Surfaces info
        for( int i = 0; i < interfaces_.size(); ++i ) {
            out << "TSURF " << interfaces_[i].name() << std::endl ;
        }

        int count = 1 ;
        // Print the TFace info
        for( int i = 0; i < surfaces_.size(); ++i ) {
            const Surface& s = surfaces_[i] ;
            out << "TFACE " << count << "  " ;
            save_type( out, s.type() ) ;
            out << " " << s.parent().name() << std::endl ;

            const KeyFacet kf = s.key_facet() ;
            out << "  " << kf.p0_ << std::endl ;
            out << "  " << kf.p1_ << std::endl ;
            out << "  " << kf.p2_ << std::endl ;

            ++count ;
        }

        int offset_layer = count ;
        // Region info + universe
        save_region( count, universe_, out ) ;
        ++count ;

        for( int i = 0; i < regions_.size(); ++i ) {
            save_region( count, regions_[i], out ) ;
            ++count ;
        }

        for( int i = 0; i < layers_.size(); ++i ) {
            save_layer( count, offset_layer, layers_[i], out ) ;
            ++count ;
        }

        out << "END" << std::endl ;

        // Now save each one of the surfaces
        for( int i = 0; i < interfaces_.size(); ++i ) {
            const BoundaryModelElement& tsurf = interfaces_[i] ;

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
            unsigned int vertex_count = 1 ;
            unsigned int offset = vertex_count ;

            std::vector< unsigned int > bstones ;
            std::vector< unsigned int > next_point ;
            std::set< unsigned int > set_end_corners ;
            
            for( uint32 j = 0; j < tsurf.nb_children(); ++j ) {
                offset = vertex_count ;

                const Surface& sp = dynamic_cast< const Surface& >( tsurf.child( j ) ) ;

                out << "TFACE" << std::endl ;
                for( uint32 k = 0; k < sp.nb_points(); ++k ) {
                    out << "VRTX " << vertex_count << " " << sp.point( k ) << std::endl ;
                    vertex_count++ ;
                }

                for( uint32 k = 0; k < sp.nb_cells(); ++k ) {
                    out << "TRGL " << sp.surf_point_id( k, 0 ) + offset << " "
                        << sp.surf_point_id( k, 1 ) + offset << " "
                        << sp.surf_point_id( k, 2 ) + offset << std::endl ;
                }

                // Info on corners and contacts
                for( uint32 k = 0; k < sp.nb_boundaries(); ++k ) {
                    const Line& cp = dynamic_cast< const Line& >( sp.boundary( k ) ) ;                  
      
                    unsigned int c = cp.model_point_id( 0 ) ;
                    unsigned int next = cp.model_point_id( 1 ) ;

                    // To be sure that we have all corners we need to ensure
                    // that all corners at the end of lines are saved too
                    set_end_corners.insert( 
                        sp.surf_point_id( cp.model_point_id( cp.nb_points()-1 ) ) + offset ) ;                    

                    /// \todo Check that this works, I am not sure for inside_border that we get both sides

                    int t = sp.facet_from_model_point_ids( c, next ) ;
                    assert( t != -1 ) ;

                    int i0 = sp.surf_point_id( t, 0 ) ;
                    int i1 = sp.surf_point_id( t, 1 ) ;
                    int i2 = sp.surf_point_id( t, 2 ) ;

                    unsigned int p0 = sp.model_point_id( i0 ) ;
                    unsigned int p1 = sp.model_point_id( i1 ) ;
                    unsigned int p2 = sp.model_point_id( i2 ) ;

                    int c_id = -1 ;
                    int next_id = -1 ;

                    if( p0 == c ) c_id = i0 ;
                    else if( p0 == next ) next_id = i0 ;

                    if( p1 == c ) c_id = i1 ;
                    else if( p1 == next ) next_id = i1 ;

                    if( p2 == c ) c_id = i2 ;
                    else if( p2 == next ) next_id = i2 ;

                    grgmesh_assert( c_id != -1 && next_id != -1 ) ;

                    bstones.push_back( c_id + offset ) ;
                    next_point.push_back( next_id + offset ) ;
                }
            }

            // Print corners and contact
            std::vector< unsigned int > end_corners( 
                set_end_corners.begin(), set_end_corners.end() ) ;
            std::vector< bool > end_corner_to_print( end_corners.size(), true ) ;
            
            for( int j = 0; j < bstones.size(); ++j ) {
                out << "BSTONE " << bstones[j] << std::endl ;

                // Determine the corners at the end of the lines that are not saved
                for( uint32 k = 0; k < end_corners.size(); k++ ) {
                    if( bstones[j] == end_corners[k] ) {
                        end_corner_to_print[k] = false ;
                        break ;
                    }
                }
            }
            // Print the corners that were at the beginning of none of the contacts
            // in this Interface
            for( uint32 j= 0; j < end_corners.size(); j++ ) {
                if( end_corner_to_print[j] ) {
                    out << "BSTONE " << end_corners[j] << std::endl ;
                }
            }
            // Print the the information to build the lines :
            // index of the point at the corner and index of the second point on the line
            for( int j = 0; j < bstones.size(); ++j ) {
                out << "BORDER " << vertex_count << " " << bstones[j] << " "
                    << next_point[j] << std::endl ;
                vertex_count++ ;
            }
            out << "END" << std::endl ;
        }
        return true ;
    }    

    /**
     *  \brief Save the surfaces of the model into an .eobj file.     
     */
    void BoundaryModel::save_as_eobj_file( const std::string& file_name ) {
        std::ofstream out ;
        out.open( file_name.c_str() );
        if( out.bad() ){
            std::cout << "Error when opening the file: " << file_name.c_str() <<std::endl ;
            return ;
        }
        out.precision( 16 ) ;
        std::vector< uint32 > offset( nb_surfaces(), 0 ) ;
        uint32 cur_offset = 0 ; 
        // Write the points once for each surface - remove non-manifold vertices
        for( unsigned int s = 0; s < nb_surfaces() ; s++ ) {
            const Surface& S = surface(s) ;
            offset[s] = cur_offset ;
            for( unsigned int p = 0; p < S.nb_points(); p++ ){
                const vec3& V = S.point(p) ;
                out << "v" 
                    << " " << V.x 
                    << " " << V.y 
                    << " " << V.z 
                    << std::endl ;
            }
            cur_offset += S.nb_points() ;
        }

        // Write the facets for a each surface 
        for( unsigned int s = 0; s < nb_surfaces() ; s++ ) {
            const Surface& S = surface(s) ;
            for( unsigned int f = 0; f < S.nb_cells(); f++ ){
                out << "f" << " " ;
                for( unsigned int v = 0; v < S.nb_points_in_facet(f); v++ ){                    
                    out << offset[s] + S.surf_point_id( f, v )+1 << " " ;            
                }
                out << std::endl ;
            }         
        }


        // Write facet attributes -- WARNING ON SUPPOSE QU'ON STOCKE DES ENTIERS
        /// \todo  Generic IO for attributes on a BoundaryModel
        {
            std::vector< std::string > names ;        
            facet_attribute_manager_.list_named_attributes(names) ;
            std::vector< AttributeStore* > attribute_stores( names.size(), nil ) ;
                    
            for( unsigned int i=0; i < names.size(); i++ ) 
            {
                // Crash if we cannot get the
                attribute_stores[i] = facet_attribute_manager_.resolve_named_attribute_store( names[i] ) ;
              
                // Output global information on the attribute
                out << "# attribute "<< names[i] << " facet "
                    << "integer" // dans graphite map à la noix entre vrai type et nom dans fichier
                    //<< attribute_stores[i]->attribute_type_id().name() 
                    << std::endl ;  
            }
            if( names.size() > 0 ) {
                uint32 count = 0 ;
                for( unsigned int s = 0; s < nb_surfaces() ; s++ ) {
                    const Surface& S = surface(s) ;
                    for( unsigned int f = 0; f < S.nb_cells(); f++ ){            
                        out << "# attrs f " << count+1 ;
                        for( uint32 j = 0; j < names.size(); j++ ) {
                            out << " " << *reinterpret_cast<uint32*>( attribute_stores[j]->data(count) ) ;
                        }
                        out << std::endl ;
                        count++ ;
                    }
                }
            }
        }

    }


/*************************************************************************************************/
/********  BoundarymodelMeasure implementation      ******************************************/
/*************************************************************************************************/
#ifdef TOTOTO
    void BoundaryModelMeasure::print_topology( std::ofstream& out ) const
    {
        bool skip_voi = true ; 

        out << "Topological information of model" << SEP << model_.name() << SEP
            << std::endl ;

        out << "Numbers of" << SEP << "Regions" << SEP << "Surfaces Parts Inside"
            << SEP << "Surface Parts Box" << SEP << "Blind Surface parts" << SEP
            << "Contact Parts Inside" << SEP << "Contact Parts Box" << SEP
            << "Corners Inside" << SEP << "Corners Box" << SEP << std::endl ;

        int s_in = nb_surface_inside() ;
        int b_in = nb_line_inside() ;
        int c_in = nb_real_corners_inside() ;

        out << SEP << model_.nb_regions() << SEP << s_in << SEP
            << model_.nb_surfaces() - s_in << SEP << nb_surface_with_free_boundary()
            << SEP << b_in << SEP << model_.nb_lines() - b_in << SEP << c_in
            << SEP << nb_real_corners() - c_in << SEP << std::endl ;

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "REGIONS" << std::endl ;
        for( uint32 i = 0; i < model_.nb_regions(); ++i ) {
            model_.region(i).print( out ) ;
        }

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;         

        out << "INTERFACES" << std::endl ;
        for( uint32 i = 0; i < model_.nb_interfaces(); ++i ) {
            if( skip_voi && model_.one_interface(i).is_on_voi() )
                continue ;
            else model_.one_interface(i).print( out ) ;
        }

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "CONTACTS" << std::endl ;
        for( uint32 i = 0; i < model_.nb_contacts(); ++i ) {
            if( skip_voi && model_.contact(i).is_on_voi() )
                continue ;
            else model_.contact(i).print( out ) ;
        }

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "SURFACE PARTS" << std::endl ;
        for( uint32 i = 0; i < model_.nb_surfaces(); ++i ) {
            if( skip_voi && model_.surface(i).is_on_voi() )
                continue ;
            else model_.surface(i).print( out ) ;
        }

        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "LINES" << std::endl ;
        for( uint32 i = 0; i < model_.nb_lines(); ++i ) {
            if( skip_voi && model_.line(i).is_on_voi() )
                continue ;
            else model_.line(i).print( out ) ;
        }
        out << std::endl ;
        BoundaryModelElement::print_categories( out ) ;

        out << "CORNERS" << std::endl ;
        for( uint32 i = 0; i < model_.nb_corners(); ++i ) {
            if( skip_voi && model_.corner(i).is_on_voi() )
                continue ;
            else model_.corner(i).print( out ) ;
        }
    }

    void BoundaryModelMeasure::print_element_info( std::ofstream& out ) const
    {
        out << model_.name() << std::endl ;

        BoundaryModelElement::print_complexity_categories( out ) ;

        for( unsigned int i = 0; i < model_.nb_regions(); ++i ) {
            model_.region(i).print_complexity( out ) ;
        }
        for( unsigned int i = 0; i < model_.nb_surfaces(); ++i ) {
            if( model_.surface(i).is_on_voi() ) continue ;
            model_.surface(i).print_complexity( out ) ;
        }
        for( unsigned int i = 0; i < model_.nb_lines(); ++i ) {
            if( model_.line(i).is_on_voi() ) continue ;
            model_.line(i).print_complexity( out ) ;
        }
        for( unsigned int i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).is_on_voi() ) continue ;
            model_.corner(i).print_complexity( out ) ;
        }
    }

    void BoundaryModelMeasure::print_type( std::ostream& out, GEOL_FEATURE t, int dim )
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

    unsigned int BoundaryModelMeasure::nb_surface_inside() const
    {
        unsigned int result = 0 ;
        for( unsigned int i = 0; i < model_.nb_surfaces(); ++i ) {
            if( model_.surface(i).type() != VOI ) ++result ;
        }
        return result ;
    }

    unsigned int BoundaryModelMeasure::nb_line_inside() const
    {
        unsigned int result = 0 ;
        for( unsigned int i = 0; i < model_.nb_lines(); ++i ) {
            if( !model_.line(i).is_on_voi() ) ++result ;
        }
        return result ;
    }

    unsigned int BoundaryModelMeasure::nb_real_corners_inside() const
    {
        unsigned int result = 0 ;
        for( int i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).is_real() && !model_.corner(i).is_on_voi() ) ++result ;
        }
        return result ;
    }

    int BoundaryModelMeasure::nb_real_corners() const
    {
        int result = 0 ;
         for( int i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).is_real() ) ++result ;
        }
        return result ;
    }

    int BoundaryModelMeasure::nb_surface_with_free_boundary() const
    {
        int result = 0 ;
        for( unsigned int i = 0; i < model_.nb_surfaces(); ++i ) {
            for( unsigned int j = 0; j < model_.surface(i).nb_boundaries(); ++j ) {
                if( model_.surface(i).boundary( j )->nb_in_boundary() == 1 ) {
                    ++result ;
                    break ;
                }
            }
        }
        return result ;
    }

    int BoundaryModelMeasure::nb_non_manifold_lines() const
    {
        int result = 0 ;
        for( unsigned int i = 0; i < model_.nb_lines(); ++i ) {
            const Line& cp = model_.line(i) ;
            if( cp.nb_in_boundary() > 1 ) ++result ;
        }
        return result ;
    }
#endif
        

/*************************************************************************************************/
/********  BoundarymodelBuilder implementation      ******************************************/
/*************************************************************************************************/
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

    /*! Determine the type of a geological feature that is at the intersection
     *  of the input features
     */
    GEOL_FEATURE BoundaryModelBuilder::determine_type(
        const std::vector< GEOL_FEATURE >& types )
    {
        if( types.size() == 0 ) return ALL ;

        // Sort and remove duplicates form the in types
        std::vector< GEOL_FEATURE > in = types ;
        std::sort( in.begin(), in.end() ) ;
        int new_size = std::unique( in.begin(), in.end() ) - in.begin() ;
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

    void BoundaryModelBuilder::set_line( 
        unsigned int id, 
        const std::vector< uint32 >& vertices )
    {
        assert( id < model_.nb_lines() ) ;
        model_.lines_[id].points_ = vertices ;

        /*for( unsigned int p = 0; p < vertices.size(); p++ ) {
            model_.lines_[id].vertices_.push_back( model_.nb_points() ) ;
            add_point( vertices[p] ) ;
        }*/
    }
    void BoundaryModelBuilder::set_line_geometry( 
        unsigned int id, const std::vector< unsigned int >& line_points ) 
    {
        model_.lines_[id].set_vertices( line_points ) ;
    }


    unsigned int BoundaryModelBuilder::create_line(
        int id ,
        const std::vector< unsigned int >& points )
    {
        if( id == -1 ) id = model_.nb_lines() ;
        assert( id == model_.nb_lines() ) ;
        model_.lines_.push_back( Line( &model_, id, points ) ) ;
        return id ;
    }

    unsigned int BoundaryModelBuilder::create_surface(
        int id,
        int parent,
        GEOL_FEATURE type )
    {
        if( id == -1 ) id = model_.nb_surfaces() ;
        assert( id == model_.nb_surfaces() ) ;
        model_.surfaces_.push_back( Surface( &model_, id, parent, type ) ) ;
        return id ;
    }

     unsigned int BoundaryModelBuilder::create_region( int id ) 
     {
        if( id == -1 ) id = model_.regions_.size() ;
        assert( id == model_.regions_.size() ) ;
        model_.regions_.push_back( BoundaryModelElement( &model_, 3, id ) ) ;
        return id ;
    }


    void BoundaryModelBuilder::remove_universe_from_regions( unsigned int id ) 
    {
        for( int i = 0; i < model_.nb_regions(); ++i ) {
            int cur_id = model_.region(i).id() ;
            assert( i == cur_id ) ;
            if( i > id ) model_.regions_[i].set_id( cur_id-1 ) ;
        }
        model_.regions_.erase( model_.regions_.begin() + id ) ;
    }

    unsigned int BoundaryModelBuilder::create_layer( const std::string& name, int id ) 
    {
        if( id == -1 ) id = model_.layers_.size() ;
        assert( id == model_.layers_.size() ) ;
        model_.layers_.push_back( BoundaryModelElement( &model_, 3, id ) ) ;
        model_.layers_[id].set_name( name ) ;
        return id ;
    }
    
    void BoundaryModelBuilder::set_corner( unsigned int corner_id,  unsigned int point_id ) 
    {
        grgmesh_debug_assert( corner_id < model_.nb_corners() ) ;
        grgmesh_debug_assert( point_id < model_.nb_points() ) ;

        model_.corners_[corner_id].set_point( point_id ) ;
    }


    unsigned int BoundaryModelBuilder::create_interface(
        const std::string& name,
        int id,
        GEOL_FEATURE type )
    {
        if( id == -1 ) id = model_.nb_interfaces() ;
        assert( id == model_.nb_interfaces() ) ;
        model_.interfaces_.push_back( BoundaryModelElement( &model_, 2, id, -1, type ) ) ;
        model_.interfaces_[id].set_name( name ) ;
        return id ;
    }

    void BoundaryModelBuilder::copy_macro_topology( const BoundaryModel* from )
    {
        model_.corners_.resize( from->nb_corners(), Corner( &model_ ) ) ;
        model_.lines_.resize( from->nb_lines(), Line( &model_ ) ) ;
        model_.surfaces_.resize( from->nb_surfaces(), Surface( &model_ ) ) ;
        model_.regions_.resize( from->nb_regions(), BoundaryModelElement( &model_ ) ) ;
        model_.layers_.resize( from->nb_layers(), BoundaryModelElement( &model_ ) ) ;
        model_.contacts_.resize( from->nb_contacts(), BoundaryModelElement( &model_ ) ) ;
        model_.interfaces_.resize( from->nb_interfaces(), BoundaryModelElement( &model_ ) ) ;
#pragma omp parallel for
        for( unsigned int i = 0; i < model_.nb_corners(); i++ ) {
            model_.corners_[i].copy_macro_topology( from->corner( i ), model_ ) ;
        }
#pragma omp parallel for
        for( unsigned int i = 0; i < model_.nb_lines(); i++ ) {
            model_.lines_[i].copy_macro_topology( from->line( i ), model_ ) ;
        }
#pragma omp parallel for
        for( unsigned int i = 0; i < model_.nb_surfaces(); i++ ) {
            model_.surfaces_[i].copy_macro_topology( from->surface( i ), model_ ) ;
        }
#pragma omp parallel for
        for( unsigned int i = 0; i < model_.nb_layers(); i++ ) {
            model_.layers_[i].copy_macro_topology( from->layer( i ), model_ ) ;
        }
#pragma omp parallel for
        for( unsigned int i = 0; i < model_.nb_regions(); i++ ) {
            model_.regions_[i].copy_macro_topology( from->region( i ), model_ ) ;
        }
#pragma omp parallel for
        for( unsigned int i = 0; i < model_.nb_contacts(); i++ ) {
            model_.contacts_[i].copy_macro_topology( from->contact( i ), model_ ) ;
        }
#pragma omp parallel for
        for( unsigned int i = 0; i < model_.nb_interfaces(); i++ ) {
            model_.interfaces_[i].copy_macro_topology( from->one_interface( i ), model_ ) ;
        }
        model_.universe_.copy_macro_topology( from->universe_, model_ ) ;
    }

    /**
     * \brief Load and build a BoundaryModel from a Gocad .ml file
     * 
     * This is pretty tricky because of the annoying not well adapted file format.
     * 
     * The correspondance between Gocad::Model3D elements and BoundaryModel
     * elements is :
     * Gocad TSurf  <-> BoundaryModel Interface
     * Gocad TFace  <-> BoundaryModel Surface
     * Gocad Region <-> BoundaryModel Region
     * Gocad Layer  <-> BoundaryModel Layer
     *
     */
    void BoundaryModelBuilder::load_file( std::istream& in )
    {
        time_t start_load, end_load ;
        time( &start_load ) ;        

        // Count the number of TSurf - Interface
        int nb_tsurf = 0 ;
        // Count the number of TFace - Surface
        int nb_tface = 0 ;

        // Counters identifying the currently read TSurf or TFace
        int tsurf_count = 0 ;
        int tface_count = 0 ;

        int current_nb_tfaces = 0 ;
        int nb_tface_in_prev_tsurf = 0 ;

        // The file contains 2 parts and is read in 2 steps
        // 1. Read model info (true)
        // 2. Read surfaces geometry (false)
        bool read_model = true ;

        // The orientation of positive Z
        // can change for each TSurf and need to be read
        int z_sign = 1 ;

        // In the .ml file - vertices are indexed TSurf by Tsurf
        // They can be duplicated inside one TSurf and betweeen TSurfs\
        
        // Indices of the points of the currently built TSurf in the model
        std::vector< unsigned int > tsurf_point_ptr ;        
        // Where the points of a TFace start in the points of the TSurf (offest)
        std::vector< int > tface_point_start ;

        // Indices of points in facets (triangles) of the currently built TFace
        std::vector< unsigned int > tface_facets ;
        // Starting and ending indices of each facet triangle in the tface_facets vector
        std::vector< unsigned int > tface_facets_ptr ;
        tface_facets_ptr.push_back( 0 ) ;

        
        // Intermediate information for contact parts building
        std::vector< Border > borders_to_build ;

        // Surfaces for which the KeyFacet orientation should be changed
        // because it does not match triangle orientations.
        std::vector< int > change_key_facet ;

        /****** File reading **********************************/
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
                    int id ;
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
                        KeyFacet( p0, p1, p2 ) ) ;
                    nb_tface++ ;
                } else if( keyword == "REGION" ) {
                    /// 3. Read Region information and create them from their name,
                    /// the surfaces on their boundary                    
                    int id ;
                    std::string name ;
                    lis >> id >> name ;

                    std::vector< std::pair< int, bool > > region_boundaries ;
                    bool end_region = false ;

                    while( !end_region ) {
                        lis.get_line() ;
                        for( unsigned int i = 0; i < 5; ++i ) {
                            int tface_id ;
                            lis >> tface_id ;
                            if( tface_id == 0 ) {
                                end_region = true ;
                                break ;
                            } else {
                                // Correction because ids begin at 1 in the file
                                tface_id = tface_id > 0 ? tface_id-1 : -tface_id-1 ;                                            
                                region_boundaries.push_back(
                                    std::pair< int, bool >( tface_id, false ) ) ;
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
                    unsigned int layer_id = create_layer( name ) ;
                    bool end_layer = false ;
                    while( !end_layer ) {
                        lis.get_line() ;
                        for( unsigned int i = 0; i < 5; ++i ) {
                            int region_id ;
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
                            std::vector< unsigned int >(
                                tsurf_point_ptr.begin() + tface_point_start.back(),
                                tsurf_point_ptr.end() ),
                            tface_facets,
                            tface_facets_ptr ) ;

                        if( !check_key_facet_orientation( tface_count-1 ) ){
                            change_key_facet.push_back( tface_count-1 ) ;
                        }

                        tface_facets.clear() ;
                        tface_facets_ptr.clear() ;
                        tface_facets_ptr.push_back( 0 ) ;

                        // End this TSurf - Interface
                        nb_tface_in_prev_tsurf += tface_point_start.size() ;
                        tsurf_point_ptr.clear() ;
                        tface_point_start.clear() ;
                    }
                } else if( keyword == "TFACE" ) {
                    // Beginning of a new TFace - Surface
                    if( tface_point_start.size() > 0 ) {
                        // End the previous TFace - Surface
                        set_surface_geometry( 
                            tface_count-1,
                            std::vector< unsigned int >(
                                tsurf_point_ptr.begin() + tface_point_start.back(),
                                tsurf_point_ptr.end() ),
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
                    tface_point_start.push_back( tsurf_point_ptr.size() ) ;

                    tface_count++ ;
                }
                /// 4. Read the surface points and facets (only triangles in Gocad Model3d files)
                else if( keyword == "VRTX" || keyword == "PVRTX" ) {
                    int id ;
                    vec3 p ;
                    lis >> id >> p ;
                    p.z *= z_sign ;
                    tsurf_point_ptr.push_back( add_point( p ) ) ;
                } else if( keyword == "PATOM" || keyword == "ATOM" ) {
                    // This keyword is used to refer to a previous VERTEX in the
                    // same TSurf
                    int id ;
                    int v_id ;
                    lis >> id >> v_id ;
                    tsurf_point_ptr.push_back( tsurf_point_ptr[v_id - 1] ) ;
                } else if( keyword == "TRGL" ) {
                    // Ids of the vertices of each triangle in the TSurf
                    int p1, p2, p3 ;
                    lis >> p1 >> p2 >> p3 ;
                    // Change to ids in the TFace
                    p1 += -tface_point_start.back()-1 ;
                    p2 += -tface_point_start.back()-1 ;
                    p3 += -tface_point_start.back()-1 ;

                    tface_facets.push_back( p1 ) ;
                    tface_facets.push_back( p2 ) ;
                    tface_facets.push_back( p3 ) ;
                    tface_facets_ptr.push_back( tface_facets.size() ) ;
                }
                /// 5. Build the corners from their position and the surface parts
                ///    containing them
                else if( keyword == "BSTONE" ) {
                    int v_id ;
                    lis >> v_id ;
                    // correction to start at 0
                    v_id-- ;

                    // Get the TFace
                    int part_id = tface_point_start.size() - 1 ;
                    for( int i = 0; i < tface_point_start.size(); ++i ) {
                        if( v_id < tface_point_start[i] ) {
                            part_id = i - 1 ;
                            break ;
                        }
                    }
                    part_id += nb_tface_in_prev_tsurf ;

                    int new_c = find_or_create_corner( tsurf_point_ptr[v_id] ) ;
                }
                /// 6. Read the Border information and store it
                else if( keyword == "BORDER" ) {
                    int id, p1, p2 ;
                    lis >> id >> p1 >> p2 ;
                    p1-- ;
                    p2-- ;

                    // Get the global corner id
                    int corner_id = find_corner( model_.point( tsurf_point_ptr[p1] ) ) ;
                    grgmesh_assert( corner_id > -1 ) ;

                    // Get the surface
                    int part_id = -1 ;
                    for( int i = 0; i < tface_point_start.size(); ++i ) {
                        if( p1 < tface_point_start[i] ) {
                            grgmesh_assert( p2 < tface_point_start[i] ) ;

                            // Get points ids in the surface
                            p1 += -tface_point_start[i - 1] ;
                            p2 += -tface_point_start[i - 1] ;

                            // i-1 is the id of the TFace in this TSurf
                            part_id = i - 1 ;

                            break ;
                        }
                    }
                    if( part_id == -1 ) {
                        // It is in the last built Tface
                        p1 += -tface_point_start[tface_point_start.size() - 1] ;
                        p2 += -tface_point_start[tface_point_start.size() - 1] ;
                        part_id = tface_point_start.size() - 1 ;
                    }
                    // The number of tfaces in previous tsurf is also to add
                    part_id += nb_tface_in_prev_tsurf ;

                    borders_to_build.push_back(
                        Border( part_id, corner_id, p1, p2 ) ) ;
                }
            }
        }

        make_points_unique() ;

        /// 7. Build the contact parts
        build_lines( borders_to_build ) ;

        /// 8. Build the contacts
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
        std::cout << "Info" << " Boundary model loading time"
            << difftime( end_load, start_load ) << " sec" << std::endl ;        
    }

    void BoundaryModelBuilder::set_surface_geometry(
        unsigned int surface_id,
        const std::vector< unsigned int >& points,
        const std::vector< unsigned int >& facets,
        const std::vector< unsigned int >& facet_ptr ) 
    {        
        if( facets.size() == 0 ) return ;

        model_.surfaces_[surface_id].set_geometry( points, facets, facet_ptr ) ;

        set_surface_adjacencies( surface_id ) ;         
    }


     /** Returns the id of the facet which first three points are those given */
    int BoundaryModelBuilder::find_key_facet( 
        uint32 surface_id, const vec3& p0, const vec3& p1,
        const vec3& p2, bool& same_sign ) const 
    {
        const Surface& surface = model_.surface( surface_id ) ;
        same_sign = false ;
    
        for( unsigned int t = 0; t < surface.nb_cells(); ++t ){            
            const vec3& pp0 = model_.point( surface.model_point_id( t, 0 ) ) ;
            const vec3& pp1 = model_.point( surface.model_point_id( t, 1 ) ) ;
            const vec3& pp2 = model_.point( surface.model_point_id( t, 2 ) ) ;
            
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
        return -1 ;
    }


    bool BoundaryModelBuilder::check_key_facet_orientation( uint32 surface_id ) const {

        Surface& surface = model_.surfaces_[surface_id] ;
        KeyFacet& key_facet = surface.key_facet_ ;

        if( key_facet.is_default() ) {
            surface.set_first_triangle_as_key() ;
            return true ;
        }
        
        vec3& p0 = key_facet.p0_ ;
        vec3& p1 = key_facet.p1_ ;
        vec3& p2 = key_facet.p2_ ;
        bool same_sign = false ;

        int t = find_key_facet( surface_id, p0,p1,p2, same_sign ) ;
        if( t == -1 ) {
            // It is because of the sign of Z that is not the same 
            p0.z *= -1 ;
            p1.z *= -1 ;
            p2.z *= -1 ;
            t = find_key_facet( surface_id, p0, p1, p2, same_sign ) ; 
        }
        grgmesh_assert( t > -1 ) ;

        return same_sign ;
   }


    /** 
     *  Compute and set the adjacencies between the facets 
     *   The adjacent facet is given for each vertex of each facet for the edge
     *   starting at this vertex.
     *   If no neighbor inside the same Surface adjacent is -1
     */
    void BoundaryModelBuilder::set_surface_adjacencies( unsigned int surface_id ) { 

        Surface& S = model_.surfaces_[surface_id] ;
        std::vector< int >& adjacent = S.adjacent_ ;
        adjacent.resize( S.facets_.size(), -1 ) ;
       
        grgmesh_assert( S.facets_.size() > 0  ) ;

        unsigned int nb_facets = S.nb_cells() ;
        unsigned int nb_points = S.nb_points() ;        
    
        // Allocate some space to store the ids of facets around each point
        std::vector< unsigned int > toto ;
        toto.reserve( 10 ) ;
        std::vector< std::vector< unsigned int > > 
            point_to_facets( nb_points, toto ) ;              

        for( unsigned int f = 0; f < nb_facets; ++f )
        {
            for( unsigned int v = 0; v < S.nb_points_in_facet( f ); v++ ) {
                point_to_facets[S.surf_point_id( f, v )].push_back( f ) ;
            }
        }
        for( unsigned int p = 0; p < nb_points; ++p ){
            std::sort( point_to_facets[p].begin(), point_to_facets[p].end() ) ;
        }

        for( unsigned int f = 0; f < nb_facets; ++f )
        {
            for( unsigned int v = 0; v < S.nb_points_in_facet(f); v++ )
            {                
                unsigned int cur = S.surf_point_id(f, v) ;
                unsigned int prev = S.surf_point_id(f, S.prev_in_facet(f,v)) ;

                const std::vector< unsigned int >& f_prev = point_to_facets[prev] ;
                const std::vector< unsigned int >& f_cur = point_to_facets[cur] ;

                std::vector< unsigned int > inter( std::min( f_prev.size(), f_cur.size() ) ) ;
                int end = std::set_intersection(
                    f_prev.begin(), f_prev.end(), f_cur.begin(), f_cur.end(), inter.begin() ) - inter.begin() ;

                if( end == 2 ) {
                    // There is one neighbor
                    unsigned int f2 = inter[0] == f ? inter[1] : inter[0] ;
                    adjacent[ S.facet_begin(f) + S.prev_in_facet(f,v) ] = f2 ;
                } else {
                    grgmesh_debug_assert( end == 1 ) ;
                }
            }
        }
    }

    void BoundaryModelBuilder::cut_surface_by_line( uint32 surface_id, uint32 line_id ) {

        Surface& S = model_.surfaces_[surface_id] ;
        const Line& L = model_.line( line_id ) ;


        //Number of segments = nb_cells 
        for( uint32 i= 0; i+1 < L.nb_points(); ++i ) {
            uint32 p0 = L.model_point_id( i ) ;
            uint32 p1 =  (i == L.nb_points()-1) ? L.model_point_id(0) : L.model_point_id( i+1 ) ;

            int f = -1 ;
            int v = -1 ;
            S.edge_from_model_point_ids(p0, p1, f, v) ;
            grgmesh_debug_assert( f != -1 && v != -1 ) ;

            int f2 = S.adjacent( f, v ) ;
            int v2 = -1 ;
            grgmesh_debug_assert( f2 != -1 ) ;
            S.edge_from_model_point_ids( p0, p1, f2, v2 ) ;
            grgmesh_debug_assert( v2 != -1 ) ;

            // Virtual cut - set adjacencies to -1 
            S.set_adjacent( f, v, -1 ) ;
            S.set_adjacent( f2, v2, -1 ) ;       
        }
        
        // Now travel on one side of the "faked" boudnary and actually duplicate
        // the point in the surface
        // Corners are not duplicated - maybe theu should be in some cases but not in general..


        // Get started in the surface - find (again) one of the edge that contains 
        // the first two points of the line
        int f = -1 ;
        int v = -1 ;
        S.oriented_edge_from_model_point_ids( L.model_point_id( 0 ), L.model_point_id( 1 ), f, v ) ;
        grgmesh_assert( f != -1 && v != -1 ) ;

        uint32 id0 = S.surf_point_id( f, v ) ;
        uint32 id1 = S.surf_point_id( f, S.next_in_facet(f,v) ) ;
        
        // Stopping criterion
        uint32 last_point = L.model_point_id( L.nb_points()-1 ) ;
        int count = 0 ;
        // On espère qu'on les récupère bien tous les points sur la ligne... A VERIFIER
        while( S.model_point_id( id1 ) != last_point ) {

            // Get the next point on the border 

            // Same algorithm than in determine_line_vertices function
            int next_f = -1 ;
            int id1_in_next = -1 ;
            int next_id1_in_next = -1 ;

            // Get the next facet and next triangle on this boundary
            S.next_on_border( f, S.facet_point_id(f, id0), S.facet_point_id(f,id1), 
                next_f, id1_in_next, next_id1_in_next ) ;
            grgmesh_assert( next_f != -1 && id1_in_next != -1 && next_id1_in_next != -1 ) ;
            
            int next_id1 = S.surf_point_id( next_f, next_id1_in_next ) ;
            

            // Duplicate the point at id1 
            // Après avoir récupéré le suivant 1 .. on doit pouvoir les deux en même 
            // temps mais là j'ai la flemme de réfléchir, faut pas que ça casse le next_on_border (jeanne)
            std::vector< uint32 > facets_around_id1 ;
            S.facets_around_point( id1, facets_around_id1, false, f ) ;

            S.points_.push_back( S.model_point_id(id1) ) ;
            int new_id1 = S.nb_points()-1 ;
            
            for( uint32 i = 0; i < facets_around_id1.size(); ++i ){
                uint32 cur_f = facets_around_id1[i] ;
                for( uint32 cur_v = 0; cur_v < S.nb_points_in_facet( cur_f ) ; cur_v++ )
                {
                    if( S.surf_point_id( cur_f, cur_v ) == id1 ) {
                        S.facets_[ S.facet_begin( cur_f ) + cur_v ] = new_id1 ;
                        break ;
                    }
                }
            }

            // Update
            f = next_f ;
            id0 = new_id1 ;
            id1 = next_id1 ;

            ++count ; // debug
        }


    }


    /**
     * When reading the file the points are duplicated between the different Surfaces
     * plus new points are added for Corners too 
     *
     * Compute the duplicates inside the points_ vector - Update the point vector - 
     * and Update the reference to points for all Corner and Surface plus for the input 
     * borders that are later used to build the model lines 
     */
    void BoundaryModelBuilder::make_points_unique()
    {
        std::vector< vec3 > new_points ;
        std::vector< uint32 > old2new ;
        compute_unique_kdtree( model_.points_, 10, old2new, new_points ) ;
       
        model_.points_.resize(0) ;
        model_.points_ = new_points ;

        for( unsigned int s = 0; s < model_.nb_surfaces(); s++ ) {
            Surface& surface = model_.surfaces_[s] ;
            for( unsigned int p = 0; p < surface.nb_points(); p++ ) {
                surface.points_[p] = old2new[surface.points_[p]] ;
            }
        }
        for( unsigned int co = 0; co < model_.nb_corners(); co++ ) {
            Corner& corner = model_.corners_[co] ;
            corner.p_ = old2new[corner.p_] ;
        }       
    }

    void BoundaryModelBuilder::set_universe(
        const std::vector< std::pair< int, bool > >& boundaries )
    {
        model_.universe_.set_name( "Universe" ) ;
        model_.universe_.set_dim( 3 ) ;

        for( unsigned int i = 0; i < boundaries.size(); ++i ) {
            assert( boundaries[i].first < model_.nb_surfaces() ) ;
            model_.universe_.add_boundary( boundaries[i].first,
                boundaries[i].second ) ;

            // If this surface have no type, set it at VOI
            model_.surfaces_[boundaries[i].first].set_type( VOI ) ;
        }
    }

    unsigned int BoundaryModelBuilder::create_region(
        const std::string& name,
        const std::vector< std::pair< int, bool > >& boundaries,
        int id )
    {
        id = create_region( id ) ;
        model_.regions_[id].set_name( name ) ;

        for( unsigned int i = 0; i < boundaries.size(); ++i ) {
            assert( boundaries[i].first < model_.nb_surfaces() ) ;
            add_region_oriented_boundary( id, boundaries[i].first,
                boundaries[i].second ) ;
        }
        return id ;
    }


    /**
     * Get the points that are on the border of the given surface, starting from 
     * point id0 in the surface towards the direction given by point id1
     *
     * WE ASSUME THAT THE STORAGE OF THE POINTS IS UNIQUE IN THE MODEL AND THAT 
     * SURFACES DO SHARE POINTS ON THEIR CONTACT LINES
     * 
     * Empties and fills the boder_vertex_model_ids with the global ids of these points
     * Returns the id of the corner at which this boundary stops
     */
    unsigned int BoundaryModelBuilder::determine_line_vertices( 
        const Surface& S, 
        unsigned int id0, 
        unsigned int id1,
        std::vector< unsigned int >& border_vertex_model_ids 
        ) const 
    {
        grgmesh_debug_assert( id0 < S.nb_points() && id1 < S.nb_points() ) ;

        border_vertex_model_ids.resize( 0 ) ;
                 
        // Starting facet that contains the two given points
        int f = S.facet_from_surface_point_ids( id0, id1 ) ;
        grgmesh_assert( f != -1 ) ;

        // Global ids at the model level 
        uint32 p0 = S.model_point_id( id0 ) ;
        uint32 p1 = S.model_point_id( id1 ) ;

        border_vertex_model_ids.push_back( p0 ) ;
        border_vertex_model_ids.push_back( p1 ) ;
            
        int p1_corner = find_corner( p1 ) ;
        while( p1_corner == -1 ) {

            int next_f = -1 ;
            int id1_in_next = -1 ;
            int next_id1_in_next = -1 ;

            // We want the next triangle that is on the boundary and share p1
            // If there is no such triangle, the third point of the current triangle is to add
            S.next_on_border( f, S.facet_point_id(f, id0), S.facet_point_id(f,id1), 
                next_f, id1_in_next, next_id1_in_next ) ;

            grgmesh_assert( next_f != -1 && id1_in_next != -1 && next_id1_in_next != -1 ) ;
            
            int next_id1 =  S.surf_point_id( next_f, next_id1_in_next ) ;

            // Update
            f = next_f ;
            id0 = id1 ;
            id1 = next_id1 ;

            p1 = S.model_point_id( next_id1 ) ;
            border_vertex_model_ids.push_back( p1 ) ;         
            p1_corner = find_corner( p1 ) ;
        }
        return p1_corner ; 
    }

    /*! Build the Lines once the storage of the points in the Model as been make unique
     * So that POINTS on lines are shared by the surfaces in contact there.
     *  from the information collected in the Border structures
     *
     *  fill their in_boundary_ vector
     */
    void BoundaryModelBuilder::build_lines( const std::vector< Border >& borders )
    {      
        std::vector< unsigned int > global_ids ;

        for( unsigned int i = 0; i < borders.size(); ++i ) {
            const Border& b = borders[i] ;

            /// 1- Build the boundary : construct the vector
            /// of vertices on the border
            const Surface& S = model_.surface( b.part_id_) ;

            unsigned int end_corner_id = determine_line_vertices( 
                S, b.p0_, b.p1_, global_ids ) ;           

            /// 2 - Check if this border already exists
            int line_id = find_or_create_line( b.corner_id_,
                end_corner_id, global_ids ) ;

            // Add the surface in which this line is
            add_line_in_boundary( line_id, b.part_id_ ) ;
        }
    }

    /** 
     * \brief Build the contacts - group of lines shared by the same interfaces
     *
     */
    void BoundaryModelBuilder::build_contacts()
    {
        for( unsigned int i = 0; i < model_.nb_lines(); ++i ) {
            // The surface part in whose boundary is the part
            std::set< int > interfaces ;
            std::vector< GEOL_FEATURE > types ;
            for( unsigned int j = 0; j < model_.lines_[i].nb_in_boundary(); ++j ) {
                unsigned int sp_id = model_.lines_[i].in_boundary_id( j ) ;
                const BoundaryModelElement& p = model_.surfaces_[sp_id].parent() ;
                interfaces.insert( p.id() ) ;
                types.push_back( p.type() ) ;
            }
            std::vector< int > toto( interfaces.begin(), interfaces.end() ) ;
            int contact_id = find_or_create_contact( toto, determine_type( types ) ) ;
            add_contact_child( contact_id, i ) ;
        }
    }

    /*! Set the surfaces delimited by each contact
     *  Warning the surfaces must be finished first
     */
    void BoundaryModelBuilder::end_contacts()
    {
        for( unsigned int i = 0; i < model_.nb_contacts(); ++i ) {
            std::set< const BoundaryModelElement* > corners ;
            for( unsigned int j = 0; j < model_.contacts_[i].nb_children(); ++j ) {
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

    /*! Set the parent contact for all the contact parts
     */
    void BoundaryModelBuilder::end_lines()
    {
        for( unsigned int i = 0; i < model_.nb_contacts(); ++i ) {
            for( unsigned int j = 0; j < model_.contacts_[i].nb_children(); ++j ) {
                unsigned int child = model_.contacts_[i].child_id( j ) ;
                model_.lines_[child].set_parent( i ) ;
                model_.lines_[child].set_type( model_.contacts_[i].type() ) ;
            }
        }
    }

    /**
     * \todo End interface function not correct, in_buondary not set
     */
    void BoundaryModelBuilder::end_interfaces()
    {
        for( unsigned int i = 0; i < model_.nb_surfaces(); ++i ) {
            int parent = model_.surfaces_[i].parent_id() ;
            add_interface_child( parent, i ) ;
        }

        for( unsigned int i = 0; i < model_.nb_contacts(); ++i ) {
            for( unsigned int j = 0; j < model_.contacts_[i].nb_in_boundary(); ++j ) {
                unsigned int b = model_.contacts_[i].in_boundary_id( j ) ;
                add_interface_boundary( b, i ) ;
            }
        }       
    }
    /**
     * \todo What do we want to have in the boundaries of layers ? surfaces or interfaces ?
     */
    void BoundaryModelBuilder::end_layers()
    {
        for( unsigned int i = 0; i < model_.nb_layers(); ++i ) {
            BoundaryModelElement& layer = model_.layers_[i] ;
            std::set< unsigned int > boundaries ;

            for( unsigned int r = 0; r < layer.nb_children(); r++ ) {
                const BoundaryModelElement& region = layer.child( r ) ;
                model_.regions_[region.id()].parent_ = i ;
                for( unsigned int sp = 0; sp < region.nb_boundaries(); sp++ ) {
                    boundaries.insert( region.boundary_id( sp ) ) ;
                }
            }

            for( std::set< unsigned int >::const_iterator it( boundaries.begin() );
                it != boundaries.end(); it++ ) {
                layer.add_boundary( *it ) ;
            }
        }
    }

    /*! Set the parent of surfaces
     */
    void BoundaryModelBuilder::end_surfaces()
    {
        for( unsigned int i = 0; i < model_.nb_lines(); ++i ) {
            for( unsigned int j = 0; j < model_.lines_[i].nb_in_boundary();
                ++j ) {
                unsigned int s_id = model_.lines_[i].in_boundary_id( j ) ;
                add_surface_boundary( s_id, i ) ;
            }
        }
        for( unsigned int i = 0; i < model_.nb_regions(); ++i ) {
            for( unsigned int j = 0; j < model_.regions_[i].nb_boundaries(); ++j ) {
                unsigned int s_id = model_.regions_[i].boundary_id( j ) ;
                add_surface_in_boundary( s_id, i ) ;
            }
        }
    }

    void BoundaryModelBuilder::end_surfaces(
        const std::vector< int >& change_orientation )
    {
        end_surfaces() ;
        for( int i = 0; i < change_orientation.size(); i++ ) {
            uint32 s_i = change_orientation[i] ; 
            
            // Change the key facet
            Surface& surface = model_.surfaces_[s_i] ;
            surface.set_first_triangle_as_key() ;

            // Change the sign of this Surface in all the regions containing it
            for( uint32 j = 0; j < surface.nb_in_boundary(); ++j ) {                
                BoundaryModelElement& region = model_.regions_[ surface.in_boundary_id( j ) ] ;              
                for( unsigned int b = 0; b < region.nb_boundaries(); ++b ){
                    if( region.boundary(b).id() == s_i ){
                        region.sides_[i] = !region.sides_[i] ;
                    }
                }
            }
        }
    }

    void BoundaryModelBuilder::end_model()
    {
        model_.nb_facets_.resize( model_.nb_surfaces()+1, 0 ) ;

        uint32 count = 0 ;
        for( uint32 i = 1; i < model_.nb_facets_.size(); ++i ) {
            count += model_.surface( i-1 ).nb_cells() ;
            model_.nb_facets_[i] = count ;
        }

        std::cout << "Model " << model_.name() <<" has " << std::endl 
            << std::setw(10) << std::left << model_.nb_points()   << " points "   << std::endl 
            << std::setw(10) << std::left << model_.nb_facets()   << " facets "   << std::endl  
            << std::setw(10) << std::left << model_.nb_regions()  << " regions "  << std::endl
            << std::setw(10) << std::left << model_.nb_surfaces() << " surfaces " << std::endl
            << std::setw(10) << std::left << model_.nb_lines()    << " lines "    << std::endl 
            << std::setw(10) << std::left << model_.nb_corners()  << " corners "  << std::endl
            << std::endl ;
    }                                                                          


    /*! Set in what surface parts are the corners
     */
    void BoundaryModelBuilder::end_corners()
    {
        for( unsigned int i = 0; i < model_.nb_lines(); ++i ) {
            unsigned int c0_id = model_.lines_[i].boundary_id( 0 ) ;
            unsigned int c1_id = model_.lines_[i].boundary_id( 1 ) ;

            add_corner_in_boundary( c0_id, i ) ;
            if( c1_id != c0_id ) add_corner_in_boundary( c1_id, i ) ;
        }
    }

    int BoundaryModelBuilder::find_or_create_contact(
        std::vector< int >& interfaces,
        GEOL_FEATURE type )
    {
        int result = find_contact( interfaces ) ;
        if( result == -1 ) {
            // Create a name for this contact
            std::string name = "contact_" ;
            for( unsigned int i = 0; i < interfaces.size(); ++i ) {
                name += model_.interfaces_[interfaces[i]].name() ;
                name += "_" ;
            }
            result = model_.nb_contacts() ;
            model_.contacts_.push_back( BoundaryModelElement( &model_, 1, result ) ) ;
            model_.contacts_[result].set_name( name ) ;
            model_.contacts_[result].set_type( type ) ;

            for( unsigned int i = 0; i < interfaces.size(); ++i ) {
                add_contact_in_boundary( result, interfaces[i] ) ;
            }
        }
        return result ;
    }

    /*! Add a surface part to the model
     *  The parent with the given name (one_interface name) MUST exist
     */
    void BoundaryModelBuilder::create_surface(
        const std::string& interface_name,
        const std::string& type,
        const KeyFacet& key )
    {
        int parent = interface_id( interface_name ) ;
        assert( parent != -1 ) ;

        int id = model_.nb_surfaces() ;
        GEOL_FEATURE t = determine_geological_type( type ) ;
        model_.surfaces_.push_back( Surface( &model_, id, parent, t ) ) ;
        model_.surfaces_[id].set_key_facet( key ) ;

        model_.interfaces_[parent].set_type( t ) ;
    }

    int BoundaryModelBuilder::find_or_create_corner( unsigned int index )
    {
        int result = find_corner( model_.point( index ) ) ;
        if( result == -1 ) {
            // Create the corner
            result = model_.nb_corners() ;
            model_.corners_.push_back( Corner( &model_, result, index ) ) ;
        }
        return result ;
    }
    
    int BoundaryModelBuilder::find_or_create_line(
        int corner0,
        int corner1,
        std::vector< uint32 >& points )
    {
        int result = find_line( corner0, corner1, points ) ;

        if( result == -1 ) {
            result = model_.nb_lines() ;
            grgmesh_assert( corner0 != -1 ) ;
            grgmesh_assert( corner1 != -1 ) ; 

        /*    if( corner1 == corner0 ) {
                // Closed contact part
                //corner1 = corner0 ;
                // Remove the last point
                grgmesh_assert( points[0] == points.back() ) ;
                points.pop_back() ;
            }*/

            /*std::vector< unsigned int > indices( points.size() ) ;
            for( unsigned int p = 0; p < points.size(); p++ ) {
                indices[p] = model_.nb_points() ;
                add_point( points[p] ) ;
            }*/

            model_.lines_.push_back(
                Line( &model_, result, corner0, corner1, points ) ) ;
        }
        return result ;
    }
  
    int BoundaryModelBuilder::interface_id( const std::string& name ) const
    {
        for( unsigned int i = 0; i < model_.nb_interfaces(); ++i ) {
            if( model_.one_interface(i).name() == name ) {
                return i ;
            }
        }
        return -1 ;
    }


    int BoundaryModelBuilder::find_corner( const vec3& p ) const
    {
        for( uint32 i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).point() == p ) return i ;
        }
        return -1 ;
    }

    int BoundaryModelBuilder::find_corner( uint32 p_id ) const 
    {
        for( uint32 i = 0; i < model_.nb_corners(); ++i ) {
            if( model_.corner(i).p_ == p_id ) return i ;
        }
        return -1 ;
    }

    int BoundaryModelBuilder::find_line(
        int corner0,
        int corner1,
        const std::vector< uint32 >& points ) const
    {
        for( unsigned int i = 0; i < model_.nb_lines(); ++i ) {
            const Line& cp = model_.line(i) ;

            // If the number of points are not the same,
            // it is not the same Line
            if( cp.nb_points() != points.size() ) continue ; 

            if( corner0 == cp.boundary_id( 0 )
                && corner1 == cp.boundary_id( 1 ) ) {
               // bool equal = true ;
                if( std::equal( points.begin(), points.end(), cp.points_.begin() ) ){
                    return i ;
                }
                /*for( unsigned int p = 0; p < cp.nb_points(); p++ ) {
                    if( cp.point( p ) != points[p] ) {
                        equal = false ; break ;
                    }
                }
                if( equal ) return i ;*/
            }

            if( corner1 == cp.boundary_id( 0 )
                && corner0 == cp.boundary_id( 1 ) ) {
                //bool equal = true ;
                if( std::equal( points.begin(), points.end(), cp.points_.rbegin() ) ){
                    return i ;
                }
                /*for( unsigned int p = 0; p < cp.nb_points(); p++ ) {
                    if( cp.point( cp.nb_points() - p - 1 ) != points[p] ) {
                        equal = false ; break ;
                    }
                }
                if( equal ) return i ;*/
            }
        }
        return -1 ;
    }

    int BoundaryModelBuilder::find_contact( const std::vector< int >& interfaces ) const
    {
        std::vector< const BoundaryModelElement* > comp( interfaces.size() ) ;
        for( unsigned int i = 0; i < interfaces.size(); ++i ) {
            comp[i] = &model_.one_interface(interfaces[i]) ;
        }

        for( unsigned int i = 0; i < model_.nb_contacts(); ++i ) {
            if( comp.size() == model_.contact(i).nb_in_boundary() ) {
                bool equal = true ;
                for( unsigned int j = 0; j < model_.contact(i).nb_in_boundary(); j++ ) {
                    if( comp[j] != &model_.contact(i).in_boundary( j ) ) {
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


    bool BoundaryModelBuilder::rebuild()
    {
        /*! \todo Valgrind finds errors !!!!!!! */
        std::vector< unsigned int > sp_to_remove ;
        std::vector< unsigned int > new_sp_id( model_.nb_surfaces() ) ;
        unsigned int offset = 0 ;
        for( unsigned int sp = 0; sp < model_.nb_surfaces(); sp++ ) {
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
        std::vector< unsigned int > cp_to_remove ;
        std::vector< unsigned int > new_cp_id( model_.nb_lines() ) ;
        for( unsigned int cp = 0; cp < model_.nb_lines(); cp++ ) {
            BoundaryModelElement& line = model_.lines_[cp] ;
            unsigned int nb_sp_removed = 0 ;
            for( unsigned int sp = 0; sp < line.nb_in_boundary(); sp++ ) {
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
        std::vector< unsigned int > s_to_remove ;
        std::vector< unsigned int > new_s_id( model_.nb_interfaces() ) ;
        for( unsigned int s = 0; s < model_.nb_interfaces(); s++ ) {
            BoundaryModelElement& surface = model_.interfaces_[s] ;
            surface.boundaries_.clear() ;
            unsigned int nb_sp_removed = 0 ;
            for( unsigned int sp = 0; sp < surface.nb_children(); sp++ ) {
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
        for( unsigned int sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            BoundaryModelElement& surface = model_.surfaces_[sp] ;
            surface.parent_ = new_s_id[surface.parent_] ;
            for( unsigned int cp = 0; cp < surface.nb_boundaries(); cp++ ) {
                surface.boundaries_[cp] =
                    new_cp_id[surface.boundaries_[cp]] ;
            }
        }
        model_.contacts_.clear() ;
        build_contacts() ;

        // Then finish the job (the order matters)
        for( unsigned int sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            model_.surfaces_[sp].boundaries_.clear() ;
        }
        offset = 0 ;
        for( unsigned int r = 0; r < model_.nb_regions(); r++ ) {
            BoundaryModelElement& region = model_.regions_[r] ;
            unsigned int nb_sp_removed = 0 ;
            for( unsigned int sp = 0; sp < region.nb_boundaries(); sp++ ) {
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

        for( unsigned int i = 0; i < model_.nb_contacts(); ++i ) {
            for( unsigned int j = 0; j < model_.contacts_[i].nb_in_boundary();
                ++j ) {
                unsigned int b = model_.contacts_[i].in_boundary_id( j ) ;
                add_interface_boundary( b, i ) ;
            }
        }

        end_lines() ;
        end_contacts() ;

        offset = 0 ;
        std::vector< unsigned int > co_to_remove ;
        std::vector< unsigned int > new_co_id( model_.nb_corners() ) ;
        for( unsigned int co = 0; co < model_.nb_corners(); co++ ) {
            BoundaryModelElement& corner = model_.corners_[co] ;
            unsigned int nb_cp_removed = 0 ;
            for( unsigned int cp = 0; cp < corner.nb_in_boundary(); cp++ ) {
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

        for( unsigned int cp = 0; cp < model_.nb_lines(); cp++ ) {
            BoundaryModelElement& line = model_.lines_[cp] ;
            line.boundaries_[0] = new_co_id[line.boundaries_[0]] ;
            line.boundaries_[1] = new_co_id[line.boundaries_[1]] ;
        }

        for( unsigned int c = 0; c < model_.nb_contacts(); c++ ) {
            BoundaryModelElement& contact = model_.contacts_[c] ;
            for( unsigned int co = 0; co < contact.nb_boundaries(); co++ ) {
                contact.boundaries_[co] = new_co_id[contact.boundaries_[co]] ;
            }
        }

        for( unsigned int l = 0; l < model_.nb_layers(); l++ ) {
            BoundaryModelElement& layer = model_.layers_[l] ;
            unsigned int nb_sp_removed = 0 ;
            for( unsigned int sp = 0; sp < layer.nb_boundaries(); sp++ ) {
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

    void BoundaryModelBuilder::update_all_ids()
    {
        for( unsigned int co = 0; co < model_.nb_corners(); co++ ) {
            model_.corners_[co].id_ = co ;
        }
        for( unsigned int cp = 0; cp < model_.nb_lines(); cp++ ) {
            model_.lines_[cp].id_ = cp ;
        }
        for( unsigned int sp = 0; sp < model_.nb_surfaces(); sp++ ) {
            model_.surfaces_[sp].id_ = sp ;
        }
        for( unsigned int c = 0; c < model_.nb_contacts(); c++ ) {
            model_.contacts_[c].id_ = c ;
        }
        for( unsigned int s = 0; s < model_.nb_interfaces(); s++ ) {
            model_.interfaces_[s].id_ = s ;
        }
        for( unsigned int r = 0; r < model_.nb_regions(); r++ ) {
            model_.regions_[r].id_ = r ;
        }
        for( unsigned int l = 0; l < model_.nb_layers(); l++ ) {
            model_.layers_[l].id_ = l ;
        }
    }
/*************************************************************************************************/
/********  End BoundarymodelBuilder implementation      ******************************************/
/*************************************************************************************************/


} // namespace 0
