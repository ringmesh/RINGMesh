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

/*! \author Jeanne Pellerin */

#include <ringmesh/boundary_model_builder.h>
#include <ringmesh/utils.h>

#include <geogram/basic/line_stream.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_preprocessing.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <set>
#include <stack>

namespace {
    using namespace RINGMesh;

    typedef BoundaryModelElement::bme_t bme_t;
    typedef BoundaryModelMeshElement BMME ;
    typedef BoundaryModelVertices::VertexInBME VBME ;    

    double read_double( GEO::LineInput& in, index_t field )
    {
        double result ;
        std::istringstream iss( in.field( field ) ) ;
        iss >> result >> std::ws ;
        return result ;
    }

    /*!
     * @brief Fill the geological info from parent or children
     *
     * @details Set the geological feature of a BME to the one of its parent
     * if it has a parent that has a geol. feature,
     * or to the one of its first child if it has one with a geol. feature.
     *
     * @note The geol feature of \b E whatever its initial value
     */
    void fill_element_geological_feature( BoundaryModelElement& E )
    {
        if( E.has_parent() && E.parent().has_geological_feature() ) {
            E.set_geological_feature( E.parent().geological_feature() ) ;
        } else if( E.nb_children() > 0 && E.child( 0 ).has_geological_feature() ) {
            E.set_geological_feature( E.child( 0 ).geological_feature() ) ;
        }

        // Paranoia : we should perhaps verify that all children
        // have the same geological features
    }

    /*! \note Copied and modified from geogram\mesh\mesh_repair.cpp
    *
    * \brief Tests whether a facet is degenerate.
    * \param[in] M the mesh that the facet belongs to
    * \param[in] f the index of the facet in \p M
    * \return true if facet \p f has duplicated vertices,
    *  false otherwise
    */
    bool facet_is_degenerate( 
        const GEO::Mesh& M,
        index_t f,
        GEO::vector<index_t>& colocated_vertices )
    {
        index_t nb_vertices = M.facets.nb_vertices( f );
        if( nb_vertices != 3 ) {
            index_t* vertices = (index_t*)alloca( nb_vertices*sizeof( index_t ) );
            for( index_t lv = 0; lv<nb_vertices; ++lv ) {
                vertices[ lv ] = colocated_vertices[ M.facets.vertex( f, lv ) ];
            }
            std::sort( vertices, vertices + nb_vertices );
            return std::unique(
                vertices, vertices + nb_vertices
                ) != vertices + nb_vertices;
        }
        index_t c1 = M.facets.corners_begin( f );
        index_t c2 = c1 + 1;
        index_t c3 = c2 + 1;
        index_t v1 = colocated_vertices[ M.facet_corners.vertex( c1 ) ];
        index_t v2 = colocated_vertices[ M.facet_corners.vertex( c2 ) ];
        index_t v3 = colocated_vertices[ M.facet_corners.vertex( c3 ) ];
        return v1 == v2 || v2 == v3 || v3 == v1;
    }

    /*! \note Copied and modified from geogram\mesh\mesh_repair.cpp
     */
    void mesh_detect_degenerate_facets(
        const GEO::Mesh& M, 
        GEO::vector<index_t>& f_is_degenerate,
        GEO::vector<index_t>& colocated_vertices
        )
    {
        f_is_degenerate.resize( M.facets.nb() );
        for( index_t f = 0; f<M.facets.nb(); ++f ) {
            f_is_degenerate[ f ] = facet_is_degenerate( M, f, colocated_vertices );
        }
    }

    /*! 
     * @brief Detect and remove degenerated facets in a Mesh
     */
    index_t detect_degenerate_facets( GEO::Mesh& M )
    {
        GEO::vector< index_t > colocated ;
        GEO::mesh_detect_colocated_vertices( M, colocated ) ;

        GEO::vector< index_t > degenerate ;
        mesh_detect_degenerate_facets( M, degenerate, colocated ) ;
        return std::count( degenerate.begin(), degenerate.end(), 1 ) ;
    }

    bool edge_is_degenerate(
        const GEO::Mesh& M,
        index_t e,
        GEO::vector<index_t>& colocated_vertices )
    {
        index_t v1 = colocated_vertices[ M.edges.vertex( e, 0 ) ];
        index_t v2 = colocated_vertices[ M.edges.vertex( e, 1 ) ];
        return v1 == v2 ;
    }

    void mesh_detect_degenerate_edges(
        const GEO::Mesh& M,
        GEO::vector<index_t>& e_is_degenerate,
        GEO::vector<index_t>& colocated_vertices
        )
    {
        e_is_degenerate.resize( M.edges.nb() );
        for( index_t e = 0; e<M.edges.nb(); ++e ) {
            e_is_degenerate[ e ] = edge_is_degenerate( M, e, colocated_vertices );
        }
    }

    /*!
    * @brief Detect and remove degenerated edges in a Mesh
    */
    index_t repair_line_mesh( GEO::Mesh& M )
    {
        GEO::vector< index_t > colocated ;
        GEO::mesh_detect_colocated_vertices( M, colocated ) ;

        GEO::vector< index_t > degenerate ;
        mesh_detect_degenerate_edges( M, degenerate, colocated ) ;
        index_t nb = std::count( degenerate.begin(), degenerate.end(), 1 ) ;
        /// We have a problem if some vertices are left isolated
        /// If we ermove them here we can kill all indices correspondances
        M.edges.delete_elements( degenerate, false ) ;
        return nb ;
    }
  
    /*************************************************************************/
    /*!
    * @brief Get the index of an Interface from its name
    *
    * @param[in] BM the model to consider
    * @param[in] name Name of the Interface
    * @return Index of the interface in the model, NO_ID if not found.
    */
    bme_t find_interface( const BoundaryModel& BM, const std::string& name )
    {
        for( index_t i = 0; i < BM.nb_interfaces(); ++i ) {
            if( BM.one_interface( i ).name() == name ) {
                return BM.one_interface( i ).bme_id() ;
            }
        }
        return bme_t() ;
    }


    /*!
    * @brief Structure used to build Line by BoundaryModelBuilderGocad
    */
    struct Border {
        Border( index_t part, index_t corner, index_t p0, index_t p1 )
            : part_id_( part ), corner_id_( corner ), p0_( p0 ), p1_( p1 )
        {}

        // Id of the Surface owning this Border
        index_t part_id_ ;

        // Id of p0 in the BoundaryModel corner vector
        index_t corner_id_ ;

        // Ids of the starting corner and second vertex on the border in the Surface
        // to which this Border belong
        index_t p0_ ;
        index_t p1_ ;
    } ;

    /*************************************************************************/

    /*!
    * @brief Utility class to sort a set of oriented triangles around a common edge
    * Used in BoundaryModelBuilderSurface.
    *
    * This code could certainly be improved.
    */
    class SortTriangleAroundEdge {
    public:
        /*!
        * @brief A triangle to sort around an edge, see SortTriangleAroundEdge
        * @details This triangle belongs to a mesh connected component identified by its index.
        */
        struct TriangleToSort {
            /*!
            * @param index Index of this TriangleToSort in SortTriangleAroundEdge
            * @param surface_index Index of the Surface
            * @param p0 point of the triangle
            * @param p1 point of the triangle
            * @param p2 point of the triangle
            */
            TriangleToSort(
                index_t index,
                index_t surface_index,
                const vec3& p0,
                const vec3& p1,
                const vec3& p2 
            ) :
                index_( index ),
                surface_index_( surface_index ),
                N_(),
                B_A_(),
                angle_( -99999 ),
                side_( false )
            {
                ringmesh_assert( p0 != p1 ) ;
                ringmesh_assert( p0 != p2 ) ;
                ringmesh_assert( p1 != p2 ) ;

                vec3 e1 = normalize( p1 - p0 ) ;
                vec3 e2 = normalize( p2 - p0 ) ;

                N_ = normalize( cross( e1, e2 ) ) ;
                ringmesh_assert( dot( N_, e1 ) < epsilon ) ;

                vec3 B = 0.5 * p1 + 0.5 * p0 ;
                vec3 p2B = p2 - B ;
                B_A_ = normalize( p2B - dot( p2B, e1 ) * e1 ) ;

                ringmesh_assert( dot( B_A_, e1 ) < epsilon ) ;
                ringmesh_assert( B_A_.length() > epsilon ) ;
            };

            bool operator<( const TriangleToSort& r ) const
            {
                return angle_ < r.angle_ ;
            }

            /// Index in SortTriangleAroundEdge
            index_t index_ ;

            /// Global index of the surface owning this triangle
            index_t surface_index_ ;

            /// Normal to the triangle - normalized vector
            vec3 N_ ;

            /// Normal to the edge p0p1 in the plane defined by the triangle - normalized
            vec3 B_A_ ;

            // Values filled by sorting function in SortTriangleAroundEdge
            double angle_ ;
            bool side_ ;
        } ;

        SortTriangleAroundEdge()
        {}

        void add_triangle(
            index_t surface_index,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2 )
        {
            triangles_.push_back(
                TriangleToSort( triangles_.size(), surface_index, p0, p1, p2 ) ) ;
        }

        /*!
        * @details Rotation around axis of an angle A of the vector V
        *
        * @param axis Oriented axis
        * @param angle Angle in RADIANS
        * @param V the vector to rotate
        */
        static vec3 rotate( const vec3& axis, double angle, const vec3& V )
        {
            vec3 q = axis ;
            if( q.length() > 0 ) {
                double s = 1.0 / q.length() ;
                q[ 0 ] *= s ;
                q[ 1 ] *= s ;
                q[ 2 ] *= s ;
            }
            q *= sinf( 0.5 * angle ) ;

            double quat[ 4 ] = { q[ 0 ], q[ 1 ], q[ 2 ], cosf( 0.5 * angle ) } ;

            double m[ 4 ][ 4 ] ;

            m[ 0 ][ 0 ] = 1 - 2.0 * ( quat[ 1 ] * quat[ 1 ] + quat[ 2 ] * quat[ 2 ] ) ;
            m[ 0 ][ 1 ] = 2.0 * ( quat[ 0 ] * quat[ 1 ] + quat[ 2 ] * quat[ 3 ] ) ;
            m[ 0 ][ 2 ] = 2.0 * ( quat[ 2 ] * quat[ 0 ] - quat[ 1 ] * quat[ 3 ] ) ;
            m[ 0 ][ 3 ] = 0.0 ;

            m[ 1 ][ 0 ] = 2.0 * ( quat[ 0 ] * quat[ 1 ] - quat[ 2 ] * quat[ 3 ] ) ;
            m[ 1 ][ 1 ] = 1 - 2.0 * ( quat[ 2 ] * quat[ 2 ] + quat[ 0 ] * quat[ 0 ] ) ;
            m[ 1 ][ 2 ] = 2.0 * ( quat[ 1 ] * quat[ 2 ] + quat[ 0 ] * quat[ 3 ] ) ;
            m[ 1 ][ 3 ] = 0.0 ;

            m[ 2 ][ 0 ] = 2.0 * ( quat[ 2 ] * quat[ 0 ] + quat[ 1 ] * quat[ 3 ] ) ;
            m[ 2 ][ 1 ] = 2.0 * ( quat[ 1 ] * quat[ 2 ] - quat[ 0 ] * quat[ 3 ] ) ;
            m[ 2 ][ 2 ] = 1 - 2.0 * ( quat[ 1 ] * quat[ 1 ] + quat[ 0 ] * quat[ 0 ] ) ;
            m[ 2 ][ 3 ] = 0.0 ;

            m[ 3 ][ 0 ] = 0.0 ;
            m[ 3 ][ 1 ] = 0.0 ;
            m[ 3 ][ 2 ] = 0.0 ;
            m[ 3 ][ 3 ] = 1.0 ;

            double x = V[ 0 ] * m[ 0 ][ 0 ] + V[ 1 ] * m[ 1 ][ 0 ] + V[ 2 ] * m[ 2 ][ 0 ] + m[ 3 ][ 0 ] ;
            double y = V[ 0 ] * m[ 0 ][ 1 ] + V[ 1 ] * m[ 1 ][ 1 ] + V[ 2 ] * m[ 2 ][ 1 ] + m[ 3 ][ 1 ] ;
            double z = V[ 0 ] * m[ 0 ][ 2 ] + V[ 1 ] * m[ 1 ][ 2 ] + V[ 2 ] * m[ 2 ][ 2 ] + m[ 3 ][ 2 ] ;
            double w = V[ 0 ] * m[ 0 ][ 3 ] + V[ 1 ] * m[ 1 ][ 3 ] + V[ 2 ] * m[ 2 ][ 3 ] + m[ 3 ][ 3 ] ;
            return vec3( x / w, y / w, z / w ) ;
        }

        void sort()
        {
            ringmesh_assert( triangles_.size() > 0 ) ;

            std::pair< index_t, bool > default_pair( index_t( -1 ), false ) ;
            sorted_triangles_.resize( 2 * triangles_.size(), default_pair ) ;

            // If there is only one Triangle to sort - nothing to do
            if( triangles_.size() == 1 ) {
                sorted_triangles_[ 0 ] = std::pair< index_t, bool >(
                    triangles_[ 0 ].surface_index_, true ) ;
                sorted_triangles_[ 1 ] = std::pair< index_t, bool >(
                    triangles_[ 0 ].surface_index_, false ) ;
                return ;
            }

            // Initialization
            // We start on the plus (true) side of the first Triangle            
            sorted_triangles_[ 0 ] = std::pair< index_t, bool >(
                triangles_[ 0 ].surface_index_, true ) ;

            // Reference vectors with wich angles will be computed
            vec3 N_ref = triangles_[ 0 ].N_ ;
            vec3 B_A_ref = triangles_[ 0 ].B_A_ ;
            vec3 Ax_ref = normalize( cross( B_A_ref, N_ref ) ) ;

            // The minus (false) side of the start triangle will the last one encountered
            triangles_[ 0 ].angle_ = 2 * M_PI ;
            triangles_[ 0 ].side_ = false ;

            for( index_t i = 1; i < triangles_.size(); ++i ) {
                TriangleToSort& cur = triangles_[ i ] ;
                // Compute the angle RADIANS between the reference and the current
                // triangle 
                double cos = dot( B_A_ref, cur.B_A_ ) ;
                // Remove invalid values
                if( cos < -1 )
                    cos = -1 ;
                else if( cos > 1 ) cos = 1 ;
                cur.angle_ = std::acos( cos ) ;
                // Put the angle between PI and 2PI if necessary
                if( dot( cross( B_A_ref, cur.B_A_ ), Ax_ref ) < 0. ) {
                    cur.angle_ = 2 * M_PI - cur.angle_ ;
                }

                // Get the side of the surface first encountered
                // when rotating in the N_ref direction
                vec3 N_rotate = rotate( Ax_ref, -cur.angle_, cur.N_ ) ;
                cur.side_ = dot( N_rotate, N_ref ) > 0 ? false : true ;
            }

            // Sort the Surfaces according to the angle
            std::sort( triangles_.begin(), triangles_.end() ) ;

            // Fill the sorted surfaces adding the side
            index_t it = 1 ;
            for( index_t i = 0; i < triangles_.size(); ++i ) {
                TriangleToSort& cur = triangles_[ i ] ;
                if( triangles_[ i ].index_ == 0 ) { // The last to add
                    ringmesh_assert( i == triangles_.size() - 1 ) ;
                    sorted_triangles_[ it ].first = cur.surface_index_ ;
                    sorted_triangles_[ it ].second = cur.side_ ;
                } else {
                    sorted_triangles_[ it ].first = cur.surface_index_ ;
                    sorted_triangles_[ it ].second = cur.side_ ;
                    sorted_triangles_[ it + 1 ].first = cur.surface_index_ ;
                    sorted_triangles_[ it + 1 ].second = !cur.side_ ;
                    it += 2 ;
                }
            }
            // All the surfaces must have been sorted
            ringmesh_assert(
                std::count( sorted_triangles_.begin(), sorted_triangles_.end(),
                default_pair ) == 0 ) ;
        }
        /*! Returns the next pair Triangle index (surface) + side
        */
        const std::pair< index_t, bool >& next(
            const std::pair< index_t, bool >& in ) const
        {
            for( index_t i = 0; i < sorted_triangles_.size(); ++i ) {
                if( sorted_triangles_[ i ] == in ) {
                    if( i == sorted_triangles_.size() - 1 )
                        return sorted_triangles_[ sorted_triangles_.size() - 2 ] ;
                    if( i == 0 ) return sorted_triangles_[ 1 ] ;

                    if( sorted_triangles_[ i + 1 ].first == sorted_triangles_[ i ].first ) {
                        // The next has the same surface id, check its sign
                        if( sorted_triangles_[ i + 1 ].second
                            != sorted_triangles_[ i ].second ) {
                            return sorted_triangles_[ i - 1 ] ;
                        } else {
                            // Sign is the same
                            return sorted_triangles_[ i + 1 ] ;
                        }
                    } else {
                        ringmesh_assert(
                            sorted_triangles_[ i - 1 ].first
                            == sorted_triangles_[ i ].first ) ;
                        if( sorted_triangles_[ i - 1 ].second
                            != sorted_triangles_[ i ].second ) {
                            return sorted_triangles_[ i + 1 ] ;
                        } else {
                            return sorted_triangles_[ i - 1 ] ;
                        }
                    }
                }
            }
            ringmesh_assert_not_reached;
        }

    private:
        std::vector< TriangleToSort > triangles_ ;
        // Pairs global triangle identifier (Surface index) and side reached
        std::vector< std::pair< index_t, bool > > sorted_triangles_ ;
    } ;

    /*!
    * @brief Mini-factory. Create an empty element of the right type
    */
    BME* new_element( BME::TYPE T, BoundaryModel* M, index_t id )
    {

        if( T == BME::CORNER ) {
            return new Corner( M, id ) ;
        } else if( T == BME::LINE ) {
            return new Line( M, id ) ;
        } else if( T == BME::SURFACE ) {
            return new Surface( M, id ) ;
        } else if( T > BME::SURFACE && T < BME::NO_TYPE ) {
            return new BoundaryModelElement( M, T, id ) ;
        } else {
            return nil ;
        }
    }


    void print_model( const BoundaryModel& model )
    {
        GEO::Logger::out( "BoundaryModel" ) << "Model " << model.name() << " has "
            << std::endl << std::setw( 10 ) << std::left << model.nb_vertices()
            << " vertices " << std::endl << std::setw( 10 ) << std::left
            << model.nb_facets() << " facets " << std::endl << std::endl
            << std::setw( 10 ) << std::left << model.nb_regions() << " regions "
            << std::endl << std::setw( 10 ) << std::left << model.nb_surfaces()
            << " surfaces " << std::endl << std::setw( 10 ) << std::left
            << model.nb_lines() << " lines " << std::endl << std::setw( 10 )
            << std::left << model.nb_corners() << " corners " << std::endl
            << std::endl << std::setw( 10 ) << std::left << model.nb_contacts()
            << " contacts " << std::endl << std::setw( 10 ) << std::left
            << model.nb_interfaces() << " interfaces " << std::endl
            << std::setw( 10 ) << std::left << model.nb_layers() << " layers "
            << std::endl << std::endl ;
    }



    /*!
    * @brief Get the index of the Corner for a given point
    * @param[in] point Geometric location to look for
    * @return NO_ID or the index of the Corner
    */
    bme_t find_corner( const BoundaryModel& BM, const vec3& point )
    {
        for( index_t i = 0; i < BM.nb_corners(); ++i ) {
            if( BM.corner( i ).vertex() == point ) {
                return bme_t( BME::CORNER, i ) ;
            }
        }
        return bme_t() ;
    }

    /*!
    * @brief Get the index of the Corner at a given model point
    * @param[in] model_point_id Index of the point in the BoudaryModel
    * @return NO_ID or the index of the Corner
    */
    bme_t find_corner( const BoundaryModel& BM, index_t model_point_id )
    {
        for( index_t i = 0; i < BM.nb_corners(); ++i ) {
            if( BM.corner( i ).model_vertex_id() == model_point_id ) {
                return bme_t( BME::CORNER, i ) ;
            }
        }
        return bme_t() ;
    }

    /*!
    * @brief Find or create a corner at given coordinates.
    *
    * @param[in] point Geometric location of the Corner
    * @return Index of the Corner
    */
    bme_t find_or_create_corner( BoundaryModelBuilder& BMB, const vec3& point )
    {
        bme_t result = find_corner( BMB.model(), point ) ;
        if( !result.is_defined() ) {
            result = BMB.create_element( BME::CORNER ) ;
            BMB.set_corner( result, point ) ;
        }
        return result ;
    }

    /*!
    * @brief Find or create a line
    *
    * @param[in] BM model to consider
    * @param[in] vertices Coordinates of the vertices of the line
    * @return Index of the Line
    */
    bme_t find_or_create_line( BoundaryModelBuilder& BMB,
                               const std::vector< vec3 >& vertices )
    {
        bme_t result ;
        for( index_t i = 0; i < BMB.model().nb_lines(); ++i ) {
            if( BMB.model().line( i ).equal( vertices ) ) {
                result = BMB.model().line( i ).bme_id() ;
            }
        }
        if( !result.is_defined() ) {
            result = BMB.create_element( BME::LINE ) ;
            BMB.set_line( result, vertices ) ;

            // Find the indices of the corner at both extremities
            // Both must be defined to have a valid LINE
            BMB.add_element_boundary(
                result, find_or_create_corner( BMB, vertices.front() ) ) ;
            BMB.add_element_boundary(
                result, find_or_create_corner( BMB, vertices.back() ) ) ;
        }
        return result ;
    }



    /*!
    * @brief Fill the boundaries of all elements of the given type
    *
    * @details If the boundary elements do not have any in_boundary
    * information, nothing is done, and model construction will eventually fail.
    */
    void fill_elements_boundaries(
        BoundaryModelBuilder& B,
        BME::TYPE type )
    {
        // We have a problem if this is called for regions
        // No way yet to know the surface orientation
        ringmesh_debug_assert( type != BME::REGION ) ;

        BME::TYPE b_type = BME::boundary_type( type ) ;
        if( b_type != BME::NO_TYPE ) {
            for( index_t i = 0; i < B.model().nb_elements( b_type ); ++i ) {
                const BME& b = B.model().element( bme_t( b_type, i ) ) ;
                for( index_t j = 0; j < b.nb_in_boundary(); ++j ) {
                    B.add_element_boundary( b.in_boundary_id( j ),
                                            bme_t( b_type, i ) ) ;
                }
            }
        }
    }

    /*!
    * @brief Fill the in_boundary vector of all elements of the given type
    *
    * @details If the in_boundary elements do not have any boundary
    * information, nothing is done, and model construction will eventually fail.
    */
    void fill_elements_in_boundaries(
        BoundaryModelBuilder& B,
        BME::TYPE type )
    {
        BME::TYPE in_b_type = BME::in_boundary_type( type ) ;
        if( in_b_type != BME::NO_TYPE ) {
            for( index_t i = 0; i < B.model().nb_elements( in_b_type ); ++i ) {
                const BME& in_b = B.element( bme_t( in_b_type, i ) ) ;
                for( index_t j = 0; j < in_b.nb_boundaries(); ++j ) {
                    B.add_element_in_boundary( in_b.boundary_id( j ),
                                               bme_t( in_b_type, i ) ) ;
                }
            }
        }
    }

    /*!
    * @brief Fill the parent of all elements of the given type
    *
    * @details If the parents do not have any child
    *  nothing is done, and model construction will eventually fail.
    */
    void fill_elements_parent(
        BoundaryModelBuilder& B,
        BME::TYPE type )
    {
        BME::TYPE p_type = BME::parent_type( type ) ;
        if( p_type != BME::NO_TYPE ) {
            for( index_t i = 0; i < B.model().nb_elements( p_type ); ++i ) {
                const BME& p = B.model().element( bme_t( p_type, i ) ) ;
                for( index_t j = 0; j < p.nb_children(); ++j ) {
                    B.set_parent( p.child_id( j ), bme_t( p_type, i ) ) ;
                }
            }
        }
    }

    /*!
    * @brief Fill the children of all elements of the given type
    *
    * @details If the children elements do not have any parent information
    * nothing is done, and model construction will eventually fail.
    */
    void fill_elements_children(
        BoundaryModelBuilder& B,
        BME::TYPE type )
    {
        BME::TYPE c_type = BME::child_type( type ) ;
        if( c_type != BME::NO_TYPE ) {
            for( index_t i = 0; i < B.model().nb_elements( c_type ); ++i ) {
                bme_t cur_child = bme_t( c_type, i ) ;
                const bme_t& parent = B.model().element( cur_child ).parent_id() ;
                if( parent.is_defined() ) {
                    B.add_child( parent, cur_child ) ;
                }
            }
        }
    }

    /*!
    * @brief Remove degenerate facets and edges from the Surface
    *        and Line of the model.
    * @pre colocated vertices have already been removed
    */
    void remove_degenerate_facet_and_edges(
        const BoundaryModel& BM, std::set< bme_t >& to_remove
        )
    {
        to_remove.clear() ;
        for( index_t i = 0; i < BM.nb_lines(); ++i ) {
            index_t nb = repair_line_mesh( BM.line( i ).mesh() ) ;
            if( nb > 0 ) {
#ifdef RINGMESH_DEBUG
                GEO::Logger::out( "BoundaryModel" )
                    << nb << " degenerated edges removed in LINE "
                    << i << std::endl ;
#endif

                // The line may be empty now - remove it from the model
                if( BM.line( i ).nb_cells() == 0 ) {
                    to_remove.insert( BM.line( i ).bme_id() ) ;
                }
            }
        }

        for( index_t i = 0; i < BM.nb_surfaces(); ++i ) {
            GEO::Mesh& M = BM.surface( i ).mesh() ;
            index_t nb = detect_degenerate_facets( M ) ;
            /// @todo Check if that cannot be simplified 
            if( nb > 0 ) {
                // If there are some degenerated facets 
                // We need to repair the model 
#ifndef RINGMESH_DEBUG
                GEO::Logger::instance()->set_quiet( true ) ;
#endif
                // Using repair function of geogram
                // Warning - This triangulates the mesh

                if( M.vertices.nb() > 0 ) {
                    // Colocated vertices must be processed before
                    // MESH_REPAIR_DUP_F 2 ;
                    GEO::MeshRepairMode mode =
                        static_cast< GEO::MeshRepairMode >( 2 ) ;
                    GEO::mesh_repair( M, mode ) ;

                    // This might create some small components - remove them
                    /// @todo How to choose the epsilon ? and the maximum number of facets ?
                    GEO::remove_small_connected_components( M, epsilon_sq, 3 ) ;

                    // This is a bit of an overkill but I am lazy
                    if( M.vertices.nb() > 0 ) {
                        GEO::mesh_repair( M, mode ) ;
                    }
                }
#ifndef RINGMESH_DEBUG
                GEO::Logger::instance()->set_quiet( false ) ;
#endif 
                if( M.vertices.nb() == 0 || M.facets.nb() == 0 ) {
                    to_remove.insert( BM.surface( i ).bme_id() ) ;
                }
                else {
                    // If the Surface has internal boundaries, we need to 
                    // re-cut the Surface along these lines
                    Surface& S = const_cast< Surface& >( BM.surface( i ) ) ;
                    std::set< index_t > cutting_lines ;
                    for( index_t l = 0; l < S.nb_boundaries(); ++l ) {
                        const Line& L = BM.line( S.boundary_id( l ).index ) ;
                        if( L.is_inside_border( S ) ) {
                            cutting_lines.insert( L.bme_id().index ) ;
                        }
                    }
                    for( std::set< index_t>::iterator it = cutting_lines.begin();
                         it != cutting_lines.end(); ++it 
                    ) {
                        S.cut_by_line( BM.line( *it ) ) ;
                    }
                }
            }
        }
    }


    /*!
    * Get the indices of the duplicated vertices that are on an inside border.
    * Only the vertex with the biggest index are added.
    * If there are more than 2 colocated vertices throws an assertion in debug mode
    */
    void vertices_on_inside_boundary(
        const BoundaryModelMeshElement& E,
        std::set< index_t >& vertices )
    {
        vertices.clear() ;
        if( E.bme_id().type == BME::CORNER ) {
            return ;
        }
        if( E.bme_id().type == BME::LINE ) {
            if( E.boundary( 0 ).is_inside_border( E ) ) {
                vertices.insert( E.nb_vertices()-1 ) ;
            }
            return ;
        }
        std::vector< const BoundaryModelMeshElement* > inside_border ;
        for( index_t i = 0; i < E.nb_boundaries(); ++i ) {
            if( E.boundary( i ).is_inside_border( E ) ) {
                inside_border.push_back(
                    dynamic_cast<const BoundaryModelMeshElement*>( &E.boundary( i ) ) ) ;
            }
        }
        if( !inside_border.empty() ) {
            // We want to get the indices of the vertices in E
            // that are colocated with those of the inside boundary
            // We assume that the model vertice are not yet computed
            GEO::NearestNeighborSearch_var kdtree =
                GEO::NearestNeighborSearch::create( 3, "BNN" ) ;
            kdtree->set_points( E.mesh().vertices.nb(),
                                E.mesh().vertices.point_ptr( 0 ) ) ;

            for( index_t i = 0; i < inside_border.size(); ++i ) {
                const GEO::Mesh& m = inside_border[ i ]->mesh() ;
                for( index_t v = 0; v < m.vertices.nb(); ++v ) {
                    index_t neighbors[ 3 ] ;
                    double sq_dist[ 3 ] ;
                    kdtree->get_nearest_neighbors(
                        3, m.vertices.point_ptr( v ), neighbors, sq_dist ) ;

                    if( sq_dist[ 2 ] < epsilon_sq ) {
                        // We have a problem there colocated points on the inside
                        // border. They are not properly processed. 
                        ringmesh_debug_assert( false ) ;
                    } else if( sq_dist[ 1 ] <epsilon_sq ) {
                        // Colocated vertices
                        vertices.insert( GEO::geo_max( neighbors[ 0 ], neighbors[ 1 ] ) ) ;
                    }
                    // Otherwise nothing to do
                }
            }
        }
    }

    void remove_colocated_element_vertices(
        BoundaryModel& BM, std::set< bme_t >& to_remove )
    {
        to_remove.clear() ;
        // For all Lines and Surfaces
        for( index_t t = BME::LINE; t < BME::REGION; ++t ) {
            BME::TYPE T = static_cast<BME::TYPE>( t ) ;

            for( index_t e = 0; e < BM.nb_elements( T ); ++e ) {
                const BMME& E = dynamic_cast< const BMME&>(
                    BM.element( bme_t( T, e ) ) ) ;

                GEO::Mesh& M = E.mesh() ;
                GEO::vector< index_t > colocated ;
                GEO::mesh_detect_colocated_vertices( M, colocated, epsilon ) ;

                // Get the vertices to delete
                std::set< index_t > inside_border ;
                vertices_on_inside_boundary( E, inside_border ) ;

                GEO::vector< index_t > to_delete( colocated.size(), 0 ) ;
                GEO::vector< index_t > old2new( colocated.size() ) ;
                index_t nb_todelete = 0 ;
                index_t cur = 0 ;
                for( index_t v = 0; v < colocated.size(); ++v ) {
                    if( colocated[ v ] == v ||
                        inside_border.find( v ) != inside_border.end()
                    ) {
                        // This point is kept 
                        // No colocated or on an inside boundary
                        old2new[ v ] = cur ;
                        ++cur ;
                    } else {
                        // The point is to remove
                        old2new[ v ] = old2new[ colocated[ v ] ] ;
                        to_delete[ v ] = 1 ;
                        nb_todelete++ ;                        
                    }
                }

                if( nb_todelete == 0 ) {
                    // Nothing to do there
                    continue ;
                } else if( nb_todelete == E.nb_vertices() ) {
                    // The complete element should be removed
                    to_remove.insert( E.bme_id() ) ;
                    continue ;
                } else {
                    // We need to update the VertexInBME at the model level
                    // if they exist. If not let's avoid this costly operation
                    // @todo This is bugged ?? I am not sure. JP
                    if( BM.vertices.is_initialized() ) {
                        // For all the vertices of this element
                        // we need to update the vertex ids in the bme_vertices_
                        // of the corresponding global vertex
                        for( index_t v = 0; v < to_delete.size(); ++v ) {
                            index_t model_id = E.model_vertex_id( v ) ;
                            const std::vector<VBME>& cur =
                                BM.vertices.bme_vertices( model_id ) ;
                            for( index_t i = 0; i < cur.size(); ++i ) {
                                if( cur[ i ] == VBME( E.bme_id(), v) ) {
                                    index_t new_id = old2new[ v ] ;
                                    BM.vertices.set_bme(
                                        model_id, i, VBME( E.bme_id(), new_id ) ) ;
                                }
                            }
                        }
                        BM.vertices.erase_invalid_vertices() ;
                    }

                    for( index_t c = 0; c < M.facet_corners.nb(); c++ ) {
                        M.facet_corners.set_vertex( c, colocated[ M.facet_corners.vertex( c ) ] );
                    }
                    for( index_t e = 0; e < M.edges.nb(); e++ ) {
                        M.edges.set_vertex( e, 0, colocated[ M.edges.vertex( e, 0 ) ] ) ;
                        M.edges.set_vertex( e, 1, colocated[ M.edges.vertex( e, 1 ) ] ) ;
                    }
                    M.vertices.delete_elements( to_delete, false ) ;
#ifdef RINGMESH_DEBUG
                    GEO::Logger::out( "BoundaryModel" )
                        << nb_todelete << " colocated vertices deleted in "
                        << E.bme_id() << std::endl ;
#endif
                }
            }
        }
    }




}

namespace RINGMesh {

    /*!
     * @brief Creates a element of the given type and add it to the correct vector
     * The BoundaryModelElement is created from its type and its index
     *
     * @param[in] type Type of the element to create
     * @return The index of the created element
     */
    bme_t BoundaryModelBuilder::create_element( BME::TYPE type )
    {
        index_t id = model_.nb_elements( type ) ;
        ringmesh_assert( id != NO_ID ) ;
        if( type >= BME::CORNER && type < BME::NO_TYPE ) {
            model_.modifiable_elements( type ).push_back(
                new_element( type, &model_, id ) ) ;
            return bme_t( type, id ) ;
        }
        else {
            ringmesh_assert_not_reached;
            return bme_t() ;
        }
    }

    /*!
     * @brief Add to the vector the elements which cannot exist if 
     *        an element in the set does not exist.
     * @details These elements are added to the set.
     *          Recursive call till nothing is added.
     *          
     * @return True if at least one element was added, otherwise false.
     */
    bool BoundaryModelBuilder::get_dependent_elements(
        std::set< bme_t >& in ) const 
    {
        index_t input_size = in.size() ;

        for( std::set< bme_t >::iterator it( in.begin() );
             it != in.end(); ++it 
           ) {
            bme_t cur = *it ;
            /// If an element has children elements - add them 
            if( BME::child_allowed( cur.type ) ) {
                index_t c = BME::child_type( cur.type ) ;
                const BME& E = model_.element( cur ) ;
                for( index_t j = 0 ; j < E.nb_children() ; ++j ) {
                    in.insert( E.child_id(j) ) ;
                }                
            }         
        }

        /// If a parent has no children anymore - add it 
        for( index_t p = BME::CONTACT; p < BME::NO_TYPE; ++p ) {
            BME::TYPE P = ( BME::TYPE ) p ;            
            for( index_t j = 0; j < model_.nb_elements( P ); ++j ) {
                bool no_child = true ;
                const BME& E = model_.element( BME::bme_t( P, j ) ) ;
                for( index_t k = 0; k < E.nb_children(); ++k ) {                    
                    if( in.count( E.child_id( k ) ) == 0 ) {
                        no_child = false ;
                        break ;
                    }
                }
                if( no_child ) {
                    in.insert( E.bme_id() ) ;
                }
            }
        }

        /// If an element is in the boundary of nothing - add it
        for( index_t t = BME::CORNER; t < BME::REGION; ++t ) {
            BME::TYPE T = ( BME::TYPE ) t ;
            for( index_t j = 0; j < model_.nb_elements( T ); ++j ) {
                bool no_incident = true ;
                const BME& E = model_.element( BME::bme_t( T, j ) ) ;
                for( index_t k = 0; k < E.nb_in_boundary(); ++k ) {
                    if( in.count( E.in_boundary_id( k ) ) == 0 ) {
                        no_incident = false ;
                        break ;
                    }
                }
                if( no_incident ) {
                    in.insert( E.bme_id() ) ;
                }
            }
        }

        if( in.size() != input_size ) {
            return get_dependent_elements( in ) ;
        }
        else {
            return false ;
        }
    }
    

    /*!
     * @brief Remove a list of elements of the model
     * @details No check is done on the consistency of this removal
     *          The elements and all references to them are removed. 
     *          All dependent elements should be in the set of elements to remove,
     *          with a prior call to get_dependent_elements function.
     * 
     * @warning NOT TESTED.
     *          The client is responsible to set the proper connectivity
     *          information between the remaining model elements.
     * 
     * @todo TEST IT
     */
    void BoundaryModelBuilder::remove_elements( 
        const std::set< bme_t >& elements )
    {
        if( elements.size() == 0 ) {
            return ;
        }

        // We need to remove elements type by type since they are 
        // stored in different vectors and since we use indices in these 
        // vectors to identify them.
        // Initialize the vector
        std::vector < std::vector < index_t > > to_erase_by_type ;
        for( index_t i = BME::CORNER; i < BME::NO_TYPE; ++i ) {
            to_erase_by_type.push_back( std::vector< index_t >( 
                model_.nb_elements( static_cast<BME::TYPE>(i) ), 0 ) ) ;
        }
        // Flag the elements to erase
        for( std::set< bme_t >::const_iterator it = elements.begin();
             it != elements.end(); ++it ) {
            bme_t cur = *it ;
            if( cur.type < BME::NO_TYPE ) {
                to_erase_by_type[ cur.type ][ cur.index ] = NO_ID ;                
            }
        }
       
        delete_elements( to_erase_by_type ) ;
        
        // Re-initialize global access to elements
        init_global_model_element_access() ;       
    }

    /*!
    * @brief Delete elements and remove all references to them in the model
    *
    * @param[in,out] to_erase For each type of element T, 
    *        store a vector of the size of model_.nb_elements(T) in which
    *        elements are flagged with NO_ID.
    *        In output it stores the mapping table between old and new indices
    *        for the elements.
    * @todo TEST IT
    */
    void BoundaryModelBuilder::delete_elements(
        std::vector< std::vector< index_t > >& to_erase )
    {        
        // Number of elements deleted for each TYPE
        std::vector< index_t > nb_removed( to_erase.size(), 0 ) ;

        /// 1. Get the mapping between old indices of the elements
        ///    and new ones (when elements to remove will actually be removed)
        for( index_t i = 0; i < to_erase.size(); ++i ) {
            for( index_t j = 0 ; j < to_erase[ i ].size(); ++j ) {
                if( to_erase[ i ][ j ] == NO_ID ) {
                    nb_removed[ i ]++ ;
                } else {
                    to_erase[ i ][ j ] = j - nb_removed[ i ] ;
                }
            }
        }

        /// 2. Effectively delete the elements 
        for( index_t i = 0; i < to_erase.size(); ++i ) {
            for( index_t j = 0 ; j < to_erase[ i ].size(); ++j ) {
                if( to_erase[ i ][ j ] == NO_ID ) {
                    BME::bme_t cur( static_cast<BME::TYPE>( i ), j ) ;
                    delete element_ptr( cur ) ;
                    set_element( cur, nil ) ;
                }
            }
            std::vector< BME* >& store = model_.modifiable_elements(
                static_cast<BME::TYPE>( i ) ) ;
            store.erase( std::remove( store.begin(), store.end(),
                static_cast<BME*>( nil ) ), store.end() ) ;
        }

        /// 3. Deal with the model vertices
        // We need to to this before changing ids of the BMEs
        // (because we need the old ids to get the new ones)
        // but after deleting the elements otherwise 
        // we are putting back the vertices of the BMMEs to delete
        for( index_t v = 0; v < model_.vertices.nb(); ++v ) {
            const std::vector<BoundaryModelVertices::VertexInBME>& cur =
                model_.vertices.bme_vertices( v ) ;
            for( index_t i = 0; i < cur.size(); ++i ) {
                bme_t id = cur[ i ].bme_id ;
                bme_t new_id( id.type, to_erase[ id.type ][ id.index ] ) ;
                model_.vertices.set_bme(
                    v, i, BoundaryModelVertices::VertexInBME( new_id, cur[ i ].v_id ) ) ;
            }
        }
        model_.vertices.erase_invalid_vertices() ;

        /// 3. Update all possible indices in remaining elements
        for( index_t i = 0; i < to_erase.size(); ++i ) {
            BME::TYPE T = static_cast<BME::TYPE>( i ) ;

            // Update all indices stored by the BME of that type 
            for( index_t j = 0; j < model_.nb_elements( T ); ++j ) {                               
                BoundaryModelElement& E = element( bme_t( T, j ) ) ;

                // Not the same than j - since we have erased some elements
                index_t old_id = E.bme_id().index ;
                ringmesh_debug_assert( to_erase[ i ][ old_id ] != NO_ID ) ;
                    
                // id_ 
                E.set_id( to_erase[ i ][ old_id ] ) ;
                // boundary_
                if( E.nb_boundaries() > 0 ) {
                    BME::TYPE B = BME::boundary_type( T ) ;
                    ringmesh_debug_assert( B < BME::NO_TYPE ) ;
                    for( index_t k = 0; k < E.nb_boundaries(); ++k ) {
                        E.set_boundary( k, bme_t(
                            B, to_erase[ B ][ E.boundary_id( k ).index ] ) ) ;
                    }
                }
                // in_boundary
                if( E.nb_in_boundary() > 0 ) {
                    BME::TYPE IB = BME::in_boundary_type( T ) ;
                    ringmesh_debug_assert( IB < BME::NO_TYPE ) ;
                    for( index_t k = 0; k < E.nb_in_boundary(); ++k ) {
                        E.set_in_boundary( k, bme_t(
                            IB, to_erase[ IB ][ E.in_boundary_id( k ).index ] ) ) ;
                    }
                }
                // parent_
                if( E.has_parent() ) {
                    BME::TYPE P = BME::parent_type( T ) ;
                    ringmesh_debug_assert( P < BME::NO_TYPE ) ;
                    E.set_parent( bme_t( P, to_erase[ P ][ E.parent_id().index ] ) ) ;
                }
                // children_ 
                if( E.nb_children() > 0 ) {
                    BME::TYPE C = BME::child_type( T ) ;
                    ringmesh_debug_assert( C < BME::NO_TYPE ) ;
                    for( index_t k = 0; k < E.nb_children(); ++k ) {
                        E.set_child( k, bme_t(
                            C, to_erase[ C ][ E.child_id( k ).index ] ) ) ;
                    }
                }
                // Clean the vectors in the element
                E.erase_invalid_element_references() ;
            }
        }
        // Do not forget the universe
        /// @todo Put the universe in the list of regions - so annoying to think of it each time
        {
            for( index_t i = 0; i < model_.universe().nb_boundaries(); ++i ) {
                model_.universe_.set_boundary( i, bme_t(
                    BME::SURFACE, to_erase[ BME::SURFACE ][
                        model_.universe().boundary_id( i ).index ] ) ) ;
            }
            model_.universe_.erase_invalid_element_references() ;
        }             
    }
    

    void BoundaryModelBuilder::resize_elements( BME::TYPE type, index_t nb )
    {
        if( type >= BME::NO_TYPE ) {
            return ;
        }
        std::vector< BME*>& store = model_.modifiable_elements( type ) ;
        store.resize( nb, nil ) ;
        for( index_t i = 0; i < nb; i++ ) {
            store[ i ] = new_element( type, &model_, i ) ;
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
        model_.universe_.set_element_type( BME::REGION ) ;
        model_.universe_.set_model( &model_ ) ;

        for( index_t i = 0; i < boundaries.size(); ++i ) {
            ringmesh_assert( boundaries[i].first < model_.nb_surfaces() ) ;
            model_.universe_.add_boundary(
                bme_t( BME::SURFACE, boundaries[i].first ),
                boundaries[i].second ) ;
        }
    }

    /*!
     * @brief Set the geometric location of a Corner
     *
     * @param[in] corner_id Index of the corner
     * @param[in] point Coordinates of the vertex
     */
    void BoundaryModelBuilder::set_corner(
        const bme_t& corner_id,
        const vec3& point )
    {
        ringmesh_assert( corner_id.index < model_.nb_corners() ) ;
        dynamic_cast< Corner* >( model_.corners_[corner_id.index] )->
            set_vertex( point, false ) ;
    }

    /*!
     * @brief Set one Line points
     *
     * @param[in] id Line index
     * @param[in] vertices Coordinates of the vertices on the line
     */
    void BoundaryModelBuilder::set_line(
        const bme_t& id,
        const std::vector< vec3 >& vertices )
    {
        ringmesh_assert( id.index < model_.nb_lines() ) ;
        dynamic_cast< Line* >( model_.lines_[ id.index ] )->
            set_vertices( vertices ) ;
    }

    /*!
     * @brief Set the points and facets for a surface
     * @details If facet_adjacencies are not given they are computed.
     *
     * @param[in] surface_id Index of the surface
     * @param[in] points Coordinates of the vertices
     * @param[in] facets Indices in the vertices vector to build facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     */
    void BoundaryModelBuilder::set_surface_geometry(
        const bme_t& surface_id,
        const std::vector< vec3 >& points,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        if( facets.size() == 0 ) {
            return ;
        }

        dynamic_cast< Surface* > (model_.surfaces_[surface_id.index] )->
            set_geometry( points, facets, facet_ptr ) ;

        set_surface_adjacencies( surface_id ) ;
    }

    /*!
     * @brief Add a point to the BoundaryModel and not to one of its elements
     * @details To use when adding the points to the model before building its elements
     */
    index_t BoundaryModelBuilder::add_unique_vertex( const vec3& p )
    {
        return model_.vertices.add_unique_vertex( p ) ;
    }

    /*!
     * @brief Set the vertex for a Corner. Store the info in the BM vertices
     *
     * @param[in] corner_id Index of the corner
     * @param[in] unique_vertex Index of the vertex in the model
     */
    void BoundaryModelBuilder::set_corner(
        const bme_t& corner_id,
        index_t unique_vertex )
    {
        ringmesh_assert( corner_id.index < model_.nb_corners() ) ;
        dynamic_cast< Corner* >( model_.corners_[corner_id.index] )->
            set_vertex( unique_vertex ) ;
    }

    /*!
     * @brief Set one Line vertices. Store the info in the BM vertices
     *
     * @param[in] id Line index
     * @param[in] unique_vertices Indices in the model of the unique vertices with which to build the Line
     */
    void BoundaryModelBuilder::set_line(
        const bme_t& id,
        const std::vector< index_t >& unique_vertices )
    {
        ringmesh_assert( id.index < model_.nb_lines() ) ;
        dynamic_cast< Line* >( model_.lines_[ id.index ] )->
            set_vertices( unique_vertices ) ;
    }

    /*!
     * @brief Set the vertices and facets for a surface
     * @details If facet_adjacencies are not given they are computed.
     *
     * @param[in] surface_id Index of the surface
     * @param[in] model_vertex_ids Indices of unique vertices in the BoundaryModel
     * @param[in] facets Indices in the vertices vector to build facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets
     */
    void BoundaryModelBuilder::set_surface_geometry(
        const bme_t& surface_id,
        const std::vector< index_t >& model_vertex_ids,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        if( facets.size() == 0 ) {
            return ;
        }
        dynamic_cast< Surface* >( model_.surfaces_[ surface_id.index ] )->
            set_geometry( model_vertex_ids, facets, facet_ptr ) ;

        set_surface_adjacencies( surface_id ) ;
    }

    /*!
     * @brief Set the facets of a surface
     *
     * @param[in] surface_id Index of the surface
     * @param[in] facets Indices of the model vertices defining the facets
     * @param[in] facet_ptr Pointer to the beginning of a facet in facets     
     */
    void BoundaryModelBuilder::set_surface_geometry(
        const bme_t& surface_id,
        const std::vector< index_t >& facets,
        const std::vector< index_t >& facet_ptr )
    {
        if( facets.size() == 0 ) {
            return ;
        }

        // Compute the vertices from the corners
        // This is quite stupid !! The real solution would be to remove
        // the vertices vector from the Surface
        std::map< index_t, index_t > old_2_new ;

        std::vector< index_t > vertices ;
        std::vector< index_t > facets_local( facets.size() ) ;
        for( index_t i = 0; i < facets.size(); ++i ) {
            index_t c = facets[i] ;
            std::map< index_t, index_t >::iterator it = old_2_new.find( c ) ;
            index_t new_corner_id = NO_ID ;

            if( it == old_2_new.end() ) {
                new_corner_id = vertices.size() ;
                old_2_new[c] = new_corner_id ;

                // Not so great to push back, but whatever
                vertices.push_back( c ) ;
            } else {
                new_corner_id = old_2_new[c] ;
            }
            facets_local[i] = new_corner_id ;
        }

        set_surface_geometry( surface_id, vertices, facets_local, facet_ptr ) ;
    }

    /*!
     * @brief Compute and set the adjacencies between the facets
     * @details The adjacent facet is given for each vertex of each facet for the edge
     * starting at this vertex.
     * If there is no neighbor inside the same Surface adjacent is set to NO_ADJACENT
     *
     * @param[in] surface_id Index of the surface
     */
    void BoundaryModelBuilder::set_surface_adjacencies(
        const bme_t& surface_id )
    {
        Surface& S = dynamic_cast< Surface& >( *model_.surfaces_[surface_id.index] );
        ringmesh_assert( S.nb_cells() > 0 ) ;

        std::vector< index_t > adjacent ;
        adjacent.resize( S.facet_end( S.nb_cells() - 1 ), Surface::NO_ADJACENT ) ;

        index_t nb_facets = S.nb_cells() ;
        index_t nb_vertices = S.nb_vertices() ;

        // Allocate some space to store the ids of facets around each vertex
        std::vector< index_t > toto ;
        toto.reserve( 10 ) ;
        std::vector< std::vector< index_t > > vertex_to_facets( nb_vertices, toto ) ;

        for( index_t f = 0; f < nb_facets; ++f ) {
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                vertex_to_facets[S.surf_vertex_id( f, v )].push_back( f ) ;
            }
        }
        for( index_t p = 0; p < nb_vertices; ++p ) {
            std::sort( vertex_to_facets[p].begin(), vertex_to_facets[p].end() ) ;
        }

        for( index_t f = 0; f < nb_facets; ++f ) {
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                index_t cur = S.surf_vertex_id( f, v ) ;
                index_t prev = S.surf_vertex_id( f, S.prev_in_facet( f, v ) ) ;

                const std::vector< index_t >& f_prev = vertex_to_facets[prev] ;
                const std::vector< index_t >& f_cur = vertex_to_facets[cur] ;

                std::vector< index_t > inter(
                    std::min( f_prev.size(), f_cur.size() ) ) ;
                index_t end = narrow_cast< index_t >(
                    std::set_intersection( f_prev.begin(), f_prev.end(),
                        f_cur.begin(), f_cur.end(), inter.begin() )
                        - inter.begin() ) ;

                if( end == 2 ) {
                    // There is one neighbor
                    index_t f2 = inter[0] == f ? inter[1] : inter[0] ;
                    adjacent[S.facet_begin( f ) + S.prev_in_facet( f, v )] = f2 ;
                } else {
                    ringmesh_debug_assert( end == 1 ) ;
                }
            }
        }
        S.set_adjacent( adjacent ) ;
    }

    /*!
     * @brief Complete missing information in BoundaryModelElements
     * boundaries - in_boundary - parent - children
     *
     * @details For all 7 types of elements, check what information is available
     * for the first one and fill the elements of the same type accordingly
     * THIS MEANS that the all the elements of the same type have been initialized with
     * the same information.
     */
    bool BoundaryModelBuilder::complete_element_connectivity()
    {
        // Lines
        if( model_.nb_lines() > 0 ) {
            if( model_.line( 0 ).nb_boundaries() == 0 ) {
                fill_elements_boundaries( *this, BME::LINE ) ;
            }
            if( model_.line( 0 ).nb_in_boundary() == 0 ) {
                fill_elements_in_boundaries( *this, BME::LINE ) ;
            }
            if( !model_.line( 0 ).parent_id().is_defined()
                && model_.nb_contacts() > 0 ) {
                fill_elements_parent( *this, BME::LINE ) ;
            }
        }

        // Corners
        if( model_.nb_corners() > 0 && model_.corner( 0 ).nb_in_boundary() == 0 ) {
            // Info from line boundaries is used here and should be available
            fill_elements_in_boundaries( *this, BME::CORNER ) ;
        }

        // Surfaces - There MUST be at least one
        if( model_.surface( 0 ).nb_boundaries() == 0 ) {
            fill_elements_boundaries( *this, BME::SURFACE ) ;
        }
        if( model_.surface( 0 ).nb_in_boundary() == 0 ) {
            fill_elements_in_boundaries( *this, BME::SURFACE ) ;
        }
        if( !model_.surface( 0 ).parent_id().is_defined() ) {
            fill_elements_parent( *this, BME::SURFACE ) ;
        }

        // Regions
        if( model_.nb_regions() > 0 ) {
            if( model_.region( 0 ).nb_boundaries() == 0 ) {
                fill_elements_boundaries( *this, BME::REGION ) ;
            }
            if( !model_.region( 0 ).parent_id().is_defined()
                && model_.nb_layers() > 0 ) {
                fill_elements_parent( *this, BME::REGION ) ;
            }
        }

        // Contacts
        if( model_.nb_contacts() > 0 && model_.contact( 0 ).nb_children() == 0 ) {
            fill_elements_children( *this, BME::CONTACT ) ;
        }

        // Interfaces
        if( model_.nb_interfaces() > 0
            && model_.one_interface( 0 ).nb_children() == 0 ) {
            fill_elements_children( *this, BME::INTERFACE ) ;
        }

        // Layers
        if( model_.nb_layers() > 0 && model_.layer( 0 ).nb_children() == 0 ) {
            fill_elements_children( *this, BME::LAYER ) ;
        }
        return true ;
    }


    /*!
     * @brief Clears and fills the model nb_elements_per_type_ vector
     * @details See global element access with BoundaryModel::element( BME::TYPE, index_t )
     */
    void BoundaryModelBuilder::init_global_model_element_access()
    {
        model_.nb_elements_per_type_.clear() ;

        index_t count = 0 ;
        model_.nb_elements_per_type_.push_back( count ) ;
        for( index_t type = BME::CORNER; type < BME::NO_TYPE; type++ ) {
            count += model_.nb_elements( (BME::TYPE) type ) ;
            model_.nb_elements_per_type_.push_back( count ) ;
        }
    }


    /*!
     * @brief This function MUST be the last function called when building a BoundaryModel
     *
     * @details Check that the model is correct and has all required information
     * Calls the complete_element_connectivity function
     * Fills nb_elements_per_type_ vector
     *
     * @return False if the model is not valid and cannot be fixed
     * otherwise returns true.
     *
     */
    bool BoundaryModelBuilder::end_model()
    {
        // The name should exist
        if( model_.name() == "" ) {
            set_model_name( "model_default_name" ) ;
        }
        
        // Get out if the model has no surface
        if( model_.nb_surfaces() == 0 ) {
            print_model( model_ ) ;
            return false ;
        }

        init_global_model_element_access() ;

        complete_element_connectivity() ;

        // Fill geological feature if missing
        for( index_t i = 0; i < model_.nb_elements( BME::ALL_TYPES ); ++i ) {
            BME& E = element( bme_t( BME::ALL_TYPES, i ) ) ;
            if( !E.has_geological_feature() ) {
                fill_element_geological_feature( E ) ;
            }
        }

        // Remove colocated vertices in each element
        std::set< bme_t > empty_elements ;
        remove_colocated_element_vertices( model_, empty_elements ) ;
        if( !empty_elements.empty() ) {
            get_dependent_elements( empty_elements ) ;
            remove_elements( empty_elements ) ;
        }

        // Basic mesh repair for surfaces and lines
        remove_degenerate_facet_and_edges( model_, empty_elements ) ;
        if( !empty_elements.empty() ) {
            get_dependent_elements( empty_elements ) ;
            remove_elements( empty_elements ) ;
        }

        // This is basic requirement ! no_colocated model vertices !
        // So remove them if there are any 
        model_.vertices.remove_colocated() ; 

        if( model_.check_model_validity() ) {
            GEO::Logger::out( "BoundaryModel" ) 
                << std::endl  
                << "Model " << model_.name() << " is valid " 
                << std::endl << std::endl ;
            print_model( model_ ) ;
            return true ;
        } else {
            GEO::Logger::out( "BoundaryModel" ) 
                << std::endl
                << "Model " << model_.name() << " is invalid "                 
                << std::endl << std::endl ;
            print_model( model_ ) ;
            return false ;
        }
    }

    /*************************************************************************/

    /*!
     * @brief Load and build a BoundaryModel from a Gocad .ml file
     *
     *  @details This is pretty tricky because of the annoying not well adapted file format.
     * The correspondance between Gocad::Model3D elements and BoundaryModel elements is :
     *  - Gocad TSurf  <-> BoundaryModel Interface
     *  - Gocad TFace  <-> BoundaryModel Surface
     *  - Gocad Region <-> BoundaryModel Region
     *  - Gocad Layer  <-> BoundaryModel Layer
     *
     * @param[in] ml_file_name Input .ml file stream
     */
    bool BoundaryModelBuilderGocad::load_ml_file( const std::string& ml_file_name )
    {
        GEO::LineInput in( ml_file_name ) ;
        if( !in.OK() ) {
            return false ;
        }

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

        /// The file contains 2 parts and is read in 2 steps
        /// 1. Read global information on model elements
        /// 2. Read surface geometries and info to build corners and contacts
        bool read_model = true ;

        // The orientation of positive Z
        // can change for each TSurf and need to be read
        int z_sign = 1 ;

        // In the .ml file - vertices are indexed TSurf by Tsurf
        // They can be duplicated inside one TSurf and betweeen TSurfs

        // Indices of the vertices of the currently built TSurf in the model
        std::vector< vec3 > tsurf_vertices ;

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

        while( !in.eof() && in.get_line() ) {
            in.get_fields() ;
            if( in.nb_fields() > 0 ) {
                if( read_model ) {
                    if( strncmp( in.field( 0 ), "name:", 5 ) == 0 ) {
                        set_model_name( &in.field( 0 )[5] ) ;
                    } else if( in.field_matches( 0, "TSURF" ) ) {
                        /// 1.1 Create Interface from its name
                        index_t f = 1 ;
                        std::ostringstream oss ;
                        do {
                            oss << in.field( f++ ) ;
                        } while( f < in.nb_fields() ) ;
                        // Create an interface and set its name
                        set_element_name( create_element( BME::INTERFACE ), oss.str() ) ;

                        nb_tsurf++ ;
                    } else if( in.field_matches( 0, "TFACE" ) ) {
                        /// 1.2 Create Surface from the name of its parent Interface
                        /// and its geological feature
                        index_t id = in.field_as_uint( 1 ) ;
                        std::string geol = in.field( 2 ) ;
                        index_t f = 3 ;
                        std::ostringstream oss ;
                        do {
                            oss << in.field( f++ ) ;
                        } while( f < in.nb_fields() ) ;
                        std::string interface_name = oss.str() ;

                        // And its key facet that give the orientation of the surface part
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 p0( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 p1( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 p2( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;

                        create_surface( interface_name, geol, p0, p1, p2 ) ;
                        nb_tface++ ;
                    } else if( in.field_matches( 0, "REGION" ) ) {
                        /// 1.3 Read Region information and create them from their name,
                        /// and the surfaces on their boundary
                        index_t id = in.field_as_uint( 1 ) ;
                        std::string name = in.field( 2 ) ;

                        std::vector< std::pair< index_t, bool > > region_boundaries ;
                        bool end_region = false ;
                        while( !end_region ) {
                            in.get_line() ;
                            in.get_fields() ;
                            for( index_t i = 0; i < 5; ++i ) {
                                int s = in.field_as_int( i ) ;
                                if( s == 0 ) {
                                    end_region = true ;
                                    break ;
                                }
                                bool side = s > 0 ;
                                if( s > 0 ) {
                                    s -= 1 ;
                                } else {
                                    s = -s - 1 ;
                                }
                                region_boundaries.push_back(
                                    std::pair< index_t, bool >( s, side ) ) ;
                            }
                        }

                        // The Universe is not a regular region
                        if( name == "Universe" ) {
                            set_universe( region_boundaries ) ;
                        } else {
                            // Create the regions and set its boundaries 
                            bme_t region_id = create_element( BME::REGION ) ;
                            set_element_name( region_id, name ) ;
                            for( index_t i = 0; i < region_boundaries.size(); ++i ) {
                                add_element_boundary( region_id,
                                    bme_t( BME::SURFACE, region_boundaries[ i ].first ),
                                    region_boundaries[ i ].second ) ;
                            }
                        }
                    } else if( in.field_matches( 0, "LAYER" ) ) {
                        /// 1.4 Build the volumetric layers from their name and
                        /// the ids of the regions they contain
                        bme_t layer_id = create_element( BME::LAYER ) ;
                        set_element_name( layer_id, in.field( 1 ) ) ;
                        bool end_layer = false ;
                        while( !end_layer ) {
                            in.get_line() ;
                            in.get_fields() ;
                            for( index_t i = 0; i < 5; ++i ) {
                                index_t region_id = in.field_as_uint( i ) ;
                                if( region_id == 0 ) {
                                    end_layer = true ;
                                    break ;
                                } else {
                                    region_id -= nb_tface + 1 ; // Remove Universe region
                                    // Correction because ids begin at 1 in the file
                                    add_child( layer_id,
                                        bme_t( BME::REGION, region_id - 1 ) ) ;
                                }
                            }
                        }
                    } else if( in.field_matches( 0, "END" ) ) {
                        // End of the high level information on the model
                        // Switch to reading the geometry of the model surfaces
                        read_model = false ;
                        continue ;
                    }
                } else {
                    if( in.field_matches( 0, "GOCAD" ) ) {
                        // This is the beginning of a new TSurf = Interface
                        tsurf_count++ ;
                    }
                    if( in.field_matches( 0, "ZPOSITIVE" ) ) {
                        if( in.field_matches( 1, "Elevation" ) ) {
                            z_sign = 1 ;
                        } else if( in.field_matches( 1, "Depth" ) ) {
                            z_sign = -1 ;
                        } else {
                            ringmesh_assert_not_reached;}
                    } else if( in.field_matches( 0, "END" ) ) {
                        // This the END of a TSurf
                        if( tsurf_count > 0 ) {
                            // End the last TFace - Surface of this TSurf
                            set_surface_geometry(
                                bme_t( BME::SURFACE, tface_count - 1 ),
                                std::vector< vec3 >(
                                    tsurf_vertices.begin() +
                                    tface_vertex_start.back(),
                                    tsurf_vertices.end() ),
                                tface_facets,
                                tface_facets_ptr ) ;

                            if( !check_key_facet_orientation( tface_count - 1 ) ) {
                                change_key_facet.push_back( tface_count - 1 ) ;
                            }

                            tface_facets.clear() ;
                            tface_facets_ptr.clear() ;
                            tface_facets_ptr.push_back( 0 ) ;

                            // End this TSurf - Interface
                            nb_tface_in_prev_tsurf += tface_vertex_start.size() ;
                            tsurf_vertices.clear() ;
                            tface_vertex_start.clear() ;
                        }
                    } else if( in.field_matches( 0, "TFACE" ) ) {
                        // Beginning of a new TFace - Surface
                        if( tface_vertex_start.size() > 0 ) {
                            // End the previous TFace - Surface  (copy from line 1180)
                            set_surface_geometry(
                                bme_t( BME::SURFACE, tface_count - 1),
                                std::vector< vec3 >(
                                    tsurf_vertices.begin() +
                                    tface_vertex_start.back(),
                                    tsurf_vertices.end() ),
                                tface_facets,
                                tface_facets_ptr ) ;

                            if( !check_key_facet_orientation( tface_count - 1 ) ) {
                                change_key_facet.push_back( tface_count - 1 ) ;
                            }

                            tface_facets.clear() ;
                            tface_facets_ptr.clear() ;
                            tface_facets_ptr.push_back( 0 ) ;
                        }

                        // Register where begin the new TFace vertices
                        tface_vertex_start.push_back( tsurf_vertices.size() ) ;

                        tface_count++ ;
                    }

                    /// 2.1 Read the surface vertices and facets (only triangles in Gocad Model3d files)
                    else if( in.field_matches( 0,
                            "VRTX" ) || in.field_matches( 0, "PVRTX" ) )
                    {
                        const vec3 p( read_double( in, 2 ), read_double( in,
                                3 ), z_sign * read_double( in, 4 ) ) ;
                        tsurf_vertices.push_back( p ) ;
                    } else if( in.field_matches( 0,
                            "PATOM" ) || in.field_matches( 0, "ATOM" ) )
                    {
                        tsurf_vertices.push_back( tsurf_vertices[ in.
                            field_as_uint(
                                2 ) - 1 ] ) ;
                    } else if( in.field_matches( 0, "TRGL" ) ) {
                        // Read ids of the vertices of each triangle in the TSurf
                        // and switch to ids in the TFace
                        tface_facets.push_back( (index_t) in.field_as_uint(
                                1 ) - tface_vertex_start.back() - 1 ) ;
                        tface_facets.push_back( (index_t) in.field_as_uint(
                                2 ) - tface_vertex_start.back() - 1 ) ;
                        tface_facets.push_back( (index_t) in.field_as_uint(
                                3 ) - tface_vertex_start.back() - 1 ) ;
                        tface_facets_ptr.push_back( tface_facets.size() ) ;
                    }

                    // 2.2 Build the corners from their position and the surface parts
                    //    containing them
                    else if( in.field_matches( 0, "BSTONE" ) ) {
                        index_t v_id = in.field_as_uint( 1 ) - 1 ;
                        if( !find_corner(model_, tsurf_vertices[v_id]).is_defined() )
                        {
                            // Create the corner
                            set_corner( create_element( BME::CORNER ), tsurf_vertices[ v_id ] ) ;
                        }
                    }

                    /// 2.3 Read the Border information and store it
                    else if( in.field_matches( 0, "BORDER" ) ) {
                        index_t p1 = in.field_as_uint( 2 ) - 1 ;
                        index_t p2 = in.field_as_uint( 3 ) - 1 ;

                        // Get the global corner id
                        bme_t corner_id =
                        find_corner(model_, tsurf_vertices[ p1 ] ) ;
                        ringmesh_assert( corner_id.is_defined() ) ;

                        // Get the surface
                        index_t part_id = NO_ID ;
                        for( index_t i = 0 ; i < tface_vertex_start.size() ; ++i ) {
                            if( p1 < tface_vertex_start[ i ] ) {
                                ringmesh_assert( p2 < tface_vertex_start[ i ] ) ;

                                // Get vertices ids in the surface
                                p1 = p1 - tface_vertex_start[ i - 1 ] ;
                                p2 = p2 - tface_vertex_start[ i - 1 ] ;

                                // i-1 is the id of the TFace in this TSurf
                                part_id = i - 1 ;
                                break ;
                            }
                        }
                        if( part_id == NO_ID ) {
                            // It is in the last built Tface
                            p1 = p1 - tface_vertex_start[ tface_vertex_start.size() - 1 ] ;
                            p2 = p2 - tface_vertex_start[ tface_vertex_start.size() - 1 ] ;

                            part_id = tface_vertex_start.size() - 1 ;
                        }

                        // The number of tfaces in previous tsurf is also to add
                        part_id += nb_tface_in_prev_tsurf ;

                        borders_to_build.push_back(
                            Border( part_id, corner_id.index, p1, p2 ) ) ;
                    }
                }
            }
        }

        // I agree that we do not need to compute the BoundaryModelVertices here
        // But perhaps the computation of Lines would be faster and safer - Jeanne

        /// 3. Build the Lines
        {
            std::vector< vec3 > line_vertices ;
            for( index_t i = 0; i < borders_to_build.size(); ++i ) {
                const Border& b = borders_to_build[i] ;

                // 1- Build the boundary : construct the vector
                // of vertices on the border
                const Surface& S = model_.surface( b.part_id_ ) ;

                bme_t end_corner_id = determine_line_vertices( S, b.p0_, b.p1_,
                    line_vertices ) ;

                // 2 - Check if this border already exists
                bme_t line_id = find_or_create_line( *this, line_vertices ) ;

                // Add the surface in which this line is
                add_element_in_boundary( line_id,
                    bme_t( BME::SURFACE, b.part_id_ ) ) ;
            }
        }

        /// 4. Build the Contacts
        build_contacts() ;

        // Modify in the Region the side of the Surface for which the key facet
        // orientation was not the same than their facet orientations
        for( index_t i = 0; i < change_key_facet.size(); i++ ) {
            const Surface& S = model_.surface( change_key_facet[i] ) ;
            for( index_t j = 0; j < S.nb_in_boundary(); ++j ) {
                BoundaryModelElement& R = element( S.in_boundary_id( j ) ) ;
                for( index_t b = 0; b < R.nb_boundaries(); ++b ) {
                    if( R.boundary_id( b ).index == change_key_facet[i] ) {
                        bool old_side = R.side( b ) ;
                        R.set_boundary( b, R.boundary_id( b ), !old_side ) ;
                    }
                }
            }
        }
        
        /// 5. Fill missing information and check model validity
        bool valid_model = end_model() ;

        time( &end_load ) ;        
        // Output of loading time only in debug mode has no meaning (JP)
        GEO::Logger::out("I/O") << "Model loading time "
            << difftime( end_load, start_load ) << " sec" << std::endl ;
        
        return valid_model ;
    }

    /*!
     * @brief Find the facet which first 3 vertices are given
     * 
     * @param[in] surface_id Index of the surface
     * @param[in] p0 First point coordinates
     * @param[in] p1 Second point coordinates
     * @param[in] p2 Third point coordinates
     * @param[out] same_sign Is true if the found facet has the same orientation than triangle p0p1p2
     * @return Index of the found facet, NO_ID if none found
     */
    index_t BoundaryModelBuilderGocad::find_key_facet(
        index_t surface_id,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        bool& same_sign ) const
    {
        const Surface& surface = model_.surface( surface_id ) ;
        same_sign = false ;

        for( index_t t = 0; t < surface.nb_cells(); ++t ) {
            const vec3& pp0 = surface.vertex( t, 0 ) ;
            const vec3& pp1 = surface.vertex( t, 1 ) ;
            const vec3& pp2 = surface.vertex( t, 2 ) ;

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
     * 
     * @param[in] surface_id Index of the surface
     * @return False if the key_facet orientation is not the same than the surface facets, else true.
     */
    bool BoundaryModelBuilderGocad::check_key_facet_orientation(
        index_t surface_id ) const
    {
        const KeyFacet& key_facet = key_facets_[surface_id] ;

        const vec3& p0 = key_facet.p0_ ;
        const vec3& p1 = key_facet.p1_ ;
        const vec3& p2 = key_facet.p2_ ;
        bool same_sign = false ;

        index_t t = find_key_facet( surface_id, p0, p1, p2, same_sign ) ;
        if( t == NO_ID ) {
            vec3 p00( p0 ) ;
            p00.z *= -1 ;
            vec3 p10( p1 ) ;
            p10.z *= -1 ;
            vec3 p20( p2 ) ;
            p20.z *= -1 ;

            // It is because of the sign of Z that is not the same
            t = find_key_facet( surface_id, p00, p10, p20, same_sign ) ;
        }
        ringmesh_assert( t != NO_ID ) ;
        return same_sign ;
    }

    /*!
     * @brief Get the points of a Line between two corners on a Surface
     *
     * @param[in] S Index of the surface
     * @param[in] id0 Index of the starting point( a corner ) in S
     * @param[in] id1 Index of the second point on the Line in S
     * @param[out] border_vertex_model_vertices Coordinates of the vertices on the Line (emptied and filled again)
     * @return Index of the Corner at which the Line ends
     */
    bme_t BoundaryModelBuilderGocad::determine_line_vertices(
        const Surface& S,
        index_t id0,
        index_t id1,
        std::vector< vec3 >& border_vertex_model_vertices ) const
    {
        ringmesh_debug_assert( id0 < S.nb_vertices() && id1 < S.nb_vertices() ) ;

        border_vertex_model_vertices.resize( 0 ) ;

        // Starting facet that contains the two given vertices
        index_t f = S.facet_from_surface_vertex_ids( id0, id1 ) ;
        ringmesh_assert( f != Surface::NO_ID ) ;

        vec3 p0 = S.vertex( id0 ) ;
        vec3 p1 = S.vertex( id1 ) ;

        border_vertex_model_vertices.push_back( p0 ) ;
        border_vertex_model_vertices.push_back( p1 ) ;

        bme_t p1_corner = find_corner( model(), p1 ) ;
        while( !p1_corner.is_defined() ) {
            index_t next_f = NO_ID ;
            index_t id1_in_next = NO_ID ;
            index_t next_id1_in_next = NO_ID ;

            // We want the next triangle that is on the boundary and share p1
            // If there is no such triangle, the third vertex of the current triangle is to add
            S.next_on_border( f, S.facet_vertex_id( f, id0 ),
                S.facet_vertex_id( f, id1 ), next_f, id1_in_next,
                next_id1_in_next ) ;

            ringmesh_assert(
                next_f != NO_ID && id1_in_next != NO_ID
                    && next_id1_in_next != NO_ID ) ;

            index_t next_id1 = S.surf_vertex_id( next_f, next_id1_in_next ) ;

            // Update
            f = next_f ;
            id0 = id1 ;
            id1 = next_id1 ;

            p1 = S.vertex( next_id1 ) ;
            border_vertex_model_vertices.push_back( p1 ) ;
            p1_corner = find_corner( model(), p1 ) ;
        }
        return p1_corner ;
    }

    /*!
     * @brief Get the points of a Line between two corners on a Surface
     *
     * WE ASSUME THAT THE MODEL VERTICES ARE AVAILABLE AND CORRECT
     *
     * @param[in] S Index of the surface
     * @param[in] id0 Index of the starting point( a corner ) in S
     * @param[in] id1 Index of the second point on the Line in S
     * @param[out] border_vertex_model_ids Indices of vertices on the Line (resized at 0 at the beginning)
     * @return Index of the Corner at which the Line ends
     */
    bme_t BoundaryModelBuilderGocad::determine_line_vertices(
        const Surface& S,
        index_t id0,
        index_t id1,
        std::vector< index_t >& border_vertex_model_ids ) const
    {
        ringmesh_debug_assert( id0 < S.nb_vertices() && id1 < S.nb_vertices() ) ;

        border_vertex_model_ids.resize( 0 ) ;

        // Starting facet that contains the two given vertices
        index_t f = S.facet_from_surface_vertex_ids( id0, id1 ) ;
        ringmesh_assert( f != Surface::NO_ID ) ;

        // Global ids at the model level
        index_t p0 = S.model_vertex_id( id0 ) ;
        index_t p1 = S.model_vertex_id( id1 ) ;

        border_vertex_model_ids.push_back( p0 ) ;
        border_vertex_model_ids.push_back( p1 ) ;

        bme_t p1_corner = find_corner( model(), p1 ) ;
        while( !p1_corner.is_defined() ) {
            index_t next_f = NO_ID ;
            index_t id1_in_next = NO_ID ;
            index_t next_id1_in_next = NO_ID ;

            // We want the next triangle that is on the boundary and share p1
            // If there is no such triangle, the third vertex of the current triangle is to add
            S.next_on_border( f, S.facet_vertex_id( f, id0 ),
                S.facet_vertex_id( f, id1 ), next_f, id1_in_next,
                next_id1_in_next ) ;

            ringmesh_assert(
                next_f != NO_ID && id1_in_next != NO_ID
                    && next_id1_in_next != NO_ID ) ;

            index_t next_id1 = S.surf_vertex_id( next_f, next_id1_in_next ) ;

            // Update
            f = next_f ;
            id0 = id1 ;
            id1 = next_id1 ;

            p1 = S.model_vertex_id( next_id1 ) ;
            border_vertex_model_ids.push_back( p1 ) ;
            p1_corner = find_corner( model(), p1 ) ;
        }
        return p1_corner ;
    }

    /*!
     * @brief Build the Contacts
     * @details One contact is a group of lines shared by the same Interfaces
     */
    void BoundaryModelBuilderGocad::build_contacts()
    {
        std::vector< std::set< bme_t > > interfaces ;
        for( index_t i = 0; i < model_.nb_lines(); ++i ) {
            const Line& L = model_.line( i ) ;
            std::set< bme_t > cur_interfaces ;
            for( index_t j = 0; j < L.nb_in_boundary(); ++j ) {
                cur_interfaces.insert( model_.element( 
                    L.in_boundary_id(j) ).parent().bme_id() ) ;
            }            
            bme_t contact_id ;
            for( index_t j = 0; j < interfaces.size(); ++j ) {
                if( cur_interfaces.size() == interfaces[j].size() && 
                    std::equal(cur_interfaces.begin(), 
                      cur_interfaces.end(), interfaces[j].begin() )
                   ) {
                    contact_id = bme_t( BME::CONTACT, j ) ;
                    break ;
                }
            }
            if( !contact_id.is_defined() ) {
                contact_id = create_element( BME::CONTACT ) ;
                ringmesh_debug_assert( contact_id.index == interfaces.size() ) ;
                interfaces.push_back( cur_interfaces ) ;
                // Create a name for this contact
                std::string name = "contact_" ;
                for( std::set< bme_t >::const_iterator it( cur_interfaces.begin() );
                     it != cur_interfaces.end(); ++it 
                    ) {
                    name += model_.element( *it ).name() ;
                    name += "_" ;
                }
                set_element_name( contact_id, name ) ;
            }
            add_child( contact_id, bme_t( BME::LINE, i ) ) ;
        }
    }

    /*!
     * @brief Add a Surface to the model
     *
     * @param[in] interface_name Name of the parent. The parent MUST exist.
     * @param[in] type Type of the Surface
     * @param[in] p0 Coordinates of the 1 point of the TFace key facet 
     * @param[in] p1 Coordinates of the 2 point of the TFace key facet 
     * @param[in] p2 Coordinates of the 3 point of the TFace key facet 
     */
    void BoundaryModelBuilderGocad::create_surface(
        const std::string& interface_name,
        const std::string& type,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2 )
    {
        bme_t parent = find_interface( model_, interface_name ) ;
        if( interface_name != "" ) {
            ringmesh_assert( parent.is_defined() ) ;
        }

        bme_t id = create_element( BME::SURFACE ) ;
        set_parent( id, parent ) ;
        set_element_geol_feature( parent, BME::determine_geological_type( type ) ) ;
        key_facets_.push_back( KeyFacet( p0, p1, p2 ) ) ;
    }

    bool BoundaryModelBuilderBM::load_file( const std::string& bm_file_name )
    {
        GEO::LineInput in( bm_file_name ) ;
        if( !in.OK() ) {
            return false ;
        }
        while( !in.eof() && in.get_line() ) {
            in.get_fields() ;
            if( in.nb_fields() > 0 ) {
                // Name of the model
                if( in.field_matches( 0, "NAME" ) ) {
                    if( in.nb_fields() > 1 ) {
                        set_model_name( in.field( 1 ) ) ;
                    }
                }
                // Number of elements of a given type
                else if( match_nb_elements( in.field( 0 ) ) != BME::NO_TYPE ) {
                    // Allocate the space
                    if( in.nb_fields() > 1 ) {
                        resize_elements( match_nb_elements( in.field( 0 ) ),
                            in.field_as_uint( 1 ) ) ;
                    }
                }

                // High-level elements
                else if( match_high_level_type( in.field( 0 ) ) ) {
                    // Read this element
                    // First line : type - id - name - geol_feature
                    if( in.nb_fields() < 4 ) {
                        GEO::Logger::err("I/O") << "Invalid line: " << in.line_number()
                            << "4 fields are expected, the type, id, name, and geological feature"
                            << std::endl ;
                        return false ;
                    }
                    BME::TYPE t = match_type( in.field( 0 ) ) ;
                    index_t id = in.field_as_uint( 1 ) ;
                    bme_t element( t, id ) ;
                    set_element_index( element ) ;
                    set_element_name( element, in.field( 2 ) ) ;
                    set_element_geol_feature( element,
                        BME::determine_geological_type( in.field( 3 ) ) ) ;
                    // Second line : indices of its children
                    in.get_line() ;
                    in.get_fields() ;
                    for( index_t c = 0; c < in.nb_fields(); c++ ) {
                        add_child( element,
                            bme_t( BME::child_type( t ),
                                in.field_as_uint( c ) ) ) ;
                    }
                }
                // Regions
                else if( match_type( in.field( 0 ) ) == BME::REGION ) {
                    // First line : type - id - name
                    if( in.nb_fields() < 3 ) {
                        GEO::Logger::err( "I/O" ) << "Invalid line: " << in.line_number()
                            << "3 fields are expected to describe a region: REGION, id, and name"
                            << std::endl ; 
                        return false ;
                    }
                    index_t id = in.field_as_uint( 1 ) ;
                    bme_t element( BME::REGION, id ) ;
                    set_element_index( element ) ;
                    set_element_name( element, in.field( 2 ) ) ;
                    // Second line : signed indices of boundaries
                    in.get_line() ;
                    in.get_fields() ;
                    for( index_t c = 0; c < in.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( in.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &in.field( c )[1], s ) ;

                        add_element_boundary( element, bme_t( BME::SURFACE, s ),
                            side ) ;
                    }
                }

                // Universe
                else if( in.field_matches( 0, "UNIVERSE" ) ) {
                    std::vector< std::pair< index_t, bool > > b_universe ;
                    // Second line: signed indices of boundaries
                    in.get_line() ;
                    in.get_fields() ;
                    for( index_t c = 0; c < in.nb_fields(); c++ ) {
                        bool side = false ;
                        if( strncmp( in.field( c ), "+", 1 ) == 0 ) {
                            side = true ;
                        }
                        index_t s ;
                        GEO::String::from_string( &in.field( c )[1], s ) ;

                        b_universe.push_back(
                            std::pair< index_t, bool >( s, side ) ) ;
                    }
                    set_universe( b_universe ) ;
                }

                // Model vertices
//                else if( in.field_matches( 0, "MODEL_VERTICES" ) ) {
//                    index_t nb_vertices = in.field_as_uint( 1 ) ;
//
//                    // Attributes
//                    in.get_line() ;
//                    in.get_fields() ;
//                    ringmesh_assert( in.field_matches( 0, "MODEL_VERTEX_ATTRIBUTES" ) ) ;
//                    index_t nb_attribs = ( in.nb_fields() - 1 ) / 2 ;
//                    std::vector< SerializedAttribute< BoundaryModel::VERTEX > >
//                    vertex_attribs( nb_attribs ) ;
//                    for( index_t i = 0; i < nb_attribs; i++ ) {
//                        vertex_attribs[ i ].bind(
//                            model_.vertex_attribute_manager(), in.field(
//                                1 + 2 * i ), in.field( 2 + 2 * i ), nb_vertices ) ;
//                    }
//                    for( index_t i = 0; i < nb_vertices; ++i ) {
//                        in.get_line() ;
//                        in.get_fields() ;
//                        add_vertex( vec3(
//                                read_double( in,
//                                    0 ), read_double( in, 1 ), read_double( in, 2 ) ) ) ;
//                        serialize_read_attributes( in, 3, i, vertex_attribs ) ;
//                    }
//                }

                // Corners
                else if( match_type( in.field( 0 ) ) == BME::CORNER ) {
                    // First line: CORNER - id - vertex id
                    if( in.nb_fields() < 5 ) {
                        GEO::Logger::err( "I/O" ) << "Invalid line: " << in.line_number()
                            << " 5 fields are expected to describe a corner: "
                            << " CORNER, index, and X, Y, Z coordinates "
                            << std::endl ;
                        return false ;
                    }
                    index_t id = in.field_as_uint( 1 ) ;
                    bme_t element( BME::CORNER, id ) ;
                    set_element_index( element ) ;
                    vec3 point( read_double( in, 2 ), read_double( in, 3 ),
                        read_double( in, 4 ) ) ;
                    set_element_vertex( element, 0, point ) ;
                }

                // Lines
                else if( match_type( in.field( 0 ) ) == BME::LINE ) {
                    index_t id = in.field_as_uint( 1 ) ;
                    bme_t cur_element( BME::LINE, id ) ;
                    Line& L = dynamic_cast< Line& >( element( cur_element ) ) ;
                    L.set_id( id ) ;

                    // Following information: vertices of the line
                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "LINE_VERTICES" ) ) ;
                    index_t nb_vertices = in.field_as_uint( 1 ) ;
                    std::vector< vec3 > vertices( nb_vertices ) ;
                    for( index_t i = 0; i < nb_vertices; i++ ) {
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 point( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;
                        vertices[i] = point ;
                    }

                    // Set the line points
                    L.set_vertices( vertices ) ;

                    // Attributes on line vertices
//                    in.get_line() ;
//                    in.get_fields() ;
//                    ringmesh_assert( in.field_matches( 0, "LINE_VERTEX_ATTRIBUTES" ) ) ;
//                    index_t nb_attribs = ( in.nb_fields() - 1 ) / 2 ;
//                    std::vector< SerializedAttribute > vertex_attribs(
//                        nb_attribs ) ;
//                    for( index_t i = 0; i < nb_attribs; i++ ) {
//                        vertex_attribs[ i ].bind(
//                            L.vertex_attribute_manager(), in.field(
//                                1 + 2 * i ), in.field( 2 + 2 * i ), nb_vertices ) ;
//                    }

                    // Read the vertices indices and attributes on vertices
//                    for( index_t i = 0; i < nb_vertices; i++ ) {
//                        in.get_line() ;
//                        in.get_fields() ;
//                        serialize_read_attributes( in, 1, i, vertex_attribs ) ;
//                    }

                    // Read attributes on line segments
//                    in.get_line() ;
//                    in.get_fields() ;
//                   ringmesh_assert( in.field_matches( 0, "LINE_SEGMENT_ATTRIBUTES" ) ) ;
//                    index_t nb_segment_attribs = ( in.nb_fields() - 1 ) / 2 ;
//                    if( nb_segment_attribs > 0 ) {
//                        std::vector< SerializedAttribute< BME::FACET > >
//                        segment_attribs( nb_segment_attribs ) ;
//                        for( index_t i = 0; i < nb_segment_attribs; i++ ) {
//                            segment_attribs[ i ].bind(
//                                L.facet_attribute_manager(), in.field(
//                                    1 + 2 * i ), in.field( 2 + 2 * i ), L.nb_cells() ) ;
//                        }
//                        for( index_t i = 0; i < L.nb_cells(); i++ ) {
//                            in.get_line() ;
//                            in.get_fields() ;
//                            serialize_read_attributes( in, 1, in.field_as_uint(
//                                    0 ), segment_attribs ) ;
//                        }
//                    }

                    // Set the corners - they can be the same
                    add_element_boundary( cur_element,
                        find_corner( model(), vertices.front() ) ) ;
                    add_element_boundary( cur_element,
                        find_corner( model(), vertices.back() ) ) ;

                    // Finally we have the in_boundary information
                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "IN_BOUNDARY" ) ) ;
                    for( index_t b = 1; b < in.nb_fields(); b++ ) {
                        L.add_in_boundary(
                            bme_t( BME::SURFACE, in.field_as_uint( b ) ) ) ;
                    }
                }

                // Surfaces
                else if( match_type( in.field( 0 ) ) == BME::SURFACE ) {
                    index_t id = in.field_as_uint( 1 ) ;
                    bme_t cur_element( BME::SURFACE, id ) ;
                    Surface& S = dynamic_cast< Surface& >( element( cur_element ) ) ;
                    S.set_id( id ) ;

                    // Read the surface vertices and their attributes
                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "SURFACE_VERTICES" ) ) ;
                    index_t nb_vertices = in.field_as_uint( 1 ) ;
                    std::vector< vec3 > vertices( nb_vertices ) ;
                    for( index_t i = 0; i < nb_vertices; i++ ) {
                        in.get_line() ;
                        in.get_fields() ;
                        vec3 point( read_double( in, 0 ), read_double( in, 1 ),
                            read_double( in, 2 ) ) ;
                        vertices[i] = point ;
                    }

//                    in.get_line() ;
//                    in.get_fields() ;
//                    ringmesh_assert( in.field_matches( 0,
//                            "SURFACE_VERTEX_ATTRIBUTES" ) ) ;
//                    index_t nb_vertex_attribs = ( in.nb_fields() - 1 ) / 2 ;
//
                    // Bind the vertex attributes
//                    std::vector< SerializedAttribute< BME::VERTEX > > vertex_attribs(
//                        nb_vertex_attribs ) ;
//                    for( index_t i = 0; i < nb_vertex_attribs; i++ ) {
//                        vertex_attribs[ i ].bind(
//                            S.vertex_attribute_manager(), in.field(
//                                1 + 2 * i ), in.field( 2 + 2 * i ), nb_vertices ) ;
//                    }

                    // Read the vertices global ids and attributes
//                    for( index_t i = 0; i < nb_vertices; i++ ) {
//                        in.get_line() ;
//                        in.get_fields() ;
//                        serialize_read_attributes( in, 1, i, vertex_attribs ) ;
//                    }

                    // Read the surface facets
                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "SURFACE_CORNERS" ) ) ;
                    index_t nb_corners = in.field_as_uint( 1 ) ;

                    in.get_line() ;
                    in.get_fields() ;
                    ringmesh_assert( in.field_matches( 0, "SURFACE_FACETS" ) ) ;
                    index_t nb_facets = in.field_as_uint( 1 ) ;

//                    in.get_line() ;
//                    in.get_fields() ;
//                    ringmesh_assert( in.field_matches( 0, "SURFACE_FACET_ATTRIBUTES" ) ) ;
//                    index_t nb_facet_attribs = ( in.nb_fields() - 1 ) / 2 ;

                    // Bind the facet attributes
//                    std::vector< SerializedAttribute< BME::FACET > > facet_attribs(
//                        nb_facet_attribs ) ;
//                    for( index_t i = 0; i < nb_facet_attribs; i++ ) {
//                        facet_attribs[ i ].bind(
//                            S.facet_attribute_manager(), in.field(
//                                1 + 2 * i ), in.field( 2 + 2 * i ), nb_facets ) ;
//                    }

                    std::vector< index_t > corners( nb_corners ) ;
                    std::vector< index_t > facet_ptr( nb_facets + 1, 0 ) ;
                    index_t count_facets = 0 ;
                    for( index_t f = 0; f < nb_facets; f++ ) {
                        in.get_line() ;
                        in.get_fields() ;
                        index_t nb_v = in.field_as_uint( 0 ) ;
                        for( index_t v = 0; v < nb_v; ++v ) {
                            corners[count_facets + v] = in.field_as_uint( v + 1 ) ;
                        }
                        count_facets += nb_v ;
                        facet_ptr[f + 1] = count_facets ;
//                        serialize_read_attributes( in, nb_v + 1, f, facet_attribs ) ;
                    }

                    S.set_geometry( vertices, corners, facet_ptr ) ;
                    set_surface_adjacencies( cur_element ) ;
                }
            }
        }
        if( !end_model() ) {
            std::cout << "Invalid BoundaryModel loaded" << std::endl ;
        }
        return true ;
    }

    BoundaryModelElement::TYPE BoundaryModelBuilderBM::match_nb_elements(
        const char* s )
    {
        // Check that the first 3 characters are NB_
        if( strncmp( s, "NB_", 3 ) != 0 ) {
            return BME::NO_TYPE ;
        } else {
            for( index_t i = BME::CORNER; i < BME::NO_TYPE; i++ ) {
                BME::TYPE type = (BME::TYPE) i ;
                if( strstr( s, BME::type_name( type ).data() ) != NULL ) {
                    return type ;
                }
            }
            return BME::NO_TYPE ;
        }
    }

    BoundaryModelElement::TYPE BoundaryModelBuilderBM::match_type( const char* s )
    {
        for( index_t i = BME::CORNER; i < BME::NO_TYPE; i++ ) {
            BME::TYPE type = (BME::TYPE) i ;
            if( strcmp( s, BME::type_name( type ).data() ) == 0 ) {
                return type ;
            }
        }
        return BME::NO_TYPE ;
    }

    /*
    * @brief Utility structure to build a BoundaryModel knowing only its surface
    * @details Store the vertices of a triangle that is on the boundary of a surface
    */
    struct BorderTriangle {
        /*!
        * @brief Constructor
        * @param s Index of the surface
        * @param f Index of the facet containing the 3 vertices
        * @param vi Indices in the BoundaryModel of the vertices defining the triangle
        *           the edge v0 - v1 is the one on the boundary
        */
        BorderTriangle(
            index_t s,
            index_t f,
            index_t v0,
            index_t v1,
            index_t v2 )
            : s_( s ), f_( f ), v0_( v0 ), v1_( v1 ), v2_( v2 )
        {}

        bool operator<( const BorderTriangle& rhs ) const
        {
            if( std::min( v0_, v1_ ) != std::min( rhs.v0_, rhs.v1_ ) ) {
                return std::min( v0_, v1_ ) < std::min( rhs.v0_, rhs.v1_ ) ;
            }
            if( std::max( v0_, v1_ ) != std::max( rhs.v0_, rhs.v1_ ) ) {
                return std::max( v0_, v1_ ) < std::max( rhs.v0_, rhs.v1_ ) ;
            }
            if( s_ != rhs.s_ ) {
                return s_ < rhs.s_ ;
            }
            if( f_ != rhs.f_ ) {
                return f_ < rhs.f_ ;
            }
            return rhs.v2_ == index_t( -1 ) ? false : v2_ < rhs.v2_ ;
        }

        bool same_edge( const BorderTriangle& rhs ) const
        {
            return std::min( v0_, v1_ ) == std::min( rhs.v0_, rhs.v1_ ) &&
                std::max( v0_, v1_ ) == std::max( rhs.v0_, rhs.v1_ ) ;
        }

        /// Indices of the points in the model. Triangle has the Surface orientation
        /// The edge v0v1 is, in surface s_, on the border.
        index_t v0_ ;
        index_t v1_ ;
        index_t v2_ ;

        // Index of the model surface
        index_t s_ ;

        // Index of the facet in the surface
        index_t f_ ;
    } ;

    /*!
    * @brief Get the BorderTriangle corresponding to the next edge on border
    * in the considered Surface
    */
    index_t get_next_border_triangle(
        const BoundaryModel& M,
        const std::vector< BorderTriangle >& BT,
        index_t from,
        bool backward = false )
    {
        const BorderTriangle& in = BT[ from ] ;
        const Surface& S = M.surface( in.s_ ) ;
        index_t NO_ID( -1 ) ;

        // Get the next edge on border in the Surface
        index_t f = in.f_ ;
        index_t f_v0 = S.facet_id_from_model( f, in.v0_ ) ;
        index_t f_v1 = S.facet_id_from_model( f, in.v1_ ) ;
        ringmesh_assert( f_v0 != NO_ID && f_v1 != NO_ID ) ;

        index_t next_f = NO_ID ;
        index_t next_f_v0 = NO_ID ;
        index_t next_f_v1 = NO_ID ;

        if( !backward ) {
            S.next_on_border( f, f_v0, f_v1, next_f, next_f_v0,
                              next_f_v1 ) ;
        } else {
            S.next_on_border( f, f_v1, f_v0, next_f, next_f_v0,
                              next_f_v1 ) ;
        }

        // Find the BorderTriangle that is correspond to this
        // It must exist and there is only one
        BorderTriangle bait( in.s_, next_f, S.model_vertex_id( next_f, next_f_v0 ),
                             S.model_vertex_id( next_f, next_f_v1 ), NO_ID ) ;

        // lower_bound returns an iterator pointing to the first element in the range [first,last)
        // which does not compare less than the given val.
        // See operator< on BorderTriangle
        index_t result = narrow_cast< index_t >( std::lower_bound( BT.begin(), BT.end(), bait ) - BT.begin() );

        ringmesh_assert( result < BT.size() ) ;
        return result ;
    }


    /*!
    * @brief Mark as visited all BorderTriangle which first edge is the same than
    * the first edge of i.
    *
    * @param[in] border_triangles Information on triangles MUST be sorted so that
    *            BorderTriangle having the same boundary edge are adjacent
    *
    */
    void visit_border_triangle_on_same_edge(
        const std::vector< BorderTriangle >& border_triangles,
        index_t i,
        std::vector< bool >& visited )
    {
        index_t j = i ;
        while( j < border_triangles.size() &&
               border_triangles[ i ].same_edge( border_triangles[ j ] ) ) {
            visited[ j ] = true ;
            j++ ;
        }
        signed_index_t k = i - 1 ;
        while( k > -1 &&
               border_triangles[ i ].same_edge( border_triangles[ k ] ) ) {
            visited[ k ] = true ;
            k-- ;
        }
    }


    /*!
    * @brief Get the indices of the Surface adjacent to the first edge of a BorderTriangle
    *
    * @param[in] border_triangles Information on triangles MUST be sorted so that
    *          BorderTriangle having the same boundary edge are adjacent
    * @param[in] i Index of the BorderTriangle
    * @param[out] adjacent_surfaces Indices of the Surface stored by the BorderTriangle sharing
    *             the first edge of i
    */
    void get_adjacent_surfaces(
        const std::vector< BorderTriangle >& border_triangles,
        index_t i,
        std::vector< index_t >& adjacent_surfaces )
    {
        adjacent_surfaces.resize( 0 ) ;

        index_t j = i ;
        while( j < border_triangles.size() &&
               border_triangles[ i ].same_edge( border_triangles[ j ] ) ) {
            adjacent_surfaces.push_back( border_triangles[ j ].s_ ) ;
            j++ ;
        }

        signed_index_t k = i - 1 ;
        while( k > -1 &&
               border_triangles[ i ].same_edge( border_triangles[ k ] ) ) {
            adjacent_surfaces.push_back( border_triangles[ k ].s_ ) ;
            k-- ;
        }

        // Sort the adjacent surfaces
        std::sort( adjacent_surfaces.begin(), adjacent_surfaces.end() ) ;

        // When the surface appear twice (the line is an internal border)
        // we keep both occurrence, otherwise this connectivity info is lost
    }


    /*!
    * @brief From a BoundaryModel in which only Surface are defined, create
    * corners, contacts and regions.
    *
    */
    bool BoundaryModelBuilderSurface::build_model( bool build_regions )
    {
        if( model_.nb_surfaces() == 0 ) {
            GEO::Logger::err( "BoundaryModel" )
                << "No surface to build the model " << std::endl ;
            return false ;
        }

        /// 1. Initialize model_ global vertices and backward information
        model_.vertices.nb() ;
        model_.vertices.bme_vertices( 0 ) ;

        /// 2.1 Get for all Surface, the triangles that have an edge
        /// on the boundary.
        std::vector< BorderTriangle > border_triangles ;
        for( index_t i = 0; i < model_.nb_surfaces(); ++i ) {
            const Surface& S = model_.surface( i ) ;
            for( index_t j = 0; j < S.nb_cells(); ++j ) {
                for( index_t v = 0; v < S.nb_vertices_in_facet( j ); ++v ) {
                    if( S.is_on_border( j, v ) ) {
                        border_triangles.push_back( BorderTriangle( i, j,
                            S.model_vertex_id( j, v ),
                            S.model_vertex_id( j, S.next_in_facet( j, v ) ),
                            S.model_vertex_id( j, S.prev_in_facet( j, v ) )
                            ) ) ;
                    }
                }
            }
        }

        /// 2.2 Sort these triangles so that triangles sharing the same edge follow one another
        std::sort( border_triangles.begin(), border_triangles.end() ) ;

        /// 3. Build the Lines and gather information to build the regions
        std::vector< SortTriangleAroundEdge > regions_info ;

        // The goal is to visit all BorderTriangle and propagate to get each Line vertices
        std::vector< bool > visited( border_triangles.size(), false ) ;
        for( index_t i = 0; i < border_triangles.size(); ++i ) {
            if( !visited[ i ] ) {
                // This is a new Line
                std::vector< index_t > vertices ;

                // Get the indices of the Surfaces around this Line
                std::vector< index_t > adjacent ;
                get_adjacent_surfaces( border_triangles, i, adjacent ) ;

                // Mark as visited the BorderTriangle around the same first edge
                visit_border_triangle_on_same_edge( border_triangles, i, visited ) ;

                // Gather information to sort triangles around the contact
                regions_info.push_back( SortTriangleAroundEdge() ) ;
                index_t j = i ;
                while( j < border_triangles.size() &&
                       border_triangles[ i ].same_edge( border_triangles[ j ] ) ) {
                    regions_info.back().add_triangle(
                        border_triangles[ j ].s_,
                        model_.vertices.unique_vertex( border_triangles[ j ].v0_ ),
                        model_.vertices.unique_vertex( border_triangles[ j ].v1_ ),
                        model_.vertices.unique_vertex( border_triangles[ j ].v2_ )
                        ) ;
                    j++ ;
                }

                // Add vertices to the Line
                vertices.push_back( border_triangles[ i ].v0_ ) ;
                vertices.push_back( border_triangles[ i ].v1_ ) ;

                // Build the contact propating forward on the border of the Surface
                // While the adjacent surfaces are the same the vertices the next edge on the
                // boundary of the Surface are added
                bool same_surfaces = true ;
                index_t next_i = get_next_border_triangle( 
                    model_, border_triangles, i ) ;
                do {
                    ringmesh_assert( next_i != NO_ID ) ;
                    if( !visited[ next_i ] ) {
                        std::vector< index_t > adjacent_next ;
                        get_adjacent_surfaces( border_triangles, next_i,
                                               adjacent_next ) ;

                        if( adjacent.size() == adjacent_next.size() &&
                            std::equal( adjacent.begin(), adjacent.end(),
                            adjacent_next.begin() )
                            ) {
                            visit_border_triangle_on_same_edge(
                                border_triangles, next_i, visited ) ;

                            // Add the next vertex
                            if( border_triangles[ next_i ].v0_ ==
                                vertices.back() ) {
                                vertices.push_back( border_triangles[ next_i ].v1_ ) ;
                            } else {
                                ringmesh_assert(
                                    border_triangles[ next_i ].v1_ == vertices.back() ) ;
                                vertices.push_back( border_triangles[ next_i ].v0_ ) ;
                            }
                        } else {
                            same_surfaces = false ;
                        }
                    } else {
                        same_surfaces = false ;
                    }
                    next_i = get_next_border_triangle( model_, border_triangles,
                                                       next_i ) ;
                } while( same_surfaces && next_i != i ) ;

                if( next_i != i ) {
                    // Propagate backward to reach the other extremity
                    same_surfaces = true ;
                    index_t prev_i =
                        get_next_border_triangle( model_, border_triangles, i,
                        true ) ;
                    do {
                        ringmesh_assert( prev_i != NO_ID && prev_i != i ) ;
                        if( !visited[ prev_i ] ) {
                            std::vector< index_t > adjacent_prev ;
                            get_adjacent_surfaces( border_triangles, prev_i,
                                                   adjacent_prev ) ;

                            if( adjacent.size() == adjacent_prev.size() &&
                                std::equal( adjacent.begin(), adjacent.end(),
                                adjacent_prev.begin() )
                                ) {
                                visit_border_triangle_on_same_edge(
                                    border_triangles, prev_i, visited ) ;

                                // Fill the Line vertices
                                if( border_triangles[ prev_i ].v0_ ==
                                    vertices.front() ) {
                                    vertices.insert(
                                        vertices.begin(),
                                        border_triangles[ prev_i ].v1_ ) ;
                                } else {
                                    ringmesh_assert(
                                        border_triangles[ prev_i ].v1_ ==
                                        vertices.front() ) ;
                                    vertices.insert(
                                        vertices.begin(),
                                        border_triangles[ prev_i ].v0_ ) ;
                                }
                            } else {
                                same_surfaces = false ;
                            }
                        } else {
                            same_surfaces = false ;
                        }
                        prev_i =
                            get_next_border_triangle( model_, border_triangles,
                            prev_i,
                            true ) ;
                    } while( same_surfaces ) ;
                }
                ringmesh_assert( vertices.size() > 1 ) ;

                // At last create the Line
                bme_t l_id = create_element( BME::LINE ) ;
                set_line( l_id, vertices ) ;
                for( index_t j = 0; j < adjacent.size(); ++j ) {
                    add_element_in_boundary(
                        l_id, bme_t( BME::SURFACE, adjacent[ j ] ) ) ;
                }

                // Find or create the corners at line extremities
                bme_t c0 = find_corner( model(), vertices.front() ) ;
                if( !c0.is_defined() ) {
                    c0 = create_element( BME::CORNER ) ;
                    set_corner( c0, vertices.front() ) ;
                }
                add_element_boundary( l_id, c0 ) ;
                bme_t c1 = find_corner( model(), vertices.back() ) ;
                if( !c1.is_defined() ) {
                    c1 = create_element( BME::CORNER ) ;
                    set_corner( c1, vertices.back() ) ;
                }
                add_element_boundary( l_id, c1 ) ;              
            }
        }

        if( build_regions ) {
            /// 4. Build the regions

            // Complete boundary information for surfaces
            // Needed to compute volumetric regions
            fill_elements_boundaries( *this, BME::SURFACE ) ;

            /// 4.1 Sort surfaces around the contacts
            for( index_t i = 0; i < regions_info.size(); ++i ) {
                regions_info[ i ].sort() ;
            }

            if( model_.nb_surfaces() == 1 ) {
                if( model_.nb_lines() != 0 ) {
                    GEO::Logger::err( "BoundaryModel" )
                        << "The unique surface provided to build the model has boundaries "
                        << std::endl ;
                    return false ;
                } else {
                    /// If there is only one surface, its inside is set to be 
                    /// the + side. No check done
                    bool inside = true ;
                    // Create the region - set the surface on its boundaries
                    bme_t cur_region_id = create_element( BME::REGION ) ;
                    add_element_boundary( cur_region_id, bme_t( BME::SURFACE, 0 ),
                                          inside ) ;

                    // Create the universe region
                    std::vector< std::pair< index_t, bool > > univ_boundaries ;
                    univ_boundaries.push_back( std::pair<index_t, bool>( 0, !inside ) ) ;
                    set_universe( univ_boundaries ) ;
                }
            } else {
                // Each side of each Surface is in one Region( +side is first )
                std::vector< index_t > surf_2_region( 2 * model_.nb_surfaces(), NO_ID ) ;

                // Start with the first Surface on its + side
                std::stack< std::pair< index_t, bool > > S ;
                S.push( std::pair< index_t, bool >( 0, true ) ) ;

                while( !S.empty() ) {
                    std::pair< index_t, bool > cur = S.top() ;
                    S.pop() ;

                    // This side is already assigned
                    if( surf_2_region[ cur.second == true ? 2 * cur.first : 2 *
                        cur.first + 1 ] != NO_ID ) {
                        continue ;
                    }

                    // Create a new region
                    bme_t cur_region_id = create_element( BME::REGION ) ;

                    std::stack< std::pair< index_t, bool > > SR ;
                    SR.push( cur ) ;
                    while( !SR.empty() ) {
                        std::pair< index_t, bool > s = SR.top() ;
                        SR.pop() ;
                        index_t s_id = s.second == true ? 2 * s.first : 2 * s.first + 1 ;

                        // This side is already assigned
                        if( surf_2_region[ s_id ] != NO_ID ) {
                            continue ;
                        }

                        // Add the surface to the current region
                        add_element_boundary( cur_region_id, bme_t( BME::SURFACE, s.first ),
                                              s.second ) ;
                        surf_2_region[ s_id ] = cur_region_id.index ;

                        // Check the other side of the surface and push it in S
                        index_t s_id_opp = !s.second == true ? 2 * s.first : 2 *
                            s.first + 1 ;
                        if( surf_2_region[ s_id_opp ] == NO_ID ) {
                            S.push( std::pair< index_t, bool >( s.first, !s.second ) ) ;
                        }

                        // For each contact, push the next oriented surface that is in the same region
                        const BoundaryModelElement& surface = model_.surface( s.first ) ;
                        for( index_t i = 0; i < surface.nb_boundaries(); ++i ) {
                            const std::pair< index_t, bool >& n =
                                regions_info[ surface.boundary_id( i ).index ].next( s ) ;
                            index_t n_id = n.second == true ? 2 * n.first : 2 *
                                n.first + 1 ;

                            if( surf_2_region[ n_id ] == NO_ID ) {
                                SR.push( n ) ;
                            }
                        }
                    }
                }

                // Check if all the surfaces were visited
                // If not, this means that there are additionnal regions included in those built
                if( std::count( surf_2_region.begin(), surf_2_region.end(), NO_ID ) != 0 ) {
                    GEO::Logger::err( "BoundaryModel" )
                        << "Small bubble regions were skipped at model building "
                        << std::endl ;
                    // Or, most probably, we have a problem before
                    ringmesh_debug_assert( false ) ;
                }

                // We need to remove from the regions_ the one corresponding
                // to the universe_, the one with the biggest volume
                double max_volume = -1. ;
                index_t universe_id = NO_ID ;
                for( index_t i = 0; i < model_.nb_regions(); ++i ) {
                    double cur_volume = BoundaryModelElementMeasure::size( &model_.region( i ) ) ;
                    if( cur_volume > max_volume ) {
                        max_volume = cur_volume ;
                        universe_id = i ;
                    }
                }

                const BoundaryModelElement& cur_region = model_.region( universe_id ) ;
                std::vector< std::pair< index_t, bool > > univ_boundaries(
                    cur_region.nb_boundaries() ) ;
                for( index_t i = 0; i < cur_region.nb_boundaries(); ++i ) {
                    univ_boundaries[ i ].first = cur_region.boundary( i ).bme_id().index ;
                    univ_boundaries[ i ].second = cur_region.side( i ) ;
                }
                set_universe( univ_boundaries ) ;

                // Erase that region
                std::set< bme_t > to_erase ;
                to_erase.insert( bme_t( BME::REGION, universe_id ) ) ;
                remove_elements( to_erase ) ;
            }
        }

        // Deliberated clear of the model vertices 
        // To force their recomputation
        model_.vertices.clear() ;

        // Finish up the model and check its validity
        return end_model() ;
    }

} // namespace
