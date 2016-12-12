/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_GEOMETRY__
#define __RINGMESH_GEOMETRY__

#include <ringmesh/basic/common.h>

#include <geogram/points/kd_tree.h>

/*!
 * @file What
 * @author Arnaud Botella
 */

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {

    bool RINGMESH_API operator==( const vec3& u, const vec3& v ) ;
    bool RINGMESH_API operator<( const vec3& u, const vec3& v ) ;
    bool RINGMESH_API operator!=( const vec3& u, const vec3& v ) ;

    /* @warning Duplicate from Geogram/basic/numeric.h */
    enum Sign {
        NEGATIVE = -1, ZERO = 0, POSITIVE = 1
    } ;
    /* @warning Duplicate from Geogram/basic/numeric.h */
    template< typename T >
    inline Sign sign( T x )
    {
        return ( x > 0 ) ? POSITIVE : ( ( x < 0 ) ? NEGATIVE : ZERO ) ;
    }

    /*!
     * See http://www.geometrictools.com/LibMathematics/Distance/Distance.html
     */
    double RINGMESH_API point_segment_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& nearest_p ) ;

    double RINGMESH_API point_triangle_distance(
        const vec3& point,
        const vec3& V0,
        const vec3& V1,
        const vec3& V2,
        vec3& closest_point,
        double& lambda0,
        double& lambda1,
        double& lambda2 ) ;

    double RINGMESH_API point_quad_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p ) ;

    double RINGMESH_API point_tetra_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p ) ;

    double RINGMESH_API point_pyramid_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        vec3& nearest_p ) ;

    double RINGMESH_API point_prism_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        vec3& nearest_p ) ;

    double RINGMESH_API point_hexa_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        const vec3& p6,
        const vec3& p7,
        vec3& nearest_p ) ;

    bool RINGMESH_API line_plane_intersection(
        const vec3& O_line,
        const vec3& D_line,
        const vec3& O_plane,
        const vec3& N_plane,
        vec3& result ) ;

    bool RINGMESH_API segment_plane_intersection(
        const vec3& seg0,
        const vec3& seg1,
        const vec3& O_plane,
        const vec3& N_plane,
        vec3& result ) ;

    bool RINGMESH_API segment_triangle_intersection(
        const vec3& seg0,
        const vec3& seg1,
        const vec3& trgl0,
        const vec3& trgl1,
        const vec3& trgl2,
        vec3& result ) ;

    bool RINGMESH_API circle_plane_intersection(
        const vec3& O_plane,
        const vec3& N_plane,
        const vec3& O_circle,
        const vec3& N_circle,
        double r,
        std::vector< vec3 >& result ) ;

    bool RINGMESH_API disk_segment_intersection(
        const vec3& p0,
        const vec3& p1,
        const vec3& O_circle,
        const vec3& N_circle,
        double r,
        vec3& result ) ;

    bool RINGMESH_API circle_triangle_intersection(
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& O_circle,
        const vec3& N_circle,
        double r,
        std::vector< vec3 >& result ) ;

    bool RINGMESH_API plane_plane_intersection(
        const vec3& O_P0,
        const vec3& N_P0,
        const vec3& O_P1,
        const vec3& N_P1,
        vec3& O_inter,
        vec3& N_inter ) ;

    bool RINGMESH_API point_inside_triangle(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        bool exact_predicates = false ) ;

    bool RINGMESH_API point_inside_quad(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        bool exact_predicates = false ) ;

    bool RINGMESH_API point_inside_tetra(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        bool exact_predicates = false ) ;

    bool RINGMESH_API point_inside_pyramid(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        bool exact_predicates = false ) ;

    bool RINGMESH_API point_inside_prism(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        bool exact_predicates = false ) ;

    bool RINGMESH_API point_inside_hexa(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        const vec3& p6,
        const vec3& p7,
        bool exact_predicates = false ) ;

    bool RINGMESH_API point_segment_projection(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& new_p ) ;

    bool RINGMESH_API tetra_barycentric_coordinates(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        double lambda[4] ) ;

    bool RINGMESH_API triangle_barycentric_coordinates(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        double lambda[3] ) ;

    void RINGMESH_API rotation_matrix_about_arbitrary_axis(
        const vec3& origin,
        const vec3& axis,
        double theta,
        bool degrees,
        GEO::Matrix< 4, double >& rot_mat ) ;

    class RINGMESH_API NNSearch {
    ringmesh_disable_copy( NNSearch ) ;
    public:
        enum MeshLocation {
            VERTICES, EDGES, FACETS, CELLS, CELL_FACETS, NB_LOCATION
        } ;
        NNSearch( const GEO::Mesh& mesh, const MeshLocation& location, bool copy =
            false ) ;
        NNSearch( const std::vector< vec3 >& vertices, bool copy = true ) ;

        ~NNSearch()
        {
            if( delete_points_ ) delete[] nn_points_ ;
        }

        /*!
         * @brief Gets the \p index_map that link all the duplicated points
         * to their first occurancy
         * @return the number of colocated vertices
         * Example:
         *     vertices = [P1, P2, P1, P3, P2, P4]
         *     index_map = [0, 1, 0, 3, 1, 5]
         *     return 2
         */
        index_t get_colocated_index_mapping(
            double epsilon,
            GEO::vector< index_t >& index_map ) const ;
        /*!
         * @brief Gets the \p index_map that link all the points
         * to a no duplicated list of index in the list of \p unique_points.
         * @return the number of colocated vertices
         * Example:
         *     vertices = [P1, P2, P1, P3, P2, P4]
         *     unique_points = [P1, P2, P3, P4]
         *     index_map = [0, 1, 0, 2, 1, 3]
         *     return 2
         */
        index_t get_colocated_index_mapping(
            double epsilon,
            GEO::vector< index_t >& index_map,
            GEO::vector< vec3 >& unique_points ) const ;
        /*!
         * Gets the closest neighbor point
         * @param[in] v the point to test
         * @param[out] distance_sq the square distance to the closest point
         * return returns the index of the closest point
         */
        index_t get_closest_neighbor( const vec3& v, double& distance_sq ) const
        {
            std::vector< index_t > result ;
            get_neighbors( v, 1, result, &distance_sq ) ;
            return result[0] ;
        }

        bool get_neighbors(
            const vec3& v,
            std::vector< index_t >& result,
            double threshold_distance ) const ;

        index_t get_neighbors(
            const vec3& v,
            index_t nb_neighbors,
            std::vector< index_t >& result,
            double* dist = nil ) const ;

        index_t nb_points() const
        {
            return nn_tree_->nb_points() ;
        }

    private:
        void build_nn_search_vertices( const GEO::Mesh& mesh, bool copy ) ;
        void build_nn_search_edges( const GEO::Mesh& mesh ) ;
        void build_nn_search_facets( const GEO::Mesh& mesh ) ;
        void build_nn_search_cells( const GEO::Mesh& mesh ) ;
        void build_nn_search_cell_facets( const GEO::Mesh& mesh ) ;
        void fill_nn_search_points( index_t index_in_nn, const vec3& center ) ;

    private:
        /// KdTree to compute the nearest neighbor search
        GEO::NearestNeighborSearch_var nn_tree_ ;
        /// Array of the points (size of 3xnumber of points)
        double* nn_points_ ;
        /*!
         * @brief Indicates if ann_points_ should ne deleted.
         * @details No need to delete nn_points_ if it is a simple pointer
         * to the mesh vertex array.
         */
        bool delete_points_ ;
    } ;

}

#endif
