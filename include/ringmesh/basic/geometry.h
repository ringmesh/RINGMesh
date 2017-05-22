/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#pragma once

#include <ringmesh/basic/common.h>

#include <geogram/points/kd_tree.h>

/*!
 * @file What
 * @author Arnaud Botella
 */

namespace GEO {
    class Mesh;
}

namespace RINGMesh {

    bool RINGMESH_API operator==( const vec3& u, const vec3& v );
    bool RINGMESH_API operator!=( const vec3& u, const vec3& v );

    /* @warning Duplicate from Geogram/basic/numeric.h */
    enum Sign {
        NEGATIVE = -1, ZERO = 0, POSITIVE = 1
    };
    /* @warning Duplicate from Geogram/basic/numeric.h */
    template< typename T >
    inline Sign sign( T x )
    {
        return ( x > 0 ) ? POSITIVE : ( ( x < 0 ) ? NEGATIVE : ZERO );
    }

    /*!
     * See http://www.geometrictools.com/LibMathematics/Distance/Distance.html
     */
    double RINGMESH_API point_segment_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& nearest_p );

    /*!
     * Computes the smallest distance between a point and a triangle
     * @param[in] point the point to test
     * @param[in] V0 the first vertex of the triangle
     * @param[in] V1 the second vertex of the triangle
     * @param[in] V2 the third vertex of the triangle
     * @param[out] closest_point the closest point on the triangle
     * @param[out] lambda0 barycentric coordinate from \p V0
     * @param[out] lambda1 barycentric coordinate from \p V1
     * @param[out] lambda2 barycentric coordinate from \p V2
     * @return the smallest distance
     */
    double RINGMESH_API point_triangle_distance(
        const vec3& point,
        const vec3& V0,
        const vec3& V1,
        const vec3& V2,
        vec3& closest_point,
        double& lambda0,
        double& lambda1,
        double& lambda2 );

    /*!
     * Computes the smallest distance between a point and a quad
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the quad
     * @param[in] p1 the second vertex of the quad
     * @param[in] p2 the third vertex of the quad
     * @param[in] p3 the fourth vertex of the quad
     * @param[out] nearest_p the closest point on the quad
     * @return the smallest distance
     */
    double RINGMESH_API point_quad_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p );

    /*!
     * Computes the distance between a point and a tetrahedron
     * @param[in] p the point
     * @param[in] p0 the first vertex of the tetrahedron
     * @param[in] p1 the second vertex of the tetrahedron
     * @param[in] p2 the third vertex of the tetrahedron
     * @param[in] p3 the fourth vertex of the tetrahedron
     * @param[out] nearest_p the nearest point on the tetrahedron
     * @return the distance between the point and the tetrahedron facets
     */
    double RINGMESH_API point_tetra_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        vec3& nearest_p );

    /*!
     * Computes the distance between a point and a pyramid
     * @param[in] p the point
     * @param[in] p0 the first vertex of the pyramid
     * @param[in] p1 the second vertex of the pyramid
     * @param[in] p2 the third vertex of the pyramid
     * @param[in] p3 the fourth vertex of the pyramid
     * @param[in] p4 the fifth vertex of the pyramid
     * @param[out] nearest_p the nearest point on the pyramid
     * @return the distance between the point and the pyramid facets
     */
    double RINGMESH_API point_pyramid_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        vec3& nearest_p );

    /*!
     * Computes the distance between a point and a prism
     * @param[in] p the point
     * @param[in] p0 the first vertex of the prism
     * @param[in] p1 the second vertex of the prism
     * @param[in] p2 the third vertex of the prism
     * @param[in] p3 the fourth vertex of the prism
     * @param[in] p4 the fifth vertex of the prism
     * @param[in] p5 the sixth vertex of the prism
     * @param[out] nearest_p the nearest point on the prism
     * @return the distance between the point and the prism facets
     */
    double RINGMESH_API point_prism_distance(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        const vec3& p4,
        const vec3& p5,
        vec3& nearest_p );

    /*!
     * Computes the distance between a point and a hexahedron
     * @param[in] p the point
     * @param[in] p0 the first vertex of the hexahedron
     * @param[in] p1 the second vertex of the hexahedron
     * @param[in] p2 the third vertex of the hexahedron
     * @param[in] p3 the fourth vertex of the hexahedron
     * @param[in] p4 the fifth vertex of the hexahedron
     * @param[in] p5 the sixth vertex of the hexahedron
     * @param[in] p6 the seventh vertex of the hexahedron
     * @param[in] p7 the heith vertex of the hexahedron
     * @param[out] nearest_p the nearest point on the hexahedron
     * @return the distance between the point and the hexahedron facets
     */
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
        vec3& nearest_p );

    /*!
     * Computes the intersection between a plane and a line
     * @param[in] O_line a point on the line
     * @param[in] D_line the direction of the plane
     * @param[in] O_plane a point on the plane
     * @param[in] N_plane the normal of the plane
     * @param[out] result the intersected point
     * @return returns true if there is an intersection
     */
    bool RINGMESH_API line_plane_intersection(
        const vec3& O_line,
        const vec3& D_line,
        const vec3& O_plane,
        const vec3& N_plane,
        vec3& result );

    /*!
     * Computes the intersection between a plane and a segment
     * @param[in] p0 the first vertex of the segment
     * @param[in] p1 the second vertex of the segment
     * @param[in] O_plane a point on the plane
     * @param[in] N_plane the normal of the plane
     * @param[out] result the intersected point
     * @return returns true if there is an intersection
     */
    bool RINGMESH_API segment_plane_intersection(
        const vec3& seg0,
        const vec3& seg1,
        const vec3& O_plane,
        const vec3& N_plane,
        vec3& result );

    /*!
     * Computes the intersection of a segment and a triangle
     * @param[in] seg0 the first vertex of the segment
     * @param[in] seg1 the second vertex of the segment
     * @param[in] trgl0 the first vertex of the triangle
     * @param[in] trgl1 the second vertex of the triangle
     * @param[in] trgl2 the third vertex of the triangle
     * @param[out] result the intersected point
     * @return true is there is an intersection
     */
    bool RINGMESH_API segment_triangle_intersection(
        const vec3& seg0,
        const vec3& seg1,
        const vec3& trgl0,
        const vec3& trgl1,
        const vec3& trgl2,
        vec3& result );

    /*!
     * Computes the intersection(s) between a circle and a plane
     * @param[in] O_plane a point on the plane
     * @param[in] N_plane the normal of the plane
     * @param[in] O_circle the center of the circle
     * @param[in] N_circle the normal of the plane supporting the circle
     * @param[in] r the radius of the circle
     * @param[out] result the intersected points
     * @return returns true if there is at least one intersection
     */
    bool RINGMESH_API circle_plane_intersection(
        const vec3& O_plane,
        const vec3& N_plane,
        const vec3& O_circle,
        const vec3& N_circle,
        double r,
        std::vector< vec3 >& result );

    /*!
     * Computes the intersection between a disk and a segment
     * @param[in] p0 the first vertex of the segment
     * @param[in] p1 the second vertex of the segment
     * @param[in] O_disk the center of the disk
     * @param[in] N_disk the normal of the plane supporting the disk
     * @param[in] r the radius of the disk
     * @param[out] result the intersected point
     * @return returns true if there is an intersection
     */
    bool RINGMESH_API disk_segment_intersection(
        const vec3& p0,
        const vec3& p1,
        const vec3& O_disk,
        const vec3& N_disk,
        double r,
        vec3& result );

    /*!
     * Computes the intersection(s) between a circle and a triangle
     * @param[in] p0 the first vertex of the triangle
     * @param[in] p1 the second vertex of the triangle
     * @param[in] p2 the third vertex of the triangle
     * @param[in] O_circle the center of the circle
     * @param[in] N_circle the normal of the plane supporting the circle
     * @param[in] r the radius of the circle
     * @param[out] result the intersected points
     * @return returns true if there is at least one intersection
     */
    bool RINGMESH_API circle_triangle_intersection(
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& O_circle,
        const vec3& N_circle,
        double r,
        std::vector< vec3 >& result );

    /*!
     * Computes the intersection between two planes
     * @param[in] O_P0 a point on the first plane
     * @param[in] N_P0 the normal of the frst plane
     * @param[in] O_P1 a point on the second plane
     * @param[in] N_P1 the normal of the second plane
     * @param[out] O_inter a point on the intersected line
     * @param[out] D_inter the direction of the intersected line
     * @return true is there is an intersection between the planes
     */
    bool RINGMESH_API plane_plane_intersection(
        const vec3& O_P0,
        const vec3& N_P0,
        const vec3& O_P1,
        const vec3& N_P1,
        vec3& O_inter,
        vec3& N_inter );

    /*!
     * Tests if a point is inside a triangle, more precisely if it is inside
     * a prism based on the triangle and its normal
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the triangle
     * @param[in] p1 the second vertex of the triangle
     * @param[in] p2 the third vertex of the triangle
     * @param[in] exact_predicates if true, the algorithm uses exact predicates
     * @return returns true if the point is inside
     */
    bool RINGMESH_API point_inside_triangle(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        bool exact_predicates = false );

    /*!
     * Tests if a point is inside a tetrahedron
     * @param[in] p the point to test
     * @param[in] p0 the first vertex of the tetrahedron
     * @param[in] p1 the second vertex of the tetrahedron
     * @param[in] p2 the third vertex of the tetrahedron
     * @param[in] p3 the fourth vertex of the tetrahedron
     * @param[in] exact_predicates if true, the algorithm uses exact predicates
     * @return returns true if the point is inside the tetrahedron
     */
    bool RINGMESH_API point_inside_tetra(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        bool exact_predicates = false );

    /*!
     * Computes the orthogonal projection of a point on a segment
     * @param[in] p the point to project
     * @param[in] p0 the first vertex of the segment
     * @param[in] p1 the second vertex of the segment
     * @param[out] new_p the projected point
     * @return returns true if the projection is possible
     */
    bool RINGMESH_API point_segment_projection(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        vec3& new_p );

    /*!
     * Computes the orthogonal projection of a point on a plane
     * @param[in] p the point to project
     * @param[in] N_plane the normal of the plane
     * @param[in] O_plane a point of the plane
     * @param[out] projected_p the projected point
     */
    void RINGMESH_API point_plane_projection(
        const vec3& p,
        const vec3& N_plane,
        const vec3& O_plane,
        vec3& projected_p );

    /*!
     * Computes barycentric coordinates of \p p
     * @param[in] p the query point
     * @param[in] p0 the first tetra vertex
     * @param[in] p1 the second tetra vertex
     * @param[in] p2 the third tetra vertex
     * @param[in] p3 the fourth tetra vertex
     * @param[out] lambda the parametric coordinates corresponding to points
     * @return false if the computation failed because of too small tetrahedron volume
     */
    bool RINGMESH_API tetra_barycentric_coordinates(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        const vec3& p3,
        double lambda[4] );

    /*!
     * Computes barycentric coordinates of \p p
     * @param[in] p the query point
     * @param[in] p0 the first triangle vertex
     * @param[in] p1 the second triangle vertex
     * @param[in] p2 the third triangle vertex
     * @param[out] lambda the parametric coordinates corresponding to points
     * @return false if the computation failed because of too small triangle area
     */
    bool RINGMESH_API triangle_barycentric_coordinates(
        const vec3& p,
        const vec3& p0,
        const vec3& p1,
        const vec3& p2,
        double lambda[3] );

    /*!
     * @brief Builds a rotational matrix about an arbitrary axis.
     *
     * Mathematical development: http://paulbourke.net/geometry/rotate/.
     *
     * @param[in] origin point in which passes the rotation axis.
     *
     * @param[in] axis vector which defines the rotation axis.
     *
     * @param[in] theta rotation angle (in radians or degrees).
     *
     * @param[in] degrees true is \p theta is in degrees, false
     * if in radians.
     *
     * @param[out] rot_mat the matrix which defines the rotation
     * of a point around the axis defined by point \p origin
     * and vector \p axis by an angle \p theta.
     * New coordinates of a point (x,y,z) are:
     * (x',y',z') = rot_mat*(x,y,z)
     */
    void RINGMESH_API rotation_matrix_about_arbitrary_axis(
        const vec3& origin,
        const vec3& axis,
        double theta,
        bool degrees,
        GEO::Matrix< 4, double >& rot_mat );

    class RINGMESH_API NNSearch {
    ringmesh_disable_copy( NNSearch );
    public:
        enum MeshLocation {
            VERTICES, EDGES, FACETS, CELLS, CELL_FACETS, NB_LOCATION
        };
        NNSearch( const GEO::Mesh& mesh, const MeshLocation& location, bool copy =
            false );
        NNSearch( const std::vector< vec3 >& vertices, bool copy = true );

        ~NNSearch()
        {
            if( delete_points_ ) delete[] nn_points_;
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
            std::vector< index_t >& index_map ) const;
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
            std::vector< index_t >& index_map,
            std::vector< vec3 >& unique_points ) const;
        /*!
         * Gets the closest neighbor point
         * @param[in] v the point to test
         * return returns the index of the closest point
         */
        index_t get_closest_neighbor( const vec3& v ) const
        {
            index_t nb_neighbors = 1;
            index_t result = get_neighbors( v, nb_neighbors ).front();
            return result;
        }

        /*!
         * Compute the neighbors of a given point, point closer than \param threshold_distance
         * @param[in] v the point to test
         * @param[in] threshold_distance distance defining the neighborhood
         * @return the point indices
         */
        std::vector< index_t > get_neighbors(
            const vec3& v,
            double threshold_distance ) const;

        /*!
         * Gets the neighboring points of a given one sorted by increasing distance
         * @param[in] v the point to test
         * @param[in] nb_neighbors the number of neighbors to return
         * @return the point indices (can be less than \p nb_neighbors
         * if there is not enough points)
         */
        std::vector< index_t > get_neighbors(
            const vec3& v,
            index_t nb_neighbors ) const;

        vec3 point( index_t v ) const
        {
            return vec3( nn_points_[3 * v], nn_points_[3 * v + 1],
                nn_points_[3 * v + 2] );
        }

        index_t nb_points() const
        {
            return nn_tree_->nb_points();
        }

    private:
        void build_nn_search_vertices( const GEO::Mesh& mesh, bool copy );
        void build_nn_search_edges( const GEO::Mesh& mesh );
        void build_nn_search_polygons( const GEO::Mesh& mesh );
        void build_nn_search_cells( const GEO::Mesh& mesh );
        void build_nn_search_cell_facets( const GEO::Mesh& mesh );
        void fill_nn_search_points( index_t index_in_nn, const vec3& center );

    private:
        /// KdTree to compute the nearest neighbor search
        GEO::NearestNeighborSearch_var nn_tree_;
        /// Array of the points (size of 3xnumber of points)
        double* nn_points_;
        /*!
         * @brief Indicates if ann_points_ should ne deleted.
         * @details No need to delete nn_points_ if it is a simple pointer
         * to the mesh vertex array.
         */
        bool delete_points_;
    };

}
