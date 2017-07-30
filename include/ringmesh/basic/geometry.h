/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
 * Applications (ASGA) All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include <array>
#include <ringmesh/basic/common.h>

/*!
 * @file Basic geometrical requests
 * @author Arnaud Botella
 */

namespace RINGMesh {

template <index_t DIMENSION>
bool operator==(const vecn<DIMENSION> &u, const vecn<DIMENSION> &v) {
    for (index_t i : range(DIMENSION)) {
        if (u[i] != v[i])
            return false;
    }
    return true;
}

template <index_t DIMENSION>
bool operator!=(const vecn<DIMENSION> &u, const vecn<DIMENSION> &v) {
    return !(u == v);
}

template <index_t DIMENSION>
bool inexact_equal(const vecn<DIMENSION> &v1, const vecn<DIMENSION> &v2,
                   double epsilon) {
    return length(v2 - v1) < epsilon;
}

double RINGMESH_API dot_perp(const vec2 &v0, const vec2 &v1);

/* @warning Duplicate from Geogram/basic/numeric.h */
enum Sign { NEGATIVE = -1, ZERO = 0, POSITIVE = 1 };
/* @warning Duplicate from Geogram/basic/numeric.h */
template <typename T> Sign sign(T x) {
    return (x > 0) ? POSITIVE : ((x < 0) ? NEGATIVE : ZERO);
}

double RINGMESH_API triangle_signed_area(const vec3 &p0, const vec3 &p1,
                                         const vec3 &p2,
                                         const vec3 &triangle_normal);

namespace Distance {
/*!
 * See http://www.geometrictools.com/LibMathematics/Distance/Distance.html
 */
std::tuple<double, vec3> RINGMESH_API point_to_segment(const vec3 &p,
                                                       const vec3 &p0,
                                                       const vec3 &p1);

std::tuple<double, vec2> RINGMESH_API point_to_segment(const vec2 &p,
                                                       const vec2 &p0,
                                                       const vec2 &p1);

/*!
 * Computes the smallest distance between a point and a triangle
 * @param[in] point the point to test
 * @param[in] V0 the first vertex of the triangle
 * @param[in] V1 the second vertex of the triangle
 * @param[in] V2 the third vertex of the triangle
 * @return a tuple containing the following elements (in this order):
 * - the smallest distance
 * - the closest point on the triangle
 */
std::tuple<double, vec3> RINGMESH_API point_to_triangle(const vec3 &point,
                                                        const vec3 &V0,
                                                        const vec3 &V1,
                                                        const vec3 &V2);

/*!
 * Computes the smallest distance between a point and a triangle
 * @param[in] point the point to test
 * @param[in] V0 the first vertex of the triangle
 * @param[in] V1 the second vertex of the triangle
 * @param[in] V2 the third vertex of the triangle
 * @return a tuple containing:
 * - the smallest distance.
 * - the closest_point the closest point on the triangle.
 */
std::tuple<double, vec2> RINGMESH_API point_to_triangle(const vec2 &point,
                                                        const vec2 &V0,
                                                        const vec2 &V1,
                                                        const vec2 &V2);

/*!
 * Computes the distance between a point and a tetrahedron
 * @param[in] p the point
 * @param[in] p0 the first vertex of the tetrahedron
 * @param[in] p1 the second vertex of the tetrahedron
 * @param[in] p2 the third vertex of the tetrahedron
 * @param[in] p3 the fourth vertex of the tetrahedron
 * @return a tuple containing:
 * - the distance between the point and the tetrahedron facets.
 * - the nearest point on the tetrahedron.
 */
std::tuple<double, vec3>
    RINGMESH_API point_to_tetra(const vec3 &p, const vec3 &p0, const vec3 &p1,
                                const vec3 &p2, const vec3 &p3);
} // namespace Distance

namespace Intersection {
/*!
 * Computes the intersection between a plane and a line
 * @param[in] O_line a point on the line
 * @param[in] D_line the direction of the plane
 * @param[in] O_plane a point on the plane
 * @param[in] N_plane the normal of the plane
 * @return returns a tuple containing a boolean (true if there is an
 * intersection) and the intersected point if any.
 */
std::tuple<bool, vec3> RINGMESH_API line_plane(const vec3 &O_line,
                                               const vec3 &D_line,
                                               const vec3 &O_plane,
                                               const vec3 &N_plane);

/*!
 * Computes the intersection between a plane and a segment
 * @param[in] p0 the first vertex of the segment
 * @param[in] p1 the second vertex of the segment
 * @param[in] O_plane a point on the plane
 * @param[in] N_plane the normal of the plane
 * @return returns a tuple containing a boolean (true if there is an
 * intersection) and the intersected point if any.
 */
std::tuple<bool, vec3> RINGMESH_API segment_plane(const vec3 &seg0,
                                                  const vec3 &seg1,
                                                  const vec3 &O_plane,
                                                  const vec3 &N_plane);

/*!
 * Computes the intersection of a segment and a triangle
 * @param[in] seg0 the first vertex of the segment
 * @param[in] seg1 the second vertex of the segment
 * @param[in] trgl0 the first vertex of the triangle
 * @param[in] trgl1 the second vertex of the triangle
 * @param[in] trgl2 the third vertex of the triangle
 * @return a tuple containing a boolean (true is there is an intersection)
 * and the intersected point if any.
 */
std::tuple<bool, vec3> RINGMESH_API segment_triangle(const vec3 &seg0,
                                                     const vec3 &seg1,
                                                     const vec3 &trgl0,
                                                     const vec3 &trgl1,
                                                     const vec3 &trgl2);

/*!
 * Computes the intersection(s) between a circle and a plane
 * @param[in] O_plane a point on the plane
 * @param[in] N_plane the normal of the plane
 * @param[in] O_circle the center of the circle
 * @param[in] N_circle the normal of the plane supporting the circle
 * @param[in] r the radius of the circle
 * @return returns a tuple containing a boolean (true if there is at least one
 * intersection) and the intersected points if any.
 */
std::tuple<bool, std::vector<vec3>>
    RINGMESH_API circle_plane(const vec3 &O_plane, const vec3 &N_plane,
                              const vec3 &O_circle, const vec3 &N_circle,
                              double r);

/*!
 * Computes the intersection between a disk and a segment
 * @param[in] p0 the first vertex of the segment
 * @param[in] p1 the second vertex of the segment
 * @param[in] O_disk the center of the disk
 * @param[in] N_disk the normal of the plane supporting the disk
 * @param[in] r the radius of the disk
 * @return returns a tuple containing a boolean (true if there is an
 * intersection) and the intersected point if any.
 */
std::tuple<bool, vec3> RINGMESH_API disk_segment(const vec3 &p0, const vec3 &p1,
                                                 const vec3 &O_disk,
                                                 const vec3 &N_disk, double r);

/*!
 * Computes the intersection(s) between a circle and a triangle
 * @param[in] p0 the first vertex of the triangle
 * @param[in] p1 the second vertex of the triangle
 * @param[in] p2 the third vertex of the triangle
 * @param[in] O_circle the center of the circle
 * @param[in] N_circle the normal of the plane supporting the circle
 * @param[in] r the radius of the circle
 * @return returns a tuple containing a boolean (true if there is at least one
 * intersection) and the intersected points if any.
 */
std::tuple<bool, std::vector<vec3>>
    RINGMESH_API circle_triangle(const vec3 &p0, const vec3 &p1, const vec3 &p2,
                                 const vec3 &O_circle, const vec3 &N_circle,
                                 double r);

/*!
 * Computes the intersection between two planes
 * @param[in] O_P0 a point on the first plane
 * @param[in] N_P0 the normal of the first plane
 * @param[in] O_P1 a point on the second plane
 * @param[in] N_P1 the normal of the second plane
 * @return a tuple containing:
 * - a boolean: true is there is an intersection between the planes.
 * - a point on the intersected line if any.
 * - the direction of the intersected line if any.
 */
std::tuple<bool, vec3, vec3> RINGMESH_API plane_plane(const vec3 &O_P0,
                                                      const vec3 &N_P0,
                                                      const vec3 &O_P1,
                                                      const vec3 &N_P1);

/*!
 * Computes the intersection between two lines
 * @param[in] O_line0 a point on the first line
 * @param[in] D_line0 the direction of the first line
 * @param[in] O_line1 a point on the second line
 * @param[in] D_line1 the direction of the second line
 * @return a tuple containing:
 * - a boolean: true if there is an intersection.
 * - the intersected point if any.
 */
std::tuple<bool, vec2> RINGMESH_API line_line(const vec2 &O_line0,
                                              const vec2 &D_line0,
                                              const vec2 &O_line1,
                                              const vec2 &D_line1);

/*!
 * Computes the intersection between two segments
 * @param[in] p0_seg0 the first vertex of the segment
 * @param[in] p1_seg0 the second vertex of the segment
 * @param[in] p0_seg1 the first vertex of the segment
 * @param[in] p1_seg1 the second vertex of the segment
 * @return a tuple containing:
 * - a boolean: true if there is an intersection.
 * - the intersected point if any.
 */
std::tuple<bool, vec2> RINGMESH_API segment_segment(const vec2 &p0_seg0,
                                                    const vec2 &p1_seg0,
                                                    const vec2 &p0_seg1,
                                                    const vec2 &p1_seg1);

/*!
 * Computes the intersection between a segment and a line
 * @param[in] p0_seg the first vertex of the segment
 * @param[in] p1_seg the second vertex of the segment
 * @param[in] O_line a point on the line
 * @param[in] D_line the direction of the line
 * @return a tuple containing:
 * - a boolean: true if there is an intersection.
 * - the intersected point if any.
 */
std::tuple<bool, vec2> RINGMESH_API segment_line(const vec2 &p0_seg,
                                                 const vec2 &p1_seg,
                                                 const vec2 &O_line,
                                                 const vec2 &D_line);
} // namespace Intersection

/*!
 * @brief Tests if a point is inside a triangle
 * @details if it is inside a prism based on the triangle and its normal
 * @param[in] p the point to test
 * @param[in] p0 the first vertex of the triangle
 * @param[in] p1 the second vertex of the triangle
 * @param[in] p2 the third vertex of the triangle
 * @param[in] exact_predicates if true, the algorithm uses exact predicates
 * @return returns true if the point is inside
 */
bool RINGMESH_API point_inside_triangle(const vec3 &p, const vec3 &p0,
                                        const vec3 &p1, const vec3 &p2,
                                        bool exact_predicates = false);

/*!
 * @brief Tests if a point is inside a triangle
 * @details if it is inside a prism based on the triangle and its normal
 * @param[in] p the point to test
 * @param[in] p0 the first vertex of the triangle
 * @param[in] p1 the second vertex of the triangle
 * @param[in] p2 the third vertex of the triangle
 * @param[in] exact_predicates if true, the algorithm uses exact predicates
 * @return returns true if the point is inside
 */
bool RINGMESH_API point_inside_triangle(const vec2 &p, const vec2 &p0,
                                        const vec2 &p1, const vec2 &p2,
                                        bool exact_predicates = false);

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
bool RINGMESH_API point_inside_tetra(const vec3 &p, const vec3 &p0,
                                     const vec3 &p1, const vec3 &p2,
                                     const vec3 &p3,
                                     bool exact_predicates = false);

/*!
 * Computes the orthogonal projection of a point on a segment
 * @param[in] p the point to project
 * @param[in] p0 the first vertex of the segment
 * @param[in] p1 the second vertex of the segment
 * @return returns a tuple containing a boolean (true if the projection is
 * possible) and the projected point if any.
 */
template <index_t DIMENSION>
std::tuple<bool, vecn<DIMENSION>>
    RINGMESH_API point_segment_projection(const vecn<DIMENSION> &p,
                                          const vecn<DIMENSION> &p0,
                                          const vecn<DIMENSION> &p1);

/*!
 * Computes the orthogonal projection of a point on a plane
 * @param[in] p the point to project
 * @param[in] N_plane the normal of the plane
 * @param[in] O_plane a point of the plane
 * @return the projected point
 */
vec3 RINGMESH_API point_plane_projection(const vec3 &p, const vec3 &N_plane,
                                         const vec3 &O_plane);

/*!
 * Computes barycentric coordinates of \p p
 * @param[in] p the query point
 * @param[in] p0 the first tetra vertex
 * @param[in] p1 the second tetra vertex
 * @param[in] p2 the third tetra vertex
 * @param[in] p3 the fourth tetra vertex
 * @return a tuple containing:
 * - a boolean (false if the computation failed because of too small tetrahedron
 * volume).
 * - the parametric coordinates corresponding to points
 */
std::tuple<bool, std::array<double, 4>>
    RINGMESH_API tetra_barycentric_coordinates(const vec3 &p, const vec3 &p0,
                                               const vec3 &p1, const vec3 &p2,
                                               const vec3 &p3);

/*!
 * Computes barycentric coordinates of \p p
 * @param[in] p the query point
 * @param[in] p0 the first triangle vertex
 * @param[in] p1 the second triangle vertex
 * @param[in] p2 the third triangle vertex
 * @return a tuple containing:
 * - a boolean (false if the computation failed because of too small triangle
 * area).
 * - the parametric coordinates corresponding to points.
 */
std::tuple<bool, std::array<double, 3>>
    RINGMESH_API triangle_barycentric_coordinates(const vec3 &p, const vec3 &p0,
                                                  const vec3 &p1,
                                                  const vec3 &p2);

/*!
 * Computes barycentric coordinates of \p p
 * @param[in] p the query point
 * @param[in] p0 the first triangle vertex
 * @param[in] p1 the second triangle vertex
 * @param[in] p2 the third triangle vertex
 * @return a tuple containing:
 * - a boolean (false if the computation failed because of too small triangle
 * area).
 * - the parametric coordinates corresponding to points.
 */
std::tuple<bool, std::array<double, 3>>
    RINGMESH_API triangle_barycentric_coordinates(const vec2 &p, const vec2 &p0,
                                                  const vec2 &p1,
                                                  const vec2 &p2);

/*!
 * @brief Builds a rotational matrix about an arbitrary axis.
 *
 * Mathematical development: http://paulbourke.net/geometry/rotate/.
 * @param[in] origin point in which passes the rotation axis.
 * @param[in] axis vector which defines the rotation axis.
 * @param[in] theta rotation angle (in radians or degrees).
 * @param[in] degrees true is \p theta is in degrees, false
 * if in radians.
 * @return the matrix which defines the rotation
 * of a point around the axis defined by point \p origin
 * and vector \p axis by an angle \p theta.
 * New coordinates of a point (x,y,z) are:
 * (x',y',z') = rot_mat*(x,y,z)
 */
GEO::Matrix<4, double> RINGMESH_API rotation_matrix_about_arbitrary_axis(
    const vec3 &origin, const vec3 &axis, double theta, bool degrees);

} // namespace RINGMesh
