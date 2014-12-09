/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef __GEOGRAM_MESH_MESH__
#define __GEOGRAM_MESH_MESH__

#include <geogram/basic/common.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/assert.h>
#include <geogram/basic/logger.h>
#include <algorithm>

/**
 * \file geogram/mesh/mesh.h
 * \brief The class to represent a surfacic or volumetric mesh
 */

namespace GEO {

    /**
     * \brief Indicates the attributes stored in a mesh and attached
     *  to the mesh elements (vertices, facets or volumes).
     * \details The set of attributes attached to a mesh is represented
     *  by a bitwise or combination of the constants.
     */
    enum MeshAttributes {
        MESH_NO_ATTRIBUTES = 0,
        MESH_VERTEX_REGION = 1,
        MESH_VERTEX_NORMAL = 2,
        MESH_VERTEX_COLOR = 4,
        MESH_FACET_REGION = 8,
        MESH_CELL_REGION = 16
    };

    /**
     * \brief Indicates the mesh elements (vertices, facets or cells)
     *  present in a mesh.
     * \details The set of elements present in a mesh is represented
     *  by a bitwise-or combination of the constants.
     */
    enum MeshElements {
        MESH_VERTICES = 1,
        MESH_FACETS = 2,
        MESH_CELLS = 4
    };


    /**
     * \brief Indicates the type of a cell in a mesh.
     * \details MESH_CONNECTOR is a special element, that
     *  represents non-conformal connections between an hex
     *  facet and two tets.
     */
    enum MeshCellType {
        MESH_TET = 0,
        MESH_HEX = 1,
        MESH_PRISM = 2,
        MESH_PYRAMID = 3,
        MESH_CONNECTOR = 4,
        MESH_NB_CELL_TYPES = 5
    };
}

/**
 * \brief Generic classes (templates) are grouped in this namespace.
 */
namespace GEOGen {

    using GEO::index_t;
    using GEO::signed_index_t;
    using GEO::coord_index_t;
    using GEO::MeshAttributes;

    template <class MESH>
    class MeshBuilder;

    template <class MESH>
    class MeshMutator;

    /**
     * \brief Base class for Mesh.
     * \details Supports the combinatorics operations and 
     *  contains some static tables. The derived class Mesh
     *  is a template, parameterized by the representation of
     *  coordinates.
     */
    class MeshBase {
    public:

        //   The following two lines are here for debugging
        // purposes, will be removed soon...
        float* debug_coords_float_;   // TODO remove it
        double* debug_coords_double_; // TODO remove it
        
        /**
         * \brief Gives the local index of a vertex in a
         *  tetrahedron from its facet and vertex local indices.
         * \param[in] lf local facet index (0,1,2 or 3)
         * \param[in] lv local vertex index in \p lf (0,1 or 2)
         * \return the local vertex index (0,1,2 or 3) of the
         * \p lv%th vertex in facet \p lf
         */
        static index_t local_tet_facet_vertex_index(index_t lf, index_t lv) {
            geo_debug_assert(lf < 4);
            geo_debug_assert(lv < 3);
            return tet_descriptor_.facet_vertex[lf][lv];
        }

        /**
         * \brief Constructs a new MeshBase.
         * \param[in] dim the dimension of the vertices
         *   (e.g., 3 for 3d)
         */
        explicit MeshBase(coord_index_t dim = 0);

        /**
         * \brief Checks whether the mesh is valid
         * \details A mesh is considered valid if the dimension has been
         * defined.
         */
        bool is_valid() const {
            return dimension_ != 0;
        }

        /**
         * \brief Resets the mesh structure
         * \param[in] keep_memory if set, memory is not
         * deallocated (may be useful for dynamic
         * meshes).
         */
        void clear(bool keep_memory = false);

        /**
         * \brief Pre-allocates the facets and corners table.
         * \details Pre-reserve memory for \p nb_facet facets
         * and \p nb_corners corners.
         * Calling this function is not mandatory (but
         * it improves performances.
         * \param[in] nb_facets number of facets to pre-allocate
         * \param[in] nb_corners number of corners to pre-allocate
         */
        void reserve_facets_and_corners(
            index_t nb_facets, index_t nb_corners
        );

        /**
         * \brief Pre-allocates the facets and corners table.
         * \details Pre-reserve memory for \p nb_facets facets and the corners
         * table for \p nb_facets considering that each facet has an average
         * number of \p avg_facet_degree corners.
         * Calling this function is not mandatory (but
         * it improves performances.
         * \param[in] nb_facets number of facets to pre-allocate
         * \param[in] avg_facet_degree an estimate of the average
         *  number of corners per facet
         */
        void reserve_facets(
            index_t nb_facets,
            index_t avg_facet_degree = 3
        ) {
            index_t nb_corners = avg_facet_degree * nb_facets;
            reserve_facets_and_corners(nb_facets, nb_corners);
        }

        /**
         * \brief Creates a chunk of tetrahedra
         * \param[in] nb_tets number of tetrahedra to create
         */
        void create_tets(index_t nb_tets);

        /**
         * \brief Begins a new facet.
         */
        void begin_facet() {
            geo_assert(!in_facet_);
            in_facet_ = true;
        }

        /**
         * \brief Adds a corner to the current facet.
         * \details The current facet needs to be created
         *  using begin_facet(). The adjacent facet of this
         *  corner is left uninitialized (-1).
         * \param[in] vindex index of the vertex to be added to
         *  the current facet.
         */
        void add_corner(index_t vindex) {
            corner_vertices_.push_back(vindex);
            corner_adjacent_facets_.push_back(-1);
        }

        /**
         * \brief Adds a corner to the current facet and sets its
         *  adjacent facet.
         * \details The current facet needs to be created
         *  using begin_facet().
         * \param[in] vindex index of the vertex to be added to the
         *  current facet.
         * \param[in] findex index of the facet opposite to this corner.
         */
        void add_corner(index_t vindex, signed_index_t findex) {
            corner_vertices_.push_back(vindex);
            corner_adjacent_facets_.push_back(findex);
        }

        /**
         * \brief Terminates the current facet.
         * \param[in] facet_region optional facet region. Ignored
         *  if the mesh does not have the MESH_FACET_REGION attribute.
         */
        void end_facet(signed_index_t facet_region = 0);

        /**
         * \brief Tests whether we are between a
         * begin_facet() / end_facet() pair.
         */
        bool in_facet() const {
            return in_facet_;
        }

        /**
         * \brief Dissociates a facet by removing the links with the
         *  adjacent facets.
         * \param[in] f index of the facet to dissociate
         */
        void dissociate_facet(index_t f);

        // TODO: fix them by splitting the facets rather than
        // dissociation.

        /**
         * \brief Fixes the degree 2 vertices in a surfacic mesh.
         * \details Degree 2 vertices that are not on the border
         *  are fixed by dissociating the facets incident to them.
         *  Degree 2 vertices are
         *  problematic for some combinatorial algorithms that
         *  require at least three incident facets per internal
         *  vertex (vertices on the boundary can be of degree 2).
         */
        void remove_degree2_vertices();

        /**
         * \brief Gets the dimension of the mesh.
         * \return the dimension of this mesh
         */
        coord_index_t dimension() const {
            return dimension_;
        }

        /**
         * \brief Gets the attributes of the mesh.
         * \return a bitwise-or combination of MeshAttributes flags.
         */
        MeshAttributes attributes() const {
            return attributes_;
        }

        /**
         * \brief Tests whether an attribute is stored in this mesh.
         * \param[in] x the attribute to test
         * \return true if this mesh has attribute \p x, false otherwise
         */
        bool has_attribute(MeshAttributes x) const {
            return (attributes_ & x) != 0;
        }

        /**
         * \brief Tests whether this mesh is known to be triangulated.
         * \details For a triangulated mesh, facets pointers are not
         *  stored, since they can be computed implicitly.
         * \retval true if the mesh is triangulated
         * \retval false otherwise
         */
        bool is_triangulated() const {
            return triangulated_;
        }

        /**
         * \brief Tests whether this mesh is known to be tetrahedralized.
         * \details For a tetrahedralized mesh, cell pointers are not
         *  stored, since they can be computed implicitly.
         * \retval true if the mesh is tetrahedralized
         * \retval false otherwise
         */
        bool is_tetrahedralized() const {
            return tetrahedralized_;
        }
        
        /**
         * \name Corners and vertices access.
         * @{
         */

        /**
         * \brief Gets the number of corners.
         * \return the number of corners of this mesh.
         */
        index_t nb_corners() const {
            return (index_t) corner_vertices_.size();
        }

        /**
         * \brief Gets a pointer to the vertex index associated
         *  with a given corner.
         * \param[in] c the index of the corner
         * \return a pointer to the vertex index associated with corner \p c
         */
        index_t* corner_vertex_index_ptr(index_t c) {
            geo_debug_assert(c < nb_corners());
            return &(corner_vertices_[c]);
        }

        /**
         * \brief Gets the number of vertices
         */
        index_t nb_vertices() const {
            return nb_vertices_;
        }

        /**
         * \brief Sets the vertex associated with a corner.
         * \param[in] c the index of the corner
         * \param[in] v the index of the vertex
         */
        void set_corner_vertex_index(index_t c, index_t v) {
            geo_debug_assert(c < nb_corners());
            geo_debug_assert(v < nb_vertices());
            corner_vertices_[c] = v;
        }

        /**
         * \brief Gets the vertex associated with a corner.
         * \param[in] c the index of the corner
         * \return the index of the vertex that corner \p c points at
         */
        index_t corner_vertex_index(index_t c) const {
            geo_debug_assert(c < nb_corners());
            return corner_vertices_[c];
        }

        /**
         * \brief Gets the facet adjacent to a corner.
         * \param[in] c index of the corner
         * \return index of the facet adjacent to corner \p c, or -1 it
         *  corner \p c is on the border
         */
        signed_index_t corner_adjacent_facet(index_t c) const {
            geo_debug_assert(c < nb_corners());
            return corner_adjacent_facets_[c];
        }

        /**
         * \brief Sets the facet adjacent to a corner
         * \param[in] c index of the corner
         * \param[in] f index of the facet, or -1 if corner \p c is on
         *  the border
         */
        void set_corner_adjacent_facet(index_t c, signed_index_t f) {
            geo_debug_assert(c < nb_corners());
            corner_adjacent_facets_[c] = f;
        }


        /**
         * \brief Sets the facet adjacent to a corner
         * \param[in] c index of the corner
         * \param[in] f index of the facet
         */
        void set_corner_adjacent_facet(index_t c, index_t f) {
            geo_debug_assert(c < nb_corners());
            corner_adjacent_facets_[c] = signed_index_t(f);
        }

        /**
         * @}
         * \name Facet access and modification
         * @{
         */

        /**
         * \brief Gets the number of facets.
         */
        index_t nb_facets() const {
            return nb_facets_;
        }

        /**
         * \brief Gets the first corner of a facet.
         * \param[in] f index of the facet
         * \return index of the first corner of facet \p f
         */
        index_t facet_begin(index_t f) const {
            geo_debug_assert(f < nb_facets());
            return triangulated_ ? f * 3 : facet_ptr_[f];
        }

        /**
         * \brief Gets one position past the last corner of a facet
         * \param[in] f index of the facet
         * \return one position past the index of the last corner of
         *  facet \p f
         */
        index_t facet_end(index_t f) const {
            geo_debug_assert(f < nb_facets());
            return triangulated_ ? (f + 1) * 3 : facet_ptr_[f + 1];
        }

        /**
         * \brief Gets the number of corners in a facet.
         * \param[in] f index of the facet
         * \return number of corners in facet \p f
         */
        index_t facet_size(index_t f) const {
            return triangulated_ ? 3 : facet_end(f) - facet_begin(f);
        }

        /**
         * \brief Gets the successor of a corner around a facet
         * \param[in] f index of the facet
         * \param[in] c index of the corner
         * \return index of the successor of corner \p c around facet
         *  \p f
         * \pre The corner \p c belongs to facet \p f
         */
        index_t next_around_facet(index_t f, index_t c) const {
            geo_debug_assert(c >= facet_begin(f) && c < facet_end(f));
            return c + 1 == facet_end(f) ? facet_begin(f) : c + 1;
        }

        /**
         * \brief Gets the predeccessor of a corner around a facet.
         * \param[in] f index of the facet
         * \param[in] c index of the corner
         * \return index of the predecessor of corner \p c around facet
         *  \p f
         * \pre The corner \p c belongs to facet \p f
         */
        index_t prev_around_facet(index_t f, index_t c) const {
            geo_debug_assert(c >= facet_begin(f) && c < facet_end(f));
            return c == facet_begin(f) ? facet_end(f) - 1 : c - 1;
        }

        /**
         * \brief Tests whether an edge exists in a facet.
         * \param[in] f index of the facet
         * \param[in] v1 index of the first extremity of the edge
         * \param[in] v2 index of the second extremity of the edge
         * \return true if the oriented edge (\p v1, \p v2) exists in facet
         *  \p f, false otherwise
         */
        inline bool has_edge(index_t f, index_t v1, index_t v2) const;

        /**
         * \brief Gets the region associated with a facet.
         * \param[in] f the facet
         * \return the index of the facet region
         * \pre the mesh needs to have the MESH_FACET_REGION attribute
         */
        signed_index_t facet_region(index_t f) const {
            geo_debug_assert(f < nb_facets());
            geo_debug_assert(has_attribute(GEO::MESH_FACET_REGION));
            return facet_regions_[f];
        }

        /**
         * \brief Sets the region associated with a facet.
         * \param[in] f index of the facet
         * \param[in] facet_region index of the facet region
         * \pre the mesh needs to have the MESH_FACET_REGION attribute
         */
        void set_facet_region(index_t f, signed_index_t facet_region) {
            geo_debug_assert(f < nb_facets());
            geo_debug_assert(has_attribute(GEO::MESH_FACET_REGION));
            facet_regions_[f] = facet_region;
        }

        /**
         * @}
         * \name Tetrahedra access and modification
         * @{
         */

        
        /**
         * \brief Gets the number of tetrahedra.
         * \return the number of tetrahedra
         */
        index_t nb_tets() const {
            geo_debug_assert(tetrahedralized_);            
            return nb_cells();
        }

        /**
         * \brief Gets a vertex of a tetrahedron.
         * \param[in] t index of the tetrahedron
         * \param[in] lv local index of the vertex in the
         *  tetrahedron (0,1,2 or 3)
         * \return the global index of tetrahedron \p t's \p lv%th vertex
         * \pre is_tetrahedralized()
         */
        index_t tet_vertex_index(index_t t, index_t lv) const {
            geo_debug_assert(tetrahedralized_);
            geo_debug_assert(t < nb_tets());
            geo_debug_assert(lv < 4);
            return cell_vertices_[4 * t + lv];
        }

        /**
         * \brief Finds the local index of a vertex in a tetrahedron.
         * \param[in] t index of the tetrahedron
         * \param[in] v global index of the vertex
         * \return the local index (0,1,2 or 3) of the vertex in
         *  tetrahedron \p t or -1 if \p t is not incident to \p v
         * \pre is_tetrahedralized()
         */
        signed_index_t find_tet_vertex(index_t t, index_t v) const {
            geo_debug_assert(tetrahedralized_);            
            geo_debug_assert(t < nb_tets());
            geo_debug_assert(v < nb_vertices());
            for(index_t lv = 0; lv < 4; lv++) {
                if(cell_vertices_[4 * t + lv] == v) {
                    return signed_index_t(lv);
                }
            }
            return -1;
        }

        /**
         * \brief Gets a vertex of a tetrahedron by local facet
         *  index and local vertex index in facet.
         * \param[in] t global index of the tetrahedron
         * \param[in] lf local facet index (0,1,2 or 3)
         * \param[in] lv local vertex index in facet (0,1 or 2)
         * \return the global index of vertex \p lv in facet \p lf of
         *  tetrahedron \p t
         * \pre is_tetrahedralized()
         */
        index_t tet_facet_vertex_index(
            index_t t, index_t lf, index_t lv
        ) const {
            geo_debug_assert(tetrahedralized_);            
            geo_debug_assert(t < nb_tets());
            geo_debug_assert(lf < 4);
            geo_debug_assert(lv < 3);
            return cell_vertices_[4 * t + local_tet_facet_vertex_index(lf,lv)];
        }

        /**
         * \brief Finds the local index of a facet in a tetrahedron
         *  by the global indices of its vertices.
         * \param[in] t index of the tetrahedron
         * \param[in] v1 global index of the first vertex
         * \param[in] v2 global index of the second vertex
         * \param[in] v3 global index of the third vertex
         * \return the local index (0,1,2 or 3) of the facet of
         *  \p t that has \p v1, \p v2, \p v3 as vertices modulo a
         *  circular permutation, or -1 if such a facet does not
         *  exist in \p t.
         * \pre is_tetrahedralized()
         */
        signed_index_t find_tet_facet(
            index_t t, index_t v1, index_t v2, index_t v3
        ) const {
            geo_debug_assert(tetrahedralized_);            
            for(index_t lf = 0; lf < 4; ++lf) {
                index_t w1 = tet_facet_vertex_index(t, lf, 0);
                index_t w2 = tet_facet_vertex_index(t, lf, 1);
                index_t w3 = tet_facet_vertex_index(t, lf, 2);
                if(
                    (v1 == w1 && v2 == w2 && v3 == w3) ||
                    (v1 == w2 && v2 == w3 && v3 == w1) ||
                    (v1 == w3 && v2 == w1 && v3 == w2)
                ) {
                    return signed_index_t(lf);
                }
            }
            return -1;
        }

        /**
         * \brief Finds the local index of the facet along
         *  which a tetrahedron is adjacent to another one.
         * \param[in] t1 index of the first tetrahedron
         * \param[in] t2 index of the second tetrahedron
         * \return lf in (0,1,2,3) such that tet_adjacent(t1,lf)==t2 or
         *  -1 if t2 is not adjacent to t1.
         * \pre is_tetrahedralized()
         */
        signed_index_t find_tet_adjacent(
            index_t t1, index_t t2
        ) const {
            geo_debug_assert(tetrahedralized_);            
            for(index_t lf = 0; lf < 4; ++lf) {
                if(tet_adjacent(t1, lf) == signed_index_t(t2)) {
                    return signed_index_t(lf);
                }
            }
            return -1;
        }

        /**
         * \brief Gets a pointer to a vertex index
         *  associated with a tetrahedron.
         * \param[in] t index of the tetrahedron
         * \param[in] lv local index of the vertex in the
         *  tetrahedron (0,1,2 or 3)
         * \return a pointer to the global index of
         *  tetrahedron \p t's \p lv%th vertex
         * \pre is_tetrahedralized()
         */
        index_t* tet_vertex_index_ptr(index_t t, index_t lv) {
            geo_debug_assert(tetrahedralized_);            
            geo_debug_assert(t < nb_tets());
            geo_debug_assert(lv < 4);
            return &(cell_vertices_[4 * t + lv]);
        }

        /**
         * \brief Gets a tetrahedron adjacent to another one
         * \param[in] t index of the tetrahedron
         * \param[in] lf local index of the facet in the
         *  tetrahedron (0,1,2 or 3). Face i is opposite to
         *  vertex i.
         * \return the global index of the tetrahedron adjacent
         *  to \p t along face \p lf, or -1 if on border
         * \pre is_tetrahedralized()
         */
        signed_index_t tet_adjacent(index_t t, index_t lf) const {
            geo_debug_assert(tetrahedralized_);            
            geo_debug_assert(t < nb_tets());
            geo_debug_assert(lf < 4);
            return cell_adjacents_[4 * t + lf];
        }

        /**
         * \brief Sets a vertex of a tetrahedron..
         * \param[in] t index of the tetrahedron
         * \param[in] lv local index of the vertex in the tetrahedron
         *  (0,1,2 or 3)
         * \param[in] v global index of the vertex
         * \pre is_tetrahedralized()
         */
        void set_tet_vertex_index(index_t t, index_t lv, index_t v) {
            geo_debug_assert(tetrahedralized_);            
            geo_debug_assert(t < nb_tets());
            geo_debug_assert(lv < 4);
            geo_debug_assert(v < nb_vertices());
            cell_vertices_[4 * t + lv] = v;
        }

        /**
         * \brief Sets a tetrahedron adjacency
         * \param[in] t index of the tetrahedron
         * \param[in] lf local index of the facet in the
         *  tetrahedron (0,1,2 or 3). Face i is opposite to
         *  vertex i.
         * \param[in] t2 global index of the tetrahedron adjacent
         *  to \p t along face \p lf, or -1 if on border
         * \pre is_tetrahedralized()
         */
        void set_tet_adjacent(index_t t, index_t lf, signed_index_t t2) {
            geo_debug_assert(tetrahedralized_);            
            geo_debug_assert(t < nb_tets());
            geo_debug_assert(lf < 4);
            geo_debug_assert(t2 < signed_index_t(nb_tets()));
            cell_adjacents_[4 * t + lf] = t2;
        }

        /**
         * \brief Sets a tetrahedron adjacency
         * \param[in] t index of the tetrahedron
         * \param[in] lf local index of the facet in the
         *  tetrahedron (0,1,2 or 3). Face i is opposite to
         *  vertex i.
         * \param[in] t2 global index of the tetrahedron adjacent
         *  to \p t along face \p lf
         * \pre is_tetrahedralized()
         */
        void set_tet_adjacent(index_t t, index_t lf, index_t t2) {
            geo_debug_assert(tetrahedralized_);            
            geo_debug_assert(t < nb_tets());
            geo_debug_assert(lf < 4);
            geo_debug_assert(t2 < nb_tets());
            cell_adjacents_[4 * t + lf] = signed_index_t(t2);
        }

        /**
         * \brief Gets the region associated with a tet.
         * \param[in] t index of the tet
         * \return the index of the tet region
         * \pre is_tetrahedralized() && has_attribute(MESH_CELL_REGION)
         */
        signed_index_t tet_region(index_t t) const {
            geo_debug_assert(tetrahedralized_);
            geo_debug_assert(has_attribute(GEO::MESH_CELL_REGION));            
            return cell_region(t);
        }

        /**
         * \brief Sets the region associated with a tet.
         * \param[in] t index of the tet
         * \param[in] tet_region index of the tet region
         * \pre is_tetrahedralized() && has_attribute(MESH_CELL_REGION)
         */
        void set_tet_region(index_t t, signed_index_t tet_region) {
            geo_debug_assert(tetrahedralized_);
            geo_debug_assert(has_attribute(GEO::MESH_CELL_REGION));
            set_cell_region(t,tet_region);
        }

        /**
         * @}
         * \name Cells access and modification
         * @{
         */

        /**
         * \brief Gets the number of tetrahedra.
         * \return the number of tetrahedra
         */
        index_t nb_cells() const {
            return nb_cells_;
        }

        /**
         * \brief Gets the type of a cell.
         * \param[in] c index of the cell
         * \return the type of the cell
         */
        GEO::MeshCellType cell_type(index_t c) const {
            geo_debug_assert(c < nb_cells());
            return
                tetrahedralized_ ?
                GEO::MESH_TET : GEO::MeshCellType(cell_types_[c]);
        }

        
        /**
         * \brief Gets the number of vertices of a cell.
         * \param[in] c index of the cell
         * \return the number of vertices of the cell
         */
        index_t cell_nb_vertices(index_t c) const {
            geo_debug_assert(c < nb_cells());
            return tetrahedralized_ ? 4 : cell_descriptor(c).nb_vertices;
        }

        /**
         * \brief Gets the index of a cell vertex by cell index 
         *  and vertex local index.
         * \param[in] c index of the cell
         * \param[in] lv local index of the vertex, 
         *  in 0 .. cell_nb_vertices(c)-1
         * \return the global index of the vertex
         */
        index_t cell_vertex_index(index_t c, index_t lv) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(lv < cell_nb_vertices(c));
            return cell_vertices_[cell_vertices_begin(c) + lv];
        }


        /**
         * \brief Gets a pointer to the vertex index associated
         *  with the first vertex of a given cell.
         * \param[in] c the index of the cell
         * \return a pointer to the vertex index associated with 
         *  the first vertex of cell \p c
         */
        index_t* cell_vertex_index_ptr(index_t c) {
            geo_debug_assert(c < nb_cells());
            return &(cell_vertices_[cell_vertices_begin(c)]);
        }

        
        /**
         * \brief Gets the number of facets of a cell.
         * \param[in] c index of the cell
         * \return the number of facets of the cell
         */
        index_t cell_nb_facets(index_t c) const {
            geo_debug_assert(c < nb_cells());
            return tetrahedralized_ ? 4 : cell_descriptor(c).nb_facets;
        }

        /**
         * \brief Gets the cell adjacents to a cell by local facet
         *  index.
         * \param[in] c index of the cell
         * \param[in] lf local index of the facet in \p c, 
         *  in 0..cell_nb_facets(c)
         * \return the index of the cell adjacent to \p c accros \p lf
         *  or -1 on the border.
         */
        signed_index_t cell_adjacent(index_t c, index_t lf) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(lf < cell_nb_facets(c));
            return cell_adjacents_[cell_adjacents_begin(c) + lf];
        }

        /**
         * \brief Sets a cell adjacency
         * \param[in] c global index of the cell
         * \param[in] lf local index of the facet in the
         *  tetrahedron (in 0..cell_nb_facets(c)). 
         * \param[in] c2 global index of the cell adjacent
         *  to \p c along face \p lf, or -1 if on border
         */
        void set_cell_adjacent(index_t c, index_t lf, signed_index_t c2) {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(lf < cell_nb_facets(c));
            geo_debug_assert(c2 < signed_index_t(nb_cells()));
            cell_adjacents_[cell_adjacents_begin(c) + lf] = c2;
        }

        /**
         * \brief Gets the number of vertices of a cell's facet
         * \param[in] c index of the cell
         * \param[in] lf local index of the facet 
         *  in \p c, in 0..cell_nb_facets(c)-1
         * \return the number of vertices in the facet \p lf of cell \p c
         */
        index_t cell_facet_nb_vertices(index_t c, index_t lf) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(lf < cell_nb_facets(c));
            return cell_descriptor(c).nb_vertices_in_facet[lf];
        }
        
        /**
         * \brief Gets a cell vertex by local facet and local vertex index.
         * \param[in] c index of the cell
         * \param[in] lf local index of the facet in \p c, 
         *   in 0..cell_nb_facets(c)-1
         * \param[in] lv local vertex index in the facet \p lf of cell \p c, 
         *   in 0..cell_facet_nb_vertices(c,lf)-1
         * \return the number of vertices in the facet \p lf of cell \p c
         */
        index_t cell_facet_vertex_index(
            index_t c, index_t lf, index_t lv
        ) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(lf < cell_nb_facets(c));
            geo_debug_assert(lv < cell_facet_nb_vertices(c,lf));
            return cell_vertex_index(
                c, cell_descriptor(c).facet_vertex[lf][lv]
            );
        }

        /**
         * \brief Tests whether two cell facets can be connected.
         * \details Two cell facets can be connected if they have the
         *  same vertices in reverse order.
         * \param[in] c1 index of the first cell
         * \param[in] f1 index of the first facet in \p c1
         * \param[in] c2 index of the second cell
         * \param[in] f2 index of the second facet in \p c2
         * \retval true if \p c1 and \p c2 can be connected by \p f1 and \p f2
         * \retval false otherwise
         */
        bool cell_facets_match(
            index_t c1, index_t f1, index_t c2, index_t f2
        ) const;
        
        /**
         * \brief Finds the local index of a vertex in a cell.
         * \param[in] c index of the cell
         * \param[in] v global index of the vertex
         * \return the local index 
         *  (in 0..cell_nb_vertices(c)-1) of the vertex in
         *  cell \p c or -1 if \p c is not incident to \p v
         */
        signed_index_t find_cell_vertex(index_t c, index_t v) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(v < nb_vertices());
            for(index_t lv=0; lv<cell_nb_vertices(c); ++lv) {
                if(cell_vertex_index(c,lv) == v) {
                    return signed_index_t(lv);
                }
            }
            return -1;
        }
        
        /**
         * \brief Finds the local index of a facet in a cell
         *  that can be connected to a facet of another cell
         * \param[in] c1 index of the cell
         * \param[in] c2 index of the other cell
         * \param[in] f2 facet of the other cell
         * \return the local index (in 0 .. cell_nb_facets(c1)) of the facet of
         *  \p c1 that has the same vertices as \p f2 in \p c2 in reverse order,
         *  modulo a circular permutation, or -1 if such a facet does not
         *  exist in \p c1.
         */
        signed_index_t find_cell_facet(
            index_t c1, index_t c2, index_t f2
        ) const {
            for(index_t f1=0; f1<cell_nb_facets(c1); ++f1) {
                if(cell_facets_match(c1,f1,c2,f2)) {
                    return signed_index_t(f1);
                }
            }
            return -1;
        }
        
        /**
         * \brief Gets the number of edges of a cell.
         * \param[in] c index of the cell
         * \return the number of edges of the cell
         */
        index_t cell_nb_edges(index_t c) const {
            geo_debug_assert(c < nb_cells());
            return cell_descriptor(c).nb_edges;
        }

        /**
         * \brief Gets the number of edges of a cell.
         * \param[in] c index of the cell
         * \param[in] le local index of the edge, in 0..cell_nb_edges(c)-1
         * \param[in] lv local index of the vertex in the edge, in {0,1}
         */
        index_t cell_edge_vertex_index(
            index_t c, index_t le, index_t lv
        ) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(le < cell_nb_edges(c));
            geo_debug_assert(lv < 2);
            return cell_vertex_index(c, cell_descriptor(c).edge_vertex[le][lv]);
        }
        
        /**
         * \brief Gets the region associated with a cell.
         * \param[in] c index of the cell
         * \return the index of the cell region
         * \pre the mesh needs to have the MESH_CELL_REGION attribute
         */
        signed_index_t cell_region(index_t c) const {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(has_attribute(GEO::MESH_CELL_REGION));
            return cell_regions_[c];
        }

        /**
         * \brief Sets the region associated with a cell.
         * \param[in] t index of the cell
         * \param[in] cell_region index of the cell region
         * \pre the mesh needs to have the MESH_CELL_REGION attribute
         */
        void set_cell_region(index_t c, signed_index_t cell_region) {
            geo_debug_assert(c < nb_cells());
            geo_debug_assert(has_attribute(GEO::MESH_CELL_REGION));
            cell_regions_[c] = cell_region;
        }


        /**
         * \brief Adds a new cell to this mesh.
         * \param[in] type the type of the cell
         * \param[in] vertices an array of vertex indices. Its size depends on 
         *  \p type (4 for MESH_TET, 8 for MESH_HEX, 6 for MESH_PRISM, 
         *  5 for MESH_PYRAMID and 4 for MESH_CONNECTOR)
         * \param[in] adjacents an array of indices of adjacent cells. 
         *  Its size depends on \p type (4 for MESH_TET, 6 for MESH_HEX, 
         *  5 for MESH_PRISM, 5 for MESH_PYRAMID and 3 for MESH_CONNECTOR), 
         *  or nil if unspecified
         * \param[in] region the cell region to be attached to this cell 
         *  (if has_attribute(MESH_CELL_REGION))
         * \return the index of the newly created cell
         */
        index_t add_cell(
            GEO::MeshCellType type,
            index_t* vertices, signed_index_t* adjacents=nil,
            signed_index_t region=0
        );

        /**
         * \brief Adds a new tetrahedron to this mesh.
         * \param[in] v1 index of the 1st vertex 
         * \param[in] v2 index of the 2nd vertex
         * \param[in] v3 index of the 3rd vertex
         * \param[in] v4 index of the 4th vertex         
         * \param[in] adj1 index of the cell adjacent across the 1st facet or -1
         * \param[in] adj2 index of the cell adjacent across the 2nd facet or -1
         * \param[in] adj3 index of the cell adjacent across the 3rd facet or -1
         * \param[in] 4dj3 index of the cell adjacent across the 4th facet or -1
         * \param[in] region the cell region to be attached to this cell 
         *  (if has_attribute(MESH_CELL_REGION))
         * \return the index of the newly created cell
         */         
        index_t add_tet(
            index_t v1, index_t v2, index_t v3, index_t v4,
            signed_index_t adj1=-1, signed_index_t adj2=-1,
            signed_index_t adj3=-1, signed_index_t adj4=-1,
            signed_index_t region=0
        );
        
        /**
         * \brief Adds a new hexahedron to this mesh.
         * \param[in] v1 index of the 1st vertex 
         * \param[in] v2 index of the 2nd vertex
         * \param[in] v3 index of the 3rd vertex
         * \param[in] v4 index of the 4th vertex         
         * \param[in] v5 index of the 5th vertex         
         * \param[in] v6 index of the 6th vertex         
         * \param[in] v7 index of the 7th vertex         
         * \param[in] v8 index of the 8th vertex         
         * \param[in] adj1 index of the cell adjacent across the 1st facet or -1
         * \param[in] adj2 index of the cell adjacent across the 2nd facet or -1
         * \param[in] adj3 index of the cell adjacent across the 3rd facet or -1
         * \param[in] adj4 index of the cell adjacent across the 4th facet or -1
         * \param[in] adj5 index of the cell adjacent across the 5th facet or -1
         * \param[in] adj6 index of the cell adjacent across the 6th facet or -1
         * \param[in] region the cell region to be attached to this cell 
         *  (if has_attribute(MESH_CELL_REGION))
         * \return the index of the newly created cell
         */         
        index_t add_hex(
            index_t v1, index_t v2, index_t v3, index_t v4,
            index_t v5, index_t v6, index_t v7, index_t v8,             
            signed_index_t adj1=-1, signed_index_t adj2=-1,
            signed_index_t adj3=-1, signed_index_t adj4=-1,
            signed_index_t adj5=-1, signed_index_t adj6=-1,
            signed_index_t region=0
        );

        /**
         * \brief Adds a new prism to this mesh.
         * \param[in] v1 index of the 1st vertex 
         * \param[in] v2 index of the 2nd vertex
         * \param[in] v3 index of the 3rd vertex
         * \param[in] v4 index of the 4th vertex         
         * \param[in] v5 index of the 5th vertex         
         * \param[in] v6 index of the 6th vertex         
         * \param[in] adj1 index of the cell adjacent across the 1st facet or -1
         * \param[in] adj2 index of the cell adjacent across the 2nd facet or -1
         * \param[in] adj3 index of the cell adjacent across the 3rd facet or -1
         * \param[in] adj4 index of the cell adjacent across the 4th facet or -1
         * \param[in] adj5 index of the cell adjacent across the 5th facet or -1
         * \param[in] region the cell region to be attached to this cell 
         *  (if has_attribute(MESH_CELL_REGION))
         * \return the index of the newly created cell
         */         
        index_t add_prism(
            index_t v1, index_t v2, index_t v3, 
            index_t v4, index_t v5, index_t v6, 
            signed_index_t adj1=-1, signed_index_t adj2=-1,
            signed_index_t adj3=-1, signed_index_t adj4=-1,
            signed_index_t adj5=-1,
            signed_index_t region=0
        );

        /**
         * \brief Adds a new pyramid to this mesh.
         * \param[in] v1 index of the 1st vertex 
         * \param[in] v2 index of the 2nd vertex
         * \param[in] v3 index of the 3rd vertex
         * \param[in] v4 index of the 4th vertex         
         * \param[in] v5 index of the 5th vertex         
         * \param[in] adj1 index of the cell adjacent across the 1st facet or -1
         * \param[in] adj2 index of the cell adjacent across the 2nd facet or -1
         * \param[in] adj3 index of the cell adjacent across the 3rd facet or -1
         * \param[in] adj4 index of the cell adjacent across the 4th facet or -1
         * \param[in] adj5 index of the cell adjacent across the 5th facet or -1
         * \param[in] region the cell region to be attached to this cell 
         *  (if has_attribute(MESH_CELL_REGION))
         * \return the index of the newly created cell
         */         
        index_t add_pyramid(
            index_t v1, index_t v2, index_t v3, 
            index_t v4, index_t v5, 
            signed_index_t adj1=-1, signed_index_t adj2=-1,
            signed_index_t adj3=-1, signed_index_t adj4=-1,
            signed_index_t adj5=-1,
            signed_index_t region=0
        );

        /**
         * \brief Adds a new connector to this mesh.
         * \details Connectors are used to represent no-conformal 
         *  connections between a hexahedron and two tetrahedra.
         * \param[in] v1 index of the 1st vertex 
         * \param[in] v2 index of the 2nd vertex
         * \param[in] v3 index of the 3rd vertex
         * \param[in] v4 index of the 4th vertex         
         * \param[in] adj1 index of the cell adjacent across the 1st facet or -1
         * \param[in] adj2 index of the cell adjacent across the 2nd facet or -1
         * \param[in] adj3 index of the cell adjacent across the 3rd facet or -1
         * \param[in] region the cell region to be attached to this cell 
         *  (if has_attribute(MESH_CELL_REGION))
         * \return the index of the newly created cell
         */         
        index_t add_connector(
            index_t v1, index_t v2, index_t v3, index_t v4, 
            signed_index_t adj1=-1, signed_index_t adj2=-1,
            signed_index_t adj3=-1,
            signed_index_t region=0
        );
        
        /**
         * @}
         * \name Global Mesh Modification
         * @{
         */

        /**
         * \brief Reconstructs facet-to-facet links.
         * \details Does not resolve non-manifold
         * configurations and orientation problems 
         *  (one may use mesh_repair() instead).
         */
        void connect_facets();

        /**
         * \brief Reconstructs tetrahedron-to-tetrahedron links.
         * \details Requires correct orientation of the tetrahedra, and
         *  does not resolve non-manifold configurations.
         * \pre is_tetrahedralized()
         */
        void connect_tets();

        /**
         * \brief Tests whether a triangular facet matches a quad facet.
         * \details Used to detect non-conformal configurations that should
         *  be resolved by a connector.
         * \param[in] c1 index of the first cell
         * \param[in] f1 index of a triangular facet in \p c1
         * \param[in] c2 index of the second cell
         * \param[in] f2 index of a quadrangular facet in \p c2
         * \retval true if the three vertices of \p f1 appear in \p f2
         *   in reverse order
         * \retval false otherwise
         */
        bool triangular_facet_matches_quad_facet(
            index_t c1, index_t lf1,
            index_t c2, index_t lf2
        ) const;
        

        /**
         * \brief Tests whether two triangular cell facets have a common edge.
         * \param[in] c1 index of the first cell
         * \param[in] f1 index of a triangular facet of \p c1
         * \param[in] c2 index of the second cell
         * \param[in] f2 index of a triangular facet of \p c2
         * \param[out] e1 index of the common edge in \p f1
         *  or index_t(-1) if no such edge exists
         * \param[out] e2 index of the common edge in \p f2
         *  or index_t(-1) if no such edge exists
         * \retval true if \p f1 and \p f2 have a common edge
         * \retval false otherwise
         */
        bool triangular_facets_have_common_edge(
            index_t c1, index_t f1,
            index_t c2, index_t f2,
            index_t& e1, index_t& e2
        ) const;
        
        /**
         * \brief Reconstructs cells-to-cells links.
         * \details Requires correct orientation of the cells, and
         *  does not resolve non-manifold configurations.
         */
        void connect_cells();

        /**
         * \brief Computes the (surfacic) border 
         *  of a tetrahedralized volumetric mesh.
         * \details The computed border replaces the facets of this mesh.
         * \pre is_tetrahedralized()
         */
        void compute_tets_boundaries();

        /**
         * \brief Computes the (surfacic) border of a volumetric mesh.
         * \details The computed border replaces the facets of this mesh.
         */
        void compute_cells_boundaries();
        
        /**
         * \brief Displays number of vertices, facets and borders.
         */
        void show_stats(const std::string& tag = "Mesh") const;

        /**
         * \brief Assigns the MESH_FACET_REGION attribute from a vector.
         * \param[in] facet_regions array of facet region indices
         * \param[in] steal_args if set, facet_regions are 'stolen' from the
         * arguments (using vector::swap).
         */
        void assign_facet_regions(
            GEO::vector<signed_index_t>& facet_regions,
            bool steal_args
        ) {
            geo_debug_assert(facet_regions.size() == nb_facets());
            if(steal_args) {
                facet_regions_.swap(facet_regions);
            } else {
                facet_regions_ = facet_regions;
            }
        }


        /** @}          
         * \name Low-level cell access and modification.
         * @{
         */

        /**
         * \brief Lookup tables that describe the combinatorics
         *  of each cell type.
         */
        struct CellDescriptor {
            /** Number of vertices */
            index_t nb_vertices;
            
            /** Number of facets */
            index_t nb_facets;

            /** Number of vertices in each facet */
            index_t nb_vertices_in_facet[6];

            /** 
             * Cell vertex index by (facet index,facet vertex index).
             */
            index_t facet_vertex[6][4];

            /**
             * Number of edges */
            index_t nb_edges;

            /**
             * Cell vertex index by (edge index, edge vertex index).
             */
            index_t edge_vertex[12][2];
        };
        
        /**
         * \brief Gets a cell descriptor by cell index.
         * \param[in] c index of the cell
         * \return a const reference to the descriptor of cell \p c
         */
        const CellDescriptor& cell_descriptor(index_t c) const {
            geo_debug_assert(c < nb_cells());
            return cell_type_to_cell_descriptor(cell_type(c));
        }

        /**
         * \brief Gets the first index of cell's elements.
         * \param[in] c index of the cell
         * \return the first index for accessing
         *  cell vertex indices, using 
         *  cell_vertices_[cell_vertices_begin(c) + lv])
         */
        index_t cell_vertices_begin(index_t c) const {
            geo_debug_assert(c < nb_cells());
            return tetrahedralized_ ? (c*4) : cell_ptr_[c];
        }

        /**
         * \brief Gets the first index of cell's elements.
         * \param[in] c index of the cell
         * \return the first index for accessing
         *  and adjacent cells cell_adjacents_[cell_adjacents_begin(c) + lf])
         */
        index_t cell_adjacents_begin(index_t c) const {
            geo_debug_assert(c < nb_cells());
            return tetrahedralized_ ? (c*4) : cell_ptr_[c];
        }
        
        /**
         * \brief Gets a cell descriptor by cell index.
         * \param[in] c index of the cell
         * \return a const reference to the descriptor of cell \p c
         */
        static CellDescriptor& cell_type_to_cell_descriptor(
            GEO::MeshCellType t
        ) {
            geo_debug_assert(index_t(t) < GEO::MESH_NB_CELL_TYPES);
            return *cell_type_to_cell_descriptor_[t];
        }
        

        /**
         * @}
         */

    protected:

        /**
         * \brief Tests whether a vector contains a non-zero value.
         * \details This function is used internally by remove_vertices() and
         * remove_facets().
         * \param[in] I a vector of signed integers
         * \retval true if \p I contains at least a non-zero value
         * \retval false otherwise
         */
        static bool has_non_zero_value(const GEO::vector<index_t>& I) {
            for(index_t i = 0; i < I.size(); i++) {
                if(I[i] != 0) {
                    return true;
                }
            }
            return false;
        }
        
        /**
         * \brief Maps a cell type to the associated cell descriptor.
         */
        static CellDescriptor*
        cell_type_to_cell_descriptor_[GEO::MESH_NB_CELL_TYPES];
        
        static CellDescriptor tet_descriptor_;
        static CellDescriptor hex_descriptor_;
        static CellDescriptor prism_descriptor_;
        static CellDescriptor pyramid_descriptor_;
        static CellDescriptor connector_descriptor_;        

        
        GEO::vector<index_t> corner_vertices_;
        GEO::vector<signed_index_t> corner_adjacent_facets_;
        GEO::vector<index_t> facet_ptr_;
        index_t nb_facets_;

        GEO::vector<index_t> cell_vertices_;
        GEO::vector<signed_index_t> cell_adjacents_;
        GEO::vector<index_t> cell_ptr_;
        GEO::vector<GEO::Numeric::int8> cell_types_;
        index_t nb_cells_;
        
        GEO::vector<signed_index_t> facet_regions_;
        GEO::vector<signed_index_t> cell_regions_;

        bool triangulated_;
        bool tetrahedralized_;
        index_t nb_vertices_;
        coord_index_t dimension_;

        bool in_facet_;
        index_t facets_reserve_;

        MeshAttributes attributes_;
    };

    /**
     * \brief Mesh stores facet connectivity of a surfacic or volumetric mesh.
     *
     * \details Facets store pointers to adjacent facets.
     * Facets are stored using the CRS (compressed row storage)
     * format, i.e. the vertices of the f-th facet are traversed
     * using:
     * \code
     * for(index_t c=M.facet_begin(f); c<M.facet_end(f); c++) {
     *    index_t v = M.corner_vertex_index(c) ;
     *    do something with M.vertex_ptr(v)
     * }
     * \endcode
     * \see MeshBuilder
     * \tparam COORD_T the numeric type used to represent point coordinates
     */
    template <class COORD_T>
    class Mesh : public MeshBase {

        /** \brief This class type */
        typedef Mesh thisclass;

    public:
        /**
         * \brief The numeric type used to represent point coordinates.
         */
        typedef COORD_T coord_t;

        // ==== [Mesh implementation] ==========================================

        /**
         * \brief Creates a new empty structure
         */
        explicit Mesh(coord_index_t dim = 0) : MeshBase(dim) {
        }

        
        /**
         * \brief Resets the mesh structure
         * \param[in] keep_memory if set, memory is not
         * deallocated (may be useful for dynamic
         * meshes).
         */
        void clear(bool keep_memory = false) {
            geo_assert(!in_facet_);
            MeshBase::clear(keep_memory);
            if(keep_memory) {
                vertices_.resize(0);
                weights_.resize(0);
            } else {
                vertices_.clear();
                weights_.clear();
            }
        }

        /**
         * \brief Pre-allocates the vertices table.
         * \details Pre-reserve memory for \p nb_vertices vertices of the
         * mesh's dimension. Calling this function is not mandatory (but
         * it improves performances).
         * \param[in] nb_vertices number of vertices to pre-allocate
         */
        void reserve_vertices(index_t nb_vertices) {
            geo_assert(dimension() != 0);
            vertices_.reserve(nb_vertices * dimension());
            update_cached_variables();
        }

        /**
         * \brief Pre-allocates the vertices table and sets mesh dimension.
         * \details Pre-reserve memory for \p nb_vertices vertices of the
         * given \p dimension.
         * Calling this function is not mandatory (but
         * it improves performances.
         * \param[in] nb_vertices number of vertices to pre-allocate
         * \param[in] dim number of coordinates of the vertices
         */
        void reserve_vertices(index_t nb_vertices, coord_index_t dim) {
            geo_assert(dim != 0);
            dimension_ = dim;
            vertices_.reserve(nb_vertices * dim);
            update_cached_variables();
        }

        /**
         * \brief Allocates and initializes the vertices table
         * \details Reserve memory for \p nb_vertices vertices of the
         * given \p dimension. All coordinates are initialized with the value
         * 0.0. It is also possible to create the vertices one by one using
         * MeshBuilder.
         * \param[in] nb_vertices number of vertices to create
         * \param[in] dim number of coordinates of the vertices
         */
        void create_vertices(index_t nb_vertices, coord_index_t dim) {
            geo_assert(dim != 0);
            dimension_ = dim;
            vertices_.assign(nb_vertices * dim, 0.0);
            update_cached_variables();
        }

        /**
         * \brief Assigns all the vertices of a mesh.
         * \param[in] vertices a const reference to a contiguous 
         *  vector of coordinates. It will be copied into this Mesh.
         * \param[in] dim number of coordinates of the vertices.
         */
        void assign_vertices(
            const GEO::vector<coord_t>& vertices, coord_index_t dim
        ) {
            geo_assert(dim != 0);
            geo_debug_assert(vertices.size()%dim == 0);
            dimension_ = dim;
            vertices_ = vertices;
            update_cached_variables();
        }
        
        /**
         * \brief Allocates and initializes the vertices weights table
         * \details Pre-reserve memory for vertices weights. The number of
         * preallocated weights is the same as the current number of
         * pre-allocated vertices. Weights are initialized with the value 1.0.
         */
        void create_vertices_weights() {
            weights_.reserve(vertices_.capacity());
            weights_.assign(vertices_.size(), 1.0);
        }

        /**
         * \brief Clears the vertices weights table
         */
        void remove_vertices_weights() {
            weights_.clear();
        }


        /**
         * \brief Changes the vertices dimension.
         * \details This operation preserves the existing coordinates. If the
         * new dimension is greater than the current dimension, new
         * coordinates are set to 0. If the new dimension is less than the
         * current dimension, old extra coordinates are lost.
         * \param[in] dim the new dimension
         */
        void set_dimension(coord_index_t dim) {
            geo_debug_assert(dim != 0);
            if(dim == dimension()) {
                return;
            }
            coord_index_t copy_dim = GEO::geo_min(dim, dimension());
            GEO::vector<coord_t> old_vertices = vertices_;
            vertices_.assign(nb_vertices() * dim, 0.0);
            for(index_t i = 0; i < nb_vertices(); i++) {
                for(index_t c = 0; c < copy_dim; c++) {
                    vertices_[dim * i + c] = old_vertices[dimension_ * i + c];
                }
            }
            dimension_ = dim;
            update_cached_variables();
        }

        /**
         * \name Vertex access and modification
         * @{
         */

        /**
         * \brief Gets a pointer to a vertex from the index of
         *  a corner incident to it.
         * \param[in] c the index of the corner
         * \return a const pointer to the coordinates of the vertex
         */
        const coord_t* corner_vertex_ptr(index_t c) const {
            return vertex_ptr(corner_vertex_index(c));
        }

        /**
         * \brief Gets the coordinates of a vertex from its index.
         * \param[in] i the index of the vertex
         * \return a const pointer to the coordinates of vertex \p i
         */
        const coord_t* vertex_ptr(index_t i) const {
            geo_debug_assert(i < nb_vertices_);
            return &(vertices_[i * dimension_]);
        }

        /**
         * \brief Gets the coordinates of a vertex from its index.
         * \param[in] i the index of the vertex
         * \return a pointer to the coordinates of vertex \p i
         */
        coord_t* vertex_ptr(index_t i) {
            geo_debug_assert(i < nb_vertices_);
            return &(vertices_[i * dimension_]);
        }

        /**
         * \brief Tests whether vertex weights are stored in this mesh.
         * \return true if vertex weights are stored, false otherwise
         */
        bool has_weights() const {
            return weights_.size() != 0;
        }

        /**
         * \brief Gets the weight associated with a given vertex.
         * \param[in] v the index of the vertex
         * \return the weight associated wigth vertex \p v
         */
        coord_t weight(index_t v) const {
            geo_debug_assert(v < nb_vertices_);
            return has_weights() ? weights_[v] : 1.0;
        }

        /**
         * \brief Gets the weight of a vertex from a corner incident to it.
         * \param[in] c the index of the corner
         * \return the weight of the vertex that corner \p c points at
         */
        coord_t corner_weight(index_t c) const {
            geo_debug_assert(c < nb_corners());
            if(!has_weights()) {
                return 1.0;
            }
            index_t v = corner_vertex_index(c);
            geo_debug_assert(v < nb_vertices_);
            return weights_[v];
        }

        /**
         * \brief Sets the weight of a vertex.
         * \param[in] i index of the vertex
         * \param[in] w weight
         */
        void set_weight(index_t i, double w) {
            geo_debug_assert(i < nb_vertices_);
            if(!has_weights() && w != 1.0) {
                create_vertices_weights();
            }
            if(has_weights()) {
                weights_[i] = coord_t(w);
            }
        }

        /**
         * \brief Computes the area of a facet.
         * \param[in] f index of the facet
         * \param[in] dim dimension that will be used to compute the area
         * \return the area of the facet, obtained by considering the
         *  \p dim first coordinates of the vertices only
         */
        double facet_area(index_t f, coord_index_t dim) const {
            geo_debug_assert(dim <= dimension());
            coord_t result = 0.0;
            // Check for empty facet, should not happen.
            if(facet_end(f) == facet_begin(f)) {
                return result;
            }
            const coord_t* p0 = corner_vertex_ptr(facet_begin(f));
            for(index_t i = facet_begin(f) + 1; i + 1 < facet_end(f); i++) {
                result += GEO::Geom::triangle_area(
                    p0, corner_vertex_ptr(i), corner_vertex_ptr(i + 1), dim
                );
            }
            return result;
        }

        /**
         * \brief Computes the area of a facet.
         * \param[in] f index of the facet
         * \return the area of facet \p f, computed with
         *  all the coordinates of the vertices.
         */
        coord_t facet_area(index_t f) const {
            return facet_area(f, dimension());
        }

        /**
         * @}
         * \name Global Mesh Modification
         * @{
         */

        /**
         * \brief Updates the cached variables.
         * \details Computes the number of vertices from the
         *  size of the vertices array.
         */
        void update_cached_variables() {
            if(dimension_ != 0) {
                nb_vertices_ = vertices_.size() / dimension_;
            } else {
                geo_assert(vertices_.size() == 0);
            }
        }

        /**
         * \brief Removes vertices using a vector of flags
         * \details Removes all vertices \c v for which the flag \c
         * to_remove[v] is non null.
         * \param[in,out]] to_remove vector of flags specifying which vertex to
         * remove. The number of flags in the vector must be the same as the
         * number of vertices.
         * \pre to_remove.size() == nb_vertices()
         * \note Vector \p to_remove is modified.
         */
        void remove_vertices(GEO::vector<index_t>& to_remove) {
            geo_assert(to_remove.size() == nb_vertices());
            if(!has_non_zero_value(to_remove)) {
                return;
            }
            GEO::vector<index_t>& old2new = to_remove;
            index_t vsize = index_t(sizeof(coord_t) * dimension());
            index_t new_nb_vertices = 0;
            for(index_t v = 0; v < nb_vertices(); v++) {
                if(old2new[v] != 0) {
                    old2new[v] = index_t(-1);
                } else {
                    if(v != new_nb_vertices) {
                        GEO::Memory::copy(
                            vertex_ptr(new_nb_vertices), vertex_ptr(v), vsize
                        );
                        if(has_weights()) {
                            set_weight(new_nb_vertices, weight(v));
                        }
                    }
                    old2new[v] = new_nb_vertices;
                    new_nb_vertices++;
                }
            }
            for(index_t c = 0; c < nb_corners(); c++) {
                corner_vertices_[c] = old2new[corner_vertices_[c]];
            }
            vertices_.resize(new_nb_vertices * dimension());
            if(has_weights()) {
                weights_.resize(new_nb_vertices);
            }
            nb_vertices_ = new_nb_vertices;
            update_cached_variables();
        }

        /**
         * \brief Removes all the vertices that have no incident facet
         *   and no incident cell.
         */
        void remove_isolated_vertices() {
            GEO::vector<index_t> to_remove(nb_vertices(), 1);
            for(index_t c = 0; c < nb_corners(); c++) {
                to_remove[corner_vertex_index(c)] = 0;
            }
            for(index_t c=0; c<nb_cells(); ++c) {
                for(index_t lv=0; lv<cell_nb_vertices(c); ++lv) {
                    to_remove[cell_vertex_index(c,lv)] = 0;
                }
            }
            remove_vertices(to_remove);
        }

        /**
         * \brief Removes facets using a vector of flags
         * \details Removes all facets \c f for which the flag \c to_remove[f]
         * is non-zero.
         * \param[in,out] to_remove vector of flags specifying which facet 
         *  to remove. The number of flags in the vector must be the same 
         *  as the number of facets.
         * \param[in] do_remove_isolated_vertices if set to true, isolated
         *  vertices are detected and removed after the facets have been
         *  removed.
         * \pre to_remove.size() == nb_facets()
         * \note Vector \p to_remove is modified.
         */
        void remove_facets(
            GEO::vector<index_t>& to_remove,
            bool do_remove_isolated_vertices=true
        ) {
            // TODO: attributes management
            geo_assert(to_remove.size() == nb_facets());
            if(!has_non_zero_value(to_remove)) {
                return;
            }
            GEO::vector<index_t>& old2new = to_remove;
            index_t new_nb_facets = 0;
            index_t new_nb_corners = 0;
            for(index_t f = 0; f < nb_facets(); f++) {
                if(old2new[f] != 0) {
                    old2new[f] = index_t(-1);
                } else {
                    old2new[f] = new_nb_facets;
                    if(!triangulated_) {
                        facet_ptr_[new_nb_facets] = new_nb_corners;
                    }
                    for(index_t c = facet_begin(f); c != facet_end(f); c++) {
                        if(c != new_nb_corners) {
                            corner_vertices_[new_nb_corners] =
                                corner_vertices_[c];
                            corner_adjacent_facets_[new_nb_corners] =
                                corner_adjacent_facets_[c];
                        }
                        new_nb_corners++;
                    }
                    new_nb_facets++;
                }
            }
            if(!triangulated_) {
                facet_ptr_[new_nb_facets] = new_nb_corners;
            }
            nb_facets_ = new_nb_facets;
            corner_vertices_.resize(new_nb_corners);
            corner_adjacent_facets_.resize(new_nb_corners);
            if(!triangulated_) {
                facet_ptr_.resize(new_nb_facets + 1);
            }
            for(index_t c = 0; c < nb_corners(); c++) {
                signed_index_t f = corner_adjacent_facets_[c];
                if(f != -1) {
                    corner_adjacent_facets_[c] = signed_index_t(old2new[f]);
                }
            }
            if(do_remove_isolated_vertices) {
                remove_isolated_vertices();
            }
        }

        /**
         * \brief Copies a mesh into this Mesh.
         * \details Facet adjacence are not computed.
         * \param[in] dim dimension of the vertices
         * \param[in] vertices coordinates of the vertices
         * \param[in] corner_vertex_index facet to vertex links
         * \param[in] facet_ptr facet pointers
         * \param[in] steal_args if set, the different vectors
         * are 'stolen' from the arguments
         * (using vector::swap).
         */
        void assign(
            coord_index_t dim,
            GEO::vector<coord_t>& vertices,
            GEO::vector<index_t> corner_vertex_index,
            GEO::vector<index_t> facet_ptr,
            bool steal_args
        ) {
            geo_assert(dim != 0);
            dimension_ = dim;
            nb_vertices_ = vertices.size() / dim;
            geo_assert(vertices.size() == dim * nb_vertices_);
            if(steal_args) {
                vertices_.swap(vertices);
                corner_vertices_.swap(corner_vertex_index);
                facet_ptr_.swap(facet_ptr);
            } else {
                vertices_ = vertices;
                corner_vertices_ = corner_vertex_index;
                facet_ptr_ = facet_ptr;
            }
            weights_.clear();
            nb_facets_ = facet_ptr_.size() - 1;
            corner_adjacent_facets_.assign(corner_vertices_.size(), -1);
            triangulated_ = false;
            update_cached_variables();
        }

        /**
         * \brief Copies a triangle mesh into this Mesh.
         * \details Facet adjacence are not computed.
         * \param[in] dim dimension of the vertices
         * \param[in] vertices coordinates of the vertices
         * \param[in] triangles facet to vertex links
         * \param[in] steal_args if set, vertices and triangles
         * are 'stolen' from the arguments
         * (using vector::swap).
         */
        void assign_triangle_mesh(
            coord_index_t dim,
            GEO::vector<coord_t>& vertices,
            GEO::vector<index_t>& triangles,
            bool steal_args
        ) {
            geo_debug_assert((triangles.size()%3)==0);
            geo_assert(dim != 0);
            dimension_ = dim;
            nb_vertices_ = vertices.size() / dim;
            geo_assert(vertices.size() == dim * nb_vertices_);
            if(steal_args) {
                vertices_.swap(vertices);
                corner_vertices_.swap(triangles);
            } else {
                vertices_ = vertices;
                corner_vertices_ = triangles;
            }
            weights_.clear();
            nb_facets_ = corner_vertices_.size() / 3;
            corner_adjacent_facets_.assign(nb_facets_ * 3, -1);
            geo_assert(corner_vertices_.size() == nb_facets_ * 3);
            facet_ptr_.clear();
            triangulated_ = true;
            cell_vertices_.clear();
            cell_adjacents_.clear();
            update_cached_variables();
        }

        /**
         * \brief Copies a tetrahedron mesh into this Mesh.
         * \details Tetrahedron adjacences are not computed.
         * \param[in] dim dimension of the vertices
         * \param[in] vertices coordinates of the vertices
         * \param[in] tets tetrahedron to vertex links
         * \param[in] steal_args if set, vertices and tets
         * are 'stolen' from the arguments
         * (using vector::swap).
         */
        void assign_tet_mesh(
            coord_index_t dim,
            GEO::vector<coord_t>& vertices,
            GEO::vector<index_t>& tets,
            bool steal_args
        ) {
            geo_assert(dim != 0);
            geo_debug_assert((tets.size()%4) == 0);
            dimension_ = dim;
            nb_vertices_ = vertices.size() / dim;
            nb_cells_ = tets.size()/4;
            geo_assert(vertices.size() == dim * nb_vertices_);
            if(steal_args) {
                vertices_.swap(vertices);
                cell_vertices_.swap(tets);
            } else {
                vertices_ = vertices;
                cell_vertices_ = tets;
            }
            cell_adjacents_.assign(cell_vertices_.size(), -1);
            weights_.clear();
            nb_facets_ = 0;
            corner_vertices_.clear();
            corner_adjacent_facets_.clear();
            facet_ptr_.clear();
            tetrahedralized_ = true;
            cell_ptr_.clear();
            cell_types_.clear();
            update_cached_variables();
        }

        /** @} */

    protected:
        GEO::vector<coord_t> vertices_;
        GEO::vector<coord_t> weights_;

    private:
        /**
         * \brief Forbids copy-construction.
         */
        Mesh(const thisclass& rhs);

        /**
         * \brief Forbids copy assignment.
         */
        thisclass& operator= (const thisclass& rhs);

        friend class MeshBuilder<thisclass>;
        friend class MeshMutator<thisclass>;
    };
}

namespace GEO {

    /**
     * \brief The mesh class used for computations (double precision).
     */
    class Mesh : public GEOGen::Mesh<double> {

        /** \brief This class type */
        typedef Mesh thisclass;

        /** \brief The base class of this class */
        typedef GEOGen::Mesh<double> baseclass;

    public:
        /**
         * \brief Constructs a new Mesh
         * \param[in] dim dimension of the vertices, or 0 for unspecified
         */
        explicit Mesh(coord_index_t dim = 0) :
            baseclass(dim) {
        }

    private:
        /**
         * \brief Forbids copy-construction.
         */
        Mesh(const thisclass& rhs);

        /**
         * \brief Forbids assignment.
         */
        thisclass& operator= (const thisclass& rhs);

        friend class GEOGen::MeshBuilder<thisclass>;
        friend class GEOGen::MeshMutator<thisclass>;
    };

    /**
     * \brief The mesh class used for display (single precision).
     */
    class SinglePrecisionMesh : public GEOGen::Mesh<float> {

        /** \brief This class type */
        typedef SinglePrecisionMesh thisclass;

        /** \brief The base class of this class */
        typedef GEOGen::Mesh<float> baseclass;

    public:
        /**
         * \brief Constructs a new SinglePrecisionMesh
         * \param[in] dim dimension of the vertices, or 0 for unspecified
         */
        explicit SinglePrecisionMesh(coord_index_t dim = 0) :
            baseclass(dim) {
        }

    private:
        /**
         * \brief Forbids copy-construction.
         */
        SinglePrecisionMesh(const thisclass& rhs);

        /**
         * \brief Forbids assignment.
         */
        thisclass& operator= (const thisclass& rhs);

        friend class GEOGen::MeshBuilder<thisclass>;
        friend class GEOGen::MeshMutator<thisclass>;
    };
}

#endif

