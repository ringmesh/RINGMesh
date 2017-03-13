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

/*! \author Francois Bonneau */

#include <ringmesh/mesh/mesh.h>

#include <stack>
#include <ringmesh/mesh/geogram_mesh.h>
#include <ringmesh/basic/algorithm.h>

namespace RINGMesh {
    MeshBase::~MeshBase()
    {
        if( vertices_nn_search_ ) delete vertices_nn_search_ ;
    }

    Mesh0D* Mesh0D::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMesh0D::type_name_static() ;
        }
        Mesh0D* mesh = Mesh0DFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "Mesh0D" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "Mesh0D" ) << "Falling back to GeogramMesh0D data structure"
                << std::endl ;

            mesh = new GeogramMesh0D ;
        }
        return mesh ;
    }

    Mesh1D* Mesh1D::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMesh1D::type_name_static() ;
        }
        Mesh1D* mesh = Mesh1DFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "Mesh1D" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "Mesh1D" ) << "Falling back to GeogramMesh1D data structure"
                << std::endl ;

            mesh = new GeogramMesh1D ;
        }
        return mesh ;
    }

    Mesh2D* Mesh2D::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMesh2D::type_name_static() ;
        }
        Mesh2D* mesh = Mesh2DFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "Mesh2D" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "Mesh2D" ) << "Falling back to GeogramMesh2D data structure"
                << std::endl ;

            mesh = new GeogramMesh2D ;
        }
        return mesh ;
    }

    void Mesh2D::next_on_border(
        index_t f,
        index_t e,
        index_t& next_f,
        index_t& next_e) const
    {
        ringmesh_assert(e < nb_facet_vertices(f));
        ringmesh_assert(is_edge_on_border(f, e));

        // Global indices in the surfaces
        index_t next_v_id = facet_vertex(f, next_facet_vertex(f, e));

        // Get the facets around the shared vertex (next_v_id) that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)
        std::vector< index_t > facets_around_next_v_id;
        index_t nb_around = facets_around_vertex(next_v_id, facets_around_next_v_id,
                                                 true, f);
        ringmesh_assert(nb_around == 1 || nb_around == 2);

        next_f = facets_around_next_v_id[0];

        if (nb_around == 2) {
            if (next_f == f) {
                next_f = facets_around_next_v_id[1];
            }
            ringmesh_assert(next_f != NO_ID);
            ringmesh_assert(is_facet_on_border(next_f));

            // Local index of next vertex in the next facet
            next_e = vertex_index_in_facet(next_f, next_v_id);
            ringmesh_assert(is_edge_on_border(next_f, next_e));
        }
        else if (nb_around == 1) {
            // next_v_id must be in two border edges of facet f
            next_e = vertex_index_in_facet(next_f, next_v_id);
            ringmesh_assert(is_edge_on_border(next_f, next_e));
        }
    }

    void Mesh2D::prev_on_border(
        index_t f,
        index_t e,
        index_t& prev_f,
        index_t& prev_e) const
    {
        ringmesh_assert(e < nb_facet_vertices(f));
        ringmesh_assert(is_edge_on_border(f, e));

        // Global indices in the surfaces
        index_t v_id = facet_vertex(f, e);

        // Get the facets around the shared vertex (v_id) that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)
        std::vector< index_t > facets_around_v_id;
        index_t nb_around = facets_around_vertex(v_id, facets_around_v_id, true,
                                                 f);
        ringmesh_assert(nb_around == 1 || nb_around == 2);

        prev_f = facets_around_v_id[0];

        if (nb_around == 2) {
            if (prev_f == f) {
                prev_f = facets_around_v_id[1];
            }
            ringmesh_assert(prev_f != NO_ID);
            ringmesh_assert(is_facet_on_border(prev_f));

            // Local index of given vertex in the prev facet
            index_t v_in_prev_f = vertex_index_in_facet(prev_f, v_id);
            // Local index of previous vertex in the prev facet
            prev_e = prev_facet_vertex(prev_f, v_in_prev_f);
            ringmesh_assert(is_edge_on_border(prev_f, prev_e));
        }
        else if (nb_around == 1) {
            // v_id must be in two border edges of facet f
            index_t v_in_next_facet = vertex_index_in_facet(prev_f, v_id);
            prev_e = prev_facet_vertex(prev_f, v_in_next_facet);
            ringmesh_assert(is_edge_on_border(prev_f, prev_e));
        }
    }

    index_t Mesh2D::facet_from_vertex_ids(index_t in0, index_t in1) const
    {
        ringmesh_assert(
            in0 < nb_vertices() && in1 < nb_vertices());

        // Another possible, probably faster, algorithm is to check if the 2 indices
        // are neighbors in facets_ and check that they are in the same facet

        // Check if the edge is in one of the facet
        for (index_t f = 0; f < nb_facets(); ++f) {
            bool found = false;
            index_t prev = facet_vertex(f,
                nb_facet_vertices(f) - 1);
            for (index_t v = 0; v < nb_facet_vertices(f); ++v) {
                index_t p = facet_vertex(f, v);
                if ((prev == in0 && p == in1) || (prev == in1 && p == in0)) {
                    found = true;
                    break;
                }
                prev = p;
            }
            if (found) {
                return f;
            }
        }
        return NO_ID;
    }

    index_t Mesh2D::vertex_index_in_facet(
        index_t facet_index,
        index_t vertex_id) const
    {
        ringmesh_assert(facet_index < nb_facets());
        for (index_t v = 0; v < nb_facet_vertices(facet_index); v++) {
            if (facet_vertex(facet_index, v)
                == vertex_id) {
                return v;
            }
        }
        return NO_ID;
    }

    index_t Mesh2D::closest_vertex_in_facet(index_t f, const vec3& v) const
    {
        index_t result = 0;
        double dist = DBL_MAX;
        for (index_t p = 0; p < nb_facet_vertices(f); p++) {
            double distance = length2(v - vertex(facet_vertex(f, p)));
            if (dist > distance) {
                dist = distance;
                result = p;
            }
        }
        return result;
    }
    index_t Mesh2D::facets_around_vertex(
        index_t surf_vertex_id,
        std::vector< index_t >& result,
        bool border_only,
        index_t f0) const
    {
        result.clear();

        index_t f = 0;
        while (f0 == NO_ID && f < nb_facets()) {
            for (index_t lv = 0; lv < nb_facet_vertices(f); lv++) {
                if (facet_vertex(f, lv) == surf_vertex_id) {
                    f0 = f;
                    break;
                }
            }
            f++;
        }

        ringmesh_assert(f0 != NO_ID);

        // Flag the visited facets
        std::vector< index_t > visited;
        visited.reserve(10);

        // Stack of the adjacent facets
        std::stack< index_t > S;
        S.push(f0);
        visited.push_back(f0);

        do {
            index_t f = S.top();
            S.pop();

            for (index_t v = 0; v < nb_facet_vertices(f); ++v) {
                if (facet_vertex(f, v) == surf_vertex_id) {
                    index_t adj_P = facet_adjacent(f, v);
                    index_t prev = prev_facet_vertex(f, v);
                    index_t adj_prev = facet_adjacent(f, prev);

                    if (adj_P != NO_ID) {
                        // The edge starting at P is not on the boundary
                        if (!contains(visited, adj_P)) {
                            S.push(adj_P);
                            visited.push_back(adj_P);
                        }
                    }
                    if (adj_prev != NO_ID) {
                        // The edge ending at P is not on the boundary
                        if (!contains(visited, adj_prev)) {
                            S.push(adj_prev);
                            visited.push_back(adj_prev);
                        }
                    }

                    if (border_only) {
                        if (adj_P == NO_ID || adj_prev == NO_ID) {
                            result.push_back(f);
                        }
                    }
                    else {
                        result.push_back(f);
                    }

                    // We are done with this facet
                    break;
                }
            }
        } while (!S.empty());

        return static_cast< index_t >(result.size());
    }
    Mesh3D* Mesh3D::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMesh3D::type_name_static() ;
        }
        Mesh3D* mesh = Mesh3DFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "Mesh3D" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "Mesh3D" ) << "Falling back to GeogramMesh3D data structure"
                << std::endl ;

            mesh = new GeogramMesh3D ;
        }
        return mesh ;
    }

    MeshAllD* MeshAllD::create_mesh( const MeshType type )
    {
        MeshType new_type = type ;
        if( new_type.empty() ) {
            new_type = GeogramMeshAllD::type_name_static() ;
        }
        MeshAllD* mesh = MeshAllDFactory::create_object( new_type ) ;
        if( !mesh ) {
            Logger::warn( "MeshAllD" )
                << "Could not create mesh data structure: " << new_type
                << std::endl ;
            Logger::warn( "MeshAllD" ) << "Falling back to GeogramMeshAllD data structure"
                << std::endl ;

            mesh = new GeogramMeshAllD ;
        }
        return mesh ;
    }


} // namespace
