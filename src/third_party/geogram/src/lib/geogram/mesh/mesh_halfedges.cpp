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

#include <geogram/mesh/mesh_halfedges.h>

namespace GEO {

    bool MeshHalfedges::move_to_next_around_vertex(Halfedge& H) const {
        geo_debug_assert(halfedge_is_valid(H));
        index_t v = mesh_.corner_vertex_index(H.corner);
        signed_index_t sf = mesh_.corner_adjacent_facet(H.corner);
        if(sf == -1) {
            return false;
        }
        index_t f = index_t(sf);
        if(
            use_facet_region_ &&
            mesh_.facet_region(H.facet) != mesh_.facet_region(f)
        ) {
            return false;
        }
        for(
            index_t c = mesh_.facet_begin(f);
            c < mesh_.facet_end(f); c++
        ) {
            index_t pc = mesh_.prev_around_facet(f, c);
            if(
                mesh_.corner_vertex_index(c) == v &&
                mesh_.corner_adjacent_facet(pc) == signed_index_t(H.facet)
            ) {
                H.corner = c;
                H.facet = f;
                return true;
            }
        }
        geo_assert_not_reached;
        return false;
    }

    bool MeshHalfedges::move_to_prev_around_vertex(Halfedge& H) const {
        geo_debug_assert(halfedge_is_valid(H));
        index_t v = mesh_.corner_vertex_index(H.corner);
        index_t pc = mesh_.prev_around_facet(H.facet, H.corner);
        signed_index_t sf = mesh_.corner_adjacent_facet(pc);
        if(sf == -1) {
            return false;
        }
        index_t f = index_t(sf);
        if(
            use_facet_region_ &&
            mesh_.facet_region(H.facet) != mesh_.facet_region(f)
        ) {
            return false;
        }
        for(index_t c = mesh_.facet_begin(f); c < mesh_.facet_end(f); c++) {
            if(
                mesh_.corner_vertex_index(c) == v &&
                mesh_.corner_adjacent_facet(c) == signed_index_t(H.facet)
            ) {
                H.corner = c;
                H.facet = f;
                return true;
            }
        }
        geo_assert_not_reached;
        return false;
    }

    void MeshHalfedges::move_to_next_around_border(Halfedge& H) const {
        geo_debug_assert(halfedge_is_valid(H));
        geo_debug_assert(halfedge_is_border(H));
        move_to_next_around_facet(H);
        index_t count = 0;
        while(move_to_next_around_vertex(H)) {
            ++count;
            geo_assert(count < 10000);
        }
    }

    void MeshHalfedges::move_to_opposite(Halfedge& H) const {
        geo_debug_assert(halfedge_is_valid(H));
        index_t v = mesh_.corner_vertex_index(
            mesh_.next_around_facet(H.facet, H.corner)
        );
        signed_index_t sf = mesh_.corner_adjacent_facet(H.corner);
        geo_assert(sf != -1);
        index_t f = index_t(sf);
        for(index_t c = mesh_.facet_begin(f); c != mesh_.facet_end(f); ++c) {
            if(mesh_.corner_vertex_index(c) == v) {
                H.facet = f;
                H.corner = c;
                return;
            }
        }
        geo_assert_not_reached;
    }
}

