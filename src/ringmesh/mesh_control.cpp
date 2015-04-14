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
 *     Ecole Nationale Sup�rieure de G�ologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/mesh_control.h>
#include <geogram/mesh/mesh_intersection.cpp>
#include <utility>      // std::swap

namespace {
void get_tet_bbox(const Mesh& M, Box& B, index_t t) {
	const double* p = M.vertices.point_ptr(M.cells.vertex(t, 0));
	for (coord_index_t coord = 0; coord < 3; ++coord) {
		B.xyz_min[coord] = p[coord];
		B.xyz_max[coord] = p[coord];
	}
	for (index_t lv = 1; lv < 4; ++lv) {
		p = M.vertices.point_ptr(M.cells.vertex(t, lv));
		for (coord_index_t coord = 0; coord < 3; ++coord) {
			B.xyz_min[coord] = geo_min(B.xyz_min[coord], p[coord]);
			B.xyz_max[coord] = geo_max(B.xyz_max[coord], p[coord]);
		}
	}
}
}

namespace RINGMesh {
DetectInter::DetectInter(const MacroMesh& mm) :
		inter_(mm.cells.nb_tet()), cur_reg_(0), cur_reg2_(0), nb_inter_(0), cur_cell_(
				0), mm_(mm) {
}

DetectInter::~DetectInter() {
}

void DetectInter::operator()(index_t idx) {
	mm_.vertices.vertex(idx);
//	mm_.cells.tet_id();

//	idx = mm_.cells.cell_index_in_mesh(idx, cur_reg_);
	if (idx > cur_cell_ || cur_reg_ != cur_reg2_) {
//		std::cout << "cur_cell " << cur_cell_ << " idx " << idx << std::endl;

// TODO simpler
//		const vec3& v0 = mm_.vertices.vertex(idx);
		vec3& v1 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 0));
		vec3& v2 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 1));
		vec3& v3 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 2));
		vec3& v4 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 2, 2));

		vec3& v1_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 0));
		vec3& v2_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 1));
		vec3& v3_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 2));
		vec3& v4_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 2, 2));

		if ((Math::point_inside_tetra(v1_2, v1, v2, v3, v4) && v1_2 != v1
				&& v1_2 != v2 && v1_2 != v3 && v1_2 != v4)
				|| (Math::point_inside_tetra(v2_2, v1, v2, v3, v4)
						&& v2_2 != v1 && v2_2 != v2 && v2_2 != v3 && v2_2 != v4)
				|| (Math::point_inside_tetra(v3_2, v1, v2, v3, v4)
						&& v3_2 != v1 && v3_2 != v2 && v3_2 != v3 && v3_2 != v4)
				|| (Math::point_inside_tetra(v4_2, v1, v2, v3, v4)
						&& v4_2 != v1 && v4_2 != v2 && v4_2 != v3 && v4_2 != v4)) {
			++nb_inter_;
			std::cout << "intersection region:cell " << cur_reg_ << ":"
					<< cur_cell_ << " " << cur_reg2_ << ":" << idx << std::endl;
		}

//		if ((Utils::point_inside_tetra(v1_2, v1, v2, v3, v4))
//				&& cur_reg_ == 1) {
//			std::cout << "intersection 1 region:cell " << cur_reg_ << ":"
//					<< cur_cell_ << " " << cur_reg2_ << ":" << idx << std::endl;
//		}
//		if ((Utils::point_inside_tetra(v2_2, v1, v2, v3, v4))
//				&& cur_reg_ == 1) {
//			std::cout << "intersection 2 region:cell " << cur_reg_ << ":"
//					<< cur_cell_ << " " << cur_reg2_ << ":" << idx << std::endl;
//		}
//		if ((Utils::point_inside_tetra(v3_2, v1, v2, v3, v4))
//				&& cur_reg_ == 1) {
//			std::cout << "intersection 3 region:cell " << cur_reg_ << ":"
//					<< cur_cell_ << " " << cur_reg2_ << ":" << idx << std::endl;
//		}
//		if ((Utils::point_inside_tetra(v4_2, v1, v2, v3, v4))
//				&& cur_reg_ == 1) {
//			std::cout << "intersection 4 region:cell " << cur_reg_ << ":"
//					<< cur_cell_ << " " << cur_reg2_ << ":" << idx << std::endl;
//		}
	}
}

index_t DetectInter::detect_mesh_intersection() {
//	for (index_t c = 0; c < mm_.cells.nb_cells(); c++) {
//	std::cout << "nb mesh " << mm_.nb_meshes() << std::endl;
	for (index_t m = 0; m < mm_.nb_meshes(); m++) {
		cur_reg_ = m;
		for (index_t c = 0; c < mm_.mesh(m).cells.nb(); c++) {
			cur_cell_ = c;

			// Get the bbox of c
			Box box;
			get_tet_bbox(mm_.mesh(m), box, c);

			for (index_t reg_idx = cur_reg_; reg_idx < mm_.nb_meshes();
					reg_idx++) {
				cur_reg2_ = reg_idx;
				mm_.tools.tet_aabb(cur_reg_).compute_bbox_cell_bbox_intersections(
						box, *this);
			}
		}
	}

	return nb_inter_;
}

}

//				index_t v1_index = mm_.mesh(cur_reg_).cells.facet_vertex(
//						cur_cell_, 0, 0);
//				index_t v2_index = mm_.mesh(cur_reg_).cells.facet_vertex(
//						cur_cell_, 0, 1);
//				index_t v3_index = mm_.mesh(cur_reg_).cells.facet_vertex(
//						cur_cell_, 0, 2);
//				index_t v4_index = mm_.mesh(cur_reg_).cells.facet_vertex(
//						cur_cell_, 2, 2);
//
//				vec3 v1 = mm_.mesh(cur_reg_).vertices.point(v1_index);
//				vec3 v2 = mm_.mesh(cur_reg_).vertices.point(v2_index);
//				vec3 v3 = mm_.mesh(cur_reg_).vertices.point(v3_index);
//				vec3 v4 = mm_.mesh(cur_reg_).vertices.point(v4_index);
//
//				double V[4][3] = { { v1.x, v1.y, v1.z }, { v2.x, v2.y, v2.z }, {
//						v3.x, v3.y, v3.z }, { v4.x, v4.y, v4.z } };
//
//				std::cout << "coord " << cur_reg_ << " " << cur_cell_ << " "
//						<< v1.x << " " << v1.y << " " << v1.z << "     " << v2.x
//						<< " " << v2.y << " " << v2.z << "     " << v3.x << " "
//						<< v3.y << " " << v3.z << "     " << v4.x << " " << v4.y
//						<< " " << v4.z << std::endl;

//*************************************
//**********  In a region   ***********
//	std::cout << "x "
//			<< mm_.mesh(0).vertices.point(
//					mm_.mesh(0).cells.facet_vertex(9, 0, 0)).x << " y "
//			<< mm_.mesh(0).vertices.point(
//					mm_.mesh(0).cells.facet_vertex(9, 0, 0)).y << " z "
//			<< mm_.mesh(0).vertices.point(
//					mm_.mesh(0).cells.facet_vertex(9, 0, 0)).z << std::endl;
//	mm_.mesh(0).vertices.point(mm_.mesh(0).cells.facet_vertex(9, 0, 0)).x = 3;
//	mm_.mesh(0).vertices.point(mm_.mesh(0).cells.facet_vertex(9, 0, 0)).y = 3;
//	mm_.mesh(0).vertices.point(mm_.mesh(0).cells.facet_vertex(9, 0, 0)).z = -3;
//*************************************

//*************************************
//*********  Between regions   ********
//	std::cout << "x "
//			<< mm_.mesh(1).vertices.point(
//					mm_.mesh(1).cells.facet_vertex(0, 0, 0)).x << " y "
//			<< mm_.mesh(1).vertices.point(
//					mm_.mesh(1).cells.facet_vertex(0, 0, 0)).y << " z "
//			<< mm_.mesh(1).vertices.point(
//					mm_.mesh(1).cells.facet_vertex(0, 0, 0)).z << std::endl;
//	mm_.mesh(1).vertices.point(mm_.mesh(1).cells.facet_vertex(0, 0, 0)).x = 9;
//	mm_.mesh(1).vertices.point(mm_.mesh(1).cells.facet_vertex(0, 0, 0)).y = 9;
//	mm_.mesh(1).vertices.point(mm_.mesh(1).cells.facet_vertex(0, 0, 0)).z = -9;
//*************************************
