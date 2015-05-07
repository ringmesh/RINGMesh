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
#include <geogram/mesh/mesh_AABB.cpp>

namespace RINGMesh {
DetectInter::DetectInter(const MacroMesh& mm) :
		cur_reg_(0), cur_reg2_(0), nb_inter_(0), cur_cell_(0), mm_(mm), inter_(
				mm_.cells.nb_cells(), false) {
}

DetectInter::~DetectInter() {
}

void DetectInter::operator()(index_t idx) {

	if (idx > cur_cell_ || cur_reg_ != cur_reg2_) {

		vec3& v1 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 0));
		vec3& v2 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 1));
		vec3& v3 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 2));
		vec3& v4 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 3));

		vec3& v1_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.vertex(idx, 0));
		vec3& v2_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.vertex(idx, 1));
		vec3& v3_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.vertex(idx, 2));
		vec3& v4_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.vertex(idx, 3));

//		// tetrahedron
		if (mm_.mesh(cur_reg_).cells.nb_vertices(cur_cell_) == 4) {
//			if (cur_reg_ == cur_reg2_) {
			if ((Math::point_inside_tetra(v1_2, v1, v2, v3, v4) && v1_2 != v1
					&& v1_2 != v2 && v1_2 != v3 && v1_2 != v4)
					|| (Math::point_inside_tetra(v2_2, v1, v2, v3, v4)
							&& v2_2 != v1 && v2_2 != v2 && v2_2 != v3
							&& v2_2 != v4)
					|| (Math::point_inside_tetra(v3_2, v1, v2, v3, v4)
							&& v3_2 != v1 && v3_2 != v2 && v3_2 != v3
							&& v3_2 != v4)
					|| (Math::point_inside_tetra(v4_2, v1, v2, v3, v4)
							&& v4_2 != v1 && v4_2 != v2 && v4_2 != v3
							&& v4_2 != v4)) {
				++nb_inter_;
				std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
						<< " contains a point of the cell " << cur_reg2_ << ":"
						<< idx << " region:cell" << std::endl;
				fill_list_intersection(idx);
				return;
			}
			// if the cell idx is a pyramid
			if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 4) {
				vec3& v5_2 = mm_.mesh(cur_reg2_).vertices.point(
						mm_.mesh(cur_reg2_).cells.vertex(idx, 4));
				if (Math::point_inside_tetra(v5_2, v1, v2, v3, v4) && v5_2 != v1
						&& v5_2 != v2 && v5_2 != v3 && v5_2 != v4) {
					++nb_inter_;
					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
							<< " contains a point of the cell " << cur_reg2_
							<< ":" << idx << " region:cell" << std::endl;
					fill_list_intersection(idx);
					return;
				}
			}
			// if the cell idx is a prism
			if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 5) {
				vec3& v6_2 = mm_.mesh(cur_reg2_).vertices.point(
						mm_.mesh(cur_reg2_).cells.vertex(idx, 5));
				if (Math::point_inside_tetra(v6_2, v1, v2, v3, v4) && v6_2 != v1
						&& v6_2 != v2 && v6_2 != v3 && v6_2 != v4) {
					++nb_inter_;
					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
							<< " contains a point of the cell " << cur_reg2_
							<< ":" << idx << " region:cell" << std::endl;
					fill_list_intersection(idx);
					return;
				}
			}
			// if the cell idx is a hexahedron
			if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) == 8) {
				vec3& v7_2 = mm_.mesh(cur_reg2_).vertices.point(
						mm_.mesh(cur_reg2_).cells.vertex(idx, 6));
				vec3& v8_2 = mm_.mesh(cur_reg2_).vertices.point(
						mm_.mesh(cur_reg2_).cells.vertex(idx, 7));
				if ((Math::point_inside_tetra(v7_2, v1, v2, v3, v4)
						&& v7_2 != v1 && v7_2 != v2 && v7_2 != v3 && v7_2 != v4)
						|| (Math::point_inside_tetra(v8_2, v1, v2, v3, v4)
								&& v8_2 != v1 && v8_2 != v2 && v8_2 != v3
								&& v8_2 != v4)) {
					++nb_inter_;
					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
							<< " contains a point of the cell " << cur_reg2_
							<< ":" << idx << " region:cell" << std::endl;
					fill_list_intersection(idx);
					return;
				}
			}
//			}

//			if (cur_reg_ != cur_reg2_) {
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 3) {
//					if (Math::point_inside_tetra(v1_2, v1, v2, v3, v4)
//							|| Math::point_inside_tetra(v2_2, v1, v2, v3, v4)
//							|| Math::point_inside_tetra(v3_2, v1, v2, v3, v4)
//							|| Math::point_inside_tetra(v4_2, v1, v2, v3, v4)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//						fill_list_intersection(idx);
//						return;
//					}
////					else if (v1_2 == v1 || v1_2 == v2 || v1_2 == v3
////							|| v1_2 == v4 || v2_2 == v1 || v2_2 == v2
////							|| v2_2 == v3 || v2_2 == v4 || v3_2 == v1
////							|| v3_2 == v2 || v3_2 == v3 || v3_2 == v4
////							|| v4_2 == v1 || v4_2 == v2 || v4_2 == v3
////							|| v4_2 == v4) {
////						// if the point of two different regions are at the same place
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////					fill_list_intersection(idx);
////						return;
////					}
//				}
//				// if the cell idx is a pyramid
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 4) {
//					vec3& v5_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 4));
//					if (Math::point_inside_tetra(v5_2, v1, v2, v3, v4)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//						fill_list_intersection(idx);
//						return;
//					}
////					else if (v5_2 == v1 || v5_2 == v2 || v5_2 == v3
////							|| v5_2 == v4) {
////						// if the point of two different regions are at the same place
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////					fill_list_intersection(idx);
////						return;
////					}
//					// if the cell idx is a prism
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 5) {
//					vec3& v6_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 5));
//					if (Math::point_inside_tetra(v6_2, v1, v2, v3, v4)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//						fill_list_intersection(idx);
//						return;
//					}
////					else if (v6_2 == v1 || v6_2 == v2 || v6_2 == v3
////							|| v6_2 == v4) {
////						// if the point of two different regions are at the same place
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////					fill_list_intersection(idx);
////						return;
////					}
//				}
//				// if the cell idx is a hexahedron
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) == 8) {
//					vec3& v7_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 6));
//					vec3& v8_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 7));
//					if (Math::point_inside_tetra(v7_2, v1, v2, v3, v4)
//							|| Math::point_inside_tetra(v8_2, v1, v2, v3, v4)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//						fill_list_intersection(idx);
//						return;
//					}
////					else if (v7_2 == v1 || v7_2 == v2 || v7_2 == v3
////							|| v7_2 == v4 || v8_2 == v1 || v8_2 == v2
////							|| v8_2 == v3 || v8_2 == v4) {
////						// if the point of two different regions are at the same place
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////					fill_list_intersection(idx);
////						return;
////					}
//				}
//			}
		}

//		// pyramid
//		if (mm_.mesh(cur_reg_).cells.nb_vertices(cur_cell_) == 5) {
//			vec3& v5 = mm_.mesh(cur_reg_).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 4));
//			if (cur_reg_ == cur_reg2_) {
//				// if the cell idx is a tetrahedron
//				if ((Math::point_inside_pyramid(v1_2, v1, v2, v3, v4, v5)
//						&& v1_2 != v1 && v1_2 != v2 && v1_2 != v3 && v1_2 != v4
//						&& v1_2 != v5)
//						|| (Math::point_inside_pyramid(v2_2, v1, v2, v3, v4, v5)
//								&& v2_2 != v1 && v2_2 != v2 && v2_2 != v3
//								&& v2_2 != v4 && v2_2 != v5)
//						|| (Math::point_inside_pyramid(v3_2, v1, v2, v3, v4, v5)
//								&& v3_2 != v1 && v3_2 != v2 && v3_2 != v3
//								&& v3_2 != v4 && v3_2 != v5)
//						|| (Math::point_inside_pyramid(v4_2, v1, v2, v3, v4, v5)
//								&& v4_2 != v1 && v4_2 != v2 && v4_2 != v3
//								&& v4_2 != v4 && v4_2 != v5)) {
//					++nb_inter_;
//					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//							<< " contains a point of the cell " << cur_reg2_
//							<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//					return;
//				}
//				// if the cell idx is a pyramid
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 4) {
//					vec3& v5_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 4));
//					if (Math::point_inside_pyramid(v5_2, v1, v2, v3, v4, v5)
//							&& v5_2 != v1 && v5_2 != v2 && v5_2 != v3
//							&& v5_2 != v4 && v5_2 != v5) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
//				}
//				// if the cell idx is a prism
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 5) {
//					vec3& v6_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 5));
//					if (Math::point_inside_pyramid(v6_2, v1, v2, v3, v4, v5)
//							&& v6_2 != v1 && v6_2 != v2 && v6_2 != v3
//							&& v6_2 != v4 && v6_2 != v5) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
//				}
//				// if the cell idx is a hexahedron
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) == 8) {
//					vec3& v7_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 6));
//					vec3& v8_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 7));
//					if ((Math::point_inside_pyramid(v7_2, v1, v2, v3, v4, v5)
//							&& v7_2 != v1 && v7_2 != v2 && v7_2 != v3
//							&& v7_2 != v4 && v7_2 != v5)
//							|| (Math::point_inside_pyramid(v8_2, v1, v2, v3, v4,
//									v5) && v8_2 != v1 && v8_2 != v2
//									&& v8_2 != v3 && v8_2 != v4 && v8_2 != v5)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
//				}
//			}
//
//			if (cur_reg_ != cur_reg2_) {
//				if (Math::point_inside_pyramid(v1_2, v1, v2, v3, v4, v5)
//						|| Math::point_inside_pyramid(v2_2, v1, v2, v3, v4, v5)
//						|| Math::point_inside_pyramid(v3_2, v1, v2, v3, v4, v5)
//						|| Math::point_inside_pyramid(v4_2, v1, v2, v3, v4,
//								v5)) {
//					++nb_inter_;
//					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//							<< " contains a point of the cell " << cur_reg2_
//							<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//					return;
//				}
////		else if (v1_2 == v1 || v1_2 == v2 || v1_2 == v3 || v1_2 == v4
////						|| v1_2 == v5 || v2_2 == v1 || v2_2 == v2 || v2_2 == v3
////						|| v2_2 == v4 || v2_2 == v5 || v3_2 == v1 || v3_2 == v2
////						|| v3_2 == v3 || v3_2 == v4 || v3_2 == v5 || v4_2 == v1
////						|| v4_2 == v2 || v4_2 == v3 || v4_2 == v4
////						|| v4_2 == v5) {
////					// if the point of two different regions are at the same place
////					++nb_inter_;
////					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////							<< " contains a point of the cell " << cur_reg2_
////							<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////					return;
////				}
//				// if the cell idx is a pyramid
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 4) {
//					vec3& v5_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 4));
//					if (Math::point_inside_pyramid(v5_2, v1, v2, v3, v4, v5)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v5_2 == v1 || v5_2 == v2 || v5_2 == v3
////							|| v5_2 == v4 || v5_2 == v5) {
////						// if the point of two different regions are at the same place
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//				// if the cell idx is a prism
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 5) {
//					vec3& v6_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 5));
//					if (Math::point_inside_pyramid(v6_2, v1, v2, v3, v4, v5)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v6_2 == v1 || v6_2 == v2 || v6_2 == v3
////							|| v6_2 == v4 || v6_2 == v5) {
////						// if the point of two different regions are at the same place
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//				// if the cell idx is a hexahedron
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) == 8) {
//					vec3& v7_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 6));
//					vec3& v8_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 7));
//					if (Math::point_inside_pyramid(v7_2, v1, v2, v3, v4, v5)
//							|| Math::point_inside_pyramid(v8_2, v1, v2, v3, v4,
//									v5)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v7_2 == v1 || v7_2 == v2 || v7_2 == v3
////							|| v7_2 == v4 || v7_2 == v5 || v8_2 == v1
////							|| v8_2 == v2 || v8_2 == v3 || v8_2 == v4
////							|| v8_2 == v5) {
////						// if the point of two different regions are at the same place
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//			}
//		}

//		// prism
//		if (mm_.mesh(cur_reg_).cells.nb_vertices(cur_cell_) == 6) {
//			vec3& v5 = mm_.mesh(cur_reg_).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 4));
//			vec3& v6 = mm_.mesh(cur_reg_).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 5));
//
//			if (cur_reg_ == cur_reg2_) {
//				if ((mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 3)
//						&& ((Math::point_inside_prism(v1_2, v1, v2, v3, v4, v5,
//								v6) && v1_2 != v1 && v1_2 != v2 && v1_2 != v3
//								&& v1_2 != v4 && v1_2 != v5 && v1_2 != v6)
//								|| (Math::point_inside_prism(v2_2, v1, v2, v3,
//										v4, v5, v6) && v2_2 != v1 && v2_2 != v2
//										&& v2_2 != v3 && v2_2 != v4
//										&& v2_2 != v5 && v2_2 != v6)
//								|| (Math::point_inside_prism(v3_2, v1, v2, v3,
//										v4, v5, v6) && v3_2 != v1 && v3_2 != v2
//										&& v3_2 != v3 && v3_2 != v4
//										&& v3_2 != v5 && v3_2 != v6)
//								|| (Math::point_inside_prism(v4_2, v1, v2, v3,
//										v4, v5, v6) && v4_2 != v1 && v4_2 != v2
//										&& v4_2 != v3 && v4_2 != v4
//										&& v4_2 != v5 && v4_2 != v6))) {
//					++nb_inter_;
//					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//							<< " contains a point of the cell " << cur_reg2_
//							<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//					return;
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 4) {
//					vec3& v5_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 4));
//					if (Math::point_inside_prism(v5_2, v1, v2, v3, v4, v5, v6)
//							&& v5_2 != v1 && v5_2 != v2 && v5_2 != v3
//							&& v5_2 != v4 && v5_2 != v5 && v5_2 != v6) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 5) {
//					vec3& v6_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 5));
//					if (Math::point_inside_prism(v6_2, v1, v2, v3, v4, v5, v6)
//							&& v6_2 != v1 && v6_2 != v2 && v6_2 != v3
//							&& v6_2 != v4 && v6_2 != v5 && v6_2 != v6) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) == 8) {
//					vec3& v7_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 6));
//					vec3& v8_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 7));
//					if ((Math::point_inside_prism(v7_2, v1, v2, v3, v4, v5, v6)
//							&& v7_2 != v1 && v7_2 != v2 && v7_2 != v3
//							&& v7_2 != v4 && v7_2 != v5 && v7_2 != v6)
//							&& (Math::point_inside_prism(v8_2, v1, v2, v3, v4,
//									v5, v6) && v8_2 != v1 && v8_2 != v2
//									&& v8_2 != v3 && v8_2 != v4 && v8_2 != v5
//									&& v8_2 != v6)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
//				}
//			}
//
//			if (cur_reg_ != cur_reg2_) {
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 3) {
//					if (Math::point_inside_prism(v1_2, v1, v2, v3, v4, v5, v6)
//							|| Math::point_inside_prism(v2_2, v1, v2, v3, v4,
//									v5, v6)
//							|| Math::point_inside_prism(v3_2, v1, v2, v3, v4,
//									v5, v6)
//							|| Math::point_inside_prism(v4_2, v1, v2, v3, v4,
//									v5, v6)) {
//
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v1_2 == v1 || v1_2 == v2 || v1_2 == v3
////							|| v1_2 == v4 || v1_2 == v5 || v1_2 == v6
////							|| v2_2 == v1 || v2_2 == v2 || v2_2 == v3
////							|| v2_2 == v4 || v2_2 == v5 || v2_2 == v6
////							|| v3_2 == v1 || v3_2 == v2 || v3_2 == v3
////							|| v3_2 == v4 || v3_2 == v5 || v3_2 == v6
////							|| v4_2 == v1 || v4_2 == v2 || v4_2 == v3
////							|| v4_2 == v4 || v4_2 == v5 || v4_2 == v6) {
////						// if a point of one region is at the same place of the point of another region
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 4) {
//					vec3& v5_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 4));
//					if (Math::point_inside_prism(v5_2, v1, v2, v3, v4, v5,
//							v6)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v5_2 == v1 || v5_2 == v2 || v5_2 == v3
////							|| v5_2 == v4 || v5_2 == v5 || v5_2 == v6) {
////						// if a point of one region is at the same place of the point of another region
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 5) {
//					vec3& v6_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 5));
//					if (Math::point_inside_prism(v6_2, v1, v2, v3, v4, v5,
//							v6)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v6_2 == v1 || v6_2 == v2 || v6_2 == v3
////							|| v6_2 == v4 || v6_2 == v5 || v6_2 == v6) {
////						// if a point of one region is at the same place of the point of another region
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) == 8) {
//					vec3& v7_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 6));
//					vec3& v8_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 7));
//					if (Math::point_inside_prism(v7_2, v1, v2, v3, v4, v5, v6)
//							|| Math::point_inside_prism(v8_2, v1, v2, v3, v4,
//									v5, v6)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v7_2 == v1 || v7_2 == v2 || v7_2 == v3
////							|| v7_2 == v4 || v7_2 == v5 || v7_2 == v6
////							|| v8_2 == v1 || v8_2 == v2 || v8_2 == v3
////							|| v8_2 == v4 || v8_2 == v5 || v8_2 == v6) {
////						// if a point of one region is at the same place of the point of another region
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//			}
//		}

//		// hexahedron
//		if (mm_.mesh(cur_reg_).cells.nb_vertices(cur_cell_) == 8) {
//			vec3& v5 = mm_.mesh(cur_reg_).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 4));
//			vec3& v6 = mm_.mesh(cur_reg_).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 5));
//			vec3& v7 = mm_.mesh(cur_reg_).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 6));
//			vec3& v8 = mm_.mesh(cur_reg_).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(cur_cell_, 7));
//
//			if (cur_reg_ == cur_reg2_) {
//				if ((Math::point_inside_hexa(v1_2, v1, v2, v3, v4, v5, v6, v7,
//						v8) && v1_2 != v1 && v1_2 != v2 && v1_2 != v3
//						&& v1_2 != v4 && v1_2 != v5 && v1_2 != v6 && v1_2 != v7
//						&& v1_2 != v8)
//						|| (Math::point_inside_hexa(v2_2, v1, v2, v3, v4, v5,
//								v6, v7, v8) && v2_2 != v1 && v2_2 != v2
//								&& v2_2 != v3 && v2_2 != v4 && v2_2 != v5
//								&& v2_2 != v6 && v2_2 != v7 && v2_2 != v8)
//						|| (Math::point_inside_hexa(v3_2, v1, v2, v3, v4, v5,
//								v6, v7, v8) && v3_2 != v1 && v3_2 != v2
//								&& v3_2 != v3 && v3_2 != v4 && v3_2 != v5
//								&& v3_2 != v6 && v3_2 != v7 && v3_2 != v8)
//						|| (Math::point_inside_hexa(v4_2, v1, v2, v3, v4, v5,
//								v6, v7, v8) && v4_2 != v1 && v4_2 != v2
//								&& v4_2 != v3 && v4_2 != v4 && v4_2 != v5
//								&& v4_2 != v6 && v4_2 != v7 && v4_2 != v8)) {
//					++nb_inter_;
//					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//							<< " contains a point of the cell " << cur_reg2_
//							<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//					return;
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 4) {
//					vec3& v5_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 4));
//					if (Math::point_inside_hexa(v5_2, v1, v2, v3, v4, v5, v6,
//							v7, v8) && v5_2 != v1 && v5_2 != v2 && v5_2 != v3
//							&& v5_2 != v4 && v5_2 != v5 && v5_2 != v6
//							&& v5_2 != v7 && v5_2 != v8) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 5) {
//					vec3& v6_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 5));
//					if (Math::point_inside_hexa(v6_2, v1, v2, v3, v4, v5, v6,
//							v7, v8) && v6_2 != v1 && v6_2 != v2 && v6_2 != v3
//							&& v6_2 != v4 && v6_2 != v5 && v6_2 != v6
//							&& v6_2 != v7 && v6_2 != v8) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) == 8) {
//					vec3& v7_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 6));
//					vec3& v8_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 7));
//					if ((Math::point_inside_hexa(v7_2, v1, v2, v3, v4, v5, v6,
//							v7, v8) && v7_2 != v1 && v7_2 != v2 && v7_2 != v3
//							&& v7_2 != v4 && v7_2 != v5 && v7_2 != v6
//							&& v7_2 != v7 && v7_2 != v8)
//							|| (Math::point_inside_hexa(v8_2, v1, v2, v3, v4,
//									v5, v6, v7, v8) && v8_2 != v1 && v8_2 != v2
//									&& v8_2 != v3 && v8_2 != v4 && v8_2 != v5
//									&& v8_2 != v6 && v8_2 != v7 && v8_2 != v8)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
//				}
//			}
//
//			if (cur_reg_ != cur_reg2_) {
//				if (Math::point_inside_hexa(v1_2, v1, v2, v3, v4, v5, v6, v7,
//						v8)
//						|| Math::point_inside_hexa(v2_2, v1, v2, v3, v4, v5, v6,
//								v7, v8)
//						|| Math::point_inside_hexa(v3_2, v1, v2, v3, v4, v5, v6,
//								v7, v8)
//						|| Math::point_inside_hexa(v4_2, v1, v2, v3, v4, v5, v6,
//								v7, v8)) {
//					++nb_inter_;
//					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//							<< " contains a point of the cell " << cur_reg2_
//							<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//					return;
//				}
////		else if (v1_2 == v1 || v1_2 == v2 || v1_2 == v3 || v1_2 == v4
////						|| v1_2 == v5 || v1_2 == v6 || v1_2 == v7 || v1_2 == v8
////						|| v2_2 == v1 || v2_2 == v2 || v2_2 == v3 || v2_2 == v4
////						|| v2_2 == v5 || v2_2 == v6 || v2_2 == v7 || v2_2 == v8
////						|| v3_2 == v1 || v3_2 == v2 || v3_2 == v3 || v3_2 == v4
////						|| v3_2 == v5 || v3_2 == v6 || v3_2 == v7 || v3_2 == v8
////						|| v4_2 == v1 || v4_2 == v2 || v4_2 == v3 || v4_2 == v4
////						|| v4_2 == v5 || v4_2 == v6 || v4_2 == v7
////						|| v4_2 == v8) {
////					// if a point of one region is at the same place of the point of another region
////					++nb_inter_;
////					std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////							<< " contains a point of the cell " << cur_reg2_
////							<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////					return;
////				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 4) {
//					vec3& v5_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 4));
//					if (Math::point_inside_hexa(v5_2, v1, v2, v3, v4, v5, v6,
//							v7, v8)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v5_2 == v1 || v5_2 == v2 || v5_2 == v3
////							|| v5_2 == v4 || v5_2 == v5 || v5_2 == v6
////							|| v5_2 == v7 || v5_2 == v8) {
////						// if a point of one region is at the same place of the point of another region
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) > 5) {
//					vec3& v6_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 5));
//					if (Math::point_inside_hexa(v6_2, v1, v2, v3, v4, v5, v6,
//							v7, v8)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v6_2 == v1 || v6_2 == v2 || v6_2 == v3
////							|| v6_2 == v4 || v6_2 == v5 || v6_2 == v6
////							|| v6_2 == v7 || v6_2 == v8) {
////						// if a point of one region is at the same place of the point of another region
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//				if (mm_.mesh(cur_reg2_).cells.nb_vertices(idx) == 8) {
//					vec3& v7_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 6));
//					vec3& v8_2 = mm_.mesh(cur_reg2_).vertices.point(
//							mm_.mesh(cur_reg2_).cells.vertex(idx, 7));
//					if (Math::point_inside_hexa(v7_2, v1, v2, v3, v4, v5, v6,
//							v7, v8)
//							&& Math::point_inside_hexa(v8_2, v1, v2, v3, v4, v5,
//									v6, v7, v8)) {
//						++nb_inter_;
//						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
//								<< " contains a point of the cell " << cur_reg2_
//								<< ":" << idx << " region:cell" << std::endl;
//		fill_list_intersection(idx);
//						return;
//					}
////		else if (v7_2 == v1 || v7_2 == v2 || v7_2 == v3
////							|| v7_2 == v4 || v7_2 == v5 || v7_2 == v6
////							|| v7_2 == v7 || v7_2 == v8 || v8_2 == v1
////							|| v8_2 == v2 || v8_2 == v3 || v8_2 == v4
////							|| v8_2 == v5 || v8_2 == v6 || v8_2 == v7
////							|| v8_2 == v8) {
////						// if a point of one region is at the same place of the point of another region
////						++nb_inter_;
////						std::cout << "the cell " << cur_reg_ << ":" << cur_cell_
////								<< " contains a point of the cell " << cur_reg2_
////								<< ":" << idx << " region:cell" << std::endl;
////		fill_list_intersection(idx);
////						return;
////					}
//				}
//			}
//		}

	}
}

index_t DetectInter::detect_mesh_intersection() {
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
				mm_.tools.cell_aabb(cur_reg2_).compute_bbox_cell_bbox_intersections(
						box, *this);
			}
		}
	}

//	//Display the elements of inter_
//	for (index_t i = 0; i < inter_.size(); i++) {
//		if (inter_[i]) {
//			std::cout << "cell " << i << std::endl;
//		}
//	}
	return nb_inter_;
}

void DetectInter::fill_list_intersection(index_t idx) {
	int reg = cur_reg_;
	int reg2 = cur_reg2_;
	index_t id = idx;
	index_t cell = cur_cell_;
	while (reg2 - 1 >= 0) {
		reg2 -= 1;
		id += mm_.mesh(reg2).cells.nb();
	}
	if (!inter_[id]) {
		inter_[id] = true;
	}

	while (reg - 1 >= 0) {
		reg -= 1;
		cell += mm_.mesh(reg).cells.nb();
	}
	if (!inter_[cell]) {
		inter_[cell] = true;
	}
}

index_t DetectInter::check_angles(float angle_min) {
	index_t nb_angle = 0;
	for (index_t m = 0; m < mm_.nb_meshes(); m++) {
		for (index_t c = 0; c < mm_.mesh(m).cells.nb(); c++) {
			index_t nb_vert = mm_.mesh(m).cells.nb_vertices(c);
			for (index_t i = 0; i < nb_vert; i++) {
				float x0 = mm_.mesh(m).vertices.point(mm_.mesh(m).cells.vertex(c, i % nb_vert)).x;
				float y0 = mm_.mesh(m).vertices.point(mm_.mesh(m).cells.vertex(c, i % nb_vert)).y;
				float z0 = mm_.mesh(m).vertices.point(mm_.mesh(m).cells.vertex(c, i % nb_vert)).z;

				float x1 = mm_.mesh(m).vertices.point(mm_.mesh(m).cells.vertex(c, (i + 1) % nb_vert)).x;
				float y1 = mm_.mesh(m).vertices.point(mm_.mesh(m).cells.vertex(c, (i + 1) % nb_vert)).y;
				float z1 = mm_.mesh(m).vertices.point(mm_.mesh(m).cells.vertex(c, (i + 1) % nb_vert)).z;

				float x2 = mm_.mesh(m).vertices.point(mm_.mesh(m).cells.vertex(c, (i + 2) % nb_vert)).x;
				float y2 = mm_.mesh(m).vertices.point(mm_.mesh(m).cells.vertex(c, (i + 2) % nb_vert)).y;
				float z2 = mm_.mesh(m).vertices.point(mm_.mesh(m).cells.vertex(c, (i + 2) % nb_vert)).z;

				float norme1 = sqrt(
						pow((x0 - x1), 2) + pow((y0 - y1), 2)
								+ pow((z0 - z1), 2));
				float norme2 = sqrt(
						pow((x2 - x1), 2) + pow((y2 - y1), 2)
								+ pow((z2 - z1), 2));

				float angle = acos(
						sqrt(
								((x0 - x1) * (x2 - x1) + (y0 - y1) * (y2 - y1)
										+ (z0 - z1) * (z2 - z1))
										/ (norme1 * norme2)));

				std::cout << i << " angle " << angle << std::endl;
				if (angle < angle_min) {
					nb_angle++;
				}
			}
			std::cout << m << " " << c << " " << nb_angle << std::endl;
		}
	}
	return nb_angle;
}

index_t DetectInter::check_volumes() {
	for (index_t m = 0; m < mm_.nb_meshes(); m++) {
		for (index_t c = 0; c < mm_.mesh(m).cells.nb(); c++) {
			if (mm_.mesh(m).cells.nb_vertices(c) == 4) {
				mm_.mesh(m).cells;
				index_t B;
				index_t h;
				index_t vol = 1 / 3 * B * h;
				//TODO: calculate the volume
			}
		}
	}

	return 0;

}
}
