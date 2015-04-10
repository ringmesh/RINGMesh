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

namespace RINGMesh {
DetectInter::DetectInter(const MacroMesh& mm, index_t nb_reg) :
		nb_reg_(nb_reg), inter_(mm.cells.nb_tet()), cur_reg_(0), cur_reg2_(0), nb_inter_(
				0), cur_cell_(0), mm_(mm), indx_(0) {

}

DetectInter::~DetectInter() {
}

void DetectInter::operator()(index_t idx) {
	indx_ = idx;
	if (idx > cur_cell_ || cur_reg_ != cur_reg2_) {

		vec3 v1 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 0));
		vec3 v2 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 1));
		vec3 v3 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 2));
		vec3 v4 = mm_.mesh(cur_reg_).vertices.point(
				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 2, 2));

		vec3 v1_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 0));
		vec3 v2_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 1));
		vec3 v3_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 2));
		vec3 v4_2 = mm_.mesh(cur_reg2_).vertices.point(
				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 2, 2));

//		double V1[4][3] = { { v3.x, v3.y, v3.z }, { v1.x, v1.y, v1.z }, { v4.x,
//				v4.y, v4.z }, { v2.x, v2.y, v2.z } };
//		double V2[4][3] = { { v1_2.x, v1_2.y, v1_2.z },
//				{ v2_2.x, v2_2.y, v2_2.z }, { v4_2.x, v4_2.y, v4_2.z }, {
//						v3_2.x, v3_2.y, v3_2.z } };

		double V3[4][3] = { { v1.x, v1.y, v1.z }, { v2.x, v2.y, v2.z }, { v3.x,
				v3.y, v3.z }, { v4.x, v4.y, v4.z } };
		double V4[4][3] = { { v1_2.x, v1_2.y, v1_2.z },
				{ v2_2.x, v2_2.y, v2_2.z }, { v3_2.x, v3_2.y, v3_2.z }, {
						v4_2.x, v4_2.y, v4_2.z } };

//		if ((V1[0] == V2[0] || V1[0] == V2[1] ||  V1[0] == V2[2] ||  V1[0] == V2[3]) &&
//				(V1[1] == V2[0] ||  V1[1] == V2[1] || V1[1] == V2[2] ||  V1[1] == V2[3]) &&
//				(V1[2] == V2[0] || V1[2] == V2[1] ||  V1[2] == V2[2] || V1[2] == V2[3]) &&
//				 V1[0] == V2[1] ||)){
		if (tet_a_tet(V3, V4)) {
			++nb_inter_;
			std::cout << "intersection region:cell " << cur_reg_ << ":"
					<< cur_cell_ << " " << cur_reg2_ << ":" << idx << std::endl;
		}
	}
}

void DetectInter::get_tet_bbox(const Mesh& M, Box& B, index_t t) {
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

void DetectInter::detect_mesh_intersection() {
	index_t nb_inter = 0;

//	double V2[4][3] = { { 20, 20, -20 }, { 0, 10, -10 }, { 10, 10, 0 }, { 6.67232, 6.65514, -3.34325 } };
//	double V3[4][3] = { { 10.1, 10.1, -20.1 }, { 13.4333, 16.7667, -15.4358 }, { 10.1, 10.1, -10.1 }, { 10.1, 20.1, -10.1 } };
//
//	if (tet_a_tet(V2, V3)) {
//		std::cout << "intersect2" << std::endl;
//	}

//	1 { 0, 10, -10 }
//	2 { 6.67232, 6.65514, -3.34325 }
//	3 { 20,20, -20 }
//	4 { 10, 10, 0 }

//	1 { 10.1, 20.1, -10.1 }
//	2 { 13.4333, 16.7667, -15.4358 }
//	3 { 10.1, 10.1, -10.1 }
//	4 { 10.1, 10.1, -20.1 }

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

//	double V[4][3] = { { 0, 0, 0 }, { 0, 1, 0 }, { 1, 0, 0 }, { 1, 1, 1 } };
//	double V2[4][3] = { { -1, -1, -1 }, { -1, 1, -1 }, { -1, 0, -1 }, { 0.5,
//			0.5, 0.5 } };
//	std::cout << "intersection (normally 1) " << tet_a_tet(V, V2) << std::endl;

//	std::cout << "nb cell " << mm_.mesh(0).cells.nb() << std::endl;
	for (index_t reg_idx = 0; reg_idx < nb_reg_; reg_idx++) {
		cur_reg_ = reg_idx;
//		std::cout << "cur_reg_ " << cur_reg_ << std::endl;
		for (index_t reg_idx2 = cur_reg_; reg_idx2 < nb_reg_; reg_idx2++) {
			cur_reg2_ = reg_idx2;
//			std::cout << "cur_reg2_ " << cur_reg2_ << std::endl;
			for (index_t c = 0; c < mm_.mesh(cur_reg_).cells.nb(); c++) {
				cur_cell_ = c;

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

				Box box;
				get_tet_bbox(mm_.mesh(cur_reg_), box, cur_cell_);

				mm_.tools.tet_aabb(cur_reg2_).compute_bbox_cell_bbox_intersections(
						box, *this);

//				std::cout << cur_reg_ << " " << cur_cell_ << "   "
//						<< mm_.mesh(cur_reg_).vertices.point(
//								mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_,
//										0, 0)).x << " "
//						<< mm_.mesh(cur_reg_).vertices.point(
//								mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_,
//										0, 0)).y << " "
//						<< mm_.mesh(cur_reg_).vertices.point(
//								mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_,
//										0, 0)).z << std::endl;

//				for (index_t c2 = 0; c2 < inter_[cur_reg2_].size(); c2++) {
//					index_t v1_index = mm_.mesh(cur_reg_).cells.facet_vertex(
//							cur_cell_, 0, 0);
//					index_t v2_index = mm_.mesh(cur_reg_).cells.facet_vertex(
//							cur_cell_, 0, 1);
//					index_t v3_index = mm_.mesh(cur_reg_).cells.facet_vertex(
//							cur_cell_, 0, 2);
//					index_t v4_index = mm_.mesh(cur_reg_).cells.facet_vertex(
//							cur_cell_, 2, 2);
//
//					vec3 v1 = mm_.mesh(cur_reg_).vertices.point(v1_index);
//					vec3 v2 = mm_.mesh(cur_reg_).vertices.point(v2_index);
//					vec3 v3 = mm_.mesh(cur_reg_).vertices.point(v3_index);
//					vec3 v4 = mm_.mesh(cur_reg_).vertices.point(v4_index);
//
//					double V[4][3] = { { v1.x, v1.y, v1.z },
//							{ v2.x, v2.y, v2.z }, { v3.x, v3.y, v3.z }, { v4.x,
//									v4.y, v4.z } };
//
//					index_t v1_2_index = mm_.mesh(cur_reg2_).cells.facet_vertex(
//							c2, 0, 0);
//					index_t v2_2_index = mm_.mesh(cur_reg2_).cells.facet_vertex(
//							c2, 0, 1);
//					index_t v3_2_index = mm_.mesh(cur_reg2_).cells.facet_vertex(
//							c2, 0, 2);
//					index_t v4_2_index = mm_.mesh(cur_reg2_).cells.facet_vertex(
//							c2, 2, 2);
//
//					vec3 v1_2 = mm_.mesh(cur_reg2_).vertices.point(v1_2_index);
//					vec3 v2_2 = mm_.mesh(cur_reg2_).vertices.point(v2_2_index);
//					vec3 v3_2 = mm_.mesh(cur_reg2_).vertices.point(v3_2_index);
//					vec3 v4_2 = mm_.mesh(cur_reg2_).vertices.point(v4_2_index);
//
//					double V2[4][3] = { { v1_2.x, v1_2.y, v1_2.z }, { v2_2.x,
//							v2_2.y, v2_2.z }, { v3_2.x, v3_2.y, v3_2.z }, {
//							v4_2.x, v4_2.y, v4_2.z } };
//
//						std::cout << "cur_cell_ " << cur_cell_ << std::endl;
//						std::cout << "coord tetra1 " << v1.x << " " << v1.y
//								<< " " << v1.z << "     " << v2.x << " " << v2.y
//								<< " " << v2.z << "     " << v3.x << " " << v3.y
//								<< " " << v3.z << "     " << v4.x << " " << v4.y
//								<< " " << v4.z << std::endl;
//
//						std::cout << "coord tetra2 " << v1_2.x << " " << v1_2.y
//								<< " " << v1_2.z << "     " << v2_2.x << " "
//								<< v2_2.y << " " << v2_2.z << "     " << v3_2.x
//								<< " " << v3_2.y << " " << v3_2.z << "     "
//								<< v4_2.x << " " << v4_2.y << " " << v4_2.z
//								<< std::endl;
//
//				}
//
//				std::cout << cur_reg_ << " " << cur_reg2_ << " " << cur_cell_
//						<< std::endl;
//			}
			}
		}
	}
	std::cout << "nb_inter " << nb_inter_ << std::endl;
}

// ----------- 3D algebraic operators -------------
#define DOT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

#define VECT(res,a,b){ \
	res[0] = a[1]*b[2]-b[1]*a[2];\
	res[1] = -a[0]*b[2]+b[0]*a[2];\
	res[2] = a[0]*b[1]-b[0]*a[1];\
}

#define SUB(res,a,b) {\
	res[0] = a[0]-b[0];\
	res[1] = a[1]-b[1];\
	res[2] = a[2]-b[2];\
}

#define SUB_DOT(a,b,c) (\
	 (a[0]-b[0])*c[0]+\
	 (a[1]-b[1])*c[1]+\
	 (a[2]-b[2])*c[2]\
)

typedef double point[3];
static point *V1, *V2;	// vertices coordinates

static int masks[4]; // for each face of the first tetrahedron stores the halfspace
// each vertex of the second tetrahedron belongs to

static double P_V1[4][3], P_V2[4][3]; // differences between the vertices of the second (first)
//  tetrahedron and the vertex 0  of the first(second) tetrahedron

static double Coord_1[4][4]; // vertices coordinates in the affine space

static double n[3];		// variable to store the normals

// FaceA ----------------------------------------------------

bool DetectInter::FaceA_1(double * Coord, int & maskEdges) {
	maskEdges = 000; // octal system

//	if (cur_reg_ == 0 && cur_cell_ == 1 && indx_ == 3) {
//		std::cout << "ok" << std::endl;
//
//		std::cout << "dot " << DOT(P_V1[0], n) << std::endl;
//		std::cout << "dot " << DOT(P_V1[1], n) << std::endl;
//		std::cout << "dot " << DOT(P_V1[2], n) << std::endl;
//		std::cout << "dot " << DOT(P_V1[3], n) << std::endl;
//	}

	if ((Coord[0] = DOT(P_V1[0], n)) >= 0) {
		maskEdges = 001;
	}
	if ((Coord[1] = DOT(P_V1[1], n)) >= 0) {
		maskEdges += 002;
	}
	if ((Coord[2] = DOT(P_V1[2], n)) >= 0) {
		maskEdges += 004;
	}
	if ((Coord[3] = DOT(P_V1[3], n)) >= 0) {
		maskEdges += 010;
	}

	return (maskEdges == 017); // if true it means that all of the vertices are out the halfspace
// defined by this face
}
// it is the same as FaceA_1, only the values V2[0]-v_ref are used only for the fourth face
// hence they do not need to be stored

bool DetectInter::FaceA_2(double * Coord, int & maskEdges) {
	maskEdges = 000;
	double * v_ref = V1[1];

	if ((Coord[0] = SUB_DOT(V2[0], v_ref, n)) >= 0) {
		maskEdges = 001;
	}
	if ((Coord[1] = SUB_DOT(V2[1], v_ref, n)) >= 0) {
		maskEdges += 002;
	}
	if ((Coord[2] = SUB_DOT(V2[2], v_ref, n)) >= 0) {
		maskEdges += 004;
	}
	if ((Coord[3] = SUB_DOT(V2[3], v_ref, n)) >= 0) {
		maskEdges += 010;
	}

	return (maskEdges == 017);
}

// FaceB --------------------------------------------------------------

inline static bool FaceB_1() {
	return ((DOT(P_V2[0] , n) >= 0) && (DOT(P_V2[1] , n) >= 0)
			&& (DOT(P_V2[2] , n) >= 0) && (DOT(P_V2[3] , n) >= 0));
}

inline static bool FaceB_2() {
	double * v_ref = V2[1];
	return (( SUB_DOT(V1[0],v_ref , n ) >= 0)
			&& ( SUB_DOT(V1[1],v_ref , n ) >= 0)
			&& ( SUB_DOT(V1[2],v_ref , n ) >= 0)
			&& ( SUB_DOT(V1[3],v_ref , n ) >= 0));
}

// EdgeA -------------------------------------------------------

inline static bool EdgeA(const int & f0, const int & f1) {
	double * coord_f0 = &Coord_1[f0][0];
	double * coord_f1 = &Coord_1[f1][0];

	int maskf0 = masks[f0];
	int maskf1 = masks[f1];

	if ((maskf0 || maskf1) != 017) { // if there is a vertex of b
		return false;	      // included in (-,-) return false
	}
	maskf0 &= (maskf0 ^ maskf1);  // exclude the vertices in (+,+)
	maskf1 &= (maskf0 ^ maskf1);

// edge 0: 0--1

	if (((maskf0 & 001) &&		// the vertex 0 of b is in (-,+)
			(maskf1 & 002)) &&		// the vertex 1 of b is in (+,-)
			(((coord_f0[1] * coord_f1[0]) - (coord_f0[0] * coord_f1[1])) >= 0)) {
		// the edge of b (0,1) intersect (-,-)
		return false;
	}
	if (((maskf0 & 002) && (maskf1 & 001))
			&& (((coord_f0[1] * coord_f1[0]) - (coord_f0[0] * coord_f1[1])) <= 0)) {
		return false;
	}
// edge 1: 0--2

	if (((maskf0 & 001) && (maskf1 & 004))
			&& (((coord_f0[2] * coord_f1[0]) - (coord_f0[0] * coord_f1[2])) >= 0)) {
		return false;
	}
	if (((maskf0 & 004) && (maskf1 & 001))
			&& (((coord_f0[2] * coord_f1[0]) - (coord_f0[0] * coord_f1[2])) <= 0)) {
		return false;
	}
// edge 2: 0--3

	if (((maskf0 & 001) && (maskf1 & 010))
			&& (((coord_f0[3] * coord_f1[0]) - (coord_f0[0] * coord_f1[3])) >= 0)) {
		return false;
	}
	if (((maskf0 & 010) && (maskf1 & 001))
			&& (((coord_f0[3] * coord_f1[0]) - (coord_f0[0] * coord_f1[3])) <= 0)) {
		return false;
	}
// edge 3: 1--2

	if (((maskf0 & 002) && (maskf1 & 004))
			&& (((coord_f0[2] * coord_f1[1]) - (coord_f0[1] * coord_f1[2])) >= 0)) {
		return false;
	}
	if (((maskf0 & 004) && (maskf1 & 002))
			&& (((coord_f0[2] * coord_f1[1]) - (coord_f0[1] * coord_f1[2])) <= 0)) {
		return false;
	}
// edge 4: 1--3

	if (((maskf0 & 002) && (maskf1 & 010))
			&& (((coord_f0[3] * coord_f1[1]) - (coord_f0[1] * coord_f1[3])) >= 0)) {
		return false;
	}
	if (((maskf0 & 010) && (maskf1 & 002))
			&& (((coord_f0[3] * coord_f1[1]) - (coord_f0[1] * coord_f1[3])) <= 0)) {
		return false;
	}
// edge 5: 2--3

	if (((maskf0 & 004) && (maskf1 & 010))
			&& (((coord_f0[3] * coord_f1[2]) - (coord_f0[2] * coord_f1[3])) >= 0)) {
		return false;
	}
	if (((maskf0 & 010) && (maskf1 & 004))
			&& (((coord_f0[3] * coord_f1[2]) - (coord_f0[2] * coord_f1[3])) <= 0)) {
		return false;
	}
	return true;// there exists a separating plane supported by the edge shared by f0 and f1
}

void side_face(std::vector<std::vector<bool> >& maskEdges, index_t i,
		index_t vect) {
	if (DOT(P_V1[i], n) >= 0) {
		maskEdges[vect].push_back(true);
	} else {
		maskEdges[vect].push_back(false);
	}
}

// main function

bool DetectInter::tet_a_tet(double V_1[4][3], double V_2[4][3]) {
	V1 = V_1;
	V2 = V_2;

	SUB(P_V1[0], V2[0], V1[0]); // V2-V1
	SUB(P_V1[1], V2[1], V1[0]);
	SUB(P_V1[2], V2[2], V1[0]);
	SUB(P_V1[3], V2[3], V1[0]);

//	if ((cur_reg_ == 0 && cur_cell_ == 6) && (cur_reg2_ == 1 && indx_ == 0)) {
////		std::cout << "V1 " << V1[0][0] << std::endl;
////		std::cout << "V1 " << V1[0][1] << std::endl;
////		std::cout << "V1 " << V1[0][2] << std::endl;
////		std::cout << "V2 " << V2[0][0] << std::endl;
////		std::cout << "V2 " << V2[0][1] << std::endl;
////		std::cout << "V2 " << V2[0][2] << std::endl;
////		std::cout << "P_V1[0]" << P_V1[0][0] << std::endl;
////		std::cout << "P_V1[1]" << P_V1[0][1] << std::endl;
////		std::cout << "P_V1[2]" << P_V1[0][2] << std::endl;
////
//		vec3 v1 = mm_.mesh(cur_reg_).vertices.point(
//				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 0));
//		vec3 v2 = mm_.mesh(cur_reg_).vertices.point(
//				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 1));
//		vec3 v3 = mm_.mesh(cur_reg_).vertices.point(
//				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 2));
//		vec3 v4 = mm_.mesh(cur_reg_).vertices.point(
//				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 2, 2));
//
//		vec3 v1_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(indx_, 0, 0));
//		vec3 v2_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(indx_, 0, 1));
//		vec3 v3_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(indx_, 0, 2));
//		vec3 v4_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(indx_, 2, 2));
//
//		std::cout << "cur_cell_ " << cur_cell_ << std::endl;
//		std::cout << "coord tetra1 " << v3.x << " " << v3.y << " " << v3.z
//				<< "     " << v2.x << " " << v2.y << " " << v2.z << "     "
//				<< v4.x << " " << v4.y << " " << v4.z << "     " << v1.x << " "
//				<< v1.y << " " << v1.z << std::endl;
//		std::cout << "coord tetra2 " << v1_2.x << " " << v1_2.y << " " << v1_2.z
//				<< "     " << v2_2.x << " " << v2_2.y << " " << v2_2.z
//				<< "     " << v4_2.x << " " << v4_2.y << " " << v4_2.z
//				<< "     " << v3_2.x << " " << v3_2.y << " " << v3_2.z
//				<< std::endl;
//	}

	static double e_v1[6][3], e_v2[6][3];   // vectors edge-oriented
	SUB(e_v1[0], V1[1], V1[0]);
//	if ((cur_reg_ == 0 && cur_cell_ == 6) && (cur_reg2_ == 1 && indx_ == 0)) {
//		std::cout << "0 " << e_v1[0][0] << std::endl;
//		std::cout << "1 " << e_v1[0][1] << std::endl;
//		std::cout << "2 " << e_v1[0][2] << std::endl;
//
//		std::cout << "0 V1 " << V1[0][0] << std::endl;
//		std::cout << "1 V1 " << V1[0][1] << std::endl;
//		std::cout << "2 V1 " << V1[0][2] << std::endl;
//
//		std::cout << "0 V1 " << V1[1][0] << std::endl;
//		std::cout << "1 V1 " << V1[1][1] << std::endl;
//		std::cout << "2 V1 " << V1[1][2] << std::endl;
//
//		std::cout << "0 V1 " << V1[2][0] << std::endl;
//		std::cout << "1 V1 " << V1[2][1] << std::endl;
//		std::cout << "2 V1 " << V1[2][2] << std::endl;
//	}
	SUB(e_v1[1], V1[2], V1[0]);

	VECT(n, e_v1[0], e_v1[1]);		// find the normal to face 0
	if ((cur_reg_ == 0 && cur_cell_ == 6) && (cur_reg2_ == 1 && indx_ == 0)) {
//		a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
		std::cout << "0 " << P_V1[0][0] << std::endl;
		std::cout << "1 " << P_V1[0][1] << std::endl;
		std::cout << "2 " << P_V1[0][2] << std::endl;
		std::cout << "0 " << n[0] << std::endl;
		std::cout << "1 " << n[1] << std::endl;
		std::cout << "2 " << n[2] << std::endl;

		std::cout << "n 1 " << n[0] + V1[0][0] << std::endl;
		std::cout << "n 2 " << n[1] + V1[0][1] << std::endl;
		std::cout << "n 3 " << n[2] + V1[0][2] << std::endl;

		std::cout << "dot 0 " << DOT(P_V1[0], n) << std::endl;
		std::cout << "dot 1 " << DOT(P_V1[1], n) << std::endl;
		std::cout << "dot 2 " << DOT(P_V1[2], n) << std::endl;
		std::cout << "dot 3 " << DOT(P_V1[3], n) << std::endl;
	}

	std::vector<std::vector<bool> > maskEdges(4);

	side_face(maskEdges, 0, 0);
	side_face(maskEdges, 1, 0);
	side_face(maskEdges, 2, 0);
	side_face(maskEdges, 3, 0);

	if ((cur_reg_ == 0 && cur_cell_ == 6) && (cur_reg2_ == 1 && indx_ == 0)) {
		std::cout << "here" << std::endl;
	}

	if (maskEdges[0][0] && maskEdges[0][1] && maskEdges[0][2]
			&& maskEdges[0][3]) {
		return false;
	}

	SUB(e_v1[2], V1[3], V1[0]);
	VECT(n, e_v1[2], e_v1[0]);

	side_face(maskEdges, 0, 1);
	side_face(maskEdges, 1, 1);
	side_face(maskEdges, 2, 1);
	side_face(maskEdges, 3, 1);

	if (maskEdges[1][0] && maskEdges[1][1] && maskEdges[1][2]
			&& maskEdges[1][3]) {
		return false;
	}

	VECT(n, e_v1[1], e_v1[2]);

	side_face(maskEdges, 0, 2);
	side_face(maskEdges, 1, 2);
	side_face(maskEdges, 2, 2);
	side_face(maskEdges, 3, 2);

	if (maskEdges[2][0] && maskEdges[2][1] && maskEdges[2][2]
			&& maskEdges[2][3]) {
		return false;
	}

	SUB(e_v1[3], V1[2], V1[1]);
	SUB(e_v1[4], V1[3], V1[1]);
	VECT(n, e_v1[4], e_v1[3]);

	side_face(maskEdges, 0, 3);
	side_face(maskEdges, 1, 3);
	side_face(maskEdges, 2, 3);
	side_face(maskEdges, 3, 3);

	if (maskEdges[3][0] && maskEdges[3][1] && maskEdges[3][2]
			&& maskEdges[3][3]) {
		return false;
	}

	if ((maskEdges[0][0] && maskEdges[1][0] && maskEdges[2][0]
			&& maskEdges[3][0])
			|| (maskEdges[0][1] && maskEdges[1][1] && maskEdges[2][1]
					&& maskEdges[3][1])
			|| (maskEdges[0][2] && maskEdges[1][2] && maskEdges[2][2]
					&& maskEdges[3][2])
			|| (maskEdges[0][3] && maskEdges[1][3] && maskEdges[2][3]
					&& maskEdges[3][3])) {
		return false;
	}

	return true;

//	if (FaceA_1(&Coord_1[0][0], masks[0])) {
//		return false;
//	}
//	SUB(e_v1[2], V1[3], V1[0]);
//	VECT(n, e_v1[2], e_v1[0]);
//
//	if (FaceA_1(&Coord_1[1][0], masks[1])) {
//		return false;
//	}
//	if (EdgeA(0, 1)) {
//		return false;
//	}
//	VECT(n, e_v1[1], e_v1[2]);
//
//	if (FaceA_1(&Coord_1[2][0], masks[2])) {
//		return false;
//	}
//	if (EdgeA(0, 2)) {
//		return false;
//	}
//	if (EdgeA(1, 2)) {
//		return false;
//	}
//
//	SUB(e_v1[4], V1[3], V1[1]);
//	SUB(e_v1[3], V1[2], V1[1]);
//
//	VECT(n, e_v1[4], e_v1[3]);
//
//	if (FaceA_2(&Coord_1[3][0], masks[3])) {
//		return false;
//	}
//	if (EdgeA(0, 3)) {
//		return false;
//	}
//
//	if (EdgeA(1, 3)) {
//		return false;
//	}
//	if (EdgeA(2, 3)) {
//		return false;
//	}
//	if ((masks[0] | masks[1] | masks[2] | masks[3]) != 017) {
//		return true;
//	}
//// from now on, if there is a separating plane it is parallel to a face of b
//
//	SUB(P_V2[0], V1[0], V2[0]);
//	SUB(P_V2[1], V1[1], V2[0]);
//	SUB(P_V2[2], V1[2], V2[0]);
//	SUB(P_V2[3], V1[3], V2[0]);
//
//	SUB(e_v2[0], V2[1], V2[0]);
//	SUB(e_v2[1], V2[2], V2[0]);
//
//	VECT(n, e_v2[0], e_v2[1]);
//	if (FaceB_1()) {
//		return false;
//	}
//	SUB(e_v2[2], V2[3], V2[0]);
//
//	VECT(n, e_v2[2], e_v2[0]);
//
//	if (FaceB_1()) {
//		return false;
//	}
//	VECT(n, e_v2[1], e_v2[2]);
//
//	if (FaceB_1()) {
//		return false;
//	}
//	SUB(e_v2[4], V2[3], V2[1]);
//	SUB(e_v2[3], V2[2], V2[1]);
//
//	VECT(n, e_v2[4], e_v2[3]);
//
//	if (FaceB_2()) {
//		return false;
//	}
//	return true;
}

#undef DOT
#undef SUB
#undef SUB_DOT
#undef VECT

}

