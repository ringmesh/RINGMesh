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
    void get_tet_bbox( const Mesh& M, Box& B, index_t t )
    {
        const double* p = M.vertices.point_ptr( M.cells.vertex( t, 0 ) ) ;
        for( coord_index_t coord = 0; coord < 3; ++coord ) {
            B.xyz_min[coord] = p[coord] ;
            B.xyz_max[coord] = p[coord] ;
        }
        for( index_t lv = 1; lv < 4; ++lv ) {
            p = M.vertices.point_ptr( M.cells.vertex( t, lv ) ) ;
            for( coord_index_t coord = 0; coord < 3; ++coord ) {
                B.xyz_min[coord] = geo_min( B.xyz_min[coord], p[coord] ) ;
                B.xyz_max[coord] = geo_max( B.xyz_max[coord], p[coord] ) ;
            }
        }
    }
}

namespace RINGMesh {
    DetectInter::DetectInter( const MacroMesh& mm )
        :
            inter_( mm.cells.nb_tet() ),
            cur_reg_( 0 ),
            cur_reg2_( 0 ),
            nb_inter_( 0 ),
            cur_cell_( 0 ),
            mm_( mm ),
            indx_( 0 )
    {

    }

    DetectInter::~DetectInter()
    {
    }

    void DetectInter::operator()( index_t idx )
    {
        indx_ = idx ;
        if( idx > cur_cell_ || cur_reg_ != cur_reg2_ ) {

            // TODO simpler
//		GEO::vecng<3, double>& v1 = mm_.mesh(cur_reg_).vertices.point(
//				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 0));
//		GEO::vecng<3, double>& v2 = mm_.mesh(cur_reg_).vertices.point(
//				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 1));
//		GEO::vecng<3, double>& v3 = mm_.mesh(cur_reg_).vertices.point(
//				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 0, 2));
//		GEO::vecng<3, double>& v4 = mm_.mesh(cur_reg_).vertices.point(
//				mm_.mesh(cur_reg_).cells.facet_vertex(cur_cell_, 2, 2));
//
//		GEO::vecng<3, double>& v1_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 0));
//		vec3 v2_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 1));
//		vec3 v3_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 2));
//		vec3 v4_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 2, 2));

            const vec3& v1 = mm_.mesh( cur_reg_ ).vertices.point(
                mm_.mesh( cur_reg_ ).cells.facet_vertex( cur_cell_, 0, 0 ) ) ;
            const vec3& v2 = mm_.mesh( cur_reg_ ).vertices.point(
                mm_.mesh( cur_reg_ ).cells.facet_vertex( cur_cell_, 0, 1 ) ) ;
            const vec3& v3 = mm_.mesh( cur_reg_ ).vertices.point(
                mm_.mesh( cur_reg_ ).cells.facet_vertex( cur_cell_, 0, 2 ) ) ;
            const vec3& v4 = mm_.mesh( cur_reg_ ).vertices.point(
                mm_.mesh( cur_reg_ ).cells.facet_vertex( cur_cell_, 2, 2 ) ) ;

            const vec3& v1_2 = mm_.mesh( cur_reg2_ ).vertices.point(
                mm_.mesh( cur_reg2_ ).cells.facet_vertex( idx, 0, 0 ) ) ;
//		vec3 v2_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 1));
//		vec3 v3_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 0, 2));
//		vec3 v4_2 = mm_.mesh(cur_reg2_).vertices.point(
//				mm_.mesh(cur_reg2_).cells.facet_vertex(idx, 2, 2));

//GEO::vecng<3, index_t> a;

            Utils::point_inside_tetra( v1_2, v1, v2, v3, v4 ) ;
//		if (Utils::point_inside_tetra(v1_2, v1, v2, v3)) {
////				|| Utils::point_inside_tetra(v2_2, v1, v2, v3)
////				|| Utils::point_inside_tetra(v3_2, v1, v2, v3)
////				|| Utils::point_inside_tetra(v4_2, v1, v2, v3)) {
//			++nb_inter_;
//			std::cout << "intersection region:cell " << cur_reg_ << ":"
//					<< cur_cell_ << " " << cur_reg2_ << ":" << idx << std::endl;
//		}
        }
    }
// bool point_inside_tetra(const GEO::vecng<unsigned int3,double> &, const GEO::vecng<unsigned int3,double> &,
// const GEO::vecng<unsigned int3,double> &, const GEO::vecng<unsigned int3,double> &, const GEO::vecng<unsigned int3,double> &)

    void DetectInter::detect_mesh_intersection()
    {
        index_t nb_inter = 0 ;

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
        // TODO for cell in mm ...
        const index_t nb_reg = mm_.nb_meshes() ;
        for( index_t reg_idx = 0; reg_idx < nb_reg; reg_idx++ ) {
            cur_reg_ = reg_idx ;
//		std::cout << "cur_reg_ " << cur_reg_ << std::endl;
            for( index_t reg_idx2 = cur_reg_; reg_idx2 < nb_reg; reg_idx2++ ) {
                cur_reg2_ = reg_idx2 ;
//			std::cout << "cur_reg2_ " << cur_reg2_ << std::endl;
                for( index_t c = 0; c < mm_.mesh( cur_reg_ ).cells.nb(); c++ ) {
                    cur_cell_ = c ;

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

                    Box box ;
                    get_tet_bbox( mm_.mesh( cur_reg_ ), box, cur_cell_ ) ;

                    mm_.tools.tet_aabb( cur_reg2_ ).compute_bbox_cell_bbox_intersections(
                        box, *this ) ;

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
        std::cout << "nb_inter " << nb_inter_ << std::endl ;
    }
}

