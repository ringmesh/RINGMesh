/*
 * mesh_control.cpp
 *
 *  Created on: 13 mars 2015
 *      Author: launoy
 */

//#include <ringmesh/mesh_control.h>
#include <geogram/mesh/mesh_intersection.cpp>

namespace RINGMesh{
	DetectInter::DetectInter( GEO::Mesh& M, GEO::vector<index_t>& inter ) :
			GEO::MeshTetsAABB(M), inter_(inter) {
			inter_.assign(0,0);//M.cells.nb(), 0);
//			for (RINGMesh::index_t i=0; i<5; i++){
//				std::cout << i << " inters " << inter_[i] << std::endl;
//			} // TODO use the constructor only once
	}

	DetectInter::~DetectInter(){}

	void DetectInter::operator() (index_t idx){
		inter_.push_back(idx);
//		inter_[idx] = 1;
//		std::cout << " idx " << idx << " " << inter_.size() << std::endl;
	}

	void DetectInter::detect_mesh_intersection(GEO::Mesh& M){
//        GEO::vector<index_t> has_intersection;
        DetectInter action(M, inter_);
        GEO::MeshTetsAABB AABB(M);

        GEO::Box box;

        GEO::vector<index_t> intersect;
        intersect.assign(bboxes_.size(), 0);

//        for (index_t i=0; i<bboxes_.size(); i++){
		box = bboxes_[1];
//			std::cout << box.xyz_max << std::endl;
		// Finds intersections between boxes
		AABB.compute_bbox_cell_bbox_intersections(box, action);
		// Check intersections between cells
		int count = 0;
		count = inter_.size();
		std::cout << count << " count " << std::endl;
		std::cout << bboxes_.size() << " bboxe size " << std::endl;

		std::cout << "coord box max x " << box.xyz_max[0] << std::endl;
		std::cout << "coord box max y " << box.xyz_max[1] << std::endl;
		std::cout << "coord box max z " << box.xyz_max[2] << std::endl;
		std::cout << "coord box min x " << box.xyz_min[0] << std::endl;
		std::cout << "coord box min y " << box.xyz_min[1] << std::endl;
		std::cout << "coord box min z " << box.xyz_min[2] << std::endl;

		std::cout << "other box" << std::endl;

		std::cout << "1 " << bboxes_[1000000].xyz_max[0] <<
				" " << bboxes_[1000000].xyz_max[1] <<
				" " << bboxes_[1000000].xyz_max[2] << std::endl;
		std::cout << "2 " << bboxes_[1000].xyz_max[0] <<
				" " << bboxes_[1000].xyz_max[1] <<
				" " << bboxes_[1000].xyz_max[2] << std::endl;
		std::cout << "3 " << bboxes_[10005].xyz_max[0] <<
				" " << bboxes_[10005].xyz_max[1] <<
				" " << bboxes_[10005].xyz_max[2] << std::endl;
		std::cout << "4 " << bboxes_[1022].xyz_max[0] <<
				" " << bboxes_[1022].xyz_max[1] <<
				" " << bboxes_[1022].xyz_max[2] << std::endl;
		std::cout << "5 " << bboxes_[7].xyz_max[0] <<
				" " << bboxes_[7].xyz_max[1] <<
				" " << bboxes_[7].xyz_max[2] << std::endl;
		std::cout << "6 " << bboxes_[20051].xyz_max[0] <<
				" " << bboxes_[20051].xyz_max[1] <<
				" " << bboxes_[20051].xyz_max[2] << std::endl;

		std::cout << "min " << std::endl;

		std::cout << "1 " << bboxes_[1000000].xyz_min[0] <<
				" " << bboxes_[1000000].xyz_min[1] <<
				" " << bboxes_[1000000].xyz_min[2] << std::endl;
		std::cout << "2 " << bboxes_[1000].xyz_min[0] <<
				" " << bboxes_[1000].xyz_min[1] <<
				" " << bboxes_[1000].xyz_min[2] << std::endl;
		std::cout << "3 " << bboxes_[10005].xyz_min[0] <<
				" " << bboxes_[10005].xyz_min[1] <<
				" " << bboxes_[10005].xyz_min[2] << std::endl;
		std::cout << "4 " << bboxes_[1022].xyz_min[0] <<
				" " << bboxes_[1022].xyz_min[1] <<
				" " << bboxes_[1022].xyz_min[2] << std::endl;
		std::cout << "5 " << bboxes_[7].xyz_min[0] <<
				" " << bboxes_[7].xyz_min[1] <<
				" " << bboxes_[7].xyz_min[2] << std::endl;
		std::cout << "6 " << bboxes_[20051].xyz_min[0] <<
				" " << bboxes_[20051].xyz_min[1] <<
				" " << bboxes_[20051].xyz_min[2] << std::endl;

//		for (RINGMesh::index_t j=0; j<inter_.size(); j++){
////				std::cout << inter_[j] << " inter " << std::endl;
//			if(inter_[j]==1){
//				count ++;
////				std::cout << count << " count " << std::endl;
//			}
//			intersect[i]=count;
//			std::cout << count << " intersection for " << i << std::endl;
//        }

//        int count = 0;
//        for (RINGMesh::index_t i=0; i<inter_.size(); i++){
//			if(inter_[i]==1)
//				count ++;
//		}
//        std::cout << "Count " << count << std::endl;
	}
}
