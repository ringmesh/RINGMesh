/*
 * mesh_control.cpp
 *
 *  Created on: 13 mars 2015
 *      Author: launoy
 */

//#include <ringmesh/mesh_control.h>
#include <geogram/mesh/mesh_intersection.cpp>

namespace RINGMesh {
DetectInter::DetectInter(GEO::Mesh& M) :
		GEO::MeshTetsAABB(M), inter_(M.cells.nb(), false) {
	std::cout << "nb cells " << M.cells.nb() << std::endl;
}

DetectInter::~DetectInter() {
}

void DetectInter::operator()(index_t idx) {
	inter_[idx] = true;
}

void DetectInter::detect_mesh_intersection(MacroMesh& mm) {
	GEO::vector<RINGMesh::DetectInter> storeMesh;
	RINGMesh::index_t ind = 1;
//	index_t nb_meshes=mm_.nb_meshes();
//	for (index_t ind=1; ind<nb_meshes; ind++){
//		RINGMesh::DetectInter macmesh(mm.mesh(ind));
////		storeMesh.push_back(macmesh);
//	}

	for (index_t i = 0; i < bboxes_.size(); i++) {
		GEO::Box box = bboxes_[i];
		compute_bbox_cell_bbox_intersections(box, *this);

		index_t nb_inter = 0;
		for (index_t j = 0; j < inter_.size(); j++) {
			if (inter_[j]) {
				nb_inter = nb_inter + 1;
				inter_[j] = false;
			}
		}
		std::cout << nb_inter << std::endl;

//	std::cout << "x "<< bboxes_[10000].xyz_max[0] << " y "<< bboxes_[10000].xyz_max[1]
//	       << " z "<< bboxes_[10000].xyz_max[2]<< std::endl;
//	std::cout << "x "<< bboxes_[10000].xyz_min[0] << " y "<< bboxes_[10000].xyz_min[1]
//		       << " z "<< bboxes_[10000].xyz_min[2]<< std::endl;
	}
	// Finds intersections between boxes
//	compute_bbox_cell_bbox_intersections(box, *this);

	int count = 0;
	count = inter_.size();
	std::cout << count << " count " << std::endl;
	std::cout << bboxes_.size() << " bboxes size " << std::endl;

	index_t nb_inter = 0;

	for (index_t i = 0; i < inter_.size(); i++) {
		if (inter_[i]) {
			nb_inter = nb_inter + 1;
		}
	}
	std::cout << "nb_inter " << nb_inter << std::endl;

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

//	for (index_t i = 0; i < bboxes_.size(); i++) {
//		std::cout << i << " max " << bboxes_[i].xyz_max[0] << " "
//				<< bboxes_[i].xyz_max[1] << " " << bboxes_[i].xyz_max[2]
//				<< std::endl;
//		std::cout << i << " min " << bboxes_[i].xyz_min[0] << " "
//				<< bboxes_[i].xyz_min[1] << " " << bboxes_[i].xyz_min[2]
//				<< std::endl;
//	}

//	std::cout << " max " << bboxes_[10].xyz_max[0] << " "
//			<< bboxes_[10].xyz_max[10] << " " << bboxes_[10].xyz_max[2]
//			<< std::endl;
//	std::cout << " min " << bboxes_[10].xyz_min[0] << " "
//			<< bboxes_[10].xyz_min[1] << " " << bboxes_[10].xyz_min[2]
//			<< std::endl;

}
}
