/*
 * main.cpp
 *
 *  Created on: 13 avr. 2015
 *      Author: launoy
 */

#include <ringmesh/io.h>
#include <ringmesh/mesh_control.h>

int main(int argc, char** argv) {
	using namespace RINGMesh;

	GEO::Logger::div(" RINGMesh - Test intersection ");
	GEO::Logger::out("TEST") << "Test for intersections" << std::endl;

	BoundaryModel boundary_model;
	RINGMeshIO::load("../data/noTouch.ml", boundary_model);

	GEO::Logger::out("Boundary Model") << boundary_model.nb_regions()
			<< " regions" << std::endl;

	MacroMesh macro_mesh(boundary_model);
	RINGMeshIO::load("../data/noTouch10.mm", macro_mesh);

	GEO::Logger::div(" Control mesh quality ");

	// In a region
//	GEO::Logger::out("Mesh control") << "x "
//			<< macro_mesh.mesh(0).vertices.point(
//					macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).x << " y "
//			<< macro_mesh.mesh(0).vertices.point(
//					macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).y << " z "
//			<< macro_mesh.mesh(0).vertices.point(
//					macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).z
//			<< std::endl;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).x = 3;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).y = 3;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).z = -3;

	// Between regions
//	GEO::Logger::out("Mesh control") << "x "
//			<< macro_mesh.mesh(1).vertices.point(
//					macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).x << " y "
//			<< macro_mesh.mesh(1).vertices.point(
//					macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).y << " z "
//			<< macro_mesh.mesh(1).vertices.point(
//					macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).z << std::endl;
//	macro_mesh.mesh(1).vertices.point(macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).x = 9;
//	macro_mesh.mesh(1).vertices.point(macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).y = 9;
//	macro_mesh.mesh(1).vertices.point(macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).z = -9;

	// Display coordinates of the four points of all the tetrahedron
//	for (index_t cur_reg_ = 0; cur_reg_ < macro_mesh.nb_meshes(); cur_reg_++) {
//		for (index_t c = 0; c < macro_mesh.mesh(cur_reg_).cells.nb(); c++) {
//
//			vec3& v1 = macro_mesh.mesh(cur_reg_).vertices.point(
//					macro_mesh.mesh(cur_reg_).cells.facet_vertex(c, 0, 0));
//			vec3& v2 = macro_mesh.mesh(cur_reg_).vertices.point(
//					macro_mesh.mesh(cur_reg_).cells.facet_vertex(c, 0, 1));
//			vec3& v3 = macro_mesh.mesh(cur_reg_).vertices.point(
//					macro_mesh.mesh(cur_reg_).cells.facet_vertex(c, 0, 2));
//			vec3& v4 = macro_mesh.mesh(cur_reg_).vertices.point(
//					macro_mesh.mesh(cur_reg_).cells.facet_vertex(c, 2, 2));
//			std::cout << "c " << c << " mesh " << cur_reg_ << "   " << v1.x << " "
//					<< v1.y << " " << v1.z << "   " << v2.x << " " << v2.y
//					<< " " << v2.z << "   " << v3.x << " " << v3.y << " "
//					<< v3.z << "   " << v4.x << " " << v4.y << " " << v4.z
//					<< std::endl;
//		}
//	}

	GEO::Logger::out("Mesh control") << "test" << std::endl;
	DetectInter inters(macro_mesh);
	GEO::Logger::out("Mesh control") << inters.detect_mesh_intersection()
			<< " intersections" << std::endl;

	return 0;
}
