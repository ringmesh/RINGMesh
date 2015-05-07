/*
 * main.cpp
 *
 *  Created on: 13 avr. 2015
 *      Author: launoy
 */

#include <ringmesh/io.h>
#include <ringmesh/mesh_control.h>
#include <ringmesh/boundary_model.h>

#include <iostream>
#include <ctime>

int main(int argc, char** argv) {
	using namespace RINGMesh;

	GEO::Logger::div(" RINGMesh - Test intersection ");
	GEO::Logger::out("TEST") << "Test for intersections" << std::endl;

	BoundaryModel boundary_model;
	RINGMeshIO::load("../data/noTouch.ml", boundary_model);

	GEO::Logger::out("Boundary Model") << boundary_model.nb_regions()
			<< " regions" << std::endl;

	MacroMesh macro_mesh(boundary_model);
	RINGMeshIO::load("../data/noTouch5.mm", macro_mesh);

	GEO::Logger::div(" Control mesh quality ");

//	// Test intersection in a region
//	GEO::Logger::out("Mesh control") << "x "
//			<< macro_mesh.mesh(0).vertices.point(
//					macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).x << " y "
//			<< macro_mesh.mesh(0).vertices.point(
//					macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).y << " z "
//			<< macro_mesh.mesh(0).vertices.point(
//					macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).z
//			<< std::endl;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).x = 1;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).y = 1;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).z = -1;

//	// Test intersection between regions
	GEO::Logger::out("Mesh control") << "x "
			<< macro_mesh.mesh(1).vertices.point(
					macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).x << " y "
			<< macro_mesh.mesh(1).vertices.point(
					macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).y << " z "
			<< macro_mesh.mesh(1).vertices.point(
					macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).z << std::endl;
	macro_mesh.mesh(1).vertices.point(macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).x = 1;
	macro_mesh.mesh(1).vertices.point(macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).y = 1;
	macro_mesh.mesh(1).vertices.point(macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).z = -1;

//  // Display coordinates of the four points of all the tetrahedron
//	for (index_t cur_reg = 0; cur_reg < mm_.nb_meshes(); cur_reg++) {
//		for (index_t c = 0; c < mm_.mesh(cur_reg).cells.nb(); c++) {
//
//			vec3& v1 = mm_.mesh(cur_reg).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(c,0));
//			vec3& v2 = mm_.mesh(cur_reg).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(c,1));
//			vec3& v3 = mm_.mesh(cur_reg).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(c,2));
//			vec3& v4 = mm_.mesh(cur_reg).vertices.point(
//					mm_.mesh(cur_reg_).cells.vertex(c,3));
//			std::cout << "c " << c << " mesh " << cur_reg << "   " << v1.x
//					<< " " << v1.y << " " << v1.z << "   " << v2.x << " "
//					<< v2.y << " " << v2.z << "   " << v3.x << " " << v3.y
//					<< " " << v3.z << "   " << v4.x << " " << v4.y << " "
//					<< v4.z << std::endl;
//		}
//	}

	// Start timer
	time_t tbegin = time(NULL); // get the current calendar time

	GEO::Logger::out("Mesh control") << "nb cell " << macro_mesh.cells.nb_cells() << std::endl;
	DetectInter inters(macro_mesh);
	GEO::Logger::out("Mesh control") << inters.detect_mesh_intersection()
			<< " intersections" << std::endl;

	// End timer
	time_t tend = time(NULL); // get the current calendar time

	// Compute execution time
	double texec = difftime(tend, tbegin); // tend-tbegin (result in second)

	GEO::Logger::out("Mesh control") << "time " << texec << std::endl;

//	double angle_min;
//	GEO::Logger::out("Mesh control") << "Enter the minimal angle wanted (rad)"
//			<< std::endl;
//	std::cin >> angle_min;
//	index_t angle = inters.check_angles(angle_min);
//	index_t count = 0;
//	for (index_t i = 0; i < angle.size(); i++) {
//		if (angle[i]) {
//			count++;
//		}
//	}
//	GEO::Logger::out("Mesh control") << "There are " << angle
//			<< " angles below the wanted value" << std::endl;

	return 0;
}
