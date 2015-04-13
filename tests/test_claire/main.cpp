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

	GEO::Logger::div(" RINGMesh - Test Claire ");
	GEO::Logger::out("TEST") << "Test for intersections" << std::endl;

	BoundaryModel boundary_model;
	RINGMeshIO::load("/home/launoy/RINGMesh/tests/data/noTouch.ml", //../data/noTouch.ml
			boundary_model);
	RINGMeshIO::save( boundary_model, "out.ml" ) ;

//	const std::string& filename = "../data/noTouch.ml";
//	std::ifstream input(filename.c_str());
//	if (!input) {
//		GEO::Logger::err("I/O") << "Cannot open file : " << std::endl;
//		std::cout << "fail " << std::endl;
//	}

	GEO::Logger::out("Boundary Model") << boundary_model.nb_regions()
			<< " regions" << std::endl;

	MacroMesh macro_mesh(boundary_model);
	RINGMeshIO::load("../data/noTouch10.mm", macro_mesh);

	GEO::Logger::div(" Control mesh quality ");

//  **********  In a region   ***********
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).x = 3;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).y = 3;
//	macro_mesh.mesh(0).vertices.point(macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).z = -3;

//	macro_mesh.mesh(0).vertices.point(macro_mesh.mesh(0).cells.vertex(5, 1)).x =
//			3;

	GEO::Logger::out("Mesh control") << "test" << std::endl;
	DetectInter inters(macro_mesh);
	GEO::Logger::out("Mesh control") << inters.detect_mesh_intersection()
			<< " intersections" << std::endl;

	return 0;
}
