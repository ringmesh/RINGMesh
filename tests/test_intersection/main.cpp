/*
 * main.cpp
 *
 *  Created on: 13 avr. 2015
 *      Author: launoy
 */

#include <ringmesh/io.h>
#include <ringmesh/mesh_control.h>

#include <ringmesh/boundary_model.h>

int main(int argc, char** argv) {
	using namespace RINGMesh;

	GEO::Logger::div(" RINGMesh - Test intersection ");
	GEO::Logger::out("TEST") << "Test for intersections" << std::endl;

	BoundaryModel boundary_model;
	RINGMeshIO::load("../data/noTouch.ml", boundary_model);

	GEO::Logger::out("Boundary Model") << boundary_model.nb_regions()
			<< " regions" << std::endl;

	MacroMesh macro_mesh(boundary_model);
	RINGMeshIO::load("../data/noTouch2Hex.mm", macro_mesh);

	GEO::Logger::div(" Control mesh quality ");

//	// In a region
//	GEO::Logger::out("Mesh control") << "x "
//			<< macro_mesh.mesh(0).vertices.point(
//					macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).x << " y "
//			<< macro_mesh.mesh(0).vertices.point(
//					macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).y << " z "
//			<< macro_mesh.mesh(0).vertices.point(
//					macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).z
//			<< std::endl;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).x = 2;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).y = 2;
//	macro_mesh.mesh(0).vertices.point(
//			macro_mesh.mesh(0).cells.facet_vertex(9, 0, 0)).z = -1;

//	//Between regions
//	GEO::Logger::out("Mesh control") << "x "
//			<< macro_mesh.mesh(1).vertices.point(
//					macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).x << " y "
//			<< macro_mesh.mesh(1).vertices.point(
//					macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).y << " z "
//			<< macro_mesh.mesh(1).vertices.point(
//					macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).z << std::endl;
//	macro_mesh.mesh(1).vertices.point(macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).x = 1;
//	macro_mesh.mesh(1).vertices.point(macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).y = 1;
//	macro_mesh.mesh(1).vertices.point(macro_mesh.mesh(1).cells.facet_vertex(0, 0, 0)).z = -1;

	GEO::Logger::out("Mesh control") << "test" << std::endl;
	GEO::Logger::out("Mesh control") << macro_mesh.cells.nb_cells() << std::endl;
	DetectInter inters(macro_mesh);
	GEO::Logger::out("Mesh control") << inters.detect_mesh_intersection()
			<< " intersections" << std::endl;

	return 0;
}
