/*
 * mesh_control.h
 *
 *  Created on: 13 mars 2015
 *      Author: launoy
 */

#ifndef MESH_CONTROL_H_
#define MESH_CONTROL_H_

#include <ringmesh/common.h>
#include <ringmesh/utils.h>
#include <ringmesh/macro_mesh.h>
#include <geogram/mesh/mesh_AABB.h>

namespace RINGMesh {
class RINGMESH_API DetectInter {
public:
	DetectInter(const MacroMesh& mm, index_t nb_reg);

	~DetectInter();

	/**
	 * \brief Operator (), changes into true the value of the given
	 * index in a list
	 * \param[in] action ACTION::operator(index_t) is
	 *  invoked for all cell that have a intersection.
	 */
	void operator()(index_t idx);

	bool mix_tetra_points(GEO::vec3& v1, GEO::vec3& v2, GEO::vec3& v3,
			GEO::vec3& v4, GEO::vec3& v5, GEO::vec3& v6, GEO::vec3& v7,
			GEO::vec3& v8);

	/**
	 * \brief Creates the box for a cell thanks to its given index.
	 * \param[M] the mesh containing the cell
	 * \param[B] the box that will be created
	 * \param[t] the index of the cell for which the box will be created
	 */
	void get_tet_bbox(const GEO::Mesh& M, GEO::Box& B, index_t t);

	/**
	 * \brief Computes all the intersections between two cells.
	 * \param[mm] macromesh in which the intersection want to be
	 * detected.
	 */
	void detect_mesh_intersection();

	/**
	 * \brief Returns true if the two tetrahedrons given intersect.
	 * \param[V] coordinates of the vertices of a cell (just done
	 * for tetrahedrons).
	 */
	bool tet_a_tet(double V_1[4][3], double V_2[4][3]);

	bool FaceA_1(double * Coord, int & maskEdges);

	bool FaceA_2(double * Coord, int & maskEdges);

private:
	std::vector<bool> inter_;
	index_t nb_reg_;
	index_t cur_reg_;
	index_t cur_reg2_;
	index_t cur_cell_;
	const MacroMesh& mm_;
	index_t nb_inter_;
	index_t indx_;
};

}

#endif /* MESH_CONTROL_H_ */
