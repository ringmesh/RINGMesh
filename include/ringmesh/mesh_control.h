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

#ifndef MESH_CONTROL_H_
#define MESH_CONTROL_H_

#include <ringmesh/common.h>
#include <ringmesh/utils.h>
#include <ringmesh/macro_mesh.h>
#include <geogram/mesh/mesh_AABB.h>

namespace RINGMesh {
class RINGMESH_API DetectInter {
public:
	DetectInter(const MacroMesh& mm);

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
	//void get_tet_bbox(const GEO::Mesh& M, GEO::Box& B, index_t t);

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
	index_t cur_reg_;
	index_t cur_reg2_;
	index_t cur_cell_;
	const MacroMesh& mm_;
	index_t nb_inter_;
	index_t indx_;
};

}

#endif /* MESH_CONTROL_H_ */
