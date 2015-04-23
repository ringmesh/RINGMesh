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
	/**
	 * DetectInter is a class whose goal is to count and display
	 * the number of intersections inside a macro mesh
	 * \brief Constructor of the class
	 * \param[mm] macro mesh in which intersections are detected
	 */
	DetectInter(const MacroMesh& mm);

	/**
	 * Destructor of the class
	 */
	~DetectInter();

	/**
	 * \brief Operator (), returns the number of the cell
	 * and its mesh that contains a point of another cell
	 * \param[in] action ACTION::operator(index_t) is
	 *  invoked for all cell that have a intersection.
	 */
	void operator()(index_t idx);

	/**
	 * \brief Computes all the intersections between two cells.
	 * \param[mm] macromesh in which the intersection are
	 * detected.
	 */
	index_t detect_mesh_intersection();

	void fill_list_intersection(index_t idx);

	/**
	 * \brief Computes the volume of each cell to know if
	 * it is high enought to be coherent
	 */
	index_t check_volumes();

private:
	index_t cur_reg_;
	index_t cur_reg2_;
	index_t cur_cell_;
	const MacroMesh& mm_;
	index_t nb_inter_;
	std::vector<bool> inter_;
};

}

#endif /* MESH_CONTROL_H_ */
