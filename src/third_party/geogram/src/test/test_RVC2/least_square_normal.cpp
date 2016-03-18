/*
 *      \V (O |R |P /A |L |I |N |E
 * (C) Bruno Levy, INRIA - ALICE, 2012-2014
 *
 *   Confidential - proprietary software
 */


#include "least_square_normal.h"
#include <geogram/numerics/matrix_util.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/process.h>
#include <geogram/basic/assert.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/algorithm.h>
#include <stack>

namespace GEO{
	LeastSquaresNormal::LeastSquaresNormal() {}

	const vec3& LeastSquaresNormal::get_normal() {
		for (index_t i = 0; i < 6; i++) {
			sym_mat_[i] = 0.0;
		}
		for (index_t i = 0; i < points_.size(); i++) {
			vec3 p = points_[i] - center_;
			sym_mat_[0] += p.x * p.x;
			sym_mat_[1] += p.x * p.y;
			sym_mat_[2] += p.y * p.y;
			sym_mat_[3] += p.x * p.z;
			sym_mat_[4] += p.y * p.z;
			sym_mat_[5] += p.z * p.z;
		}

		MatrixUtil::semi_definite_symmetric_eigen(sym_mat_, 3, eig_vec_,
				eig_val_);

		N_ = vec3(eig_vec_[6], eig_vec_[7], eig_vec_[8]);
		if (Numeric::is_nan(N_.x) || Numeric::is_nan(N_.y) ||
				Numeric::is_nan(N_.z)) {
			Logger::warn("Co3Ne") << "NAN detected!" << std::endl;
		}
		return N_;
	}
}
