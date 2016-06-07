/*
 *      \V (O |R |P /A |L |I |N |E
 * (C) Bruno Levy, INRIA - ALICE, 2012-2014
 *
 *   Confidential - proprietary software
 */

#include <vorpalib/basic/common.h>
#include <geogram/basic/numeric.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/memory.h>
#include <geogram/mesh/index.h>
#include <geogram/delaunay/delaunay_nn.h>
#include <vector>

namespace GEO {
    class Mesh;


	/**
	 * \brief Computes the normal to the best approximating point
	 *  of a set of points.
	 * \code
	 * LeastSquaresNormal LSN ;
	 * LSN.begin() ;
	 *    LSN.add_point(P1) ;
	 *    ...
	 *    LSN.add_point(Pn) ;
	 * LSN.end() ;
	 * vec3 N = LSN.get_normal() ;
	 * \endcode
	 */
	class LeastSquaresNormal {
		public:
			/**
			 * \brief Creates a new LeastSquaresNormal
			 */
			LeastSquaresNormal();

			/**
			 * \brief Begins a least squares normal estimation.
			 */
			inline void begin() {
				points_.resize(0);
				center_ = vec3(0.0, 0.0, 0.0);
			}

			/**
			 * \brief Terminates a least squares normal estimation.
			 */
			inline void end() {
				if (nb_points() != 0) {
					center_ = (1.0 / double(nb_points())) * center_;
				}
			}

			/**
			 * \brief Adds a point to the current least squares normal
			 *  estimation.
			 * \param[in] p the current point
			 */
			inline void add_point(const vec3& p) {
				points_.push_back(p);
				center_ += p;
			}

			/**
			 * \brief Gets the estimated normal.
			 * \return the estimated normal
			 */
			const vec3& get_normal(); 

			/**
			 * \brief Gets the number of points
			 * \return the number of points used to estimate the normal
			 */
			inline index_t nb_points() const { return points_.size(); }

			/**
			 * \brief Gets a point by its index.
			 * \param[in] i index of the point
			 * \return a const reference to the point
			 */
			inline const vec3& point(index_t i) const { return points_[i]; }

			/**
			 * \brief Gets the center of the points.
			 * \return the center of the points
			 */
			inline const vec3& center() const { return center_; }

		private:
			vector<vec3> points_;
			vec3 center_;
			vec3 N_;
			double sym_mat_[6];
			double eig_vec_[9];
			double eig_val_[3];
	};
}
