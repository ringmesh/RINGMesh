/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram/points/nn_search.h>
#include <geogram/points/kd_tree.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>
#include <geogram/third_party/ANN/ANN.h>

namespace {

    using namespace GEO;

    /**
     * \brief Implementation of NearestNeighborSearch using the ANN library.
     */
    class NearestNeighborSearch_ANN : public NearestNeighborSearch {
    public:
        /**
         * \brief Constructs a new NearestNeighborSearch_ANN.
         * \param[in] dim dimension of the points
         */
        NearestNeighborSearch_ANN(
            coord_index_t dim
        ) :
            NearestNeighborSearch(dim),
            ann_tree_(nil) {
        }

        virtual void set_points(index_t nb_points, const double* points) {
            set_points(nb_points, points, dimension());
        }

        virtual bool stride_supported() const {
            return true;
        }

        virtual void set_points(
            index_t nb_points, const double* points, index_t stride
        ) {
            nb_points_ = nb_points;
            points_ = points;
            stride_ = stride;

            // Patched ANN so that we no longer need
            // to generate an array of pointers to
            // the points, See ANN.h
#ifdef ANN_CONTIGUOUS_POINT_ARRAY
            delete ann_tree_;
            ann_tree_ = new ANNkd_tree(
                ANNpointArray(points_, stride_),
                int(nb_points), 
                int(dimension())
            );
#else
            delete ann_tree_;
            ann_tree_ = nil;
            ann_points_.resize(nb_points);
            for(index_t i = 0; i < nb_points; i++) {
                ann_points_[i] = const_cast<double*>(points) + stride_ * i;
            }
            ann_tree_ = new ANNkd_tree(
                &ann_points_[0], int(nb_points), int(dimension())
            );
#endif
        }

        virtual void get_nearest_neighbors(
            index_t nb_neighbors,
            const double* query_point,
            index_t* neighbors,
            double* neighbors_sq_dist
        ) const {
            ann_tree_->annkSearch(
                const_cast<double*>(query_point),
                int(nb_neighbors), (ANNidxArray) neighbors, neighbors_sq_dist,
                (exact_ ? 0.0 : 0.1)
            );
        }

    protected:
        /**
         * \brief NearestNeighborSearch_ANN destructor
         */
        virtual ~NearestNeighborSearch_ANN() {
            delete ann_tree_;
            ann_tree_ = nil;
        }

    private:
#ifndef ANN_CONTIGUOUS_POINT_ARRAY
        std::vector<ANNcoord*> ann_points_;
#endif
        ANNkd_tree * ann_tree_;
    };
}

/****************************************************************************/

namespace GEO {

    NearestNeighborSearch::NearestNeighborSearch(
        coord_index_t dimension
    ) :
        dimension_(dimension),
        nb_points_(0),
        stride_(0),
        points_(nil),
        exact_(true) {
    }

    void NearestNeighborSearch::get_nearest_neighbors(
        index_t nb_neighbors,
        index_t query_point,
        index_t* neighbors,
        double* neighbors_sq_dist
    ) const {
        get_nearest_neighbors(
            nb_neighbors, 
            point_ptr(query_point), 
            neighbors, 
            neighbors_sq_dist
        );
    }

    void NearestNeighborSearch::set_points(
        index_t nb_points, const double* points
    ) {
        nb_points_ = nb_points;
        points_ = points;
        stride_ = dimension_;
    }

    bool NearestNeighborSearch::stride_supported() const {
        return false;
    }

    void NearestNeighborSearch::set_points(
        index_t nb_points, const double* points, index_t stride
    ) {
        if(stride == index_t(dimension())) {
            set_points(nb_points, points);
            return;
        }
        geo_assert(stride_supported());
        nb_points_ = nb_points;
        points_ = points;
        stride_ = stride;
    }

    void NearestNeighborSearch::set_exact(bool x) {
        exact_ = x;
    }

    NearestNeighborSearch::~NearestNeighborSearch() {
    }

    NearestNeighborSearch* NearestNeighborSearch::create(
        coord_index_t dimension, const std::string& name_in
    ) {
        geo_register_NearestNeighborSearch_creator(
            NearestNeighborSearch_ANN, "ANN"
        );
        geo_register_NearestNeighborSearch_creator(
            KdTree, "BNN"
        );

        std::string name = name_in;
        if(name == "default") {
            name = CmdLine::get_arg("algo:nn_search");
        }

        NearestNeighborSearch* nns =
            NearestNeighborSearchFactory::create_object(name, dimension);
        if(nns != nil) {
            return nns;
        }

        Logger::warn("NNSearch")
            << "Could not create NNSearch algorithm: " << name
            << std::endl
            << "Falling back to ANN"
            << std::endl;

        return new NearestNeighborSearch_ANN(dimension);
    }
}

