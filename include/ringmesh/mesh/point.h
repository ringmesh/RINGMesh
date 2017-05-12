/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#pragma once

#include <ringmesh/basic/common.h>

/*!
 * @file ringmesh/geomodel.h
 * @author Antoine Mazuyer and Pierre Anquez
 */
namespace RINGMesh {

    class RINGMESH_API Point {
        static const index_t MAX_DIM = 3;

    public:
        Point() = delete;

        double operator[]( index_t i ) const
        {
            ringmesh_assert( i < dim_ );
            return data_[i];
        }

        double& operator[]( index_t i )
        {
            ringmesh_assert( i < dim_ );
            return data_[i];
        }

        double* data() const
        {
            return data_;
        }
    protected:
        Point( index_t dim )
            : dim_( dim )
        {
        }
    private:
        const double data_[MAX_DIM];
        index_t dim_;
    };

    class RINGMESH_API Point1D: public Point {

    public:
        Point1D()
            : Point( 1 )
        {
            data_[0] = 0.;
        }

        Point1D( double x )
            : dim_( 1 )
        {
            data_[0] = x;
        }
        double x() const
        {
            return data_[0];
        }
        double& x()
        {
            return data_[0];
        }
    };

    class RINGMESH_API Point2D: public Point {
    public:
        Point2D()
            : Point( 2 )
        {
            data_[0] = 0.;
            data_[1] = 0.;
        }
        Point2D( double x, double y )
            : dim_( 2 )
        {
            data_[0] = x;
            data_[1] = y;
        }
        double x() const
        {
            return data_[0];
        }
        double y() const
        {
            return data_[1];
        }
        double& x()
        {
            return data_[0];
        }
        double& y()
        {
            return data_[1];
        }
    };

    class RINGMESH_API Point3D: public Point {
    public:
        Point3D()
            : Point( 3 )
        {
            data_[0] = 0.;
            data_[1] = 0.;
            data_[2] = 0.;

        }
        Point3D( double x, double y, double z )
            : dim_( 3 )
        {
            data_[0] = x;
            data_[1] = y;
            data_[2] = z;

        }
        double x() const
        {
            return data_[0];
        }
        double y() const
        {
            return data_[1];
        }
        double z() const
        {
            return data_[2];
        }

        double& x()
        {
            return data_[0];
        }
        double& y()
        {
            return data_[1];
        }
        double& z()
        {
            return data_[2];
        }
    };
}
