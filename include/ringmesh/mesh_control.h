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
#include <geogram/mesh/mesh_AABB.h>

//namespace GEO {
//class MeshTetsAABB ;
////class StoreIntersections;
//}

namespace RINGMesh {
    class RINGMESH_API DetectInter: public GEO::MeshTetsAABB {
    public:
        DetectInter( GEO::Mesh& M) ; //: GEO::MeshTetsAABB(M), inter_(inter){};

        ~DetectInter() ; // TODO check if need to be virtual

        void operator()( index_t idx ) ;

        void detect_mesh_intersection() ;

    protected:
        GEO::vector< bool > inter_ ;
    } ;

}

#include "mesh_control.hpp"

#endif /* MESH_CONTROL_H_ */
