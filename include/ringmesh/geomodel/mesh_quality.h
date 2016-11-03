#ifndef __RINGMESH_MESH_QUALITY__
#define __RINGMESH_MESH_QUALITY__

#include <ringmesh/basic/common.h>

/*!
 * @author Benjamin Chauvin
 */

namespace RINGMesh {
    class GeoModel ;
}

namespace RINGMesh {

    enum MeshQualityMode {
        INSPHERE_RADIUS_BY_CIRCUMSPHERE_RADIUS,
        INSPHERE_RADIUS_BY_MAX_EDGE_LENGTH,
        VOLUME_BY_SUM_SQUARE_EDGE,
        MIN_SOLID_ANGLE
    } ;
    void RINGMESH_API compute_prop_tet_mesh_quality(
        MeshQualityMode mesh_qual_mode,
        const GeoModel& geo_model ) ;
}

#endif
