#ifndef __RINGMESH_MESH_QUALITY__
#define __RINGMESH_MESH_QUALITY__

#include <ringmesh/basic/common.h>

/*!
 * @author Benjamin Chauvin
 * This code is inspired from
 * http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_quality/tet_mesh_quality.html
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

    /*!
     * @brief Computes and stores mesh quality in the GeoModel.
     *
     * Enables to have a metrics on the quality of the 3D mesh. This is
     * important for instance for numerical processes. For now all the qualities
     * are based on the regular tetrahedron: a good element here is an element
     * near equilaterality (this definition may change in function of the
     * requirements of the numerical process).
     * The quality is between 0 and 1. 0 corresponds to a bad tetrahedron, and
     * 1 to a good tetrahedron (equilaterality).
     * For more information about the mesh quality, see
     * <a href="http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_quality/tet_mesh_quality.html">
     * TET_MESH_QUALITY Interactive Program for Tet Mesh Quality</a>
     *
     * @param[in] mesh_qual_mode mesh quality to compute.
     * @param[in,out] geo_model GeoModel in which the mesh quality is performed.
     * The quality is stored on the cells of each Region.
     *
     * @warning The GeoModel must have at least one region. All the regions
     * must be meshed by simplexes (tetrahedra).
     */
    void RINGMESH_API compute_prop_tet_mesh_quality(
        MeshQualityMode mesh_qual_mode,
        const GeoModel& geo_model ) ;
}

#endif
