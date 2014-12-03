/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#include <grgmesh/grgmesh_assert.h>

#include <geogram/basic/assert.h>

namespace GRGMesh {

    static void grgmesh_abort()
    {
        GEO::geo_abort() ;
    }

    void GRGMESH_API grgmesh_assertion_failed(
        const std::string& condition_string,
        const std::string& file,
        int line )
    {
        GEO::geo_assertion_failed( condition_string, file, line ) ;
    }

    void GRGMESH_API grgmesh_should_not_have_reached(
        const std::string& file,
        int line )
    {
        GEO::geo_should_not_have_reached( file, line ) ;
    }

}

