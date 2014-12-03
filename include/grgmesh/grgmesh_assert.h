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

#ifndef __GRGMESH_ASSERT__
#define __GRGMESH_ASSERT__

#include <grgmesh/common.h>
#include <string>

namespace GRGMesh {

    static void grgmesh_abort() ;

    void GRGMESH_API grgmesh_assertion_failed(
        const std::string& condition_string,
        const std::string& file,
        int line ) ;

    void GRGMESH_API grgmesh_should_not_have_reached(
        const std::string& file,
        int line ) ;

}


#define grgmesh_assert(x) {                                        \
    if(!(x)) {                                                 \
        ::GRGMesh::grgmesh_assertion_failed(#x,__FILE__, __LINE__) ;   \
    }                                                          \
} 


#define grgmesh_assert_not_reached {                               \
    ::GRGMesh::grgmesh_should_not_have_reached(__FILE__, __LINE__) ;   \
}

#ifdef GRGMESH_DEBUG
#define grgmesh_debug_assert(x) grgmesh_assert(x)
#else
#define grgmesh_debug_assert(x)
#endif

#endif
