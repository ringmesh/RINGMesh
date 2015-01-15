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

#ifndef __GRGMESH_IO__
#define __GRGMESH_IO__

#include <grgmesh/common.h>
#include <grgmesh/macro_mesh.h>
#include <grgmesh/boundary_model.h>


#define MAX_FILENAME 512
#define READ_SIZE 8192
class string ;

namespace GRGMesh {
    class BoundaryModel ;
}

namespace GRGMesh {
    namespace GRGMeshIO {
        bool GRGMESH_API load_BoundaryModel_from_Model3D(
            const std::string& filename,
            BoundaryModel& model ) ;
        bool GRGMESH_API save_BoundaryModel(
            BoundaryModel& model,
            const std::string& filename ) ;


        bool GRGMESH_API save_macro_mesh(
            const MacroMesh& mm,
            const std::string& filename ) ;
        bool GRGMESH_API load_macro_mesh(
            MacroMesh& mm,
            const std::string& mesh_file) ;
    }
}
#endif
