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

#include <grgmesh/boundary_model.h>
#include <grgmesh/io.h>
#include <iostream>
#include <fstream>

namespace GRGMesh {

    bool GRGMeshIO::load_BoundaryModel_from_Model3D(
        const std::string& filename,
        BoundaryModel& model )
    {
        BoundaryModelBuilder builder( model ) ;
        builder.load_file( filename ) ;
        return true ;
    }

}
