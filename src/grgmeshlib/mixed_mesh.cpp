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

#include <grgmeshlib/mixed_mesh.h>

namespace GRGMesh {

    CellType MixedMesh::cell_type( uint64 c, uint64& c_index ) const
    {
        uint8 result = 0 ;
        uint64 size = 0 ;
        for( ; result < 7; result++ ) {
            size += cells_[result].size() ;
            if( c > size ) break ;
        }
        c_index = c - size ;
        return CellType( result ) ;
    }

}

