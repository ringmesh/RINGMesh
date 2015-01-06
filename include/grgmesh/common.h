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
 

#ifndef __GRGMESH_COMMON__
#define __GRGMESH_COMMON__

#if defined(_WIN32)
#    ifndef WIN32
#        define WIN32
#    endif
#endif

#ifdef WIN32
#   ifdef GRGMESH_EXPORTS
#        define GRGMESH_API __declspec( dllexport )
#    else
#        define GRGMESH_API __declspec( dllimport )
#    endif
#else
#   define GRGMESH_API
#endif

#ifndef NDEBUG
#   define GRGMESH_DEBUG
#endif

#ifdef WIN32
#   pragma warning( disable: 4267 )
#endif

#include <grgmesh/types.h>
#include <grgmesh/grgmesh_assert.h>

#endif

