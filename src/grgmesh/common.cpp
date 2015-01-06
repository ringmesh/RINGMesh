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
 

#ifdef _MSC_VER

#pragma section(".CRT$XCU",read)
#define INITIALIZER(f) \
   static void __cdecl f(void); \
   __declspec(allocate(".CRT$XCU")) void (__cdecl*f##_)(void) = f; \
   static void __cdecl f(void)

#elif defined(__GNUC__)

#define INITIALIZER(f) \
   static void f(void) __attribute__((constructor)); \
   static void f(void)

#endif

#include <geogram/basic/common.h>

static bool inited = false ;
INITIALIZER(initialize)
{
    if( !inited ) {
        inited = true ;
        GEO::initialize() ;
        GEO::CmdLine::import_arg_group( "standard" ) ;
        GEO::CmdLine::import_arg_group( "algo" ) ;
        GEO::CmdLine::import_arg_group( "opt" ) ;
//        GEO::CmdLine::set_arg( "sys:assert", "abort" ) ;
        GEO::CmdLine::set_arg( "sys:FPE", false ) ;
    }
}


