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

#include <grgmeshlib/grgmesh_assert.h>

#include <iostream>
#include <cstdlib>

namespace GRGMesh {

    static void grgmesh_abort()
    {
#ifdef WIN32
        // Under windows, rather than calling abort(),
        // we trigger a seg fault by deferencing the null pointer,
        // because abort() is more difficult to see in the debugger.
        // DebugBreak(); // Windows API function to trigger a breakpoint in the debugger.
        *((int*)0) = 0xbadbeef ;
#else
        abort() ;
#endif
    }

    void GRGMESH_API grgmesh_assertion_failed(
        const std::string& condition_string,
        const std::string& file,
        int line )
    {
        std::cerr << "Assertion failed: " << condition_string << std::endl ;
        std::cerr << "File: " << file << std::endl ;
        std::cerr << "Line: " << line << std::endl ;
        grgmesh_abort() ;
    }

    void GRGMESH_API grgmesh_range_assertion_failed(
        double value,
        double min_value,
        double max_value,
        const std::string& file,
        int line )
    {
        std::cerr << "Range assertion failed: " << value << " in " << "[ "
            << min_value << " ... " << max_value << " ]" << std::endl ;
        std::cerr << "File: " << file << std::endl ;
        std::cerr << "Line: " << line << std::endl ;
        grgmesh_abort() ;
    }

    void GRGMESH_API grgmesh_should_not_have_reached(
        const std::string& file,
        int line )
    {
        std::cerr << "Control should not have reached this point:" << std::endl ;
        std::cerr << "File: " << file << std::endl ;
        std::cerr << "Line: " << line << std::endl ;
        grgmesh_abort() ;
    }

}

