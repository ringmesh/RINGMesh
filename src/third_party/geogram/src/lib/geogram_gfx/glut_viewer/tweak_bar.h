/*
 *    _____   _       _   _   _____        _     _   _   _____   _          __  _____   _____
 *   /  ___| | |     | | | | |_   _|      | |   / / | | | ____| | |        / / | ____| |  _  \
 *   | |     | |     | | | |   | |        | |  / /  | | | |__   | |  __   / /  | |__   | |_| |
 *   | |  _  | |     | | | |   | |        | | / /   | | |  __|  | | /  | / /   |  __|  |  _  /
 *   | |_| | | |___  | |_| |   | |        | |/ /    | | | |___  | |/   |/ /    | |___  | | \ \
 *   \_____/ |_____| \_____/   |_|        |___/     |_| |_____| |___/|___/     |_____| |_|  \_\
 *
 *  Version 1.0
 *  Bruno Levy, August 2009
 *  INRIA, Project ALICE
 *
 * Helper C++ functions for AntTweakBar
 *
 */

#ifndef __GLUT_VIEWER_TWEAK_BAR_
#define __GLUT_VIEWER_TWEAK_BAR_

#include <AntTweakBar.h>
#include <GL/gl.h>
#include <string>

/****************************************************************************/

inline bool TwAddVarRW(TwBar* bar, const char* name, bool& var, const char* def = 0) {
    return TwAddVarRW(bar, name, TW_TYPE_BOOLCPP, &var, def);
}

inline bool TwAddVarRW(TwBar* bar, const char* name, GLboolean& var, const char* def = 0) {
    return TwAddVarRW(bar, name, TW_TYPE_BOOL8, &var, def);
}

inline bool TwAddVarRW(TwBar* bar, const char* name, int& var, const char* def = 0) {
    return TwAddVarRW(bar, name, TW_TYPE_INT32, &var, def);
}

inline bool TwAddVarRW(TwBar* bar, const char* name, unsigned int& var, const char* def = 0) {
    return TwAddVarRW(bar, name, TW_TYPE_UINT32, &var, def);
}

inline bool TwAddVarRW(TwBar* bar, const char* name, float& var, const char* def = 0) {
    return TwAddVarRW(bar, name, TW_TYPE_FLOAT, &var, def);
}

inline bool TwAddVarRW(TwBar* bar, const char* name, double& var, const char* def = 0) {
    return TwAddVarRW(bar, name, TW_TYPE_DOUBLE, &var, def);
}

inline bool TwAddVarRW(TwBar* bar, const char* name, std::string& var, const char* def = 0) {
    return TwAddVarRW(bar, name, TW_TYPE_STDSTRING, &var, def);
}

/****************************************************************************/

inline bool TwAddVarRO(TwBar* bar, const char* name, const bool& var, const char* def = 0) {
    return TwAddVarRO(bar, name, TW_TYPE_BOOLCPP, &var, def);
}

inline bool TwAddVarRO(TwBar* bar, const char* name, const GLboolean& var, const char* def = 0) {
    return TwAddVarRO(bar, name, TW_TYPE_BOOL8, &var, def);
}

inline bool TwAddVarRO(TwBar* bar, const char* name, const int& var, const char* def = 0) {
    return TwAddVarRO(bar, name, TW_TYPE_INT32, &var, def);
}

inline bool TwAddVarRO(TwBar* bar, const char* name, const unsigned int& var, const char* def = 0) {
    return TwAddVarRO(bar, name, TW_TYPE_UINT32, &var, def);
}

inline bool TwAddVarRO(TwBar* bar, const char* name, const float& var, const char* def = 0) {
    return TwAddVarRO(bar, name, TW_TYPE_FLOAT, &var, def);
}

inline bool TwAddVarRO(TwBar* bar, const char* name, const double& var, const char* def = 0) {
    return TwAddVarRO(bar, name, TW_TYPE_DOUBLE, &var, def);
}

inline bool TwAddVarRO(TwBar* bar, const char* name, const std::string& var, const char* def = 0) {
    return TwAddVarRO(bar, name, TW_TYPE_STDSTRING, &var, def);
}

/****************************************************************************/

/**
 * 'in' contains a ','-separated list of enum values. The numeric value of the
 *  i-th one is supposed to be equal to i (starting from 0).
 */
inline TwType TwDefineEnum(const char* name, const std::string& in) {
    std::vector<TwEnumVal> vals;
    std::vector<std::string> fields;
    // Split the 'in' string.
    char separator = ',';
    int length = (int) in.length();
    int start = 0;
    while(start < length) {
        int end = (int) in.find(separator, start);
        if(end < 0) {
            end = length;
        }
        if(end - start > 0) {
            while(start < end && in[start] == ' ') {
                start++;
            }
            while(end != length && in[end] == ' ' && end > start) {
                end--;
            }
            if(end - start > 0) {
                fields.push_back(in.substr(start, end - start));
            }
        }
        start = end + 1;
    }
    // Create the enum description for AntTweakBar
    vals.resize(fields.size());
    for(unsigned int i = 0; i < fields.size(); i++) {
        vals[i].Value = i;
        vals[i].Label = &(fields[i][0]);
    }
    // Note: AntTweakBar copies the definition into
    // its own data structure, so we can safely destroy
    // the local variables 'vals' and 'fields'.
    return TwDefineEnum(name, &(vals[0]), (unsigned int) vals.size());
}

#endif

