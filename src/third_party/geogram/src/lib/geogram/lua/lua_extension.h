/*
 *  Copyright (c) 2012-2016, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef GEOGRAM_LUA_LUA_EXTENSION

#include <geogram/basic/common.h>
#include <geogram/basic/assert.h>
#include <geogram/basic/numeric.h>
#include <geogram/basic/string.h>
#include <geogram/basic/memory.h>

extern "C" {
#include <geogram/third_party/lua/lua.h>    
#include <geogram/third_party/lua/lauxlib.h>
#include <geogram/third_party/lua/lualib.h>
}

/**
 * \file geogram/lua/lua_extension.h
 * \brief Utilities to write lua bindings.
 */

namespace GEO {


    typedef int (*lua_isfunc)(lua_State* L, int idx);

    int my_lua_isboolean(lua_State* L, int idx) {
	return lua_isboolean(L,idx);
    }

    int my_lua_ispositiveinteger(lua_State* L, int idx) {
	if(!lua_isinteger(L,idx)) {
	    return 0;
	}
	lua_Integer x = lua_tointeger(L,idx);
	return (x > 0);
    }
    
    inline void lua_check_type(
	lua_State* L, int idx, lua_isfunc test = nil, int expected_type=LUA_TNIL
    ) {
	int type = lua_type(L, idx);
	if(expected_type != LUA_TNIL && type != expected_type) {
	    std::string error =
		std::string("Argument ") + String::to_string(idx) +
		": expected " + lua_typename(L,expected_type) + " but got "
		+ lua_typename(L,type);
	    throw(error);
	}
	if(test != nil && !test(L,idx)) {
	    std::string error = std::string("Argument ") +
		": conversion error";
	    throw(error);
	}
    }

    
    template <class T> struct lua_to {
	lua_to(lua_State* L, int idx) {
	    geo_argused(L);
	    geo_argused(idx);
	    geo_assert_not_reached;
	}
    };

    template <class T> class lua_to_x {
      public:
	lua_to_x(lua_State* L, int idx, lua_isfunc test=nil, int lua_type=LUA_TNIL) {
	    lua_check_type(L,idx,test,lua_type);
	}
	operator T() const {
	    return x_;
	}
      protected:
	T x_;
    };
    
    template<> struct lua_to<int> : public lua_to_x<int> {
      lua_to(lua_State* L, int idx) : lua_to_x(L,idx,lua_isinteger) {
	    x_ = int(lua_tointeger(L,idx));
	}
    };

    template<> struct lua_to<index_t> : public lua_to_x<index_t> {
      lua_to(lua_State* L, int idx) : lua_to_x(L,idx,my_lua_ispositiveinteger) {
	    x_ = index_t(lua_tointeger(L,idx));
	}
    };
    
    template<> struct lua_to<float> : public lua_to_x<float> {
      lua_to(lua_State* L, int idx) : lua_to_x(L,idx,lua_isnumber) {
	    x_ = float(lua_tonumber(L,idx));
	}
    };
    
    template<> struct lua_to<double> : public lua_to_x<double> {
      lua_to(lua_State* L, int idx) : lua_to_x(L,idx,lua_isnumber) {
	    x_ = double(lua_tonumber(L,idx));
	}
    };

    template<> struct lua_to<bool> : public lua_to_x<bool> {
      lua_to(lua_State* L, int idx) : lua_to_x(L,idx,my_lua_isboolean) {
	    x_ = (lua_toboolean(L,idx) != 0);
	}
    };

    template<> struct lua_to<const char*> : public lua_to_x<const char*> {
      lua_to(lua_State* L, int idx) : lua_to_x(L,idx,lua_isstring) {
	    x_ = lua_tostring(L,idx);
	}
    };

    template<> struct lua_to<const std::string&> : public lua_to_x<std::string> {
      lua_to(lua_State* L, int idx) : lua_to_x(L,idx,lua_isstring) {
	    x_ = lua_tostring(L,idx);
	}
      operator const std::string&() {
	  return x_;
      }
    };

    /**********************************************************************/
    
    template<class T> void lua_push(lua_State* L, T x) {
	geo_argused(L);
	geo_argused(x);
	geo_assert_not_reached;
    }

    template<> void lua_push(lua_State* L, int x) {
	lua_pushinteger(L,lua_Integer(x));
    }

    template<> void lua_push(lua_State* L, index_t x) {
	lua_pushinteger(L,lua_Integer(x));
    }

    template<> void lua_push(lua_State* L, float x) {
	lua_pushnumber(L,lua_Number(x));
    }

    template<> void lua_push(lua_State* L, double x) {
	lua_pushnumber(L,lua_Number(x));
    }

    template<> void lua_push(lua_State* L, bool x) {
	lua_pushboolean(L,x?1:0);
    }
    
    template<> void lua_push(lua_State* L, const char* x) {
	lua_pushstring(L,x);
    }

    template<> void lua_push(lua_State* L, const std::string& x) {
	lua_pushstring(L,x.c_str());
    }
    
    /**********************************************************************/


    #define LUA_WRAP_BEGIN try {

    #define LUA_WRAP_END } catch(std::string err) {			       \
	return luaL_error(L,(std::string(__FUNCTION__) + ": " + err).c_str()); \
    } 
	
    template <class R> int lua_wrap(lua_State* L, R (*fptr)(void)) {
	R retval = fptr();
	lua_push(L,retval);
	return 0;
    }
    
    template <class R, class T1> int lua_wrap(lua_State* L, R (*fptr)(T1)) {
	LUA_WRAP_BEGIN
        T1 arg1 = lua_to<T1>(L,1);
	R retval = fptr(arg1);
        lua_push(L,retval);
        LUA_WRAP_END
	return 1;
    }

    template <class R, class T1, class T2> int lua_wrap(lua_State* L, R (*fptr)(T1,T2)) {
        LUA_WRAP_BEGIN
        T1 arg1 = lua_to<T1>(L,1);
        T2 arg2 = lua_to<T2>(L,2);	    
        R retval = fptr(arg1,arg2);
        lua_push(L,retval);
        LUA_WRAP_END
        return 1;
    }

}

#endif
