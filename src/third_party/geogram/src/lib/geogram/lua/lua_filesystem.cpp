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

#include <geogram/lua/lua_filesystem.h>
#include <geogram/lua/lua_extension.h>
#include <geogram/basic/file_system.h>


namespace {
    using namespace GEO;
    
    int lua_filesystem_is_file(lua_State* L) {
	return lua_wrap(L, FileSystem::is_file);
    }

    int lua_filesystem_is_directory(lua_State* L) {
	return lua_wrap(L, FileSystem::is_directory);
    }

    int lua_filesystem_create_directory(lua_State* L) {
	return lua_wrap(L, FileSystem::create_directory);
    }

    int lua_filesystem_delete_directory(lua_State* L) {
	return lua_wrap(L, FileSystem::delete_directory);
    }

    int lua_filesystem_delete_file(lua_State* L) {
	return lua_wrap(L, FileSystem::delete_file);
    }

    int lua_filesystem_get_directory_entries(lua_State* L) {
	if(lua_gettop(L) != 1) {
	    return luaL_error(
		L, "'FileSystem.get_directory_entries()' invalid number of arguments"
	    );
	}
	if(!lua_isstring(L,1)) {
	    return luaL_error(
		L, "'FileSystem.get_directory_entries()' argument should be a string"
	    );
	}
	const char* path = lua_tostring(L,1);
	std::vector<std::string> result;
	FileSystem::get_directory_entries(path, result);
	lua_newtable(L);
	for(size_t i=0; i<result.size(); ++i) {
	    lua_pushstring(L,result[i].c_str());
	    lua_seti(L,-2,lua_Integer(i+1));
	}
	return 1;
    }

    int lua_filesystem_get_current_working_directory(lua_State* L) {
	return lua_wrap(L, FileSystem::get_current_working_directory);
    }

    int lua_filesystem_set_current_working_directory(lua_State* L) {
	return lua_wrap(L, FileSystem::set_current_working_directory);
    }

    int lua_filesystem_rename_file(lua_State* L) {
	return lua_wrap(L, FileSystem::rename_file);
    }

    int lua_filesystem_get_time_stamp(lua_State* L) {
	return lua_wrap(L, FileSystem::get_time_stamp);
    }
}

#define DECLARE_FS_FUNC(F)                   \
    lua_pushliteral(L,#F);                   \
    lua_pushcfunction(L,lua_filesystem_##F); \
    lua_settable(L,1)

namespace GEO {
    void init_lua_filesystem(lua_State* L) {
	lua_newtable(L);
	DECLARE_FS_FUNC(is_file);
	DECLARE_FS_FUNC(is_directory);
	DECLARE_FS_FUNC(create_directory);
	DECLARE_FS_FUNC(delete_directory);
	DECLARE_FS_FUNC(delete_file);
	DECLARE_FS_FUNC(get_directory_entries);
	DECLARE_FS_FUNC(get_current_working_directory);
	DECLARE_FS_FUNC(set_current_working_directory);
	DECLARE_FS_FUNC(rename_file);
	DECLARE_FS_FUNC(get_time_stamp);

	lua_setglobal(L, "FileSystem");
    }
}

