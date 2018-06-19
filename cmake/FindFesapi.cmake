set (FESAPI_INCLUDE_DIR
    ${FESAPI_INSTALL_PREFIX}/include/
    ${FESAPI_INSTALL_PREFIX}/include/fesapi
    CACHE PATH "Directories containing the FESAPI header files" FORCE )

IF (UNIX)
	set (FESAPI_LIBRARY
		"${FESAPI_INSTALL_PREFIX}/lib/libFesapiCppUnderDev.so"
		CACHE PATH "Directories containing the FESAPI lib file" FORCE )
ELSEIF (WIN32)
	set (FESAPI_LIBRARY
		"${FESAPI_INSTALL_PREFIX}/lib/FesapiCppUnderDev.0.12.0.0.lib"
		CACHE PATH "Directories containing the FESAPI lib file" FORCE )
ENDIF (UNIX)

include(FindPackageHandleStandardArgs)

mark_as_advanced(FESAPI_INCLUDE_DIR FESAPI_LIBRARY)

if(NOT TARGET fesapi)

    add_library(fesapi UNKNOWN IMPORTED)

    # Interface include directory
    set_target_properties(fesapi PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${FESAPI_INCLUDE_DIR}"
    )

    # Link to library file
    set_target_properties(fesapi PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${FESAPI_LIBRARY}"
    )

endif()
