@echo OFF
setlocal ENABLEDELAYEDEXPANSION 
rem This file is to be called with the following arguments 
rem "cmake_generator"
rem By default it is "Visual Studio 10 Win64"


echo.
echo ============= Check CMake ============
echo.

echo Using cmake --version
if %errorlevel%==0 (
    echo Found CMake
) else (
    echo Error: CMake not found, please install it - see http://www.cmake.org/ - and ensure that is in your PATH
	pause
    exit /B 1
)

rem Set CMake GENERATOR by default it is Visual Studio 2010 in 64 bits
if [%1] == [] (
    set GENERATOR="Visual Studio 10 Win64"
) else (
    set GENERATOR=%1
)
echo Using cmake generator %GENERATOR%
set cmake_generator_options=-G %GENERATOR%


echo.
echo ============= Check RINGMesh platform ============
echo.

rem Get the Geogram and RINGMesh platform from CMake GENERATOR
rem The available platforms can be found in RINGMesh\src\third_party\geogram\cmake\platforms
rem For CMake > 3, the year of the compiler should be provided : "Visual Studio 10 2010 Win64"
if %GENERATOR%=="Visual Studio 10 2010 Win64" (
	set opsys=Win64-vs2010
) else if %GENERATOR%=="Visual Studio 10 Win64" (
	set opsys=Win64-vs2010
) else if %GENERATOR%=="Visual Studio 11 2012 Win64" (
	set opsys=Win64-vs2012
) else if %GENERATOR%=="Visual Studio 11 Win64" (
	set opsys=Win64-vs2012
) else if %GENERATOR%=="Visual Studio 12 2013 Win64" (
	set opsys=Win64-vs2013
) else if %GENERATOR%=="Visual Studio 12 Win64" (
	set opsys=Win64-vs2013
) else (
	echo The platform %GENERATOR% is not supported by RINGMesh
	pause
    exit /B 1
)

rem Set Windows as geogram platform
set geoplatform=Win-vs-dynamic-generic

rem Checking the geogram platform existence
cd %~dp0
if not exist "src\third_party\geogram\cmake\platforms\%opsys%" (
    echo Error: Geogram does not suppport the given platform: %opsys%
	pause
    exit /B 1
)

rem Set Visual Variables - TODO check this is really useful, I have doubts (Jeanne)
call "src\third_party\geogram\cmake\platforms\%opsys%\setvars.bat" 


echo.
echo ============= Check Visual Studio compilers ============
echo.


if "%CMAKE_VS_GENERATOR_TOOLSET%" neq "" (
    echo Using cmake generator toolset %CMAKE_VS_GENERATOR_TOOLSET%
    set cmake_generator_options=%cmake_generator_options% -T "%CMAKE_VS_GENERATOR_TOOLSET%"
)


echo.
echo ============= GeoGram ============
echo.


echo set(GEOGRAM_WITH_TETGEN TRUE) > src/third_party/geogram/CMakeOptions.txt
echo set(GEOGRAM_WITH_MEDIT FALSE) >> src/third_party/geogram/CMakeOptions.txt
echo set(GEOGRAM_WITH_GRAPHICS TRUE) >> src/third_party/geogram/CMakeOptions.txt


echo.
echo ============= Configuring and compiling Geogram  ============
echo.


if not exist build\geogram\%opsys% (
    mkdir build\geogram\%opsys%
)

cd build\geogram\%opsys%
cmake ..\..\..\src\third_party\geogram %cmake_debug_options% %cmake_generator_options% -DVORPALINE_PLATFORM:STRING=%geoplatform% 
cmake --build . --config Release
cmake --build . --config Debug

cd %~dp0


echo.
echo ============= Configuring RINGMesh ============
echo.


if not exist build\ringmesh\%opsys% (
    mkdir build\ringmesh\%opsys%
)

cd build\ringmesh\%opsys%
cmake ..\..\.. %cmake_debug_options% %cmake_generator_options% -DPLATFORM:STRING=%opsys%


cd %~dp0


echo.
echo ============== RINGMesh build configured ==================
echo.
echo To build RINGMesh:
echo - go to build/ringmesh/%opsys%
echo - run 'cmake --build . --config Release'
echo.
echo Note: local configuration can be specified in CMakeOptions.txt
echo See CMakeOptions.txt.sample for an example
echo You'll need to re-run configure.bat if you create or modify CMakeOptions.txt
echo.


set opsys=

pause
