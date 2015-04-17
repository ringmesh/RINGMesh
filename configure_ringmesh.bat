@echo OFF
setlocal ENABLEDELAYEDEXPANSION 
rem This file is to be called with the following arguments 
rem "cmake_generator" "path_to_cmake" "geogram_platform"
rem Default arguments are  "Visual Studio 10 Win64" "C:\Program Files (x86)\CMake 2.8"  "Win64-vs2010"

echo.
echo ============= Checking CMake path ============
echo.

rem Set the path to CMake. 
rem CMAKE3 is in by default in "C:\Program Files (x86)\CMake"
IF [%2] == [] (
    set CMAKE_PATH="C:\Program Files (x86)\CMake 2.8"
) else (
    set CMAKE_PATH=%2
)
set CMAKE_EXECUTABLE=%CMAKE_PATH%\bin\cmake


%CMAKE_EXECUTABLE% --version
if %errorlevel% == 0 (
    echo Found CMake
) else (
    echo Error: CMake not found, please install it - see http://www.cmake.org/
	pause
    exit /B 1
)

rem Set CMake GENERATOR by default it is Visual Studio 2010 in 64 bits
IF [%1] == [] (
    set GENERATOR="Visual Studio 10 Win64"
) else (
    set GENERATOR=%1
)
echo Using cmake generator %GENERATOR%
set cmake_generator_options=-G %GENERATOR%


echo.
echo ============= Check RINGMesh platform ============
echo.

rem Set Geogram and RINGMesh platforms 
rem WARNING: It should match the GENERATOR set for CMake
rem The available platforms can be found in RINGMesh\src\third_party\geogram\cmake\platforms
IF [%3] == [] (
    set opsys=Win64-vs2010
) else (
    set opsys=%3
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
rem Define WIN32, since it seems it is not well defined with Visual Studio 2013
rem echo ADD_DEFINITIONS(-DWIN32")" >> src/third_party/geogram/CMakeOptions.txt


echo.
echo ============= Configuring and compiling Geogram  ============
echo.


if not exist build\geogram\%opsys% (
    mkdir build\geogram\%opsys%
)

cd build\geogram\%opsys%
%CMAKE_EXECUTABLE% ..\..\..\src\third_party\geogram %cmake_debug_options% %cmake_generator_options% -DVORPALINE_PLATFORM:STRING=%geoplatform% 
rem Launching compiling does not work with CMake 3, the question is WHY ?
%CMAKE_EXECUTABLE% --build . --config Release
%CMAKE_EXECUTABLE% --build . --config Debug

cd %~dp0


echo.
echo ============= Configuring RINGMesh ============
echo.


if not exist build\ringmesh\%opsys% (
    mkdir build\ringmesh\%opsys%
)

cd build\ringmesh\%opsys%
%CMAKE_EXECUTABLE% ..\..\.. %cmake_debug_options% %cmake_generator_options% -DGEOGRAM_PLATFORM:STRING=%opsys%


cd %~dp0



echo.
echo ============== RINGMesh build configured ==================
echo.
echo To build RINGMesh:
echo - go to build/ringmesh/%opsys%
echo - run 'cmake --build . --config=Release(or Debug) [--target=target_to_build]'
echo.
echo Note: local configuration can be specified in CMakeOptions.txt
echo See CMakeOptions.txt.sample for an example
echo You'll need to re-run configure.bat if you create or modify CMakeOptions.txt
echo.


set opsys=

pause
