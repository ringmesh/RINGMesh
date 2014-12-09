@echo OFF
setlocal ENABLEDELAYEDEXPANSION 

IF [%1] == [] (
    set GENERATOR="Visual Studio 10 Win64"
) else (
    set GENERATOR=%1
)

IF [%2] == [] (
    set CMAKE_PATH="C:\Program Files (x86)\CMake 2.8"
) else (
    set CMAKE_PATH=%2
)
set CMAKE_EXECUTABLE=%CMAKE_PATH%\bin\cmake

rem Possible values for the generator (-G) are as follow. Ues the one suitable for you
rem
rem 	"Borland Makefiles"
rem 	"MSYS Makefiles"
rem 	"MinGW Makefiles"
rem 	"NMake Makefiles"
rem 	"Unix Makefiles"
rem 	"Visual Studio 6"
rem 	"Visual Studio 7"
rem 	"Visual Studio 7 .NET 2003"
rem 	"Visual Studio 8 2005"
rem 	"Visual Studio 8 2005 Win64"
rem 	"Visual Studio 9 2008"
rem 	"Visual Studio 9 2008 Win64"
rem 	"Visual Studio 10 Win64"
rem 	"Watcom WMake"
rem 	"CodeBlocks - MinGW Makefiles"
rem 	"CodeBlocks - Unix Makefiles"
rem 	"Eclipse CDT4 - MinGW Makefiles"
rem 	"Eclipse CDT4 - NMake Makefiles"
rem 	"Eclipse CDT4 - Unix Makefiles"
rem *****************************************************************************



set opsys=%3

:: Checking for CMake

echo.
echo ============= Checking for CMake ============
echo.

%CMAKE_EXECUTABLE% --version
if %errorlevel% == 0 (
    echo Found CMake
) else (
    echo Error: CMake not found, please install it - see http://www.cmake.org/
    exit /B 1
)

:: Checking the current OS

if "%opsys%" == "" (
    if "%PROCESSOR_ARCHITECTURE%" == "x86" (
        set opsys=Win32-vs2012
    ) else if "%PROCESSOR_ARCHITECTURE%" == "AMD64" (
        set opsys=Win64-vs2012
    ) else (
        echo Error: OS not supported
        exit /B 1
    )
)

if not exist "src\third_party\geogram\cmake\platforms\%opsys%" (
    echo Error: unsupported platform: %opsys%
    exit /B 1
)


:: Import the platform specific configuration

echo.
echo ============= Checking for Visual Studio ============
echo.


call "src\third_party\geogram\cmake\platforms\%opsys%\setvars.bat" 

echo Using cmake generator %GENERATOR%
set cmake_generator_options=-G %GENERATOR%

if "%CMAKE_VS_GENERATOR_TOOLSET%" neq "" (
    echo Using cmake generator toolset %CMAKE_VS_GENERATOR_TOOLSET%
    set cmake_generator_options=%cmake_generator_options% -T "%CMAKE_VS_GENERATOR_TOOLSET%"
)


:: Generate build tree


echo.
echo ============= GeoGram ============
echo.


echo.
echo ============= Creating build system for %opsys% ============
echo.


if not exist build\geogram\%opsys% (
    mkdir build\geogram\%opsys%
)

cd build\geogram\%opsys%
%CMAKE_EXECUTABLE% ..\..\..\src\third_party\geogram %cmake_debug_options% %cmake_generator_options% -DVORPALINE_PLATFORM:STRING=%opsys%
%CMAKE_EXECUTABLE% --build . --config Release
%CMAKE_EXECUTABLE% --build . --config Debug

cd %~dp0

:: Generate build tree

echo.
echo ============= GRGMesh ============
echo.


echo.
echo ============= Creating build system for %opsys% ============
echo.


if not exist build\grgmesh\%opsys% (
    mkdir build\grgmesh\%opsys%
)

cd build\grgmesh\%opsys%
%CMAKE_EXECUTABLE% ..\..\.. %cmake_debug_options% %cmake_generator_options% -DGEOGRAM_PLATFORM:STRING=%opsys%


cd %~dp0



echo.
echo ============== GRGMesh build configured ==================
echo.
echo To build GRGMesh:
echo - go to build/%opsys%
echo - run 'cmake --build . --config=Release(or Debug) [--target=target_to_build]'
echo.
echo Note: local configuration can be specified in CMakeOptions.txt
echo See CMakeOptions.txt.sample for an example
echo You'll need to re-run configure.bat if you create or modify CMakeOptions.txt
echo.

:: Clear globals

set opsys=

pause
