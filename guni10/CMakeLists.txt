###
#  @file CMakeLists.txt
#  @license
#    Copyright (c) 2013-2014
#    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao 
#
#    This file is part of Uni10.
#
#    Uni10 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Uni10 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
#  @endlicense
#  @brief Main specification file for CMake
#  @author Ying-Jer Kao
#  @date 2014-05-06
#  @since 0.1.0
###


######################################################################
### CMAKE VERSION
######################################################################
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


######################################################################
### DISABLING IN-SOURCE BUILD
### (DO NOT MODIFY)
######################################################################

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR 
    "
     In-source builds not allowed.\n
     Please\n
       1) make a build directory by:
            mkdir -p path/to/my/build/directory\n 
       2) run CMake in that build directory:
            cd path/to/my/build/directory
            cmake path/to/my/src/directory 
         
     You have to first remove the file \"CMakeCache.txt\" by:
         rm CMakeCache.txt
    "
  )
endif()

######################################################################
# Version information
######################################################################
set(UNI10_VERSION_MAJOR 0)
set(UNI10_VERSION_MINOR 9)
set(UNI10_VERSION_PATCH 0)
set(UNI10_VERSION COUNT 3)
set(UNI10_VERSION
  ${UNI10_VERSION_MAJOR}.${UNI10_VERSION_MINOR}.${UNI10_VERSION_PATCH}
)



 
######################################################################
### COMMAND LINE OPTIONS
######################################################################
#option(BUILD_STATIC_LIBS "Build the static library" ON)
#option(BUILD_SHARED_LIBS "Build the shared library" ON)
#option(BUILD_TESTS "Build test programs" OFF)
option(CUDA "Build using Nvidia CUDA for GPU library" OFF)
option(BUILD_EXAMPLES "Build Example Codes" ON)
option(BUILD_PYTHON_WRAPPER "Build Python wrapper" ON)
option(BUILD_WITH_MKL "Build Uni10 with MKL" OFF)
option(BUILD_WITH_INTEL_COMPILERS "Build Uni10 with Intel Compilers" OFF)
if (BUILD_WITH_MKL)
  option(MKL_SDL "Link to a single MKL dynamic libary." ON)
  option(MKL_MLT "Use multi-hreading libary. [Default]" ON)
  mark_as_advanced(MKL_SDL MKL_MLT)
endif()
if (BUILD_WITH_INTEL_COMPILERS)
  message("Build with Intel Compilers")
  message("Finding Intel Compilers....")
  find_program(CMAKE_C_COMPILER NAMES icc)
  find_program(CMAKE_CXX_COMPILER NAMES icpc)
  find_program(CMAKE_AR NAMES xiar)
  find_program(CMAKE_LINKER NAMES xild)
 
  if (NOT (CMAKE_C_COMPILER AND
           CMAKE_CXX_COMPILER AND
           CMAKE_AR AND
           CMAKE_LINKER ))
    message("Cannot find Intel compilers.")
    message("You may need to run `. /opt/intel/bin/compilervars.sh intel64/ia32'")
    message("Falling back to system compilers.")
    find_program(CMAKE_C_COMPILER NAMES cc)
    find_program(CMAKE_CXX_COMPILER NAMES cpp)
    find_program(CMAKE_AR NAMES ar)
    find_program(CMAKE_LINKER NAMES ld)
  endif ()
endif ()


######################################################################
### PROJECT
######################################################################

project(uni10)

######################################################################
### CMAKE INSTALL OPTIONS 
######################################################################
if(UNIX) 
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX "/usr/local/uni10" 
		CACHE PATH "UNI10 install prefix" 
   		FORCE
	)
  endif()
endif()
if(WIN32)
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "C:\\Program Files\\uni10" 
                CACHE PATH "UNI10 install prefix" 
                FORCE
        )
  endif()
endif()
######################################################################
### CMAKE BUILD TYPE
######################################################################

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)





######################################################################
### Find BLAS and LAPACK
######################################################################
if( NOT (DEFINED BLAS_LIBRARIES AND DEFINED LAPACK_LIBRARIES))
  if (BUILD_WITH_MKL)
     find_package( MKL REQUIRED)
     ADD_DEFINITIONS("-DMKL")
     
     set(LAPACK_LIBRARIES ${MKL_MINIMAL_LIBRARY})
     message("MKL: ${MKL_MINIMAL_LIBRARY}")
  else()
    find_package( LAPACK REQUIRED)
    message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
  endif()
else()
  set(LAPACK_LIBRARIES  ${BLAS_LIBRARIES}  ${LAPACK_LIBRARIES})
  message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
endif()

if(CUDA) 
  find_package(CUDA)
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_20,code=sm_20)
  cuda_compile(VECADD_O vecadd.cu)
endif()
######################################################################
### FLAGS
######################################################################
if(UNIX ) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()
if(MKL_MLT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if( UNIX OR MINGW OR CYGWIN) 
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x ")
    elseif (WINDOWS) 
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qstd=c++0x")
    endif()
endif()

if (MKL_MLT) 
  set(CMAKE_EXE_LINKER_FLAGS "-pthread")
endif()

######################################################################
### PATHS
######################################################################

#set the default path for built executables to the "bin" directory
if (UNIX OR MINGW OR CYGWIN)
	set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
endif()

#set the default path for built libraries to the "lib" directory
if (UNIX OR MINGW OR CYGWIN)
  set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
endif()

# path to include directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

######################################################################
### ADD SUBDIRECTORIES
######################################################################

add_subdirectory(src)
if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
if (BUILD_PYTHON_WRAPPER)
  add_subdirectory(Swig)
endif()

######################################################################
### ADD LIBRARY
######################################################################
set(uni10-objects
$<TARGET_OBJECTS:uni10-datatype>
$<TARGET_OBJECTS:uni10-data-structure>
$<TARGET_OBJECTS:uni10-tools>
$<TARGET_OBJECTS:uni10-tensor-network>
$<TARGET_OBJECTS:uni10-numeric>
)



if ("${CMAKE_GENERATOR}" STREQUAL "Xcode")
# Add Dummy file to Xcode Project to build the library
add_custom_command(
PRE_BUILD
OUTPUT ${CMAKE_BINARY_DIR}/src/dummy.cpp
COMMAND touch ${CMAKE_BINARY_DIR}/src/dummy.cpp
)
add_library(uni10 SHARED ${uni10-objects} ${CMAKE_BINARY_DIR}/src/dummy.cpp)
add_library(uni10-static STATIC ${uni10-objects} ${CMAKE_BINARY_DIR}/src/dummy.cpp)
else()
add_library(uni10 SHARED ${uni10-objects} )
add_library(uni10-static STATIC ${uni10-objects})
endif()
SET_TARGET_PROPERTIES(uni10-static PROPERTIES OUTPUT_NAME "uni10")
if(APPLE)
  set(MACOSX_RPATH 1)
  SET_TARGET_PROPERTIES(uni10 PROPERTIES VERSION ${UNI10_VERSION} SOVERSION ${UNI10_VERSION_MAJOR} MACOSX_RPATH 1)
else()
  SET_TARGET_PROPERTIES(uni10 PROPERTIES VERSION ${UNI10_VERSION} SOVERSION ${UNI10_VERSION_MAJOR})
endif()
target_link_libraries(uni10 ${LAPACK_LIBRARIES})
target_link_libraries(uni10-static ${LAPACK_LIBRARIES})
######################################################################
### RPATH SETTINGS
######################################################################
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
######################################################################
### INSTALL
######################################################################
install(TARGETS uni10 uni10-static
        LIBRARY 
        DESTINATION lib 
        COMPONENT libraries
        ARCHIVE
        DESTINATION lib 
        COMPONENT libraries)



######################################################################
### INSTALL HEADER FILES
######################################################################

install(DIRECTORY src/
	DESTINATION include 
        COMPONENT headers
	FILES_MATCHING PATTERN 		"*.h*"	
)

######################################################################
### INSTALL EXAMPLES
######################################################################
if(BUILD_EXAMPLES)
  install(DIRECTORY examples/
	DESTINATION examples/
	COMPONENT examples
	PATTERN "CMakeLists*" EXCLUDE
  )
endif()
if (BUILD_PYTHON_WRAPPER)
  install(DIRECTORY Swig/examples/
	DESTINATION examples/python/
        COMPONENT python_examples
  )
endif()
######################################################################
### INSTALL DOCUMENTATION
######################################################################

install(DIRECTORY doc/
	DESTINATION doc 
        COMPONENT documentation
)
install(FILES INSTALL.txt README.txt GPL.txt LGPL.txt release_notes.txt 
        DESTINATION .
        COMPONENT common)

set(libname "${CMAKE_SHARED_LIBRARY_PREFIX}uni10.${UNI10_VERSION}${CMAKE_SHARED_LIBRARY_SUFFIX}")



#if (UNIX)
  #### Make Symlink to /usr/local/lib ####
#  install(CODE " EXECUTE_PROCESS(
#     COMMAND ln -sf ${CMAKE_INSTALL_PREFIX}/lib/${libname} /usr/local/lib/${CMAKE_SHARED_LIBRARY_PREFIX}uni10.${UNI10_VERSION_MAJOR}${CMAKE_SHARED_LIBRARY_SUFFIX})")

#  install(CODE " EXECUTE_PROCESS(
#     COMMAND ln -sf ${CMAKE_INSTALL_PREFIX}/lib/${libname} /usr/local/lib/${CMAKE_SHARED_LIBRARY_PREFIX}uni10${CMAKE_SHARED_LIBRARY_SUFFIX})")
#endif()
######################################################################
# Information Messages
######################################################################
message(STATUS "")
message(STATUS "------------------------------------------------")
message(STATUS "  Uni10: The Universal Tensor Network Library")
message(STATUS "------------------------------------------------")
message(STATUS "")
message(STATUS " Version: ${UNI10_VERSION}")
message(STATUS " Generator: ${CMAKE_GENERATOR}")
message(STATUS " Build Target: ${CMAKE_SYSTEM_NAME}")
message(STATUS " Installation Prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS " CXX Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS " CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS " BLAS and LAPACK Libraries: ${LAPACK_LIBRARIES}")
if(BUILD_EXAMPLES)
  message(STATUS " Build Examples: YES")
else()
  message(STATUS " Build Examples: NO")
endif()
if(BUILD_PYTHON_WRAPPER)
  message(STATUS " Build Python Wrapper: YES")
  message(STATUS " Python Excutable  : ${PYTHON_EXECUTABLE}")
  message(STATUS " Python Headers    : ${PYTHON_INCLUDE_DIR}")
  message(STATUS " Python Library    : ${PYTHON_LIBRARY}")
  message(STATUS " Python Module Dir : ${PYTHON_SITE_PACKAGES}")
  message(STATUS " Swig Version      : ${SWIG_VERSION}")
  message(STATUS " Swig Directory    : ${SWIG_DIR}")
else()
  message(STATUS " Build Python Wrapper: NO")
endif()
message(STATUS "")
message(STATUS "")

######################################################################
### CPACK 
######################################################################
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Uni10: the Universal Tensor Network Library")
SET(CPACK_PACKAGE_VENDOR "Uni10")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LGPL.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "${UNI10_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${UNI10_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${UNI10_VERSION_PATCH}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "uni10 ${UNI10_VERSION_MAJOR}.${UNI10_VERSION_MINOR}")
#IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
#  SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
#  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
#  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} My Famous Project")
#  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.my-project-home-page.org")
#  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.my-personal-home-page.com")
#  SET(CPACK_NSIS_CONTACT "me@my-personal-home-page.com")
#  SET(CPACK_NSIS_MODIFY_PATH ON)
#ELSE(WIN32 AND NOT UNIX)
#  SET(CPACK_STRIP_FILES "bin/MyExecutable")
#  SET(CPACK_SOURCE_STRIP_FILES "")
#ENDIF(WIN32 AND NOT UNIX)
#SET(CPACK_PACKAGE_EXECUTABLES "MyExecutable" "My Executable")
INCLUDE(CPack)

