###
#  @file CMakeLists.txt
#  @license
#    Copyright (c) 2013-2014
#    National Taiwan University
#    National Tsing-Hua University
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
if(APPLE)
  # OS X requires 2.8.12 for rpath
  cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
else()
  # For RHEL/CentOS 7
  cmake_minimum_required(VERSION 2.8.11 FATAL_ERROR)
endif()

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
set(UNI10_VERSION_MAJOR 1)
set(UNI10_VERSION_MINOR 0)
set(UNI10_VERSION_PATCH 0)
set(UNI10_VERSION COUNT 6)
set(UNI10_VERSION
  ${UNI10_VERSION_MAJOR}.${UNI10_VERSION_MINOR}.${UNI10_VERSION_PATCH}
)




######################################################################
### COMMAND LINE OPTIONS
######################################################################
#option(BUILD_STATIC_LIBS "Build the static library" ON)
#option(BUILD_SHARED_LIBS "Build the shared library" ON)
#option(BUILD_TESTS "Build test programs" OFF)
option(BUILD_CUDA_SUPPORT "Build using Nvidia CUDA for GPU library" OFF)
option(BUILD_PACKAGES "Build Source/Binary Packages" ON)
option(BUILD_EXAMPLES "Build Example Codes" ON)
option(BUILD_PYTHON_WRAPPER "Build Python wrapper" OFF)
option(BUILD_WITH_MKL "Build Uni10 with MKL" OFF)
option(BUILD_WITH_INTEL_COMPILERS "Build Uni10 with Intel Compilers" OFF)
option(BUILD_ARPACK_SUPPORT "Build the arpack wrapper" OFF)
option(BUILD_DOC "Build API docuemntation" OFF)
option(BUILD_HDF5_SUPPORT "Build HDF5" OFF)

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
    if(WIN32)
      if ("$CMAKE_SYSTEM_PROCESSOR" STREQUAL "x86")
        set(BLAS_LIBRARIES "${CMAKE_CURRENT_SOURCE_DIR}/dep/x86/libblas.lib")
        set(LAPACK_LIBRARIES ${BLAS_LIBRARIES} "${CMAKE_CURRENT_SOURCE_DIR}/x86/liblapack.lib")
      else()
        set(BLAS_LIBRARIES "${CMAKE_CURRENT_SOURCE_DIR}/dep/x64/libblas.lib")
        set(LAPACK_LIBRARIES ${BLAS_LIBRARIES} "${CMAKE_CURRENT_SOURCE_DIR}/dep/x64/liblapack.lib")
      endif()
    else()
      find_package( LAPACK REQUIRED)
    endif()
    message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
  endif()
else()
  set(LAPACK_LIBRARIES  ${BLAS_LIBRARIES}  ${LAPACK_LIBRARIES})
  message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
endif()

######################################################################
### Find CUDA
######################################################################
if(BUILD_CUDA_SUPPORT)
  find_package(CUDA REQUIRED)
  list(APPEND CUDA_NVCC_FLAGS "-arch=sm_20;-std=c++11;-O2;-DCUDA_SUPPORT") #-DCUDA_SUPPORT to control the macro define
  SET(CUDA_PROPAGATE_HOST_FLAGS OFF)
  string(REPLACE "libcusparse" "libcusolver" CUDA_cusolver_LIBRARY ${CUDA_cusparse_LIBRARY})#Workaround
  #set(CUDA_VERBOSE_BUILD true)

  #FIND_PATH( CULA_INCLUDE_DIR "cula.h"
  #         PATH_SUFFIXES "cula/include" )
  #FIND_LIBRARY( CULA_LIBRARY
  #            NAMES "cula_lapack"
  #            PATH_SUFFIXES "cula/lib64" )
  #MESSAGE(STATUS "CULA_INCLUDE_DIR = ${CULA_INCLUDE_DIR}")
endif()



######################################################################
### Find ARPACK
######################################################################
IF( BUILD_ARPACK_SUPPORT )
 IF( BUILD_CUDA_SUPPORT)
   message( FATAL_ERROR "CANNOT BUILD ARPACK and CUDA SIMULTANEOUSLY")
 ELSE()
  if( NOT (DEFINED ARPACK_LIBRARIES) )
    find_package( ARPACK REQUIRED )
    set(ARPACK_LIBRARIES ${ARPACK_LIBRARY})
  else()
    set(ARPACK_LIBRARIES  ${ARPACK_LIBRARIES})
    message( STATUS "ARPACK found: ${ARPACK_LIBRARIES}")
  endif()
 ENDIF()
ENDIF()
######################################################################
### Find HDF5 Library and Include dirs
######################################################################
IF(BUILD_HDF5_SUPPORT)
  ADD_DEFINITIONS("-DHDF5")
  if( NOT (DEFINED HDF5_C_LIBRARIES AND DEFINED HDF5_CXX_LIBRARIES AND DEFINED HDF5_INCLUDE_DIRS))
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
    set(HDF5_LIBs ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    find_package(ZLIB)
    set( HDF5_LIBs ${HDF5_LIBs} ${ZLIB_LIBRARIES})
    #TODO: Need win32
  else()
    set(HDF5_LIBs ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    message( STATUS "HDF5 found: ${HDF5_LIBRARIES}")
  endif()
  include_directories(${HDF5_INCLUDE_DIRS})
ENDIF()
######################################################################
### FLAGS
######################################################################
if(UNIX )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()
if(MKL_MLT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
if( UNIX OR MINGW OR CYGWIN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
elseif (WINDOWS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qstd=c++11")
endif()

if (MKL_MLT)
  set(CMAKE_EXE_LINKER_FLAGS "-pthread")
endif()

######################################################################
### PATHS
######################################################################

#set the default path for built executables to the "bin" directory
if (UNIX OR MINGW OR CYGWIN)
  set(EXECUTABLE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

#set the default path for built libraries to the "lib" directory
#if (UNIX OR MINGW OR CYGWIN)
#  set(LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#endif()

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
  add_subdirectory(Swig pyUni10 )
endif()
if (NOT WIN32)
 add_subdirectory(gtest-1.7.0)
 add_subdirectory(test)
endif()
######################################################################
### ADD LIBRARY
######################################################################
set(uni10-objects
$<TARGET_OBJECTS:uni10-datatype>
$<TARGET_OBJECTS:uni10-data-structure>
$<TARGET_OBJECTS:uni10-tools>
$<TARGET_OBJECTS:uni10-tensor-network>
)
IF( BUILD_ARPACK_SUPPORT )
  set(uni10-objects ${uni10-objects}
  $<TARGET_OBJECTS:uni10-arpack>
  )
ENDIF()
IF(BUILD_HDF5_SUPPORT)
  set(uni10-objects ${uni10-objects}
  $<TARGET_OBJECTS:uni10-hdf5io>
  )
ENDIF()

IF(BUILD_CUDA_SUPPORT)
#### Build CUDA
# set(CUDA_SEPARABLE_COMPILATION ON)
#find_library(CULA_LIBRARY libcula_lapack /usr/local/cula/lib64 /opt/local/cula/lib64)
 CUDA_ADD_LIBRARY(uni10gpu SHARED ${uni10-objects}
    src/uni10/numeric/lapack/lib/uni10_dgemm.cu
    src/uni10/numeric/lapack/lib/uni10_lapack_gpu.cu
    src/uni10/tools/lib/uni10_tools_gpu.cu )
 CUDA_ADD_LIBRARY(uni10gpu-static STATIC ${uni10-objects}
    src/uni10/numeric/lapack/lib/uni10_dgemm.cu
    src/uni10/numeric/lapack/lib/uni10_lapack_gpu.cu
    src/uni10/tools/lib/uni10_tools_gpu.cu )
 SET_TARGET_PROPERTIES(uni10gpu-static PROPERTIES OUTPUT_NAME "uni10gpu")
 if(APPLE)
  set(MACOSX_RPATH 1)
  SET_TARGET_PROPERTIES(uni10gpu PROPERTIES VERSION ${UNI10_VERSION} SOVERSION ${UNI10_VERSION_MAJOR} MACOSX_RPATH 1)
 else()
  SET_TARGET_PROPERTIES(uni10gpu PROPERTIES VERSION ${UNI10_VERSION} SOVERSION ${UNI10_VERSION_MAJOR})
 endif()
 CUDA_ADD_CUBLAS_TO_TARGET(uni10gpu)
 CUDA_ADD_CUBLAS_TO_TARGET(uni10gpu-static)
 target_link_libraries(uni10gpu ${CUDA_cusolver_LIBRARY})
 target_link_libraries(uni10gpu-static ${CUDA_cusolver_LIBRARY})
 #target_link_libraries(uni10gpu ${CULA_LIBRARY})
 #target_link_libraries(uni10gpu-static ${CULA_LIBRARY})
 target_link_libraries(uni10gpu ${LAPACK_LIBRARIES})
 target_link_libraries(uni10gpu-static ${LAPACK_LIBRARIES})
 IF(BUILD_HDF5_SUPPORT)
   target_link_libraries(uni10gpu ${HDF5_LIBs})
   target_link_libraries(uni10gpu-static ${HDF5_LIBs})
 ENDIF()
ELSE()
 set(uni10-objects ${uni10-objects} $<TARGET_OBJECTS:uni10-numeric>)
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
 IF( BUILD_ARPACK_SUPPORT )
  target_link_libraries(uni10 ${ARPACK_LIBRARIES})
  target_link_libraries(uni10-static ${ARPACK_LIBRARIES})
 ENDIF()
 IF(BUILD_HDF5_SUPPORT)
   target_link_libraries(uni10 ${HDF5_LIBs})
   target_link_libraries(uni10-static ${HDF5_LIBs})
 ENDIF()
ENDIF()
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
### Build Documents
######################################################################
IF(BUILD_DOC)
  FIND_PACKAGE(Doxygen REQUIRED)
  if(DOXYGEN_FOUND)
    configure_file(Uni10_Doxyfile.in Uni10Doxyfile)
    add_custom_target(doc ALL
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Uni10Doxyfile
     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
     COMMENT "Generating API documentation with Doxygen" VERBATIM
     )
  endif(DOXYGEN_FOUND)
ENDIF()




######################################################################
### INSTALL
######################################################################
IF(BUILD_CUDA_SUPPORT)
### GPU
install(TARGETS uni10gpu uni10gpu-static
        LIBRARY
        DESTINATION lib
        COMPONENT libraries
        ARCHIVE
        DESTINATION lib
        COMPONENT libraries)
ELSE()
### CPU
install(TARGETS uni10 uni10-static
        LIBRARY
        DESTINATION lib
        COMPONENT libraries
        ARCHIVE
        DESTINATION lib
        COMPONENT libraries)
ENDIF()

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
if(DOXYGEN_FOUND)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/
	DESTINATION doc
        COMPONENT documentation
)
endif()
install(FILES INSTALL README GPL LGPL ChangeLog
        DESTINATION .
        COMPONENT common)
######################################################################
### INSTALL TESTS
######################################################################
install(DIRECTORY test/
        DESTINATION test/
        COMPONENT test
	PATTERN "CMakeLists*" EXCLUDE
)



set(libname "${CMAKE_SHARED_LIBRARY_PREFIX}uni10.${UNI10_VERSION}${CMAKE_SHARED_LIBRARY_SUFFIX}")
IF(BUILD_CUDA_SUPPORT)
set(libname "${CMAKE_SHARED_LIBRARY_PREFIX}uni10gpu.${UNI10_VERSION}${CMAKE_SHARED_LIBRARY_SUFFIX}")
ENDIF()



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
message(STATUS " Build Target: ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS " Installation Prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS " CXX Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS " CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS " BLAS and LAPACK Libraries: ${LAPACK_LIBRARIES}")
message(STATUS " ARPACK Libraries: ${ARPACK_LIBRARIES}")
message(STATUS " HDF5 Libraries: ${HDF5_LIBs}")
message(STATUS " HDF5 Includes: ${HDF5_INCLUDE_DIRS}")
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
  message(STATUS " Swig Version      : ${SWIG_VERSION}")
  message(STATUS " Swig Directory    : ${SWIG_DIR}")
else()
  message(STATUS " Build Python Wrapper: NO")
endif()
if(BUILD_CUDA_SUPPORT)
  message( STATUS " CUDA Version: ${CUDA_VERSION_STRING}")
  message( STATUS " CUDA Toolkit Root: ${CUDA_TOOLKIT_ROOT_DIR}")
  message( STATUS " CuSolver library: ${CUDA_cusolver_LIBRARY}")
  #message( STATUS " CULA library: ${CULA_LIBRARY}")
endif()
if(BUILD_DOC)
  message( STATUS " Doxygen: ${DOXYGEN_EXECUTABLE}")
endif()
message(STATUS "")
message(STATUS "")
message(STATUS "")

if (BUILD_PACKAGES)
include("Package.cmake")
include(CPack)
endif()
