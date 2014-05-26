#  Copyright Uni10 2013 - 2014.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
######################################################################
### CMAKE VERSION
######################################################################
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

######################################################################
# Version information
######################################################################
set(UNI10_VERSION_MAJOR 0)
set(UNI10_VERSION_MINOR 1)
set(UNI10_VERSION_PATCH 0)
set(UNI10_VERSION
  ${UNI10_VERSION_MAJOR}.${UNI10_VERSION_MINOR}.${UNI10_VERSION_PATCH}
)



 
######################################################################
### COMMAND LINE OPTIONS
######################################################################
#option(BUILD_STATIC_LIBS "Build the static library" ON)
#option(BUILD_SHARED_LIBS "Build the shared library" ON)
#option(BUILD_TESTS "Build test programs" OFF)
option(USE_MKL "Use the Intel Math Kernal Library." OFF)
option(INTEL_CC "Use the Intel compiler." OFF)

if (INTEL_CC)
  find_program(CMAKE_C_COMPILER NAMES icc)
  find_program(CMAKE_CXX_COMPILER NAMES icpc)
  find_program(CMAKE_AR NAMES xiar)
  find_program(CMAKE_LINKER NAMES xild)
 
  if (CMAKE_C_COMPILER MATCHES CMAKE_C_COMPILER-NOTFOUND OR
    CMAKE_CXX_COMPILER MATCHES CMAKE_CXX_COMPILER-NOTFOUND OR
    CMAKE_AR MATCHES CMAKE_AR-NOTFOUND OR
    CMAKE_LINKER MATCHES CMAKE_LINKER-NOTFOUND)
    message(FATAL_ERROR "Cannot find Intel compiler.  You may need to run `. /opt/intel/bin/compilervars.sh intel64'")
  endif ()
endif (INTEL_CC)

######################################################################
### PROJECT
######################################################################

project(uni10)


######################################################################
### CMAKE INSTALL OPTIONS 
######################################################################
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX "/usr/local/uni10" 
		CACHE PATH "UNI10 install prefix" 
   		FORCE
	)
endif()
######################################################################
### CMAKE BUILD TYPE
######################################################################

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)



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
### Find BLAS and LAPACK
######################################################################

if(USE_MKL)
  find_package( MKL REQUIRED)
  ADD_DEFINITIONS("-DMKL")
  #message( STATUS "MKL found: ${MKL_MINIMAL_LIBRARY}")
else ()
#  find_package( BLAS REQUIRED)
#  message( STATUS "BLAS found: ${BLAS_LIBRARIES}" )
  find_package( LAPACK REQUIRED)
  message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
endif()


######################################################################
### FLAGS
######################################################################
if(UNIX) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif(UNIX)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif()

######################################################################
### PATHS
######################################################################

#set the default path for built executables to the "bin" directory
if (UNIX)
	set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
endif (UNIX)

#set the default path for built libraries to the "lib" directory
if (UNIX)
  set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
endif (UNIX)

# path to include directory
if (UNIX)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
endif (UNIX)


######################################################################
### ADD SUBDIRECTORIES
######################################################################

add_subdirectory(src)
add_subdirectory(examples)
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
add_library(uni10 SHARED ${uni10-objects} )
add_library(uni10-static STATIC ${uni10-objects})
SET_TARGET_PROPERTIES(uni10-static PROPERTIES OUTPUT_NAME "uni10")
if(USE_MKL)
set(LAPACK_LIBRARIES ${MKL_MINIMAL_LIBRARY})
message("MKL: ${MKL_MINIMAL_LIBRARY}")
endif()
target_link_libraries(uni10 ${LAPACK_LIBRARIES})
target_link_libraries(uni10-static ${LAPACK_LIBRARIES})

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
	FILES_MATCHING PATTERN 		"*.h*"	
)

######################################################################
### INSTALL EXAMPLES
######################################################################

install(DIRECTORY examples/
	DESTINATION examples/
	PATTERN "CMakeLists*" EXCLUDE
)
######################################################################
### INSTALL DOCUMENTATION
######################################################################

install(DIRECTORY doc/
	DESTINATION doc )
install(FILES INSTALL.txt README.txt DESTINATION .)
#include(CPACK)
