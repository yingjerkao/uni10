#  Copyright Uni10 2013 - 2014.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
######################################################################
### CMAKE VERSION
######################################################################
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)


######################################################################
# Version information
######################################################################
set(UNI10_VERSION_MAJOR 0)
set(UNI10_VERSION_MINOR 0)
set(UNI10_VERSION_PATCH 0)
set(UNI10_VERSION
  ${UNI10_VERSION_MAJOR}.${UNI10_VERSION_MINOR}.${UNI10_VERSION_PATCH}
)


######################################################################
### CMAKE INSTALL OPTIONS 
######################################################################

if (UNIX)
	set(CMAKE_INSTALL_PREFIX 
	  "/opt/uni10/Library/uni10-${UNI10_VERSION}" 
		CACHE PATH "UNI10 install prefix" 
  	FORCE
	)
endif (UNIX)

if (APPLE)
	set(CMAKE_PREFIX_PATH 
		/opt/local 
		${CMAKE_PREFIX_PATH}
	)
endif (APPLE)


######################################################################
### PROJECT
######################################################################

project(uni10)


######################################################################
### CMAKE BUILD TYPE
######################################################################

set(CMAKE_BUILD_TYPE Release)


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
### Check for BOOST
######################################################################

find_package( Boost 1.36.0 )
          if(Boost_FOUND)
             include_directories(${Boost_INCLUDE_DIRS})
          endif()


######################################################################
### FLAGS
######################################################################
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


######################################################################
### ADD HEADER FILES
######################################################################

file(	COPY 					src/
			DESTINATION 	include
			FILES_MATCHING
				PATTERN			"*.h*"
)


######################################################################
### INSTALL HEADER FILES
######################################################################

install(	DIRECTORY 			src/ 
					DESTINATION 		include 
					FILES_MATCHING 
							PATTERN 		"*.hpp"
)


