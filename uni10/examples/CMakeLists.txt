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

cmake_minimum_required(VERSION 2.8.12)


######################################################################
### BUILD SHARED LIBRARY
######################################################################
include_directories(${CMAKE_SOURCE_DIR})
set(EXAMPLES egQ1 egQ2 egB1 egB2 egN1 egM1 egM2 egU1 egU2 egU3)
foreach( EXAMPLE ${EXAMPLES} )
    add_executable(${EXAMPLE} ${EXAMPLE}.cpp)
    target_link_libraries(${EXAMPLE}  ${LAPACK_LIBRIARIES} uni10-static)
    #target_link_libraries(${EXAMPLE}  uni10-static)
endforeach()
######################################################################
### BUILD EXAMPLES
######################################################################
install(TARGETS ${EXAMPLES} DESTINATION examples/bin COMPONENT examples)
install(FILES egN1_network DESTINATION examples/bin COMPONENT examples)

if(WIN32)
  if ("$CMAKE_SYSTEM_PROCESSOR" STREQUAL "x86")
    file(GLOB DEP_LIBS "../dep/x86/lib*")
    foreach(DLIB ${DEP_LIBS})
	  install(FILES ${DLIB} DESTINATION examples/bin COMPONENT examples  )
    endforeach()
  else()
    file(GLOB DEP_LIBS "../dep/x64/lib*")
    foreach(DLIB ${DEP_LIBS})
	  install(FILES ${DLIB} DESTINATION examples/bin COMPONENT examples  )
    endforeach()
  endif()
endif()
