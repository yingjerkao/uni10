###
#  @file CMakeLists.txt
#  @license
#    Copyright (c) 2013-2014
#    National Taiwan University
#    National Tsing-Hua University
#
#    This file is part of Uni10, the Universal Tensor Network Library.
#
#    Uni10 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Uni10 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
#  @endlicense
#  @brief Specification file for CMake
#  @author Ying-Jer Kao
#  @date 2014-05-06
#  @since 0.1.0
###
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)


######################################################################
### LIST OF FILES
######################################################################
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

################################
# Unit Tests
################################

SET(CMAKE_CXX_FLAGS "-O3 -std=c++11")
set(test_sources testSource/testQnum.cpp)
# Add test cpp file
add_executable( runUnitTests ${test_sources})
# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests pthread gtest gtest_main uni10)
add_test( runUnitTests runUnitTests )


