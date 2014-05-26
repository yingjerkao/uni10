###
#  @file CMakeLists.txt
#  @license
#    Copyright (c) 2013-2014
#    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao 
#
#    This file is part of Uni10, the Universal Tensor Network Library.
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
#  Copyright Uni10  2014.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)


######################################################################
### LIST OF FILES
######################################################################

set(numeric_lib_sources 
  uni10_lapack.cpp 
)

######################################################################
### BUILD SHARED LIBRARY
######################################################################

add_library(uni10-numeric OBJECT ${numeric_lib_sources})


######################################################################
### INSTALL
######################################################################

#install(TARGETS uni10-numeric DESTINATION lib COMPONENT libraries)


