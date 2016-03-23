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

set(data_structure_lib_sources
  Block.cpp
  BlockReal.cpp
  BlockComplex.cpp
  BlockTools.cpp
  Bond.cpp
)

######################################################################
### BUILD SHARED LIBRARY
######################################################################

add_library(uni10-data-structure OBJECT ${data_structure_lib_sources})


######################################################################
### INSTALL
######################################################################

#install(TARGETS uni10-data-structure DESTINATION lib COMPONENT libraries)


