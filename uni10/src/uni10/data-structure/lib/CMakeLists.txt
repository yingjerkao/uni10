#  Copyright Uni10  2014.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)


######################################################################
### LIST OF FILES
######################################################################

set(data_structure_lib_sources 
  Block.cpp
  Bond.cpp
)

######################################################################
### BUILD SHARED LIBRARY
######################################################################

add_library(uni10-data-structure SHARED ${data_structure_lib_sources})


######################################################################
### INSTALL
######################################################################

install(TARGETS uni10-data-structure DESTINATION lib COMPONENT libraries)


