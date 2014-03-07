#  Copyright Uni10  2014.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)


######################################################################
### LIST OF FILES
######################################################################

set(tools_lib_sources 
  uni10_tools.cpp 
)

######################################################################
### BUILD SHARED LIBRARY
######################################################################

add_library(uni10-tools OBJECT ${tools_lib_sources})


######################################################################
### INSTALL
######################################################################

#install(TARGETS uni10-tools DESTINATION lib COMPONENT libraries)


