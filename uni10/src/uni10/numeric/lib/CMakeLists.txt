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

add_library(uni10-numeric SHARED ${numeric_lib_sources})


######################################################################
### INSTALL
######################################################################

install(TARGETS uni10-numeric DESTINATION lib COMPONENT libraries)


