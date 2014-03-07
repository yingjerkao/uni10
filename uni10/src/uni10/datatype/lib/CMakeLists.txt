#  Copyright Uni10  2014.
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)


######################################################################
### LIST OF FILES
######################################################################

set(datatype_lib_sources 
  Qnum.cpp 
)

######################################################################
### BUILD SHARED LIBRARY
######################################################################

add_library(uni10-datatype OBJECT ${datatype_lib_sources})


######################################################################
### INSTALL
######################################################################

#install(TARGETS uni10-datatype DESTINATION lib COMPONENT libraries)


