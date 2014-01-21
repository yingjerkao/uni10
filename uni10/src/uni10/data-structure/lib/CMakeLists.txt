#  Copyright Uni10  2014.


######################################################################
### LIST OF FILES
######################################################################

set(data_structure_lib_sources 
  Block.cpp       
  Matrix.cpp
)

######################################################################
### BUILD SHARED LIBRARY
######################################################################

add_library(uni10-data-structure SHARED ${data_structure_lib_sources})


######################################################################
### INSTALL
######################################################################

install(TARGETS uni10-data-structure DESTINATION lib COMPONENT libraries)


