#  Copyright Uni10  2014.


######################################################################
### LIST OF FILES
######################################################################

set(datatype_lib_sources 
  quantum-state.cpp 
  fermionic-quantum-state.cpp 
)

######################################################################
### BUILD SHARED LIBRARY
######################################################################

add_library(uni10-datatype SHARED ${datatype_lib_sources})


######################################################################
### INSTALL
######################################################################

install(TARGETS uni10-datatype DESTINATION lib COMPONENT libraries)


