#  Copyright Uni10  2014.


######################################################################
### LIST OF FILES
######################################################################

set(tensor_network_lib_sources 
  Network.cpp
  UniTensor.cpp
  UniTensor_elem.cpp
  Matrix.cpp
  operate.cpp 
)


######################################################################
### BUILD SHARED LIBRARY
######################################################################

add_library(uni10-tensor-network SHARED ${tensor_network_lib_sources})


######################################################################
### INSTALL
######################################################################

install(TARGETS uni10-tensor-network DESTINATION lib COMPONENT libraries)


