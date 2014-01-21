#  Copyright Uni10  2014.


######################################################################
### LIST OF FILES
######################################################################

set(tensor_network_lib_sources 
  Bond.cpp
  Network.cpp
  SyTensor.cpp
  SyTensor_elem.cpp
  myTools.cpp
  operate.cpp 
)


######################################################################
### BUILD SHARED LIBRARY
######################################################################

add_library(uni10-tensor-network SHARED ${uni10_tensor_network_lib_sources})


######################################################################
### INSTALL
######################################################################

install(TARGETS uni10-tensor-network DESTINATION lib COMPONENT libraries)


