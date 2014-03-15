#  Copyright Uni10  2014.
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)


######################################################################
### BUILD SHARED LIBRARY
######################################################################
include_directories(${CMAKE_SOURCE_DIR})
  
add_executable(asd.e Ascend.cpp)
add_executable(egQ1.e egQ1.cpp)
target_link_libraries(asd.e ${LAPACK_LIBRARIES} uni10)
target_link_libraries(egQ1.e ${LAPACK_LIBRARIES} uni10)

######################################################################
### BUILD EXAMPLES
######################################################################
install(TARGETS asd.e egQ1.e DESTINATION examples COMPONENT executables)


