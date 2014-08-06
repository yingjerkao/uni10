#  Copyright Uni10  2014.
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)


######################################################################
### BUILD SHARED LIBRARY
######################################################################
include_directories(${CMAKE_SOURCE_DIR})
<<<<<<< HEAD
set(EXAMPLES egQ1 egQ2 egB1 egB2 egM1 egM2 egU1 egU2 egU3)
=======
set(EXAMPLES egQ1 egQ2 egB1 egB2 egN1 egM1 egM2 egU1 egU2 egU3)
>>>>>>> eb8ff2a4ac324948ab953f0306a14dc5c5d03a77
foreach( EXAMPLE ${EXAMPLES} )
    add_executable(${EXAMPLE} ${EXAMPLE}.cpp)
    target_link_libraries(${EXAMPLE} ${LAPACK_LIBRARIES} uni10)
endforeach()
######################################################################
### BUILD EXAMPLES
######################################################################
install(TARGETS ${EXAMPLES} DESTINATION examples/bin COMPONENT executables)
