#  Copyright Uni10  2014.
cmake_minimum_required(VERSION 2.8.12)


######################################################################
### BUILD SHARED LIBRARY
######################################################################
include_directories(${CMAKE_SOURCE_DIR})
set(EXAMPLES egQ1 egQ2 egB1 egB2 egN1 egM1 egM2 egU1 egU2 egU3)
foreach( EXAMPLE ${EXAMPLES} )
    add_executable(${EXAMPLE} ${EXAMPLE}.cpp)
    target_link_libraries(${EXAMPLE}  ${LAPACK_LIBRIARIES} uni10-static)
    #target_link_libraries(${EXAMPLE}  uni10-static)
endforeach()
######################################################################
### BUILD EXAMPLES
######################################################################
install(TARGETS ${EXAMPLES} DESTINATION examples/bin COMPONENT executables)
install(FILES egN1_network DESTINATION examples/bin COMPONENT executables)

if(WIN32)
  if ("$CMAKE_SYSTEM_PROCESSOR" STREQUAL "x86")
    file(GLOB DEP_LIBS "../dep/x86/lib*")
    foreach(DLIB ${DEP_LIBS})
	  install(FILES ${DLIB} DESTINATION examples/bin COMPONENT executables  )
    endforeach()
  else()
    file(GLOB DEP_LIBS "../dep/x64/lib*")
    foreach(DLIB ${DEP_LIBS})
	  install(FILES ${DLIB} DESTINATION examples/bin COMPONENT executables  )
    endforeach()
  endif()
endif()
