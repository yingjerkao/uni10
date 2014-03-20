######################################################################
### Generate Python Bindings
######################################################################
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
set(MODULE_NAME "uni10")
set(INTERFACE_FILE "uni10.i")
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")
FILE(GLOB_RECURSE headers ${CMAKE_CURRENT_SOURCE_DIR}/../src *.h *.hpp EXCEPT uni10.hpp) 
message("Headers:${headers}") 
FILE(GLOB_RECURSE codes ${CMAKE_CURRENT_SOURCE_DIR}../src *.ccp ) 
message("codes:${codes}") 
SET_SOURCE_FILES_PROPERTIES(uni10.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(uni10.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(uni10 python uni10.i codes)
SWIG_LINK_LIBRARIES(uni10 ${PYTHON_LIBRARIES})

