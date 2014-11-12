# This is a CMake example for Python
cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})
FIND_PACKAGE(PythonInterp)
set(_CURRENT_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})

### Workaround for non-system Python distributions

if (NOT DEFINED PYTHON_INCLUDE_DIR OR NOT DEFINED PYTHON_LIBRARY)
  if(UNIX)
    if (NOT ${PYTHON_EXECUTABLE} STREQUAL "/usr/bin/python")
        ### Homebrew work around 
       if (${PYTHON_EXECUTABLE} STREQUAL "/usr/local/bin/python" AND APPLE) 
          set(BREW_PREFIX /usr/local/Cellar/python/PYTHON_VERSION_STRING/Frmaeworks/Python.frameworks/Versions/Current) ## Homebrew Python
          message(WARNING "Using Homebrew Python.") 
          find_path(PYTHON_INCLUDE_DIR NAMES Python.h
               PATHS ${BREW_PREFIX}/Headers
               DOC "Location of Python Header files"
               NO_DEFAULT_PATH )
         if (PYTHON_INCLUDE_DIR)
           message(STATUS "Found Python.h in ${PYTHON_INCLUDE_DIR}" )
         else ()
           message(FATAL_ERROR "Python.h not found! Please set PYTHON_INCLUDE_DIR manually.")
         endif()
         find_library(PYTHON_LIBRARY
                    NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}python${_CURRENT_VERSION}${CMAKE_SHARED_LIBRARY_SUFFIX}
                          ${CMAKE_SHARED_LIBRARY_PREFIX}python${CMAKE_SHARED_LIBRARY_SUFFIX}
                    PATHS ${BREW_PREFIX}/lib/
                    DOC "Location of Python Libraries"
               NO_DEFAULT_PATH )
         if (PYTHON_LIBRARY)
           message(STATUS "Found PythonLibs: ${PYTHON_LIBRARY}")
         else ()
           message(FATAL_ERROR "PythonLibs not found! Please set PYTHON_LIBRARY manully.")
         endif()
         if (PYTHON_LIBRARY AND PYTHON_INCLUDE_DIR)
            set(PYTHONLIBS_FOUND)
         endif()
       else() ### Not Homebrew, other Pythons
         STRING(REGEX REPLACE "/bin/python.*$" "" PYTHON_ROOT "${PYTHON_EXECUTABLE}")
         find_path(PYTHON_INCLUDE_DIR NAMES Python.h 
               PATHS ${PYTHON_ROOT}/include/ 
               DOC "Location of Python Header files"
               PATH_SUFFIXES 
                    python${_CURRENT_VERSION}
               NO_DEFAULT_PATH )
         if (PYTHON_INCLUDE_DIR) 
           message(STATUS "Found Python.h in ${PYTHON_INCLUDE_DIR}" )
         else ()
           message(FATAL_ERROR "Python.h not found! Please set PYTHON_INCLUDE_DIR manually.")
         endif()
         find_library(PYTHON_LIBRARY 
                    NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}python${_CURRENT_VERSION}${CMAKE_SHARED_LIBRARY_SUFFIX}
                          ${CMAKE_SHARED_LIBRARY_PREFIX}python${CMAKE_SHARED_LIBRARY_SUFFIX} 
                    PATHS ${PYTHON_ROOT}/lib/
                    DOC "Location of Python Libraries"
                    PATH_SUFFIXES 
                        python${_CURRENT_VERSION}
               NO_DEFAULT_PATH )
         if (PYTHON_LIBRARY) 
           message(STATUS "Found PythonLibs: ${PYTHON_LIBRARY}")
         else ()
           message(FATAL_ERROR "PythonLibs not found! Please set PYTHON_LIBRARY manully.")
         endif()
         if (PYTHON_LIBRARY AND PYTHON_INCLUDE_DIR) 
            set(PYTHONLIBS_FOUND)
         endif()
       endif()
    else()
     FIND_PACKAGE(PythonLibs)
    endif()
  endif(UNIX)
  if (WIN32)
    STRING(REGEX REPLACE "/bin/python.*$" "" PYTHON_ROOT "${PYTHON_EXECUTABLE}")
    find_path(PYTHON_INCLUDE_DIR NAMES Python.h 
          PATHS ${PYTHON_ROOT}/include/ 
          DOC "Location of Python Header files"
          PATH_SUFFIXES 
               python${_CURRENT_VERSION}
          NO_DEFAULT_PATH )
    if (PYTHON_INCLUDE_DIR)
      message(STATUS "Found Python.h in ${PYTHON_INCLUDE_DIR}" )
    else ()
      message(FATAL_ERROR "Python.h not found! Please set PYTHON_INCLUDE_DIR manually.")
    endif()
    find_library(PYTHON_LIBRARY 
               NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}python${_CURRENT_VERSION}${CMAKE_SHARED_LIBRARY_SUFFIX}
                     ${CMAKE_SHARED_LIBRARY_PREFIX}python${CMAKE_SHARED_LIBRARY_SUFFIX} 
               PATHS ${PYTHON_ROOT}/lib/
               DOC "Location of Python Libraries"
               PATH_SUFFIXES 
                   python${_CURRENT_VERSION}
               NO_DEFAULT_PATH )
    if (PYTHON_LIBRARY)
       message(STATUS "Found PythonLibs: ${PYTHON_LIBRARY}")
    else ()
       message(FATAL_ERROR "PythonLibs not found! Please set PYTHON_LIBRARY manully.")
    endif()
    if (PYTHON_LIBRARY AND PYTHON_INCLUDE_DIR)
       set(PYTHONLIBS_FOUND)
    endif()
  else()
    FIND_PACKAGE(PythonLibs)
  endif(WIN32) 
endif()

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(pyUni10.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(pyUni10.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(pyUni10 python pyUni10.i )
SWIG_LINK_LIBRARIES(pyUni10 ${PYTHON_LIBRARY} uni10)

set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR})
if(UNIX) 
execute_process(COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

if(WIN32)
execute_process(COMMAND python.exe -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "PYTHON_SITE_PACKAGES:${PYTHON_SITE_PACKAGES}")
endif()

install(TARGETS _pyUni10 DESTINATION ${PYTHON_SITE_PACKAGES})

install(FILES ${CMAKE_SWIG_OUTDIR}/pyUni10.py DESTINATION ${PYTHON_SITE_PACKAGES})

## Workaround for multiple versions of python on OS X
if (APPLE)
  string(REGEX MATCH "libpython.*.dylib$" LIBPYTHON ${PYTHON_LIBRARY})
  install(CODE "execute_process(COMMAND install_name_tool -change ${LIBPYTHON} ${PYTHON_LIBRARY} ${PYTHON_SITE_PACKAGES}/_pyUni10.so)")
endif()
