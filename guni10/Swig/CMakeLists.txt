###
#  @file CMakeLists.txt
#  @license
#    Copyright (c) 2013-2014
#    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao 
#
#    This file is part of Uni10.
#
#    Uni10 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Uni10 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
#  @endlicense
#  @brief SWIG specification file for CMake 
#  @author Ying-Jer Kao
#  @date 2014-10-20
#  @since 0.9.0
###

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
    STRING(REGEX REPLACE "python.*$" "" PYTHON_ROOT "${PYTHON_EXECUTABLE}")
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
               NAMES python${_CURRENT_VERSION}.lib
                     python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.lib
                     python.lib        
               PATHS ${PYTHON_ROOT}/libs/
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
SWIG_LINK_LIBRARIES(pyUni10 ${PYTHON_LIBRARY} uni10-static)

set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR})
if(UNIX) 
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
 OUTPUT_VARIABLE _SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

if(WIN32)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
 OUTPUT_VARIABLE _SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
set(PYTHON_SITE_PACKAGES "${_SITE_PACKAGES}" CACHE PATH "Python Module Direcotry")
if( "${PYTHON_SITE_PACKAGES}" STREQUAL "" )
  message(FATAL_ERROR "Cannot find Python site package directory. Set  PYTHON_SITE_PACKAGES manully.")
endif()
install(TARGETS _pyUni10 DESTINATION ${PYTHON_SITE_PACKAGES}/pyUni10 COMPONENT pyUni10 )

install(FILES ${CMAKE_SWIG_OUTDIR}/pyUni10.py 
        DESTINATION ${PYTHON_SITE_PACKAGES}/pyUni10 COMPONENT pyUni10)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
        DESTINATION ${PYTHON_SITE_PACKAGES}/pyUni10 COMPONENT pyUni10)
if (WIN32)
  if ("$CMAKE_SYSTEM_PROCESSOR" STREQUAL "x86")
    file(GLOB DEP_LIBS "../dep/x86/lib*")
    foreach(DLIB ${DEP_LIBS})
       install(FILES ${DLIB} DESTINATION ${PYTHON_SITE_PACKAGES}/pyUni10
           COMPONENT pyUni10  )
    endforeach()
  else()
    file(GLOB DEP_LIBS "../dep/x64/lib*")
    foreach(DLIB ${DEP_LIBS})
       install(FILES ${DLIB} DESTINATION ${PYTHON_SITE_PACKAGES}/pyUni10
           COMPONENT pyUni10 )
    endforeach()
  endif()
endif()
## Workaround for multiple versions of python on OS X
if (APPLE)
  string(REGEX MATCH "libpython.*.dylib$" LIBPYTHON ${PYTHON_LIBRARY})
  install(CODE "execute_process(COMMAND install_name_tool -change ${LIBPYTHON} ${PYTHON_LIBRARY} ${PYTHON_SITE_PACKAGES}/_pyUni10.so)")
endif()
