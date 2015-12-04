###
#  @file CMakeLists.txt
#  @license
#    Copyright (c) 2013-2014
#    Natioanl Taiwan University
#    Natioanl Tsing-Hua University
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


FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})
FIND_PACKAGE(PythonInterp)
set(_CURRENT_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})

### Workaround for non-system Python distributions

if (NOT DEFINED PYTHON_INCLUDE_DIR OR NOT DEFINED PYTHON_LIBRARY)
  if(UNIX)
      if (NOT ${PYTHON_EXECUTABLE} MATCHES "^(/usr)?/bin/python")
        ### Homebrew work around 
        if(APPLE)
          if (${PYTHON_EXECUTABLE} STREQUAL "/usr/local/bin/python") 
            set(BREW_PREFIX 
     /usr/local/Cellar/python/PYTHON_VERSION_STRING/Frmaeworks/Python.frameworks/Versions/Current) 
           ## Homebrew Python
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
         endif()
       endif(APPLE)
       if (NOT PYTHONLIBS_FOUND)
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
if ( PYTHON_VERSION_STRING AND PYTHONLIBS_VERSION_STRING ) 
   if ( NOT PYTHON_VERSION_STRING VERSION_EQUAL 
PYTHONLIBS_VERSION_STRING ) 
     message ( FATAL_ERROR 
             "Version mismatch between python interpreter and libraries" 
) 
   endif ()
endif ( )

set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(pyUni10.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(pyUni10.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(pyUni10 python pyUni10.i )
SWIG_LINK_LIBRARIES(pyUni10 ${PYTHON_LIBRARY} uni10-static)
#SET_SOURCE_FILES_PROPERTIES(${swig_generated_file_fullname} PROPERTIES COMPILE_FLAGS "-Wl,rpath,@loader_path/../../../${PYTHON_LIBRARY}")
if (WIN32)
  if ("$CMAKE_SYSTEM_PROCESSOR" STREQUAL "x86")
    file(GLOB DEP_LIBS "../dep/x86/lib*")
    foreach(DLIB ${DEP_LIBS})
       innstall(FILES ${DLIB} DESTINATION pyUni10
           COMPONENT pyUni10  )
    endforeach()
  else()
    file(GLOB DEP_LIBS "../dep/x64/lib*")
    foreach(DLIB ${DEP_LIBS})
       install(FILES ${DLIB} DESTINATION pyUni10
           COMPONENT pyUni10 )
    endforeach()
  endif()
endif()

configure_file(setup.py.in setup.py)
configure_file(__init__.py.in __init__.py)

#add_custom_command(TARGET _pyUni10
#                   POST_BUILD
#                   COMMAND ${PYTHON_EXECUTABLE} ARGS setup.py install
#                   WORKING_DIRECTORY ${CMAKE_SWIG_OUTDIR})
if ( APPLE)
  string(REGEX MATCH "libpython.*.dylib$" LIBPYTHON ${PYTHON_LIBRARY})
  install(CODE "execute_process( COMMAND install_name_tool -change ${LIBPYTHON} @loader_path/../../../${LIBPYTHON} ${CMAKE_SWIG_OUTDIR}/_pyUni10.so)"
         COMPONENT pyUni10)
endif()

install(TARGETS _pyUni10 DESTINATION pyUni10 COMPONENT pyUni10 )

install(FILES ${CMAKE_SWIG_OUTDIR}/pyUni10.py 
        DESTINATION pyUni10 COMPONENT pyUni10)

install(FILES ${CMAKE_SWIG_OUTDIR}/__init__.py 
        DESTINATION pyUni10 COMPONENT pyUni10)
install(FILES ${CMAKE_SWIG_OUTDIR}/setup.py 
        DESTINATION pyUni10 COMPONENT pyUni10)




