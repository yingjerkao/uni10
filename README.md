Uni10 Documentation             {#mainpage}
===================



Introduction                    {#introduction}
============

  Uni10 is an open-source C++ library designed for the development of
tensor network algorithms. Programming tensor network algorithms is 
tedious and  prone to errors.  The task of keeping track of tensor
indices while performing contraction of a complicated tensor network 
can be daunting. It is desirable to have a platform that provides
 bookkeeping capability and optimization.

  This software distinguishes itself from  other available software
solutions by providing the following advantages:

  * Fully implemented in C++.
 
  * Aimed toward applications in tensor network algorithms. 
 
  * Provides basic tensor operations with an easy-to-use interface.
 
  * Provides a Network class to process and store the  details of the 
    graphical representations of the networks.
 
  * Implements a heuristic algorithm to search for an optimal pairwise 
    contraction order based on the available computation and memory
    resources. 
 
  * Provides a collection of Python wrappers which interact with the 
    compiled C++ library to take advantage of  the Python language 
    for better code readability and faster prototyping,  without 
    sacrificing the speed. 
 
  * Provides behind-the-scene optimization and acceleration.



Copyright and Changes
---------------------

  See GPL and LGPL for copyright conditions.

  See [Release Notes](\ref ChangeLog) for release notes and changes.

  
Installation
============

Requirements
------------
  * <a href="http://cmake.org/" target="_blank">cmake</a> version > 2.8.12
  * C++ compiler 
  * BLAS and LAPACK libraries and header files

Build 
-----
To build Un10, follow the following steps:

  1. Create a build directory
  
  2. Use Cmake to generate makefile
  
  3. Build library and exmamples

  4. Install library and examples (May require root access)

Examples
--------

Using system c++, blas and lapack

    > mkdir build
    > cd build
    > cmake </path/to/uni10/>
    > make
    > sudo make install 

The installation path defaults to `/usr/local/uni10`.

To override the default path, use `CMAKE_INSTALL_PREFIX` :

    > cmake -DCMAKE_INSTALL_PREFIX=</installation_path> </path/to/uni10/>

To use MKL and Intel compiler:

    > cmake -DBUILD_WITH_MKL=on -DBUILD_WITH_INTEL_COMPILER=on </path/to/uni10/>

If cmake failes to find blas and lapack, specify the libraries by

    > cmake -DBLAS_LIBRARIES=</path/to/blas> -DLAPACK_LIBRARIES=</path/to/lapack> </path/to/uni10/>


Developers
----------

  * Yun-Da Hsieh (National Taiwan University)

  * Ying-Jer Kao (National Taiwan University)

  * Pochung Chen (National Tsing-Hua University)

  * Tama Ma (Singapore National University)

  * Sukhbinder Singh (Macquarie University)


Help and Bug Reports
--------------------





Known issues
------------


