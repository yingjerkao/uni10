Release Notes and Change Log {#ChangeLog}
============================

version 0.9.1
=============

Released on 2015-01-21
This is a maintenance release with bug fixes

ChangeLog
---------

  * Fix doc install problem when no doxygen is found 

version 0.9.0
=============

Released on 2015-01-01 
This release adds Python wrapper pyUni10, API name chaneges, and fixes bugs 

ChangeLog 
---------

  * Add support rank-0 tensor as scalar

  * Can now assign contraction order to a network

  * Add methods combine(), degeneray() and Qlist() to class Bond

  * Add method combineBond() to class UniTensor

  * Add uni10_lapack_wrapper.h to wrap blas/lapack calls to fortran blas/lapack libraries

  * Modify uni10_lapack.cpp to accmmonodate different blas/lapack header files

  * Add new method contract() for UniTensor

  * Fix invalid type conversion bug in Network.cpp

  * Add UniTensor::assign function

  * Updated Copyright files

  * Add Python wrapper

  * Add Python examples

  * Change function name outer to otimes

  * Add Matrix::resize

  * Add Matrix::identity()
  
  * Name change
    + Matrix::diagonalize -> eigh

  * Add methods 
    + UniTensor::getElem();
    + UniTensor::setElem(double*);
    Name Change
    + UniTensor::addLabel -> setLabel
    + UniTensor::addRawElem -> setRawElem

  * Add  functions
    + otimes(const Matrix&, const Matrix&);
    + UniTensor::putBlock(Matrix& mat);
    + UniTensor::getBlock();
  
  * Add Qnum::hash();

  * Add support for all-out-bond tensors

  * Fix bugs in Matrix::trace(), vectorSum()

version 0.1.0b
==============
Released on 2014-05-20 
Beta release of Uni10 with bug fixes


version 0.1.0a
==============
Released on 2014-03-15
Alpha release of Uni10

