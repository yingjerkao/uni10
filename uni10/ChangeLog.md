Release Notes and Change Log 
============================

version 1.0.0
=============

Released on 2016-06-06
This release adds Complex data type, HDF5 io, and connector to ARPACK.
This release also prepares for a future release with cuda support. 
Bug fixes. 

ChangeLog
---------
  * Added HOSVD

  * Added LQ, QR, RQ, QL

  * Added HDF5 support

  * Added Complex datatype support

version 0.9.1
=============

Released on 2015-01-21
This is a maintenance release with bug fixes

ChangeLog
---------

  * Fixed doc install problem when no doxygen is found 

version 0.9.0
=============

Released on 2015-01-01 
This release adds Python wrapper pyUni10, API name chaneges, and fixes bugs 

ChangeLog 
---------

  * Added support rank-0 tensor as scalar

  * Can now assign contraction order to a network

  * Added methods combine(), degeneray() and Qlist() to class Bond

  * Added method combineBond() to class UniTensor

  * Added uni10_lapack_wrapper.h to wrap blas/lapack calls to fortran blas/lapack libraries

  * Modify uni10_lapack.cpp to accmmonodate different blas/lapack header files

  * Added new method contract() for UniTensor

  * Fix invalid type conversion bug in Network.cpp

  * Added UniTensor::assign function

  * Updated Copyright files

  * Added Python wrapper

  * Added Python examples

  * Change function name outer to otimes

  * Added Matrix::resize

  * Added Matrix::identity()
  
  * Name change
    + Matrix::diagonalize -> eigh

  * Added methods 
    + UniTensor::getElem();
    + UniTensor::setElem(double*);
    Name Change
    + UniTensor::addLabel -> setLabel
    + UniTensor::addRawElem -> setRawElem

  * Added  functions
    + otimes(const Matrix&, const Matrix&);
    + UniTensor::putBlock(Matrix& mat);
    + UniTensor::getBlock();
  
  * Added Qnum::hash();

  * Added support for all-out-bond tensors

  * Fix bugs in Matrix::trace(), vectorSum()

version 0.1.0b
==============
Released on 2014-05-20 
Beta release of Uni10 with bug fixes


version 0.1.0a
==============
Released on 2014-03-15
Alpha release of Uni10

