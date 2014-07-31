/****************************************************************************
*  @file uni10_lapack_wrapper.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao
*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief C wrapper functions for fortran Blas and Lapack libraries
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef UNI10_LAPACK_WRAPPER_H
#define UNI10_LAPACK_WRAPPER_H
#include <stdint.h>
extern "C" {
// BLAS functions
void dgemm_(const char *transa, const char *transb, const int32_t *m, const int32_t *n, const int32_t *k,
           const double *alpha, const double *a, const int32_t *lda, const double *b, const int32_t *ldb,
           const double *beta, double *c, const int32_t *ldc);
double dasum_(const int32_t *n, const double *x, const int32_t *incx);

void daxpy_(const int32_t *n, const double *alpha, const double *x, const int32_t *incx, double *y, const int32_t *incy);

void dscal_(const int32_t *n, const double *a, double *x, const int32_t *incx);
double dnrm2_(const int32_t *n, const double *x, const int32_t *incx);
// LAPACK functions
void dgesvd_( const char* jobu, const char* jobvt, const int32_t* m,
              const int32_t* n, double* a, const int32_t* lda, double* s,
              double* u, const int32_t* ldu, double* vt, const int32_t* ldvt,
              double* work, const int32_t* lwork, int32_t* info );
void dsyev_( const char* jobz, const char* uplo, const int32_t* n, double* a,
             const int32_t* lda, double* w, double* work, const int32_t* lwork,
             int32_t* info );

}
// Wrappers for BLAS and LAPACK functions used in uni10_lapack.cpp
inline void dgemm(const char *transa, const char *transb, const int32_t *m, const int32_t *n, const int32_t *k,
           const double *alpha, const double *a, const int32_t *lda, const double *b, const int32_t *ldb,
           const double *beta, double *c, const int32_t *ldc)
{
  dgemm_(transa, transb, m, n, k,alpha, a, lda, b, ldb, beta, c, ldc);
}

inline double dasum(const int32_t *n, const double *x, const int32_t *incx)
{ return dasum_(n, x, incx); }

inline void daxpy(const int32_t *n, const double *alpha, const double *x, const int32_t *incx, double *y, const int32_t *incy)
{ daxpy_(n, alpha, x, incx, y, incy); }

inline double dnrm2(const int32_t *n, const double *x, const int32_t *incx)
{ return dnrm2_(n, x, incx); }

inline void dscal(const int32_t *n, const double *a, double *x, const int32_t *incx)
{   dscal_(n, a, x, incx);}

inline void dsyev( const char* jobz, const char* uplo, const int32_t* n, double* a,
             const int32_t* lda, double* w, double* work, const int32_t* lwork,
             int32_t* info )
{ dsyev_(  jobz,  uplo,  n,  a, lda, w,  work,  lwork, info ); }

inline void dgesvd( const char* jobu, const char* jobvt, const int32_t* m,
              const int32_t* n, double* a, const int32_t* lda, double* s,
              double* u, const int32_t* ldu, double* vt, const int32_t* ldvt,
              double* work, const int32_t* lwork, int32_t* info )
{ dgesvd_( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info ); }
#endif
