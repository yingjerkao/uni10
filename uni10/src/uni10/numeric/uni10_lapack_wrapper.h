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
#include <stdexcept>
#include <cstdint>
#include <complex>
#include <iostream>
extern "C" {
// BLAS functions
void dgemm_(const char *transa, const char *transb, const int32_t *m, const int32_t *n, const int32_t *k,
           const double *alpha, const double *a, const int32_t *lda, const double *b, const int32_t *ldb,
           const double *beta, double *c, const int32_t *ldc);
void zgemm_(const char *transa, const char *transb, const int32_t *m, const int32_t *n, const int32_t *k,
           const std::complex<double> *alpha, const std::complex<double> *a, const int32_t *lda, const std::complex<double> *b, const int32_t *ldb,
           const std::complex<double> *beta, std::complex<double> *c, const int32_t *ldc);

double dasum_(const int32_t *n, const double *x, const int32_t *incx);

void daxpy_(const int32_t *n, const double *alpha, const double *x, const int32_t *incx, double *y, const int32_t *incy);
void zaxpy_(const int32_t *n, const std::complex<double> *alpha, const std::complex<double> *x, const int32_t *incx, std::complex<double> *y, const int32_t *incy);

void dscal_(const int32_t *n, const double *a, double *x, const int32_t *incx);
void zscal_(const int32_t *n, const std::complex<double> *a, std::complex<double> *x, const int32_t *incx);
void zdscal_(const int32_t *n, const double *a, std::complex<double> *x, const int32_t *incx);

double dnrm2_(const int32_t *n, const double *x, const int32_t *incx);
double dznrm2_(const int32_t *n, const std::complex<double> *x, const int32_t *incx);

void dgemv_(const char *trans, const int32_t *m, const int32_t *n, const double *alpha, const double *a, const int32_t *lda, const double *x,
           const int32_t *incx, const double *beta, const double *y, const int32_t *incy);
void zgemv_(const char *trans, const int32_t *m, const int32_t *n, const std::complex<double> *alpha, const std::complex<double> *a, const int32_t *lda,
            const std::complex<double> *x,const int32_t *incx, const std::complex<double> *beta, const std::complex<double> *y, const int32_t *incy);

double ddot_(const int32_t *n, const double *x, const int32_t *incx, const double *y, const int32_t *incy);
void zdotc_(std::complex<double>* res, const int32_t *n, const std::complex<double> *x, const int32_t *incx, const std::complex<double> *y, const int32_t *incy);

// LAPACK functions
void dgesvd_( const char* jobu, const char* jobvt, const int32_t* m,
              const int32_t* n, double* a, const int32_t* lda, double* s,
              double* u, const int32_t* ldu, double* vt, const int32_t* ldvt,
              double* work, const int32_t* lwork, int32_t* info );
void zgesvd_( const char* jobu, const char* jobvt, const int32_t* m,
              const int32_t* n, std::complex<double>* a, const int32_t* lda, double* s,
              std::complex<double>* u, const int32_t* ldu, std::complex<double>* vt, const int32_t* ldvt,
              std::complex<double>* work, const int32_t* lwork, double* rwork, int32_t* info );
void dsyev_( const char* jobz, const char* uplo, const int32_t* n, double* a,
             const int32_t* lda, double* w, double* work, const int32_t* lwork,
             int32_t* info );
void zheev_( const char* jobz, const char* uplo, const int32_t* n, std::complex<double>* a,
             const int32_t* lda, double* w, std::complex<double>* work, const int32_t* lwork,
             const double* rwork, int32_t* info );
void zgeev_( const char* jobvl, const char* jobvr, const int32_t* n, const std::complex<double>* a,
    const int32_t* lda, const std::complex<double>* w, const std::complex<double> *vl, const int32_t *ldvl,
    const std::complex<double> *vr, const int32_t *ldvr, const std::complex<double> *work, const int32_t* lwork,
    const double *rwork, int32_t* info );

void dstev_( const char* jobz, const int32_t* n, const double* d, const double* e, const double* z,
             const int32_t* ldaz, const double* work, int32_t* info );

void dgetrf_( const int32_t *m, const int32_t *n, const double *a,  const int32_t *lda, const int32_t *ipiv, int32_t* info );
void zgetrf_( const int32_t *m, const int32_t *n, const std::complex<double> *a,  const int32_t *lda, const int32_t *ipiv, int32_t* info );
void dgetri_( const int32_t *n, const double *a,  const int32_t *lda, const int32_t *ipiv, const double* work, const int32_t* lwork, int32_t* info );
void zgetri_( const int32_t *n, const std::complex<double> *a, const int32_t *lda, const int32_t *ipiv, const std::complex<double> *work, const int32_t *lwork, int32_t *info );

/*===========================  real  qr rq lq ql  =================================*/
void dgelqf_( const int32_t* m, const int32_t* n, double* a, const int32_t* lda, double* tau, double* work, const int32_t* lwork, int32_t* info );
void dorglq_( const int32_t* m, const int32_t* n, const int32_t* k, double* a, const int32_t* lda, const double* tau, double* work, const int32_t* lwork, int32_t* info );

void dgeqlf_( const int32_t* m, const int32_t* n, double* a, const int32_t* lda, double* tau, double* work, const int32_t* lwork, int32_t* info );
void dorgql_( const int32_t* m, const int32_t* n, const int32_t* k, double* a, const int32_t* lda, const double* tau, double* work, const int32_t* lwork, int32_t* info );

void dgeqrf_( const int32_t* m, const int32_t* n, double* a, const int32_t* lda, double* tau, double* work, const int32_t* lwork, int32_t* info );
void dorgqr_( const int32_t* m, const int32_t* n, const int32_t* k, double* a, const int32_t* lda, const double* tau, double* work, const int32_t* lwork, int32_t* info );

void dgerqf_( const int32_t* m, const int32_t* n, double* a, const int32_t* lda, double* tau, double* work, const int32_t* lwork, int32_t* info );
void dorgrq_( const int32_t* m, const int32_t* n, const int32_t* k, double* a, const int32_t* lda, const double* tau, double* work, const int32_t* lwork, int32_t* info );
/*=================================================================================*/
/*===========================  complex  qr rq lq ql  ==============================*/
void zgeqrf_( const int32_t* m, const int32_t* n, std::complex<double>* a, const int32_t* lda, std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info );
void zungqr_( const int32_t* m, const int32_t* n, const int32_t* k, std::complex<double>* a, const int32_t* lda, const std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info );

void zgerqf_( const int32_t* m, const int32_t* n, std::complex<double>* a, const int32_t* lda, std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info );
void zungrq_( const int32_t* m, const int32_t* n, const int32_t* k, std::complex<double>* a, const int32_t* lda, const std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info );

void zgelqf_( const int32_t* m, const int32_t* n, std::complex<double>* a, const int32_t* lda, std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info );
void zunglq_( const int32_t* m, const int32_t* n, const int32_t* k, std::complex<double>* a, const int32_t* lda, const std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info );

void zgeqlf_( const int32_t* m, const int32_t* n, std::complex<double>* a, const int32_t* lda, std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info );
void zungql_( const int32_t* m, const int32_t* n, const int32_t* k, std::complex<double>* a, const int32_t* lda, const std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info );

/*=================================================================================*/
}


// Wrappers for BLAS and LAPACK functions used in uni10_lapack.cpp
inline void dgemm(const char *transa, const char *transb, const int32_t *m, const int32_t *n, const int32_t *k,
           const double *alpha, const double *a, const int32_t *lda, const double *b, const int32_t *ldb,
           const double *beta, double *c, const int32_t *ldc)
{
  dgemm_(transa, transb, m, n, k,alpha, a, lda, b, ldb, beta, c, ldc);
}

inline void zgemm(const char *transa, const char *transb, const int32_t *m, const int32_t *n, const int32_t *k,
           const std::complex<double> *alpha, const std::complex<double> *a, const int32_t *lda, const std::complex<double> *b, const int32_t *ldb,
           const std::complex<double> *beta, std::complex<double> *c, const int32_t *ldc)
{
  zgemm_(transa, transb, m, n, k,alpha, a, lda, b, ldb, beta, c, ldc);
}

inline double dasum(const int32_t *n, const double *x, const int32_t *incx)
{ return dasum_(n, x, incx); }

inline void daxpy(const int32_t *n, const double *alpha, const double *x, const int32_t *incx, double *y, const int32_t *incy)
{
  daxpy_(n, alpha, x, incx, y, incy);
}

inline void zaxpy(const int32_t *n, const std::complex<double> *alpha, const std::complex<double> *x, const int32_t *incx, std::complex<double> *y, const int32_t *incy)
{
  zaxpy_(n, alpha, x, incx, y, incy);
}

inline double dnrm2(const int32_t *n, const double *x, const int32_t *incx)
{ return dnrm2_(n, x, incx); }

inline double dznrm2(const int32_t *n, const std::complex<double> *x, const int32_t *incx)
{
  return dznrm2_(n, x, incx);
}

inline void dscal(const int32_t *n, const double *a, double *x, const int32_t *incx)
{
  dscal_(n, a, x, incx);
}
inline void zscal(const int32_t *n, const std::complex<double> *a, std::complex<double> *x, const int32_t *incx)
{
  zscal_(n, a, x, incx);
}

inline void zdscal(const int32_t *n, const double *a, std::complex<double> *x, const int32_t *incx)
{
  zdscal_(n, a, x, incx);
}
inline void dsyev( const char* jobz, const char* uplo, const int32_t* n, double* a,
             const int32_t* lda, double* w, double* work, const int32_t* lwork,
             int32_t* info )
{ dsyev_(  jobz,  uplo,  n,  a, lda, w,  work,  lwork, info ); }

inline void zheev( const char* jobz, const char* uplo, const int32_t* n, std::complex<double>* a,
             const int32_t* lda, double* w, std::complex<double>* work, const int32_t* lwork,
             const double* rwork, int32_t* info )
{ zheev_(  jobz,  uplo,  n,  a, lda, w,  work,  lwork, rwork, info ); }

inline void zgeev( const char* jobvl, const char* jobvr, const int32_t* n, const std::complex<double>* a,
    const int32_t* lda, const std::complex<double>* w, const std::complex<double> *vl, const int32_t *ldvl,
    const std::complex<double> *vr, const int32_t *ldvr, const std::complex<double> *work, const int32_t* lwork,
    const double *rwork, int32_t* info )
{
  zgeev_(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

inline void dgesvd( const char* jobu, const char* jobvt, const int32_t* m,
              const int32_t* n, double* a, const int32_t* lda, double* s,
              double* u, const int32_t* ldu, double* vt, const int32_t* ldvt,
              double* work, const int32_t* lwork, int32_t* info )
{
  dgesvd_( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info );
}

inline void zgesvd( const char* jobu, const char* jobvt, const int32_t* m,
              const int32_t* n, std::complex<double>* a, const int32_t* lda, double* s,
              std::complex<double>* u, const int32_t* ldu, std::complex<double>* vt, const int32_t* ldvt,
              std::complex<double>* work, const int32_t* lwork, double* rwork, int32_t* info )
{
  zgesvd_( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info );
}

inline void dgemv(const char *trans, const int32_t *m, const int32_t *n, const double *alpha, const double *a, const int32_t *lda, const double *x,
           const int32_t *incx, const double *beta, const double *y, const int32_t *incy)
{
  dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline void zgemv(const char *trans, const int32_t *m, const int32_t *n, const std::complex<double> *alpha, const std::complex<double> *a, const int32_t *lda,
                  const std::complex<double> *x,const int32_t *incx, const std::complex<double> *beta, const std::complex<double> *y, const int32_t *incy)
{
    zgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

inline double ddot( const int32_t *n, const double *x, const int32_t *incx, const double *y, const int32_t *incy)
{
  return ddot_(n, x, incx, y, incy);
}

inline void zdotc(std::complex<double>* res, const int32_t *n, const std::complex<double> *x, const int32_t *incx, const std::complex<double> *y, const int32_t *incy)
{   
    zdotc_(res, n, x, incx, y, incy);
}

inline void dstev( const char* jobz, const int32_t* n, const double* d, const double* e, const double* z,
             const int32_t* ldaz, const double* work, int32_t* info )
{
  dstev_( jobz, n, d, e, z, ldaz, work, info );
}

inline void dgetrf( const int32_t *m, const int32_t *n, const double *a,  const int32_t *lda, const int32_t *ipiv, int32_t* info )
{
  dgetrf_( m, n, a, lda, ipiv, info );
}
inline void zgetrf( const int32_t *m, const int32_t *n, const std::complex<double> *a,  const int32_t *lda, const int32_t *ipiv, int32_t* info )
{
  zgetrf_(m, n, a, lda, ipiv, info);
}

inline void dgetri( const int32_t *n, const double *a,  const int32_t *lda, const int32_t *ipiv, const double* work, const int32_t* lwork, int32_t* info )
{
  dgetri_(n, a, lda, ipiv, work, lwork, info);
}

inline void zgetri( const int32_t *n, const std::complex<double> *a, const int32_t *lda, const int32_t *ipiv, const std::complex<double> *work, const int32_t *lwork, int32_t *info )
{
  zgetri_(n, a, lda, ipiv, work, lwork, info);
}

/*=================================  qr rq lq ql  =================================*/
inline void dgelqf( const int32_t* m, const int32_t* n, double* a, const int32_t* lda, double* tau, double* work, const int32_t* lwork, int32_t* info )
{
  dgelqf_(m, n, a, lda, tau, work, lwork, info );
}
inline void dorglq( const int32_t* m, const int32_t* n, const int32_t* k, double* a, const int32_t* lda, const double* tau, double* work, const int32_t* lwork, int32_t* info )
{
  dorglq_(m, n, k, a, lda, tau, work, lwork, info );
}
inline void dgeqlf( const int32_t* m, const int32_t* n, double* a, const int32_t* lda, double* tau, double* work, const int32_t* lwork, int32_t* info )
{
  dgeqlf_(m, n, a, lda, tau, work, lwork, info );
}
inline void dorgql( const int32_t* m, const int32_t* n, const int32_t* k, double* a, const int32_t* lda, const double* tau, double* work, const int32_t* lwork, int32_t* info )
{
  dorgql_(m, n, k, a, lda, tau, work, lwork, info );
}
inline void dgeqrf( const int32_t* m, const int32_t* n, double* a, const int32_t* lda, double* tau, double* work, const int32_t* lwork, int32_t* info )
{
  dgeqrf_(m, n, a, lda, tau, work, lwork, info );
}
inline void dorgqr( const int32_t* m, const int32_t* n, const int32_t* k, double* a, const int32_t* lda, const double* tau, double* work, const int32_t* lwork, int32_t* info )
{
  dorgqr_(m, n, k, a, lda, tau, work, lwork, info );
}
inline void dgerqf( const int32_t* m, const int32_t* n, double* a, const int32_t* lda, double* tau, double* work, const int32_t* lwork, int32_t* info )
{
  dgerqf_(m, n, a, lda, tau, work, lwork, info );
}
inline void dorgrq(const int32_t* m, const int32_t* n, const int32_t* k, double* a, const int32_t* lda, const double* tau, double* work, const int32_t* lwork, int32_t* info )
{
  dorgrq_(m, n, k, a, lda, tau, work, lwork, info );
}
/*=================================================================================*/

/*=================================  qr rq lq ql  =================================*/
inline void zgeqrf( const int32_t* m, const int32_t* n, std::complex<double>* a, const int32_t* lda, std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info )
{
  zgeqrf_(m, n, a, lda, tau, work, lwork, info );
}
inline void zungqr( const int32_t* m, const int32_t* n, const int32_t* k, std::complex<double>* a, const int32_t* lda, const std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info )
{
  zungqr_(m, n, k, a, lda, tau, work, lwork, info );
}
inline void zgerqf( const int32_t* m, const int32_t* n, std::complex<double>* a, const int32_t* lda, std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info )
{
  zgerqf_(m, n, a, lda, tau, work, lwork, info );
}
inline void zungrq( const int32_t* m, const int32_t* n, const int32_t* k, std::complex<double>* a, const int32_t* lda, const std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info )
{
  zungrq_(m, n, k, a, lda, tau, work, lwork, info );
}
inline void zgelqf( const int32_t* m, const int32_t* n, std::complex<double>* a, const int32_t* lda, std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info )
{
  zgelqf_(m, n, a, lda, tau, work, lwork, info );
}
inline void zunglq( const int32_t* m, const int32_t* n, const int32_t* k, std::complex<double>* a, const int32_t* lda, const std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info )
{
  zunglq_(m, n, k, a, lda, tau, work, lwork, info );
}
inline void zgeqlf( const int32_t* m, const int32_t* n, std::complex<double>* a, const int32_t* lda, std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info )
{
  zgeqlf_(m, n, a, lda, tau, work, lwork, info );
}
inline void zungql( const int32_t* m, const int32_t* n, const int32_t* k, std::complex<double>* a, const int32_t* lda, const std::complex<double>* tau, std::complex<double>* work, const int32_t* lwork, int32_t* info )
{
  zungql_(m, n, k, a, lda, tau, work, lwork, info );
}
/*=================================================================================*/
#endif
