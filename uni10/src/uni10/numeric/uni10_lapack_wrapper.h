#ifndef UNI10_LAPACK_WRAPPER_H
#define UNI10_LAPACK_WRAPPER_H
#include <stdint.h>
extern "C" {
// BLAS functions
void dgemm_(const char *transa, const char *transb, const int32_t *m, const int32_t *n, const int32_t *k,
           const double *alpha, const double *a, const int32_t *lda, const double *b, const int32_t *ldb,
           const double *beta, double *c, const int32_t *ldc);
void    daxpy_(const int32_t *n, const double *alpha, const double *x, const int32_t *incx, double *y, const int32_t *incy);

void    dscal_(const int32_t *n, const double *a, double *x, const int32_t *incx);
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

inline void    daxpy(const int32_t *n, const double *alpha, const double *x, const int32_t *incx, double *y, const int32_t *incy)
{ daxpy_(n, alpha, x, incx, y, incy); }

inline void    dscal(const int32_t *n, const double *a, double *x, const int32_t *incx)
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
