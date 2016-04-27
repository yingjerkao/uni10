/****************************************************************************
*  @file uni10_lapack_cpu.cpp
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University
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
*  @brief Implementation file for the BLAS and LAPACK wrappers
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <complex>
#ifdef MKL
  #define MKL_Complex8 std::complex<float>
  #define MKL_Complex16 std::complex<double>
  #include "mkl.h"
#else
  #include <uni10/numeric/lapack/uni10_lapack_wrapper.h>
#endif
#include <string.h>
#include <uni10/numeric/lapack/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
#include <iostream>
namespace uni10{
void matrixMul(double* A, double* B, int M, int N, int K, double* C, bool ongpuA, bool ongpuB, bool ongpuC){
	double alpha = 1, beta = 0;
	dgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);
}

void diagRowMul(double* mat, double* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){
	for(size_t i = 0; i < M; i++)
		vectorScal(diag[i], &(mat[i * N]), N, false);
}

void diagColMul(double *mat, double* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){
  for(size_t i = 0; i < M; i++){
    size_t ridx = i * N;
    for(size_t j = 0; j < N; j++)
      mat[ridx + j] *= diag[j];
  }
}

void vectorAdd(double* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu){	// Y = Y + X
	double a = 1.0;
	int inc = 1;
	int64_t left = N;
	size_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		daxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
		offset += chunk;
		left -= INT_MAX;
	}
}
void vectorScal(double a, double* X, size_t N, bool ongpu){
	int inc = 1;
	int64_t left = N;
	size_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		dscal(&chunk, &a, X + offset, &inc);
		offset += chunk;
		left -= INT_MAX;
	}
}

void vectorMul(double* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu){ // Y = Y * X, element-wise multiplication;
  for(size_t i = 0; i < N; i++)
    Y[i] *= X[i];
}

void vectorExp(double a, double* X, size_t N, bool ongpu){
	for(size_t i = 0; i < N; i++)
		X[i] = std::exp(a * X[i]);
}

/*Generate a set of row vectors which form a othonormal basis
 *For the incoming matrix "elem", the number of row <= the number of column, M <= N
 */
void orthoRandomize(double* elem, int M, int N, bool ongpu){
	int eleNum = M*N;
	double *random = (double*)malloc(eleNum * sizeof(double));
	elemRand(random, M * N, false);
	int min = M < N ? M : N;
	double *S = (double*)malloc(min*sizeof(double));
	if(M <= N){
		double *U = (double*)malloc(M * min * sizeof(double));
		matrixSVD(random, M, N, U, S, elem, false);
		free(U);
	}
	else{
		double *VT = (double*)malloc(min * N * sizeof(double));
		matrixSVD(random, M, N, elem, S, VT, false);
		free(VT);
	}
	free(random);
	free(S);
}

void eigDecompose(double* Kij_ori, int N, std::complex<double>* Eig, std::complex<double>* EigVec, bool ongpu){
  std::complex<double> *Kij = (std::complex<double>*) malloc(N * N * sizeof(std::complex<double>));
  elemCast(Kij, Kij_ori, N * N, ongpu, ongpu);
  eigDecompose(Kij, N, Eig, EigVec, ongpu);
  free(Kij);
}

void eigSyDecompose(double* Kij, int N, double* Eig, double* EigVec, bool ongpu){
	memcpy(EigVec, Kij, N * N * sizeof(double));
	int ldA = N;
	int lwork = -1;
	double worktest;
	int info;
	dsyev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, &worktest, &lwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dsyev': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
	lwork = (int)worktest;
	double* work= (double*)malloc(sizeof(double)*lwork);
	dsyev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, work, &lwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dsyev': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
	free(work);
}
// lapack is builded by fortran which is load by column, so we use 
// dorgqr -> lq
// dorglq -> qr
// dorgrq -> ql 
// dorgql -> rq
void matrixQR(double* Mij_ori, int M, int N, double* Q, double* R, bool ongpu){
  assert(M >= N);
  double* Mij = (double*)malloc(N*M*sizeof(double));
  memcpy(Mij, Mij_ori, N*M*sizeof(double));
  double* tau = (double*)malloc(M*sizeof(double));
  int lda = N;
  int lwork = -1;
  double worktestdge;
  double worktestdor;
  int info;
  int K = N;
  dgelqf(&N, &M, Mij, &lda, tau, &worktestdge, &lwork, &info);
  dorglq(&N, &M, &K, Mij, &lda, tau, &worktestdor, &lwork, &info);
  lwork = (int)worktestdge;
  double* workdge = (double*)malloc(lwork*sizeof(double));
  dgelqf(&N, &M, Mij, &lda, tau, workdge, &lwork, &info);
  //getQ
  lwork = (int)worktestdor;
  double* workdor = (double*)malloc(lwork*sizeof(double));
  dorglq(&N, &M, &K, Mij, &lda, tau, workdor, &lwork, &info);
  memcpy(Q, Mij, N*M*sizeof(double));
  //getR
  double alpha = 1, beta = 0;
  dgemm((char*)"N", (char*)"T", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

  free(Mij);
  free(tau);
  free(workdge);
  free(workdor);
}

void matrixRQ(double* Mij_ori, int M, int N, double* Q, double* R, bool ongpu){

  assert(N >= M);
  double* Mij = (double*)malloc(M*N*sizeof(double));
  memcpy(Mij, Mij_ori, M*N*sizeof(double));
  double* tau = (double*)malloc(M*sizeof(double));
  int lda = N;
  int lwork = -1;
  double worktestdge;
  double worktestdor;
  int info;
  int K = M;
  dgeqlf(&N, &M, Mij, &lda, tau, &worktestdge, &lwork, &info);
  dorgql(&N, &M, &K, Mij, &lda, tau, &worktestdor, &lwork, &info);
  lwork = (int)worktestdge;
  double* workdge = (double*)malloc(lwork*sizeof(double));
  dgeqlf(&N, &M, Mij, &lda, tau, workdge, &lwork, &info);
  free(workdge);
  ///getQ
  lwork = (int)worktestdor;
  double* workdor = (double*)malloc(lwork*sizeof(double));
  dorgql(&N, &M, &K, Mij, &lda, tau, workdor, &lwork, &info);
  memcpy(Q, Mij, N*M*sizeof(double));
  //getR
  double alpha = 1, beta = 0;
  dgemm((char*)"T", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, R, &M);

  free(Mij);
  free(tau);
  free(workdge);
  free(workdor);

}

void matrixLQ(double* Mij_ori, int M, int N, double* Q, double* L, bool ongpu){

  assert(N >= M);
  double* Mij = (double*)malloc(M*N*sizeof(double));
  memcpy(Mij, Mij_ori, M*N*sizeof(double));
  double* tau = (double*)malloc(M*sizeof(double));
  int lda = N;
  int lwork = -1;
  double worktestdge;
  double worktestdor;
  int info;
  int K = M;
  dgeqrf(&N, &M, Mij, &lda, tau, &worktestdge, &lwork, &info);
  dorgqr(&N, &M, &K, Mij, &lda, tau, &worktestdor, &lwork, &info);
  lwork = (int)worktestdge;
  double* workdge = (double*)malloc(lwork*sizeof(double));
  dgeqrf(&N, &M, Mij, &lda, tau, workdge, &lwork, &info);
  //getQ
  lwork = (int)worktestdor;
  double* workdor = (double*)malloc(lwork*sizeof(double));
  dorgqr(&N, &M, &K, Mij, &lda, tau, workdor, &lwork, &info);
  memcpy(Q, Mij, N*M*sizeof(double));
  //getR
  double alpha = 1, beta = 0;
  dgemm((char*)"T", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, L, &M);

  free(Mij);
  free(tau);
  free(workdge);
  free(workdor);
}

void matrixQL(double* Mij_ori, int M, int N, double* Q, double* R, bool ongpu){
  assert(M >= N);
  double* Mij = (double*)malloc(N*M*sizeof(double));
  memcpy(Mij, Mij_ori, N*M*sizeof(double));
  double* tau = (double*)malloc(M*sizeof(double));
  int lda = N;
  int lwork = -1;
  double worktestdge;
  double worktestdor;
  int info;
  int K = N;
  dgerqf(&N, &M, Mij, &lda, tau, &worktestdge, &lwork, &info);
  dorgrq(&N, &M, &K, Mij, &lda, tau, &worktestdor, &lwork, &info);
  lwork = (int)worktestdge;
  double* workdge = (double*)malloc(lwork*sizeof(double));
  dgerqf(&N, &M, Mij, &lda, tau, workdge, &lwork, &info);
  //getQ
  lwork = (int)worktestdor;
  double* workdor = (double*)malloc(lwork*sizeof(double));
  dorgrq(&N, &M, &K, Mij, &lda, tau, workdor, &lwork, &info);
  memcpy(Q, Mij, N*M*sizeof(double));
  //getR
  double alpha = 1, beta = 0;
  dgemm((char*)"N", (char*)"T", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

  free(Mij);
  free(tau);
  free(workdge);
  free(workdor);
}

void matrixSVD(double* Mij_ori, int M, int N, double* U, double* S, double* vT, bool ongpu){
	double* Mij = (double*)malloc(M * N * sizeof(double));
	memcpy(Mij, Mij_ori, M * N * sizeof(double));
	int min = std::min(M, N);
	int ldA = N, ldu = N, ldvT = min;
	int lwork = -1;
	double worktest;
	int info;
	dgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, &worktest, &lwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgesvd': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
	lwork = (int)worktest;
	double *work = (double*)malloc(lwork*sizeof(double));
	dgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, work, &lwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgesvd': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
	free(work);
	free(Mij);
}

void matrixInv(double* A, int N, bool diag, bool ongpu){
  if(diag){
    for(int i = 0; i < N; i++)
      A[i] = A[i] == 0 ? 0 : 1.0/A[i];
    return;
  }
  int *ipiv = (int*)malloc((N+1)*sizeof(int));
  int info;
  dgetrf(&N, &N, A, &N, ipiv, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgetrf': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  int lwork = -1;
  double worktest;
  dgetri(&N, A, &N, ipiv, &worktest, &lwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgetri': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  lwork = (int)worktest;
  double *work = (double*)malloc(lwork * sizeof(double));
  dgetri(&N, A, &N, ipiv, work, &lwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgetri': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  free(ipiv);
  free(work);
}


void setTranspose(double* A, size_t M, size_t N, double* AT, bool ongpu, bool ongpuT){
  for(size_t i = 0; i < M; i++)
    for(size_t j = 0; j < N; j++)
      AT[j * M + i] = A[i * N + j];
}

void setTranspose(double* A, size_t M, size_t N, bool ongpu){
  size_t memsize = M * N * sizeof(double);
  double *AT = (double*)malloc(memsize);
  setTranspose(A, M, N, AT, ongpu, ongpu);
  memcpy(A, AT, memsize);
  free(AT);
}

void setCTranspose(double* A, size_t M, size_t N, double* AT, bool ongpu, bool ongpuT){
  setTranspose(A, M, N, AT, ongpu, ongpuT);
}
void setCTranspose(double* A, size_t M, size_t N, bool ongpu){
  setTranspose(A, M, N, ongpu);
}

void setIdentity(double* elem, size_t M, size_t N, bool ongpu){
	size_t min;
	if(M < N)	min = M;
	else		min = N;
	memset(elem, 0, M * N * sizeof(double));
	for(size_t i = 0; i < min; i++)
		elem[i * N + i] = 1;
}

double vectorSum(double* X, size_t N, int inc, bool ongpu){
  double sum = 0;
  size_t idx = 0;
  for(size_t i = 0; i < N; i++){
    sum += X[idx];
    idx += inc;
  }
	return sum;
}

double vectorNorm(double* X, size_t N, int inc, bool ongpu){
	double norm2 = 0;
	double tmp = 0;
	int64_t left = N;
	size_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		tmp = dnrm2(&chunk, X + offset, &inc);
		norm2 += tmp * tmp;
		offset += chunk;
		left -= INT_MAX;
	}
	return sqrt(norm2);
}

bool lanczosEV(double* A, double* psi, size_t dim, size_t& max_iter, double err_tol, double& eigVal, double* eigVec, bool ongpu){
  int N = dim;
  const int min_iter = 2;
  const double beta_err = 1E-15;
  if(!(max_iter > min_iter)){
    std::ostringstream err;
    err<<"Maximum iteration number should be set greater than 2.";
    throw std::runtime_error(exception_msg(err.str()));
  }
  double a = 1;
  double alpha;
  double beta = 1;
  int inc = 1;
  size_t M = max_iter;
  double *Vm = (double*)malloc((M + 1) * N * sizeof(double));
  double *As = (double*)malloc(M * sizeof(double));
  double *Bs = (double*)malloc(M * sizeof(double));
  double *d = (double*)malloc(M * sizeof(double));
  double *e = (double*)malloc(M * sizeof(double));
  int it = 0;
  memcpy(Vm, psi, N * sizeof(double));
  vectorScal(1 / vectorNorm(psi, N, 1, false), Vm, N, false);
  memset(&Vm[(it+1) * N], 0, N * sizeof(double));
  memset(As, 0, M * sizeof(double));
  memset(Bs, 0, M * sizeof(double));
  double e_diff = 1;
  double e0_old = 0;
  bool converged = false;
  while((((e_diff > err_tol) && it < max_iter) || it < min_iter) && beta > beta_err){
    double minus_beta = -beta;
	  dgemv((char*)"T", &N, &N, &a, A, &N, &Vm[it * N], &inc, &minus_beta, &Vm[(it+1) * N], &inc);
    alpha = ddot(&N, &Vm[it*N], &inc, &Vm[(it+1) * N], &inc);
    double minus_alpha = -alpha;
    daxpy(&N, &minus_alpha, &Vm[it * N], &inc, &Vm[(it+1) * N], &inc);

    beta = vectorNorm(&Vm[(it+1) * N], N, 1, false);
		if(it < max_iter - 1)
			memcpy(&Vm[(it + 2) * N], &Vm[it * N], N * sizeof(double));
    As[it] = alpha;
    if(beta > beta_err){
      vectorScal(1/beta, &Vm[(it+1) * N], N, false);
      if(it < max_iter - 1)
        Bs[it] = beta;
    }
    else
      converged = true;
    it++;
    if(it > 1){
      double* z = NULL;//(double*)malloc(it * it * sizeof(double));
      double* work = NULL;//(double*)malloc(4 * it * sizeof(double));
      int info;
      memcpy(d, As, it * sizeof(double));
      memcpy(e, Bs, it * sizeof(double));
      dstev((char*)"N", &it, d, e, z, &it, work, &info);
      if(info != 0){
        std::ostringstream err;
        err<<"Error in Lapack function 'dstev': Lapack INFO = "<<info;
        throw std::runtime_error(exception_msg(err.str()));
      }
      double base = std::abs(d[0]) > 1 ? std::abs(d[0]) : 1;
      e_diff = std::abs(d[0] - e0_old) / base;
      e0_old = d[0];
      if(e_diff <= err_tol)
        converged = true;
    }
  }
  if(it > 1){
    memcpy(d, As, it * sizeof(double));
    memcpy(e, Bs, it * sizeof(double));
    double* z = (double*)malloc(it * it * sizeof(double));
    double* work = (double*)malloc(4 * it * sizeof(double));
    int info;
    dstev((char*)"V", &it, d, e, z, &it, work, &info);
    if(info != 0){
      std::ostringstream err;
      err<<"Error in Lapack function 'dstev': Lapack INFO = "<<info;
      throw std::runtime_error(exception_msg(err.str()));
    }
    memset(eigVec, 0, N * sizeof(double));

    for(int k = 0; k < it; k++){
      daxpy(&N, &z[k], &Vm[k * N], &inc, eigVec, &inc);
    }
    max_iter = it;
    eigVal = d[0];
    free(z), free(work);
  }
  else{
    max_iter = 1;
    eigVal = 0;
  }
  free(Vm), free(As), free(Bs), free(d), free(e);
  return converged;
}

bool lanczosEVL(double* A, double* psi, size_t dim, size_t& max_iter, double err_tol, double& eigVal, double* eigVec, bool ongpu){
  int N = dim;
  const int min_iter = 2;
  const double beta_err = 1E-15;
  if(!(max_iter > min_iter)){
    std::ostringstream err;
    err<<"Maximum iteration number should be set greater than 2.";
    throw std::runtime_error(exception_msg(err.str()));
  }
  double a = 1;
  double alpha;
  double beta = 1;
  int inc = 1;
  size_t M = max_iter;
  double *Vm = (double*)malloc((M + 1) * N * sizeof(double));
  double *As = (double*)malloc(M * sizeof(double));
  double *Bs = (double*)malloc(M * sizeof(double));
  double *d = (double*)malloc(M * sizeof(double));
  double *e = (double*)malloc(M * sizeof(double));
  int it = 0;
  memcpy(Vm, psi, N * sizeof(double));
  vectorScal(1 / vectorNorm(psi, N, 1, false), Vm, N, false);
  memset(&Vm[(it+1) * N], 0, N * sizeof(double));
  memset(As, 0, M * sizeof(double));
  memset(Bs, 0, M * sizeof(double));
  double e_diff = 1;
  double e0_old = 0;
  bool converged = false;
  while((((e_diff > err_tol) && it < max_iter) || it < min_iter) && beta > beta_err){
    double minus_beta = -beta;
    dgemv((char*)"T", &N, &N, &a, A, &N, &Vm[it * N], &inc, &minus_beta, &Vm[(it+1) * N], &inc);
    alpha = ddot(&N, &Vm[it*N], &inc, &Vm[(it+1) * N], &inc);
    double minus_alpha = -alpha;
    daxpy(&N, &minus_alpha, &Vm[it * N], &inc, &Vm[(it+1) * N], &inc);

    beta = vectorNorm(&Vm[(it+1) * N], N, 1, false);
    if(it < max_iter - 1)
      memcpy(&Vm[(it + 2) * N], &Vm[it * N], N * sizeof(double));
    As[it] = alpha;
    if(beta > beta_err){
      vectorScal(1/beta, &Vm[(it+1) * N], N, false);
      if(it < max_iter - 1)
	Bs[it] = beta;
    }
    else
      converged = true;
    it++;
    if(it > 1){
      double* z = NULL;//(double*)malloc(it * it * sizeof(double));
      double* work = NULL;//(double*)malloc(4 * it * sizeof(double));
      int info;
      memcpy(d, As, it * sizeof(double));
      memcpy(e, Bs, it * sizeof(double));
      dstev((char*)"N", &it, d, e, z, &it, work, &info);
      if(info != 0){
	std::ostringstream err;
	err<<"Error in Lapack function 'dstev': Lapack INFO = "<<info;
	throw std::runtime_error(exception_msg(err.str()));
      }
      double base = std::abs(d[0]) > 1 ? std::abs(d[0]) : 1;
      e_diff = std::abs(d[0] - e0_old) / base;
      e0_old = d[0];
      if(e_diff <= err_tol)
	converged = true;
    }
  }
  if(it > 1){
    memcpy(d, As, it * sizeof(double));
    memcpy(e, Bs, it * sizeof(double));
    double* z = (double*)malloc(it * it * sizeof(double));
    double* work = (double*)malloc(4 * it * sizeof(double));
    int info;
    dstev((char*)"V", &it, d, e, z, &it, work, &info);
    if(info != 0){
      std::ostringstream err;
      err<<"Error in Lapack function 'dstev': Lapack INFO = "<<info;
      throw std::runtime_error(exception_msg(err.str()));
    }
    memset(eigVec, 0, N * sizeof(double));

    for(int k = 0; k < it; k++){
      daxpy(&N, &z[k], &Vm[k * N], &inc, eigVec, &inc);
    }
    max_iter = it;
    eigVal = d[0];
    free(z), free(work);
  }
  else{
    max_iter = 1;
    eigVal = 0;
  }
  free(Vm), free(As), free(Bs), free(d), free(e);
  return converged;
}

/***** Complex version *****/
void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, double *S, std::complex<double>* vT, bool ongpu){
	std::complex<double>* Mij = (std::complex<double>*)malloc(M * N * sizeof(std::complex<double>));
	memcpy(Mij, Mij_ori, M * N * sizeof(std::complex<double>));
	int min = std::min(M, N);
	int ldA = N, ldu = N, ldvT = min;
	int lwork = -1;
  std::complex<double> worktest;
	int info;
  double *rwork = (double*) malloc(std::max(1, 5 * min) * sizeof(double));
	zgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, &worktest, &lwork, rwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgesvd': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
	lwork = (int)(worktest.real());
	std::complex<double> *work = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
	zgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, work, &lwork, rwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'zgesvd': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  free(rwork);
	free(work);
	free(Mij);
}
void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, std::complex<double>* S_ori, std::complex<double>* vT, bool ongpu){
	int min = std::min(M, N);
  double* S = (double*)malloc(min * sizeof(double));
  matrixSVD(Mij_ori, M, N, U, S, vT, ongpu);
  elemCast(S_ori, S, min, false, false);
  free(S);
}

void matrixInv(std::complex<double>* A, int N, bool diag, bool ongpu){
  if(diag){
    for(int i = 0; i < N; i++)
      A[i] = std::abs(A[i]) == 0 ? 0.0 : 1.0/A[i];
    return;
  }
  int *ipiv = (int*)malloc((N+1) * sizeof(int));
  int info;
  zgetrf(&N, &N, A, &N, ipiv, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgetrf': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
	int lwork = -1;
  std::complex<double> worktest;
  zgetri(&N, A, &N, ipiv, &worktest, &lwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgetri': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
	lwork = (int)(worktest.real());
  std::complex<double> *work = (std::complex<double>*)malloc(lwork * sizeof(std::complex<double>));
  zgetri(&N, A, &N, ipiv, work, &lwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgetri': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  free(ipiv);
  free(work);
}

double vectorNorm(std::complex<double>* X, size_t N, int inc, bool ongpu){
	double norm2 = 0;
	double tmp = 0;
	int64_t left = N;
	size_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		tmp = dznrm2(&chunk, X + offset, &inc);
		norm2 += tmp * tmp;
		offset += chunk;
		left -= INT_MAX;
	}
	return sqrt(norm2);
}

std::complex<double> vectorSum(std::complex<double>* X, size_t N, int inc, bool ongpu){
  std::complex<double> sum = 0.0;
  size_t idx = 0;
  for(size_t i = 0; i < N; i++){
    sum += X[idx];
    idx += inc;
  }
	return sum;
}
void matrixMul(std::complex<double>* A, std::complex<double>* B, int M, int N, int K, std::complex<double>* C, bool ongpuA, bool ongpuB, bool ongpuC){
  std::complex<double> alpha = 1.0, beta = 0.0;
	zgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);
}

void vectorAdd(std::complex<double>* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu){	// Y = Y + X
  for(size_t i = 0; i < N; i++)
    Y[i] += X[i];
}
void vectorAdd(std::complex<double>* Y, std::complex<double>* X, size_t N, bool y_ongpu, bool x_ongpu){	// Y = Y + X
  std::complex<double> a = 1.0;
	int inc = 1;
	int64_t left = N;
	size_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		zaxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
		offset += chunk;
		left -= INT_MAX;
	}
}

void vectorScal(double a, std::complex<double>* X, size_t N, bool ongpu){
	int inc = 1;
	int64_t left = N;
	size_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		zdscal(&chunk, &a, X + offset, &inc);
		offset += chunk;
		left -= INT_MAX;
	}
}

void vectorScal(const std::complex<double>& a, std::complex<double>* X, size_t N, bool ongpu){
	int inc = 1;
	int64_t left = N;
	size_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		zscal(&chunk, &a, X + offset, &inc);
		offset += chunk;
		left -= INT_MAX;
	}
}
void vectorMul(std::complex<double>* Y, std::complex<double>* X, size_t N, bool y_ongpu, bool x_ongpu){ // Y = Y * X, element-wise multiplication;
  for(size_t i = 0; i < N; i++)
    Y[i] *= X[i];
}

void diagRowMul(std::complex<double>* mat, std::complex<double>* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){
	for(size_t i = 0; i < M; i++)
		vectorScal(diag[i], &(mat[i * N]), N, false);
}

void diagColMul(std::complex<double> *mat, std::complex<double>* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){
  for(size_t i = 0; i < M; i++){
    size_t ridx = i * N;
    for(size_t j = 0; j < N; j++)
      mat[ridx + j] *= diag[j];
  }
}

void vectorExp(double a, std::complex<double>* X, size_t N, bool ongpu){
	for(size_t i = 0; i < N; i++)
		X[i] = std::exp(a * X[i]);
}

void vectorExp(const std::complex<double>& a, std::complex<double>* X, size_t N, bool ongpu){
	for(size_t i = 0; i < N; i++)
		X[i] = std::exp(a * X[i]);
}

void orthoRandomize(std::complex<double> *elem, int M, int N, bool ongpu){
	int eleNum = M*N;
  std::complex<double> *random = (std::complex<double>*)malloc(eleNum * sizeof(std::complex<double>));
	elemRand(random, M * N, false);
	int min = M < N ? M : N;
	double *S = (double*)malloc(min*sizeof(double));
	if(M <= N){
    std::complex<double> *U = (std::complex<double>*)malloc(M * min * sizeof(std::complex<double>));
		matrixSVD(random, M, N, U, S, elem, false);
		free(U);
	}
	else{
		std::complex<double> *VT = (std::complex<double>*)malloc(min * N * sizeof(std::complex<double>));
		matrixSVD(random, M, N, elem, S, VT, false);
		free(VT);
	}
	free(random);
	free(S);
}
void setTranspose(std::complex<double>* A, size_t M, size_t N, std::complex<double>* AT, bool ongpu, bool ongpuT){
	for(size_t i = 0; i < M; i++)
		for(size_t j = 0; j < N; j++)
			AT[j * M + i] = A[i * N + j];
}
void setTranspose(std::complex<double>* A, size_t M, size_t N, bool ongpu){
  size_t memsize = M * N * sizeof(std::complex<double>);
  std::complex<double> *AT = (std::complex<double>*)malloc(memsize);
  setTranspose(A, M, N, AT, ongpu, ongpu);
  memcpy(A, AT, memsize);
  free(AT);
}

void setCTranspose(std::complex<double>* A, size_t M, size_t N, std::complex<double> *AT, bool ongpu, bool ongpuT){
	for(size_t i = 0; i < M; i++)
		for(size_t j = 0; j < N; j++)
			AT[j * M + i] = std::conj(A[i * N + j]);
}
void setCTranspose(std::complex<double>* A, size_t M, size_t N, bool ongpu){
  size_t memsize = M * N * sizeof(std::complex<double>);
  std::complex<double> *AT = (std::complex<double>*)malloc(memsize);
  setCTranspose(A, M, N, AT, ongpu, ongpu);
  memcpy(A, AT, memsize);
  free(AT);
}

void eigDecompose(std::complex<double>* Kij, int N, std::complex<double>* Eig, std::complex<double>* EigVec, bool ongpu){
  size_t memsize = N * N * sizeof(std::complex<double>);
  std::complex<double> *A = (std::complex<double>*) malloc(memsize);
  memcpy(A, Kij, memsize);
  int ldA = N;
  int ldvl = 1;
  int ldvr = N;
  int lwork = -1;
  double *rwork = (double*) malloc(2 * N * sizeof(double));
  std::complex<double> worktest;
  int info;
  zgeev((char*)"N", (char*)"V", &N, A, &ldA, Eig, NULL, &ldvl, EigVec, &ldvr, &worktest, &lwork, rwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'zgeev': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  lwork = (int)worktest.real();
  std::complex<double>* work = (std::complex<double>*)malloc(sizeof(std::complex<double>)*lwork);
  zgeev((char*)"N", (char*)"V", &N, A, &ldA, Eig, NULL, &ldvl, EigVec, &ldvr, work, &lwork, rwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'zgeev': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  free(work);
  free(rwork);
  free(A);
}

void eigSyDecompose(std::complex<double>* Kij, int N, double* Eig, std::complex<double>* EigVec, bool ongpu){
  //eigDecompose(Kij, N, Eig, EigVec, ongpu);
  memcpy(EigVec, Kij, N * N * sizeof(std::complex<double>));
  int ldA = N;
  int lwork = -1;
  std::complex<double> worktest;
  double* rwork = (double*) malloc((3*N+1) * sizeof(double));
  int info;
  zheev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, &worktest, &lwork, rwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'zheev': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  lwork = (int)worktest.real();
  std::complex<double>* work= (std::complex<double>*)malloc(sizeof(std::complex<double>)*lwork);
  zheev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, work, &lwork, rwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'zheev': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  free(work);
  free(rwork);
}

void setConjugate(std::complex<double> *A, size_t N, bool ongpu){
  for(size_t i = 0; i < N; i++)
    A[i] = std::conj(A[i]);
}

void setIdentity(std::complex<double>* elem, size_t M, size_t N, bool ongpu){
	size_t min;
	if(M < N)	min = M;
	else		min = N;
	memset(elem, 0, M * N * sizeof(std::complex<double>));
	for(size_t i = 0; i < min; i++)
		elem[i * N + i] = 1.0;
}

bool lanczosEV(std::complex<double>* A, std::complex<double>* psi, size_t dim, size_t& max_iter, double err_tol, double& eigVal, std::complex<double>* eigVec, bool ongpu){
  int N = dim;
  const int min_iter = 2;
  const double beta_err = 1E-15;
  if(!(max_iter > min_iter)){
    std::ostringstream err;
    err<<"Maximum iteration number should be set greater than 2.";
    throw std::runtime_error(exception_msg(err.str()));
  }
  std::complex<double> a = 1.0;
  std::complex<double> alpha_tmp=0.;
  double alpha;
  double beta = 1;
  int inc = 1;
  size_t M = max_iter;
  std::complex<double> *Vm = (std::complex<double>*)malloc((M + 1) * N * sizeof(std::complex<double>));
  double *As = (double*)malloc(M * sizeof(double));
  double *Bs = (double*)malloc(M * sizeof(double));
  double *d = (double*)malloc(M * sizeof(double));
  double *e = (double*)malloc(M * sizeof(double));
  int it = 0;
  memcpy(Vm, psi, N * sizeof(std::complex<double>));
  vectorScal(1 / vectorNorm(psi, N, 1, false), Vm, N, false);
  memset(&Vm[(it+1) * N], 0, N * sizeof(std::complex<double>));
  memset(As, 0, M * sizeof(double));
  memset(Bs, 0, M * sizeof(double));
  double e_diff = 1;
  double e0_old = 0;
  bool converged = false;
  while((((e_diff > err_tol) && it < max_iter) || it < min_iter) && beta > beta_err){
    std::complex<double> minus_beta = -beta;
	  zgemv((char*)"T", &N, &N, &a, A, &N, &Vm[it * N], &inc, &minus_beta, &Vm[(it+1) * N], &inc);
    zdotc(&alpha_tmp, &N, &Vm[it*N], &inc, &Vm[(it+1) * N], &inc);
    alpha=alpha_tmp.real();
    std::complex<double> minus_alpha = -alpha;
    zaxpy(&N, &minus_alpha, &Vm[it * N], &inc, &Vm[(it+1) * N], &inc);
    beta = vectorNorm(&Vm[(it+1) * N], N, 1, false);
		if(it < max_iter - 1)
			memcpy(&Vm[(it + 2) * N], &Vm[it * N], N * sizeof(std::complex<double>));
    As[it] = alpha;
    if(beta > beta_err){
      vectorScal(1/beta, &Vm[(it+1) * N], N, false);
      if(it < max_iter - 1)
        Bs[it] = beta;
    }
    else
      converged = true;
    it++;
    if(it > 1){
      double* z = NULL;//(double*)malloc(it * it * sizeof(double));
      double* work = NULL;//(double*)malloc(4 * it * sizeof(double));
      int info;
      memcpy(d, As, it * sizeof(double));
      memcpy(e, Bs, it * sizeof(double));
      dstev((char*)"N", &it, d, e, z, &it, work, &info);
      if(info != 0){
        std::ostringstream err;
        err<<"Error in Lapack function 'dstev': Lapack INFO = "<<info;
        throw std::runtime_error(exception_msg(err.str()));
      }
      double base = std::abs(d[0]) > 1 ? std::abs(d[0]) : 1;
      e_diff = std::abs(d[0] - e0_old) / base;
      e0_old = d[0];
      if(e_diff <= err_tol)
        converged = true;
    }
  }
  if(it > 1){
    memcpy(d, As, it * sizeof(double));
    memcpy(e, Bs, it * sizeof(double));
    double* z = (double*)malloc(it * it * sizeof(double));
    double* work = (double*)malloc(4 * it * sizeof(double));
    int info;
    dstev((char*)"V", &it, d, e, z, &it, work, &info);
    if(info != 0){
      std::ostringstream err;
      err<<"Error in Lapack function 'dstev': Lapack INFO = "<<info;
      throw std::runtime_error(exception_msg(err.str()));
    }
    memset(eigVec, 0, N * sizeof(std::complex<double>));
    std::complex<double> cz;
    for(int k = 0; k < it; k++){
      cz = z[k];
      zaxpy(&N, &cz, &Vm[k * N], &inc, eigVec, &inc);
    }
    max_iter = it;
    eigVal = d[0];
    free(z), free(work);
  }
  else{
    max_iter = 1;
    eigVal = 0;
  }
  free(Vm), free(As), free(Bs), free(d), free(e);
  return converged;
}

void matrixQR(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* R, bool ongpu){
  std::complex<double>* Mij = (std::complex<double>*)malloc(N*M*sizeof(std::complex<double>));
  memcpy(Mij, Mij_ori, N*M*sizeof(std::complex<double>));
  std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
  int lda = N;
  int lwork = -1;
  std::complex<double> worktestzge;
  std::complex<double> worktestzun;
  int info;
  int K = N;
  zgelqf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
  zunglq(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
  lwork = (int)worktestzge.real();
  std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
  zgelqf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
  //getQ
  lwork = (int)worktestzun.real();
  std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
  zunglq(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
  memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
  //getR
  std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
  zgemm((char*)"N", (char*)"C", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

  free(Mij);
  free(tau);
  free(workzge);
  free(workzun);
}

void matrixRQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* R, bool ongpu){

  std::complex<double>* Mij = (std::complex<double>*)malloc(M*N*sizeof(std::complex<double>));
  memcpy(Mij, Mij_ori, M*N*sizeof(std::complex<double>));
  std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
  int lda = N;
  int lwork = -1;
  std::complex<double> worktestzge;
  std::complex<double> worktestzun;
  int info;
  int K = M;
  zgeqlf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
  zungql(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
  lwork = (int)worktestzge.real();
  std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
  zgeqlf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
  //getQ
  lwork = (int)worktestzun.real();
  std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
  zungql(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
  memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
  //getR
  std::complex<double> alpha (1.0, 0.0), beta (0.0, 0.0);
  zgemm((char*)"C", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, R, &M);

  free(Mij);
  free(tau);
  free(workzge);
  free(workzun);

}

void matrixLQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* L, bool ongpu){

  std::complex<double>* Mij = (std::complex<double>*)malloc(M*N*sizeof(std::complex<double>));
  memcpy(Mij, Mij_ori, M*N*sizeof(std::complex<double>));
  std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
  int lda = N;
  int lwork = -1;
  std::complex<double> worktestzge;
  std::complex<double> worktestzun;
  int info;
  int K = M;
  zgeqrf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
  zungqr(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
  lwork = (int)worktestzge.real();
  std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
  zgeqrf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
  //getQ
  lwork = (int)worktestzun.real();
  std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
  zungqr(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
  memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
  //getR
  std::complex<double> alpha (1.0, 0.0), beta (0.0, 0.0);
  zgemm((char*)"C", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, L, &M);

  free(Mij);
  free(tau);
  free(workzge);
  free(workzun);
}

void matrixQL(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* L, bool ongpu){
  assert(M >= N);
  std::complex<double>* Mij = (std::complex<double>*)malloc(N*M*sizeof(std::complex<double>));
  memcpy(Mij, Mij_ori, N*M*sizeof(std::complex<double>));
  std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
  int lda = N;
  int lwork = -1;
  std::complex<double> worktestzge;
  std::complex<double> worktestzun;
  int info;
  int K = N;
  zgerqf(&N, &M, Mij, &lda, tau, &worktestzge, &lwork, &info);
  zungrq(&N, &M, &K, Mij, &lda, tau, &worktestzun, &lwork, &info);
  lwork = (int)worktestzge.real();
  std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
  zgerqf(&N, &M, Mij, &lda, tau, workzge, &lwork, &info);
  //getQ
  lwork = (int)worktestzun.real();
  std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
  zungrq(&N, &M, &K, Mij, &lda, tau, workzun, &lwork, &info);
  memcpy(Q, Mij, N*M*sizeof(std::complex<double>));
  //getR
  std::complex<double> alpha (1.0, 0.0), beta (1.0, 1.0);
  zgemm((char*)"N", (char*)"C", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, L, &N);

  free(Mij);
  free(tau);
  free(workzge);
  free(workzun);
}

};	/* namespace uni10 */
