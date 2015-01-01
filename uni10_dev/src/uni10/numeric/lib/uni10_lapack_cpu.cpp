/****************************************************************************
*  @file uni10_lapack_wrapper.cpp
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
*  @brief Implementation file for the BLAS and LAPACK wrappers
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifdef MKL
  #include "mkl.h"
#else
  #include <uni10/numeric/uni10_lapack_wrapper.h>
#endif
#include <string.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
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

void syDiag(double* Kij, int N, double* Eig, double* EigVec, bool ongpu){
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
  int *ipiv = (int*)malloc(N+1 * sizeof(int));
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

void setTranspose(double* A, size_t M, size_t N, double* AT, bool ongpu){
	for(size_t i = 0; i < M; i++)
		for(size_t j = 0; j < N; j++)
			AT[j * M + i] = A[i * N + j];
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
      double* z = (double*)malloc(it * it * sizeof(double));
      double* work = (double*)malloc(4 * it * sizeof(double));
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
void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, std::complex<double>* S, std::complex<double>* vT, bool ongpu){
	std::complex<double>* Mij = (std::complex<double>*)malloc(M * N * sizeof(std::complex<double>));
	memcpy(Mij, Mij_ori, M * N * sizeof(std::complex<double>));
	int min = std::min(M, N);
	int ldA = N, ldu = N, ldvT = min;
	int lwork = -1;
  std::complex<double> worktest;
	int info;
  double* _S = (double*)malloc(min * sizeof(double));
  double *rwork = (double*) malloc(std::max(1, 5 * min) * sizeof(double));
	zgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, _S, vT, &ldu, U, &ldvT, &worktest, &lwork, rwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'dgesvd': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
	lwork = (int)(worktest.real());
	std::complex<double> *work = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
	zgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, _S, vT, &ldu, U, &ldvT, work, &lwork, rwork, &info);
  if(info != 0){
    std::ostringstream err;
    err<<"Error in Lapack function 'zgesvd': Lapack INFO = "<<info;
    throw std::runtime_error(exception_msg(err.str()));
  }
  for(int i = 0; i < min; i++)
    S[i] = _S[i];
  free(rwork);
  free(_S);
	free(work);
	free(Mij);
}

void matrixInv(std::complex<double>* A, int N, bool diag, bool ongpu){
  if(diag){
    for(int i = 0; i < N; i++)
      A[i] = std::abs(A[i]) == 0 ? 0.0 : 1.0/A[i];
    return;
  }
  int *ipiv = (int*)malloc(N+1 * sizeof(int));
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

void vectorScal(std::complex<double> a, std::complex<double>* X, size_t N, bool ongpu){
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

};	/* namespace uni10 */
