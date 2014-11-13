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
#include "stdlib.h"
#ifdef MKL
  #include "mkl.h"
#else
  #include <uni10/numeric/uni10_lapack_wrapper.h>
#endif
#include <string.h>
#include <stdio.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
namespace uni10{
void matrixMul(double* A, double* B, int M, int N, int K, double* C, bool ongpuA, bool ongpuB, bool ongpuC){
	double alpha = 1, beta = 0;
	dgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);
}

void diagMM(double* diag, double* mat, size_t M, size_t N, bool diag_ongpu, bool mat_ongpu){
	for(size_t i = 0; i < M; i++)
		vectorScal(diag[i], &(mat[i * N]), N, false);
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

void vectorExp(double a, double* X, size_t N, bool ongpu){
	for(size_t i = 0; i < N; i++)
		X[i] = exp(a * X[i]);

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
	assert(info == 0);
	lwork = (int)worktest;
	double* work= (double*)malloc(sizeof(double)*lwork);
	dsyev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, work, &lwork, &info);
	assert(info == 0);
	free(work);
}

void matrixSVD(double* Mij_ori, int M, int N, double* U, double* S, double* vT, bool ongpu){
	//Mij = U * S * VT
	double* Mij = (double*)malloc(M * N * sizeof(double));
	memcpy(Mij, Mij_ori, M * N * sizeof(double));
	int min = M < N ? M : N;	//min = min(M,N)
	int ldA = N, ldu = N, ldvT = min;
	//int lwork = 12 * N;
	int lwork = -1;
	double worktest;
	int info;
	dgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, &worktest, &lwork, &info);
	assert(info == 0);
	lwork = (int)worktest;
	double *work = (double*)malloc(lwork*sizeof(double));
	dgesvd((char*)"S", (char*)"S", &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, work, &lwork, &info);
	assert(info == 0);
	free(work);
	free(Mij);
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
	int64_t left = N;
	size_t offset = 0;
	int chunk;
	while(left > 0){
		if(left > INT_MAX)
			chunk = INT_MAX;
		else
			chunk = left;
		sum += dasum(&chunk, X + offset, &inc);
		offset += chunk;
		left -= INT_MAX;
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

bool lanczosEV(double* A, double* psi, size_t dim, int& max_iter, double err_tol, double& eigVal, double* eigVec, bool ongpu){
  int N = dim;
  const int min_iter = 2;
  const double beta_err = 1E-15;
  if(max_iter > N)
    max_iter = N;
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
  while(((e_diff > err_tol && it < max_iter) || it < min_iter) && beta > beta_err){
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
    it++;
    if(it > 1){
      double* z = (double*)malloc(it * it * sizeof(double));
      double* work = (double*)malloc(4 * it * sizeof(double));
      int info;
      memcpy(d, As, it * sizeof(double));
      memcpy(e, Bs, it * sizeof(double));
      dstev((char*)"N", &it, d, e, z, &it, work, &info);
      assert(info == 0);
      double base = fabs(d[0]) > 1 ? fabs(d[0]) : 1;
      e_diff = fabs(d[0] - e0_old) / base;
      e0_old = d[0];
    }
  }
  if(it > 1){
    memcpy(d, As, it * sizeof(double));
    memcpy(e, Bs, it * sizeof(double));
    double* z = (double*)malloc(it * it * sizeof(double));
    double* work = (double*)malloc(4 * it * sizeof(double));
    int info;
    dstev((char*)"V", &it, d, e, z, &it, work, &info);
    assert(info == 0);
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
}
};	/* namespace uni10 */
