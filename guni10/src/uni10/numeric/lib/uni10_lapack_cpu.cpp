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

};	/* namespace uni10 */
