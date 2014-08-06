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
void matrixMul(double* A, double* B, int M, int N, int K, double* C){
	double alpha = 1, beta = 0;
	dgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);
}

void vectorAdd(double* X, double* Y, size_t N){
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
void vectorScal(double a, double* X, size_t N){
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

void matvecMul(double* A, double* X, int M, int N, double* Y){
    double alpha = 1, beta = 0;
    int incx=1,incy=1;
    
    dgemv((char*)"T", &N, &M, &alpha, A, &N, X, &incx, &beta, Y, &incy);
    
}
/*Generate a set of row vectors which form a othonormal basis
 *For the incoming matrix "elem", the number of row <= the number of column, M <= N
 */
void orthoRandomize(double* elem, int M, int N){
	int eleNum = M*N;
	double *random = (double*)malloc(eleNum * sizeof(double));
	elemRand(random, M * N, 0);
	assert(M <= N);
	int min = M; //min = min(M,N)
	int ldA = M, ldu = M, ldvT = min;
	double *S = (double*)malloc(min*sizeof(double));
	double *u = (double*)malloc(ldu*M*sizeof(double));
	//int lwork = 16*N;
	int lwork = -1;
	double worktest;
	int info;
	dgesvd((char*)"N", (char*)"S", &M, &N, random, &ldA, S, u, &ldu, elem, &ldvT, &worktest, &lwork, &info);
	assert(info == 0);
	lwork = (int)worktest;
	//tmpT = u*S*vT
	double *work = (double*)malloc(lwork*sizeof(double));
	dgesvd((char*)"N", (char*)"S", &M, &N, random, &ldA, S, u, &ldu, elem, &ldvT, work, &lwork, &info);
	//reshape from Fortran format to C format
	memcpy(random, elem, eleNum * sizeof(double));
	setTranspose(random, N, M, elem, 0);
	free(random);
	free(S);
	free(u);
	free(work);
}

void syDiag(double* Kij, int N, double* Eig, double* EigVec){
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

void matrixSVD(double* Mij_ori, int M, int N, double* U, double* S, double* vT){
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

};	/* namespace uni10 */
