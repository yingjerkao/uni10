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
#include <cula.h>
#include <cublas.h>
namespace uni10{
bool CULAINIT = false;
const size_t GPU_OPERATE_MEM = GPU_GLOBAL_MEM / 3;
void culaInit(){
	if(!CULAINIT)
		culaInitialize();
	CULAINIT = true;
}

void matrixMul(double* A, double* B, int M, int N, int K, double* C, bool ongpuA, bool ongpuB, bool ongpuC){
	double alpha = 1, beta = 0;
	if(ongpuA && ongpuB && ongpuC){
		cublasDgemm('N', 'N', N, M, K, alpha, B, N, A, K, beta, C, N);
	}
	else{
		mmtype types[] = {MM_DDD, MM_DDH, MM_DHD, MM_DHH, MM_HDD, MM_HDH, MM_HHD, MM_HHH};
		int mm_idx = 0;
		int p, q;
		if(!ongpuA)
			mm_idx |= 4;
		if(!ongpuB)
			mm_idx |= 2;
		if(!ongpuC)
			mm_idx |= 1;
		printf("mm_idx = %d\n", mm_idx);
		printf("M = %u, K = %u, N = %u\n", M, K, N);
		mmtype mm_t = types[mm_idx];
		size_t elemPool = GPU_OPERATE_MEM / sizeof(double);
		size_t min_chunk_size = 8;
		int KM_min_ratio = 4;
		if(mm_t == MM_DDH){
			p = ((M * N) + elemPool - 1) / elemPool;
			q = 1;
		}
		else if(mm_t == MM_DHD){
			if(K * N < elemPool){	//allocate K * N
				p = 1;
				q = 1;
			}
			else{	// allocate K * qN + M * qN;
				if(K / M < KM_min_ratio)
					p = (KM_min_ratio * M + K - 1) / K;
				if(M / p < min_chunk_size)
					p = M / min_chunk_size;
				int pM = M / p;
				q = ((K + pM) * N + elemPool - 1) / elemPool;
			}
		}
		else if(mm_t == MM_HDD){
			p = (M * K + elemPool - 1) / elemPool;
			q = 1;
		}
		else if(mm_t == MM_DHH){
			if(K * N + M * N < elemPool){
				p = 1;
				q = 1;
			}
			else{	// The same as MM_DHD
				if(K / M < KM_min_ratio)
					p = (KM_min_ratio * M + K - 1) / K;
				if(M / p < min_chunk_size)
					p = M / min_chunk_size;
				int pM = M / p;
				q = ((K + pM) * N + elemPool - 1) / elemPool;
			}
		}
		else if(mm_t == MM_HDH){
			q = 1;
			p = (M * (K + N) + elemPool - 1) / elemPool;
		}
		else if(mm_t == MM_HHD){
			if((K * N + min_chunk_size * K) < elemPool){
				q = 1;
				size_t elem_left = elemPool - K * N;
				p = (M * K + elem_left - 1) / elem_left;
			}
			else{
				size_t elem_left = elemPool - min_chunk_size * K;
				if(K / M < KM_min_ratio)
					p = (KM_min_ratio * M + K - 1) / K;
				if(M / p < min_chunk_size)
					p = M / min_chunk_size;
				int pM = M / p;
				q = ((K + pM) * N + elem_left - 1) / elem_left;
				int qN = N / q;
				elem_left = elemPool - (K * qN + M * qN);
				p = (M * K + elem_left - 1) / elem_left;
			}
		}
		else{	// MM_HHH
			if((K * N + M * N + min_chunk_size * K) < elemPool){
				q = 1;
				size_t elem_left = elemPool - (K * N + M * N);
				p = (M * K + elem_left - 1) / elem_left;
			}
			else{	// The same as MM_HHD
				size_t elem_left = elemPool - min_chunk_size * K;
				if(K / M < KM_min_ratio)
					p = (KM_min_ratio * M + K - 1) / K;
				if(M / p < min_chunk_size)
					p = M / min_chunk_size;
				int pM = M / p;
				q = ((K + pM) * N + elem_left - 1) / elem_left;
				int qN = N / q;
				elem_left = elemPool - (K * qN + M * qN);
				p = (M * K + elem_left - 1) / elem_left;
			}
		}
		printf("p = %d, q = %d, mm_t = %d\n", p, q, mm_idx);
		uni10Dgemm(p, q, M, N, K, A, B, C, mm_t);
	}	
}

__global__ void _diagMM(double* diag, double* mat, size_t M, size_t N){
	size_t idx = blockIdx.y * BLOCKMAX * THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	double scalar = diag[idx / N];
	if(idx < M * N)
		mat[idx] *= scalar;
}

void diagMM(double* diag, double* mat, size_t M, size_t N, bool diag_ongpu, bool mat_ongpu){
	double* d_elem = diag;
	size_t d_memsize = M * sizeof(double);
	if(mat_ongpu){
		if(!diag_ongpu){
			d_elem = (double*)elemAllocForce(d_memsize, true);
			elemCopy(d_elem, diag, d_memsize, true, diag_ongpu);
		}
		size_t blockNum = (M * N + THREADMAX - 1) / THREADMAX;
		dim3 gridSize(blockNum % BLOCKMAX, (blockNum + BLOCKMAX - 1) / BLOCKMAX);
		_diagMM<<<gridSize, THREADMAX>>>(d_elem, mat, M, N);
		if(!diag_ongpu)
			elemFree(d_elem, d_memsize, true);
	}
	else{
		if(diag_ongpu){
			d_elem = (double*)malloc(d_memsize);
			elemCopy(d_elem, diag, d_memsize, false, diag_ongpu);
		}
		for(size_t i = 0; i < M; i++)
			vectorScal(d_elem[i], &(mat[i * N]), N, mat_ongpu);
		if(diag_ongpu)
			free(d_elem);
	}

}

void vectorAdd(double* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu){	// Y = X + Y
	double a = 1.0;
	int inc = 1;
	if(y_ongpu){
		if(x_ongpu)
			cublasDaxpy(N, a, X, inc, Y, inc);
		else{
			size_t memsize = N * sizeof(double);
			double* elem = (double*)elemAllocForce(memsize, true);
			elemCopy(elem, X, memsize, true, false);
			cublasDaxpy(N, a, elem, inc, Y, inc);
			elemFree(elem, memsize, true);
		}
	}
	else{
		double *elem;
		size_t memsize = N * sizeof(double);
		if(x_ongpu){
			double* elem = (double*)elemAllocForce(memsize, false);
			elemCopy(elem, X, memsize, false, true);
		}
		else
			elem = X;
		int64_t left = N;
		size_t offset = 0;
		int chunk;
		while(left > 0){
			if(left > INT_MAX)
				chunk = INT_MAX;
			else
				chunk = left;
			daxpy(&chunk, &a, elem + offset, &inc, Y + offset, &inc);
			offset += chunk;
			left -= INT_MAX;
		}
		if(x_ongpu)
			elemFree(elem, memsize, false);
	}
}
void vectorScal(double a, double* X, size_t N, bool ongpu){
	int inc = 1;
	if(ongpu)
		cublasDscal(N, a, X, inc);
	else{
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
}

__global__ void _vectorExp(double a, double* X, size_t N){
	size_t idx = blockIdx.y * BLOCKMAX * THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	if(idx < N)
		X[idx] = exp(a * X[idx]);
}

void vectorExp(double a, double* X, size_t N, bool ongpu){
	if(ongpu){
		size_t blockNum = (N + THREADMAX - 1) / THREADMAX;
		dim3 gridSize(blockNum % BLOCKMAX, (blockNum + BLOCKMAX - 1) / BLOCKMAX);
		_vectorExp<<<gridSize, THREADMAX>>>(a, X, N);
	}
	else
		for(size_t i = 0; i < N; i++)
			X[i] = exp(a * X[i]);
}

/*Generate a set of row vectors which form a othonormal basis
 *For the incoming matrix "elem", the number of row <= the number of column, M <= N
 */
void orthoRandomize(double* elem, int M, int N, bool ongpu){
	int eleNum = M*N;
	double *random = (double*)elemAllocForce(eleNum * sizeof(double), ongpu);
	elemRand(random, M * N, ongpu);
	int min = M < N ? M : N;
	double *S = (double*)elemAllocForce(min*sizeof(double), ongpu);
	if(M <= N){
		double *U = (double*)elemAllocForce(M * min * sizeof(double), ongpu);
		matrixSVD(random, M, N, U, S, elem, ongpu);
		elemFree(U, M * min * sizeof(double), ongpu);
	}
	else{
		double *VT = (double*)elemAllocForce(min * N * sizeof(double), ongpu);
		matrixSVD(random, M, N, elem, S, VT, ongpu);
		elemFree(VT, min * N * sizeof(double), ongpu);
	}
	elemFree(random, eleNum * sizeof(double), ongpu);
	elemFree(S, min * sizeof(double), ongpu);
}

void syDiag(double* Kij, int N, double* Eig, double* EigVec, bool ongpu){
	elemCopy(EigVec, Kij, N * N * sizeof(double), ongpu, ongpu);
	int ldA = N;
	if(ongpu){
		culaInit();
		culaDeviceDsyev('V', 'U', N, EigVec, ldA, Eig);	
	}
	else{
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
}

void matrixSVD(double* Mij_ori, int M, int N, double* U, double* S, double* vT, bool ongpu){
	//Mij = U * S * VT
	int min = M < N ? M : N;	//min = min(M,N)
	int ldA = N, ldu = N, ldvT = min;
	if(ongpu){
		size_t memsize = M * N * sizeof(double);
		//assert(cudaMalloc((void**)&Mij, memsize) == cudaSuccess);
		double* Mij = (double*)elemAllocForce(memsize, ongpu);
		assert(cudaMemcpy(Mij, Mij_ori, memsize, cudaMemcpyDeviceToDevice) == cudaSuccess);
		culaInit();
		culaDeviceDgesvd('S', 'S', N, M, Mij, ldA, S, vT, ldu, U, ldvT);
		//cudaFree(Mij);
		elemFree(Mij, memsize, ongpu);
	}
	else{
		double* Mij = (double*)malloc(M * N * sizeof(double));
		memcpy(Mij, Mij_ori, M * N * sizeof(double));
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
}

__global__ void _transpose(double* A, size_t M, size_t N, double* AT){
	size_t y = blockIdx.y * blockDim.y + threadIdx.y;
	size_t x = blockIdx.x * blockDim.x + threadIdx.x;
	if(y < M && x < N)
		AT[x * M + y] = A[y * N + x];
}

void setTranspose(double* A, size_t M, size_t N, double* AT, bool ongpu){
	if(ongpu){
		int thread = 32;
		size_t blockXNum = (N + thread - 1) / thread;
		size_t blockYNum = (M + thread - 1) / thread;
		dim3 blockSize(thread, thread);
		dim3 gridSize(blockXNum, blockYNum);
		_transpose<<<gridSize, blockSize>>>(A, M, N, AT);
	}
	else{
		for(size_t i = 0; i < M; i++)
			for(size_t j = 0; j < N; j++)
				AT[j * M + i] = A[i * N + j];
	}
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
	if(ongpu){
		return cublasDasum(N, X, inc);
	}
	else{
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
}

double vectorNorm(double* X, size_t N, int inc, bool ongpu){
	if(ongpu){
		return cublasDnrm2(N, X, inc);
	}
	else{
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
}

};	/* namespace uni10 */
