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
#include <cula.h>
#include <cublas.h>
namespace uni10{
bool CULAINIT = false;
const size_t GPU_OPERATE_MEM = UNI10_GPU_GLOBAL_MEM / 3;
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

__global__ void _diagRowMul(double* mat, double* diag, size_t M, size_t N){
	size_t idx = blockIdx.y * UNI10_BLOCKMAX * UNI10_THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	double scalar = diag[idx / N];
	if(idx < M * N)
		mat[idx] *= scalar;
}

void diagRowMul(double* mat, double* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){
	double* d_elem = diag;
	size_t d_memsize = M * sizeof(double);
	if(mat_ongpu){
		if(!diag_ongpu){
			d_elem = (double*)elemAllocForce(d_memsize, true);
			elemCopy(d_elem, diag, d_memsize, true, diag_ongpu);
		}
		size_t blockNum = (M * N + UNI10_THREADMAX - 1) / UNI10_THREADMAX;
		dim3 gridSize(blockNum % UNI10_BLOCKMAX, (blockNum + UNI10_BLOCKMAX - 1) / UNI10_BLOCKMAX);
		_diagRowMul<<<gridSize, UNI10_THREADMAX>>>(mat, d_elem, M, N);
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
void diagColMul(double *mat, double* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){
    bool = GPU_READY = false;
    assert(GPU_READY);
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
void vectorMul(double* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu){ // Y = Y * X, element-wise multiplication;
  bool = GPU_READY = false;
  assert(GPU_READY);
}

__global__ void _vectorExp(double a, double* X, size_t N){
	size_t idx = blockIdx.y * UNI10_BLOCKMAX * UNI10_THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	if(idx < N)
		X[idx] = std::exp(a * X[idx]);
}

void vectorExp(double a, double* X, size_t N, bool ongpu){
	if(ongpu){
		size_t blockNum = (N + UNI10_THREADMAX - 1) / UNI10_THREADMAX;
		dim3 gridSize(blockNum % UNI10_BLOCKMAX, (blockNum + UNI10_BLOCKMAX - 1) / UNI10_BLOCKMAX);
		_vectorExp<<<gridSize, UNI10_THREADMAX>>>(a, X, N);
	}
	else
		for(size_t i = 0; i < N; i++)
			X[i] = std::exp(a * X[i]);
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
		assert(culaDeviceDsyev('V', 'U', N, EigVec, ldA, Eig) == culaNoError);
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
		assert(culaDeviceDgesvd('S', 'S', N, M, Mij, ldA, S, vT, ldu, U, ldvT) == culaNoError);
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

void matrixInv(double* A, int N, bool diag bool ongpu){
	if(ongpu){
    bool = GPU_READY = false;
    assert(GPU_READY);
	}
  else{
    if(diag){
      for(int i = 0; i < N; i++)
        A[i] = A[i] == 0 ? 0 : 1/A[i];
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
__global__ void _identity(double* mat, size_t elemNum, size_t col){
	size_t idx = blockIdx.y * UNI10_BLOCKMAX * UNI10_THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	if(idx < elemNum)
		mat[idx * col + idx] = 1;
}

void setIdentity(double* elem, size_t M, size_t N, bool ongpu){
	size_t min;
	if(M < N)	min = M;
	else		min = N;
	size_t blockNum = (min + UNI10_THREADMAX - 1) / UNI10_THREADMAX;
	dim3 gridSize(blockNum % UNI10_BLOCKMAX, (blockNum + UNI10_BLOCKMAX - 1) / UNI10_BLOCKMAX);
	_identity<<<gridSize, UNI10_THREADMAX>>>(elem, min, N);
}


double vectorSum(double* X, size_t N, int inc, bool ongpu){
	if(ongpu){
    bool = GPU_READY = false;
    assert(GPU_READY);
		return cublasDasum(N, X, inc);
	}
	else{
    double sum = 0;
    size_t idx = 0;
    for(size_t i = 0; i < N; i++){
      sum += X[idx];
      idx += inc;
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

void lanczosEV(double* A, double* psi, size_t dim, int& max_iter, double err_tol, double& eigVal, double* eigVec, bool ongpu){
	int N = dim;
	const int min_iter = 2;
	const double beta_err = 1E-15;
	if(max_iter > N)
		max_iter = N;
	assert(max_iter > min_iter);
	double a = 1;
	double alpha;
	double beta = 1;
	int inc = 1;
	size_t M = max_iter;
	double *Vm = (double*)elemAllocForce((M + 1) * N * sizeof(double), ongpu);
	double *As = (double*)elemAllocForce(M * sizeof(double), ongpu);
	double *Bs = (double*)elemAllocForce(M * sizeof(double), ongpu);
	double *d = (double*)elemAllocForce(M * sizeof(double), ongpu);
	double *e = (double*)elemAllocForce(M * sizeof(double), ongpu);
	int it = 0;
	elemCopy(Vm, psi, N * sizeof(double), ongpu, ongpu);
	//memcpy(Vm, psi, N * sizeof(double));
	vectorScal(1 / vectorNorm(psi, N, 1, ongpu), Vm, N, ongpu);
	elemBzero(&Vm[(it+1) * N], N * sizeof(double), ongpu);
	elemBzero(As, M * sizeof(double), ongpu);
	elemBzero(Bs, M * sizeof(double), ongpu);
	double e_diff = 1;
	double e0_old = 0;
	while(((e_diff > err_tol && it < max_iter) || it < min_iter) && beta > beta_err){
		// q1 = Vm[it*N], v = Vm[(it+1) * N], q0 = v
		double minus_beta = -beta;
		//v = A * q1 - beta * q0 = A * q1 - beta * v
		if(ongpu){
			cublasDgemv('T', N, N, a, A, N, &Vm[it * N], inc, minus_beta, &Vm[(it+1) * N], inc);
			alpha = cublasDdot(N, &Vm[it*N], inc, &Vm[(it+1) * N], inc);
			double minus_alpha = -alpha;
			cublasDaxpy(N, minus_alpha, &Vm[it*N], inc, &Vm[(it+1) * N], inc);
		}
		else{
			dgemv((char*)"T", &N, &N, &a, A, &N, &Vm[it * N], &inc, &minus_beta, &Vm[(it+1) * N], &inc);
			alpha = ddot(&N, &Vm[it*N], &inc, &Vm[(it+1) * N], &inc);
			double minus_alpha = -alpha;
			daxpy(&N, &minus_alpha, &Vm[it * N], &inc, &Vm[(it+1) * N], &inc);
		}
		beta = vectorNorm(&Vm[(it+1) * N], N, 1, ongpu);
		if(it < max_iter - 1)
			elemCopy(&Vm[(it + 2) * N], &Vm[it * N], N * sizeof(double), ongpu, ongpu);
		setElemAt(it, alpha, As, ongpu);
		if(beta > beta_err){
			vectorScal(1/beta, &Vm[(it+1) * N], N, ongpu);
			if(it < max_iter - 1)
				setElemAt(it, beta, Bs, ongpu);
		}
		it++;
		if(it > 1){
			double *work = NULL;
			double *z = NULL;
			elemCopy(d, As, it * sizeof(double), ongpu, ongpu);
			elemCopy(e, Bs, it * sizeof(double), ongpu, ongpu);
			if(ongpu){
				culaInit();
				assert(culaDeviceDsteqr('N', it, d, e, z, it) == culaNoError);
			}
			else{
				int info;
				dstev((char*)"N", &it, d, e, z, &it, work, &info);
				assert(info == 0);
			}
			double ev = getElemAt(0, d, ongpu);
			double base = std::abs(ev) > 1 ? std::abs(ev) : 1;
			e_diff = std::abs(ev - e0_old) / base;
			e0_old = ev;
		}
	}
	if(it > 1){
		elemCopy(d, As, it * sizeof(double), ongpu, ongpu);
		elemCopy(e, Bs, it * sizeof(double), ongpu, ongpu);
		double* z = (double*)elemAllocForce(it * it * sizeof(double), ongpu);
		size_t lwork = 4 * it * sizeof(double);
		double* work = (double*)elemAllocForce(lwork, ongpu);
		if(ongpu){
			assert(culaDeviceDsteqr('I', it, d, e, z, it) == culaNoError);
		}
		else{
			int info;
			dstev((char*)"V", &it, d, e, z, &it, work, &info);
			assert(info == 0);
		}
		elemBzero(eigVec, N * sizeof(double), ongpu);

		if(ongpu){
			double *z_H = (double*)malloc(it * sizeof(double));
			elemCopy(z_H, z, it * sizeof(double), false, ongpu);
			for(int k = 0; k < it; k++)
				cublasDaxpy(N, z_H[k], &Vm[k * N], inc, eigVec, inc);
		}
		else{
			for(int k = 0; k < it; k++)
				daxpy(&N, &z[k], &Vm[k * N], &inc, eigVec, &inc);
		}
		max_iter = it;
		eigVal = getElemAt(0, d, ongpu);
		elemFree(z, it * it * sizeof(double), ongpu);
		elemFree(work, lwork, ongpu);
	}
	else{
		max_iter = 1;
		eigVal = 0;
	}
	elemFree(Vm, (M + 1) * N * sizeof(double), ongpu);
	elemFree(As, M * sizeof(double), ongpu);
	elemFree(Bs, M * sizeof(double), ongpu);
	elemFree(d, M * sizeof(double), ongpu);
	elemFree(e, M * sizeof(double), ongpu);
}

};	/* namespace uni10 */
