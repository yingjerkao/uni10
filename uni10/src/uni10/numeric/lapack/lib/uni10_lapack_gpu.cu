/****************************************************************************
*  @file uni10_lapack_gpu.cpp
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
#ifdef MKL
  #include "mkl.h"
#else
  #include <uni10/numeric/lapack/uni10_lapack_wrapper.h>
#endif
#include <string.h>
#include <uni10/numeric/lapack/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
namespace uni10{

const size_t GPU_OPERATE_MEM = UNI10_GPU_GLOBAL_MEM / 3;

bool IN_MEM(size_t memsize);

bool IN_MEM(size_t memsize){
  size_t free_db, total_db;
  cudaError_t cuda_status;
  cuda_status = cudaMemGetInfo(&free_db, &total_db); 
  if(free_db < memsize && cuda_status == cudaSuccess)
    return false;
  return true;
}

void matrixMul(double* A, double* B, int M, int N, int K, double* C, bool ongpuA, bool ongpuB, bool ongpuC){

  double alpha = 1, beta = 0;
  cublasStatus_t status;
  cublasHandle_t handle;
  status = cublasCreate(&handle);
  if(ongpuA && ongpuB && ongpuC){
    status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, M, K, &alpha, B, N, A, K, &beta, C, N);
    assert(status == CUBLAS_STATUS_SUCCESS);
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
    //printf("mm_idx = %d\n", mm_idx);
    //printf("M = %u, K = %u, N = %u\n", M, K, N);
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
    //printf("p = %d, q = %d, mm_t = %d\n", p, q, mm_idx);
    uni10Dgemm(p, q, M, N, K, A, B, C, mm_t);
  }
}

void vectorAdd(double* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}// Y = Y + X
void vectorScal(double a, double* X, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}// X = a * X

void vectorMul(double* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}// Y = Y * X, element-wise multiplication;
double vectorSum(double* X, size_t N, int inc, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

double vectorNorm(double* X, size_t N, int inc, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void vectorExp(double a, double* X, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void diagRowMul(double* mat, double* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){
	
  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void diagColMul(double* mat, double* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
/*Generate a set of row vectors which form a othonormal basis
 *For the incoming matrix "elem", the number of row <= the number of column, M <= N
 */
void orthoRandomize(double* elem, int M, int N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void eigDecompose(double* Kij, int N, std::complex<double>* Eig, std::complex<double> *EigVec, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void eigSyDecompose(double* Kij, int N, double* Eig, double* EigVec, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixSVD(double* Mij_ori, int M, int N, double* U, double* S, double* vT, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixInv(double* A, int N, bool diag, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void _transposeCPU(double* A, size_t M, size_t N, double* AT){

  for(size_t i = 0; i < M; i++)
    for(size_t j = 0; j < N; j++)
      AT[j * M + i] = A[i * N + j];

}

__global__ void _transposeGPU(double* A, size_t M, size_t N, double* AT){

  size_t y = blockIdx.y * blockDim.y + threadIdx.y;
  size_t x = blockIdx.x * blockDim.x + threadIdx.x;
  if(y < M && x < N)
    AT[x * M + y] = A[y * N + x];

}

void setTranspose(double* A, size_t M, size_t N, double* AT, bool ongpu, bool ongpuT){

  if(ongpu && ongpuT){
    int thread = 32;
    size_t blockXNum = (N + thread - 1) / thread;
    size_t blockYNum = (M + thread - 1) / thread;
    dim3 blockSize(thread, thread);
    dim3 gridSize(blockXNum, blockYNum);
    _transposeGPU<<<gridSize, blockSize>>>(A, M, N, AT);
  }
  else if((!ongpu) && (!ongpuT)){
    _transposeCPU(A, M, N, AT);
  }
  else{
    size_t memsize = M*N*sizeof(double);
    double* h_A = (double*)malloc(memsize);
    cudaMemcpy(h_A, A, memsize, cudaMemcpyDeviceToHost);
    _transposeCPU(h_A, M, N, AT);
    free(h_A);
  }

}

void setTranspose(double* A, size_t M, size_t N, bool ongpu){

  size_t memsize = M * N * sizeof(double);
  double *AT;
  if(ongpu && IN_MEM(memsize)){
    cudaError_t cuflag = cudaMalloc(&AT, memsize);
    assert(cuflag == cudaSuccess);
    setTranspose(A, M, N, AT, ongpu, true);
    elemCopy(A, AT, memsize, ongpu, true);
    cudaFree(AT);
  }else{
    AT = (double*)malloc(memsize);
    setTranspose(A, M, N, AT, ongpu, false);
    elemCopy(A, AT, memsize, ongpu, false);
    free(AT);
  }

}

void setCTranspose(double* A, size_t M, size_t N, double *AT, bool ongpu, bool ongpuT){
  // conj = trans in real 
  setTranspose(A, M, N, AT, ongpu, ongpuT);
}

void setCTranspose(double* A, size_t M, size_t N, bool ongpu){
  // conj = trans in real 
  setTranspose(A, M, N, ongpu);
}

__global__ void _identity(double* mat, size_t elemNum, size_t col){
	size_t idx = blockIdx.y * UNI10_BLOCKMAX * UNI10_THREADMAX +  blockIdx.x * blockDim.x + threadIdx.x;
	if(idx < elemNum)
		mat[idx * col + idx] = 1;
}

void setIdentity(double* elem, size_t M, size_t N, bool ongpu){

  size_t min = std::min(M, N);
  if(ongpu){
    cudaMemset(elem, 0, M * N * sizeof(double));
    size_t blockNum = (min + UNI10_THREADMAX - 1) / UNI10_THREADMAX;
    dim3 gridSize(blockNum % UNI10_BLOCKMAX, (blockNum + UNI10_BLOCKMAX - 1) / UNI10_BLOCKMAX);
    _identity<<<gridSize, UNI10_THREADMAX>>>(elem, min, N);
  }else{
    memset(elem, 0, M * N * sizeof(double));
    for(size_t i = 0; i < min; i++)
      elem[i * N + i] = 1;
  }

}

void reseapeElem(double* elem, size_t* transOffset){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

bool lanczosEV(double* A, double* psi, size_t dim, size_t& max_iter, double err_tol, double& eigVal, double* eigVec, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixQR(double* Mij_ori, int M, int N, double* Q, double* R, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixRQ(double* Mij_ori, int M, int N, double* Q, double* R, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixQL(double* Mij_ori, int M, int N, double* Q, double* L, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixLQ(double* Mij_ori, int M, int N, double* Q, double* L, bool ongpu){

  assert(N >= M);
  if(ongpu){

    cusolverDnHandle_t cusolverHandle = NULL;
    cublasHandle_t cublasHandle = NULL;
    cusolverDnCreate(&cusolverHandle);
    cublasCreate(&cublasHandle);
    // elem copy
    size_t memsize = M * N * sizeof(double);
    double* Mij = NULL;
    cudaError_t cuflag = cudaMalloc(&Mij, memsize);
    assert(cuflag == cudaSuccess);
    cuflag = cudaMemcpy(Mij, Mij_ori, memsize, cudaMemcpyDeviceToDevice);
    assert(cuflag == cudaSuccess);
    // cuda info
    int* info = NULL;
    cuflag = cudaMalloc(&info, sizeof(int));
    assert(cuflag == cudaSuccess);
    cuflag = cudaMemset(info, 0, sizeof(int));
    assert(cuflag == cudaSuccess);
    // cuda workdge
    int lwork = 0;
    double* workdge = NULL;
    double* tau = NULL;
    int lda = N;
    int ldc = N;
    //int K = M;
    cusolverStatus_t cusolverflag = cusolverDnDgeqrf_bufferSize(cusolverHandle, N, M, (double*)Mij, lda, &lwork);
    assert(cusolverflag == CUSOLVER_STATUS_SUCCESS);
    cuflag = cudaMalloc(&workdge, sizeof(double)*lwork);
    assert(cuflag == cudaSuccess);
    cuflag = cudaMalloc((void**)&tau, sizeof(double)*M);
    assert(cuflag == cudaSuccess);

    cusolverflag = cusolverDnDgeqrf(cusolverHandle, N, M, Mij, lda, tau, workdge, lwork, info);
    assert(cusolverflag == CUSOLVER_STATUS_SUCCESS);
    // error info to host	
    int h_info = 0;
    cuflag = cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);

    double* elemI;
    cuflag = cudaMalloc(&elemI, N*M*sizeof(double));
    assert(cuflag == cudaSuccess);
    setIdentity(elemI, M, N, true);

    cusolverflag = cusolverDnDormqr(cusolverHandle, CUBLAS_SIDE_LEFT, CUBLAS_OP_N, N, M, M, Mij, lda, tau, elemI, ldc, workdge, lwork, info);
    assert(cusolverflag == CUSOLVER_STATUS_SUCCESS);

    cuflag = cudaMemcpy(Q, elemI, memsize, cudaMemcpyDeviceToDevice);
    assert(cuflag == cudaSuccess);

    cuflag = cudaDeviceSynchronize();
    assert(cuflag == cudaSuccess);

    double alpha = 1, beta = 0;
    cublasStatus_t cublasflag = cublasDgemm(cublasHandle, CUBLAS_OP_T, CUBLAS_OP_N, M, M, N, &alpha, elemI, N, Mij_ori, N, &beta, L, M);
    assert(cublasflag == CUBLAS_STATUS_SUCCESS);
    
    cudaFree(elemI);
    cudaFree(Mij);
    cudaFree(tau);
    cudaFree(workdge);
    cudaFree(info);

  }
  else{

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

}

/***** Complex version *****/

void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, double *S, std::complex<double>* vT, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, std::complex<double>* S, std::complex<double>* vT, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
void matrixInv(std::complex<double>* A, int N, bool diag, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

std::complex<double> vectorSum(std::complex<double>* X, size_t N, int inc, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
double vectorNorm(std::complex<double>* X, size_t N, int inc, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
void matrixMul(std::complex<double>* A, std::complex<double>* B, int M, int N, int K, std::complex<double>* C, bool ongpuA, bool ongpuB, bool ongpuC){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
void vectorAdd(std::complex<double>* Y, double* X, size_t N, bool y_ongpu, bool x_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}// Y = Y + X
void vectorAdd(std::complex<double>* Y, std::complex<double>* X, size_t N, bool y_ongpu, bool x_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}// Y = Y + X
void vectorScal(double a, std::complex<double>* X, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));
	
}// X = a * X
void vectorScal(const std::complex<double>& a, std::complex<double>* X, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}// X = a * X

void vectorMul(std::complex<double>* Y, std::complex<double>* X, size_t N, bool y_ongpu, bool x_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

} // Y = Y * X, element-wise multiplication;
void diagRowMul(std::complex<double>* mat, std::complex<double>* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void diagColMul(std::complex<double>* mat, std::complex<double>* diag, size_t M, size_t N, bool mat_ongpu, bool diag_ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void vectorExp(double a, std::complex<double>* X, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));
	
}

void vectorExp(const std::complex<double>& a, std::complex<double>* X, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void orthoRandomize(std::complex<double>* elem, int M, int N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
void setTranspose(std::complex<double>* A, size_t M, size_t N, std::complex<double>* AT, bool ongpu, bool ongpuT){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
void setTranspose(std::complex<double>* A, size_t M, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void setCTranspose(std::complex<double>* A, size_t M, size_t N, std::complex<double>* AT, bool ongpu, bool ongpuT){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void setCTranspose(std::complex<double>* A, size_t M, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
void eigDecompose(std::complex<double>* Kij, int N, std::complex<double>* Eig, std::complex<double> *EigVec, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void eigSyDecompose(std::complex<double>* Kij, int N, double* Eig, std::complex<double>* EigVec, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
void setConjugate(std::complex<double> *A, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void setIdentity(std::complex<double>* elem, size_t M, size_t N, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

bool lanczosEV(std::complex<double>* A, std::complex<double>* psi, size_t dim, size_t& max_iter, double err_tol, double& eigVal, std::complex<double>* eigVec, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

bool lanczosEVL(std::complex<double>* A, std::complex<double>* psi, size_t dim, size_t& max_iter, double err_tol, double& eigVal, std::complex<double>* eigVec, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixQR(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* R, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}
void matrixRQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* R, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixQL(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* L, bool ongpu){

  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

void matrixLQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* L, bool ongpu){
	
  std::ostringstream err;
  err<<"GPU version is not ready !!!!";
  throw std::runtime_error(exception_msg(err.str()));

}

};	/* namespace uni10 */
