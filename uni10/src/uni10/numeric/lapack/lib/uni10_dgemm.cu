#include <uni10/numeric/lapack/uni10_lapack.h>
#include "cuda_runtime.h"
#include "cublas_v2.h"
//#include "cusolverDn.h"

namespace uni10{
void getRows(int M, int N, int start, int span, double* iA, double* fA, mmtype how);
void getCols(int M, int N, int start, int span, double* iA, double* fA, mmtype how);
void putBack(int startR, int startC, int spanR, int spanC, int M, int N, double *subA, double *A, mmtype how);
void uni10Dgemm(int p, int q, int M, int N, int K, double* A, double* B, double* C, mmtype how){
	int pM = M / p;
	int qN = N / q;
	int rows, cols;
	rows = pM;// + (M % p);
	cols = qN;// + (N % q);
	double *subA;
	double *subB;
	double *subC;
	double alpha = 1.0;
	double beta = 0.0;
	cublasStatus_t status;
	cublasHandle_t handle;
	status = cublasCreate(&handle);
	if(how & 4)
		assert(cudaMalloc((void**)&subA, rows * K * sizeof(double)) == cudaSuccess);
	if(how & 2 || q > 1)
		assert(cudaMalloc((void**)&subB, K * cols * sizeof(double)) == cudaSuccess);
	if(how & 1 || q > 1)
		assert(cudaMalloc((void**)&subC, rows * cols * sizeof(double)) == cudaSuccess);
	int rChunkNum = (M + pM - 1) / pM;
	int cChunkNum = (N + qN - 1) / qN;
	for(int i = 0; i < rChunkNum; i++){
		for(int j = 0; j < cChunkNum; j++){
			rows = pM; // + (M % p) * ((i + 1) / p);
			cols = qN; // + (N % q) * ((j + 1) / q);
			if(i == rChunkNum - 1 && (M % pM) > 0)
				rows = M % pM ;
			if(j == cChunkNum - 1 && (N % qN) > 0)
				cols = N % qN;
			if(how & 4)
				getRows(M, K, i * pM, rows, A, subA, how);
			else
				subA = A + i * pM * K;
			if(how & 2 || q > 1)
				getCols(K, N, j * qN, cols, B, subB, how);
			else
				subB = B;
			if(how & 1 || q > 1){
				status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, cols, rows, K, &alpha, subB, cols, subA, K, &beta, subC, cols);
				putBack(i * pM, j * qN, rows, cols, M, N, subC, C, how);
			}
			else{
				subC = C + i * pM * N;
				status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, cols, rows, K, &alpha, subB, cols, subA, K, &beta, subC, cols);
			}
			//cublasStatus_t status = cublasGetError();
			assert(status == CUBLAS_STATUS_SUCCESS);
		}
	}
	if(how & 4)
		cudaFree(subA);
	if(how & 2 || q > 1)
		cudaFree(subB);
	if(how & 1 || q > 1)
		cudaFree(subC);
}

void putBack(int startR, int startC, int spanR, int spanC, int M, int N, double *subC, double *C, mmtype how){
	if(how & 1){
		double* tmp = (double*)malloc(spanR * spanC * sizeof(double));
		assert(cudaMemcpy(tmp, subC, sizeof(double) * spanR * spanC, cudaMemcpyDeviceToHost) == cudaSuccess);
		for(int pos = 0; pos < spanR; pos++)
			memcpy(C + (pos + startR) * N + startC, tmp + pos * spanC, spanC * sizeof(double));
		free(tmp);
	}
	else{
		for(int pos = 0; pos < spanR; pos++)
			assert(cudaMemcpy(C + (pos + startR) * N + startC, subC + pos * spanC, spanC * sizeof(double), cudaMemcpyDeviceToDevice) == cudaSuccess);
	}
}

void getRows(int M, int N, int start, int span, double* iA, double* fA, mmtype how){
	assert(start + span <= M);
	if(how & 4){
		assert(cudaMemcpy(fA, iA + start * N, sizeof(double) * N * span, cudaMemcpyHostToDevice) == cudaSuccess);
	}
}

void getCols(int M, int N, int start, int span, double* iB, double* fB, mmtype how){
	assert(start + span <= N);
	if(how & 2){
		double* tmp = (double*)malloc(M * span * sizeof(double));
		for(int i = 0; i < M; i++)
			memcpy(tmp + i * span, iB + i * N + start, span * sizeof(double));
		assert(cudaMemcpy(fB, tmp, sizeof(double) * M * span, cudaMemcpyHostToDevice) == cudaSuccess);
		free(tmp);
	}
	else{
		for(int i = 0; i < M; i++)
			cudaMemcpy(fB + i * span, iB + i * N + start, span * sizeof(double), cudaMemcpyDeviceToDevice);
	}
}
};	/* namespace uni10 */
