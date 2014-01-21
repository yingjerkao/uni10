#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <boost/random.hpp>
using namespace boost;

// Device code
__global__ void VecAdd(float* A, float* B, float* C, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N)
        C[i] = A[i] + B[i];
}
            
// Host code
int main()
{
    int N = 256;
    size_t size = N * sizeof(float);

    // Allocate input vectors h_A and h_B in host memory
    float* h_A = (float*)malloc(size);
    float* h_B = (float*)malloc(size);
    float* h_C = (float*)malloc(size);
    float* hret_C = (float*)malloc(size);

    // Initialize input vectors
	for(int i = 0; i < N; i++){
		h_A[i] = 1;
		h_B[i] = 1;
		hret_C[i] = h_A[i] + h_B[i];
	}

    // Allocate vectors in device memory
    float* d_A;
    assert(cudaMalloc(&d_A, size) == cudaSuccess);
    float* d_B;
    cudaMalloc(&d_B, size);
    float* d_C;
    cudaMalloc(&d_C, size);

    // Copy vectors from host memory to device memory
    cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);
	cudaError_t cuflag;
	assert(cuflag == cudaSuccess);
	

    // Invoke kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    VecAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N);

    // Copy result from device memory to host memory
    // h_C contains the result in host memory
    cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);
	for(int i = 0; i < N; i++){
		printf("%f, %f\n", h_C[i], hret_C[i]);
	}

    // Free device memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
            
    // Free host memory

}
