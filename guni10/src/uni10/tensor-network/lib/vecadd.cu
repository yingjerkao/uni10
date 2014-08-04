#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <uni10/tensor-network/Matrix.h>

// Device code
namespace uni10{
__global__ void VecAdd(double* A, double* B, double* C, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N)
        C[i] = A[i] + B[i];
}

void vecAdd(double* h_A, double* h_B, int N, double* h_C){
    size_t size = N * sizeof(double);

    // Allocate vectors in device memory
    double* d_A;
    assert(cudaMalloc(&d_A, size) == cudaSuccess);
    double* d_B;
    cudaMalloc(&d_B, size);
    double* d_C;
    cudaMalloc(&d_C, size);

    // Copy vectors from host memory to device memory
    cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);
	

    // Invoke kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    VecAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N);

    // Copy result from device memory to host memory
    // h_C contains the result in host memory
    cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
            
    // Free host memory
}
}
