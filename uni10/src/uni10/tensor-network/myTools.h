#include <stdint.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
void* myMalloc(void* ptr, size_t memsize, int& status);
void* myMemcpy(void* des, const void* src, size_t memsize, int des_st, int src_st);
void myFree(void* ptr, int status);
void membzero(void* ptr, size_t memsize, int status);
void randomNums(double* elem, int N, int status);
void myTranspose(double* A, int M, int N, double* AT, int status);
void myEye(double* elem, int M, int N, int status);
