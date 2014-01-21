#pragma once
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

#ifndef USEGPU
#define USEGPU 0
#endif
#define FERMIONIC 1

#define ELEM_PER_THREAD 8
const int ONGPU = 2048;
void* myMalloc(size_t memsize, int& status);
void* myMemcpy(void* des, const void* src, size_t memsize, int des_st, int src_st);
void myFree(void* ptr, size_t memsize, int status);
void membzero(void* ptr, size_t memsize, int status);
void randomNums(double* elem, int N, int status);
void myTranspose(double* A, int M, int N, double* AT, int status);
void myEye(double* elem, int M, int N, int status);
void reshapeElem(int bondNum, double* des, int des_metasz, int* des_meta, double* src, int src_metasz, int* src_meta, int allocThread, int rsp_infosz, int* rsp_info);

