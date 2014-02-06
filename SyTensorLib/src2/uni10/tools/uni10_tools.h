#ifndef UNI10_TOOLS_H
#define UNI10_TOOLS_H
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

void* myMalloc(void* ptr, size_t memsize, int& status);
void* myMemcpy(void* des, const void* src, size_t memsize, int des_st, int src_st);
void myFree(void* ptr, size_t memsize, int status);
void membzero(void* ptr, size_t memsize, int status);
void randomNums(double* elem, int N, int status);
#endif /* UNI10_TOOLS_H */
