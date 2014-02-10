#ifndef UNI10_TOOLS_H
#define UNI10_TOOLS_H
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <uni10/data-structure/uni10_struct.h>
namespace uni10{
void* myMalloc(void* ptr, size_t memsize, int& status);
void* myMemcpy(void* des, const void* src, size_t memsize, int des_st, int src_st);
void myFree(void* ptr, size_t memsize, int status);
void membzero(void* ptr, size_t memsize, int status);
void randomNums(double* elem, int N, int status);
std::vector<_Swap> recSwap(int* ord, int n, int* ordF);
std::vector<_Swap> recSwap(int* _ord, int n);	//Given the reshape order out to in. 
};	/* namespace uni10 */	

#endif /* UNI10_TOOLS_H */
