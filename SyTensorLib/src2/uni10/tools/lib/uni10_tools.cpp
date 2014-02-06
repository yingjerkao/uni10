#include <boost/random.hpp>
#include <uni10/tools/uni10_tools.h>
using namespace boost;

mt19937 lapack_rng(777);
uniform_01<mt19937> lapack_uni01_sampler(lapack_rng);

size_t MEM_USAGE = 0;

void* myMalloc(void* ptr, size_t memsize, int& status){
	ptr = malloc(memsize);
	assert(ptr != NULL);
	MEM_USAGE += memsize;
	return ptr;
}

void* myMemcpy(void* des, const void* src, size_t memsize, int des_st, int src_st){
	return memcpy(des, src, memsize);
}

void myFree(void* ptr, size_t memsize, int status){
	free(ptr);
	MEM_USAGE -= memsize;
	ptr = NULL;
}
void membzero(void* ptr, size_t memsize, int status){
	memset(ptr, 0, memsize);
}

void randomNums(double* elem, int N, int status){
	for(int i = 0; i < N; i++)
		elem[i] = lapack_uni01_sampler();
}

