#include <boost/random.hpp>
#include <uni10/tools/uni10_tools.h>
using namespace boost;

namespace uni10{
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

std::vector<_Swap> recSwap(int* _ord, int n){	//Given the reshape order out to in. 
	int ordF[n];
	for(int i = 0; i < n; i++)
		ordF[i] = i;
	return recSwap(_ord, n, ordF);
}
std::vector<_Swap> recSwap(int* _ord, int n, int* ordF){	//Given the reshape order out to in. 
	int* ord = (int*)malloc(sizeof(int) * n);
	memcpy(ord, _ord, sizeof(int) * n);
	std::vector<_Swap> swaps;
	_Swap sg; 
	int tmp;
	for(int i = 0; i < n - 1; i++)
		for(int j = 0; j < n - i - 1; j++)
			if(ord[j] > ord[j + 1]){
				sg.b1 = ordF[ord[j + 1]]; 
				sg.b2 = ordF[ord[j]];
				tmp = ord[j];
				ord[j] = ord[j + 1];
				ord[j + 1] = tmp;
				swaps.push_back(sg);
			}
	free(ord);
	return swaps;
}

};	/* namespace uni10 */	
