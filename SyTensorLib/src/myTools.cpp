#include <boost/random.hpp>
using namespace boost;

mt19937 lapack_rng(777);
uniform_01<mt19937> lapack_uni01_sampler(lapack_rng);

int64_t GPU_MEM_MAX = 0;
int64_t GPU_MEM_USAGE = 0;

void* myMalloc(void* ptr, size_t memsize, int& status){
	ptr = realloc(ptr, memsize);
	assert(ptr != NULL);
	return ptr;
}

void* myMemcpy(void* des, const void* src, size_t memsize, int des_st, int src_st){
	return memcpy(des, src, memsize);
}

void myFree(void* ptr, int status){
	free(ptr);
}
void membzero(void* ptr, size_t memsize, int status){
	memset(ptr, 0, memsize);
}

void randomNums(double* elem, int N, int status){
	for(int i = 0; i < N; i++)
		elem[i] = lapack_uni01_sampler();
}

void myTranspose(double* A, int M, int N, double* AT, int status){
	for(int i = 0; i < M; i++)
		for(int j = 0; j < N; j++)
			AT[j * M + i] = A[i * N + j];
}

void myEye(double* elem, int M, int N, int status){
	int min;
	if(M < N)	min = M;
	else		min = N;
	memset(elem, 0, M * N * sizeof(double));
	for(int i = 0; i < min; i++)
		elem[i * N + i] = 1;
}

