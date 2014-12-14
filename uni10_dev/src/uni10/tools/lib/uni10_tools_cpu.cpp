#include <uni10/tools/uni10_tools.h>
#include <string.h>
#include <stdio.h>

namespace uni10{

size_t MEM_USAGE = 0;

void* elemAlloc(size_t memsize, bool& ongpu){
	void* ptr = NULL;
	ptr = malloc(memsize);
  if(ptr == NULL){
    std::ostringstream err;
    err<<"Fails in allocating memory.";
    throw std::runtime_error(exception_msg(err.str()));
  }
	MEM_USAGE += memsize;
	ongpu = false;
	return ptr;
}

void* elemAllocForce(size_t memsize, bool ongpu){
	void* ptr = NULL;
	ptr = malloc(memsize);
  if(ptr == NULL){
    std::ostringstream err;
    err<<"Fails in allocating memory.";
    throw std::runtime_error(exception_msg(err.str()));
  }
	MEM_USAGE += memsize;
	return ptr;
}

void* elemCopy(void* des, const void* src, size_t memsize, bool des_ongpu, bool src_ongpu){
	return memcpy(des, src, memsize);
}

void elemFree(void* ptr, size_t memsize, bool ongpu){
	free(ptr);
	MEM_USAGE -= memsize;
	ptr = NULL;
}
void elemBzero(void* ptr, size_t memsize, bool ongpu){
	memset(ptr, 0, memsize);
}

void elemRand(double* elem, size_t N, bool ongpu){
	for(size_t i = 0; i < N; i++)
		elem[i] = ((double)rand()) / RAND_MAX; //lapack_uni01_sampler();
}

void setDiag(double* elem, double* diag_elem, size_t m, size_t n, size_t diag_n, bool ongpu, bool diag_ongpu){
	int min = m < n ? m : n;
	min = min < diag_n ? min : diag_n;
	for(size_t i = 0; i < min; i++)
		elem[i * n + i] = diag_elem[i];
}
void getDiag(double* elem, double* diag_elem, size_t m, size_t n, size_t diag_n, bool ongpu, bool diag_ongpu){
	int min = m < n ? m : n;
	min = min < diag_n ? min : diag_n;
	for(size_t i = 0; i < min; i++)
		diag_elem[i] = elem[i * n + i];
}
void* mvGPU(void* elem, size_t memsize, bool& ongpu){
	ongpu = false;
	return elem;
}
void* mvCPU(void* elem, size_t memsize, bool& ongpu){
	ongpu = false;
	return elem;
}
void syncMem(void** elemA, void** elemB, size_t memsizeA, size_t memsizeB, bool& ongpuA, bool& ongpuB){
	ongpuA = false;
	ongpuB = false;
}

void shrinkWithoutFree(size_t memsize, bool ongpu){
	MEM_USAGE -= memsize;
}

void reshapeElem(double* oldElem, int bondNum, size_t elemNum, size_t* offset, double* newElem){
  std::ostringstream err;
  err<<"Fatal error(code = T1). Please contact the developer of the uni10 library.";
  throw std::runtime_error(exception_msg(err.str()));
}

double getElemAt(size_t idx, double* elem, bool ongpu){
	return elem[idx];
}

void setElemAt(size_t idx, double val, double* elem, bool ongpu){
	elem[idx] = val;
}

};	/* namespace uni10 */
