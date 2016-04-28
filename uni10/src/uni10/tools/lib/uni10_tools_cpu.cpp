#include <uni10/tools/uni10_tools.h>
#include <string.h>

namespace uni10{

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
	size_t min = m < n ? m : n;
	min = min < diag_n ? min : diag_n;
	for(size_t i = 0; i < min; i++)
		elem[i * n + i] = diag_elem[i];
}
void getDiag(double* elem, double* diag_elem, size_t m, size_t n, size_t diag_n, bool ongpu, bool diag_ongpu){
	size_t min = m < n ? m : n;
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
    

double  elemMax(double* elem, size_t elemNum, bool ongpu){
    
    if (ongpu) {
        // GPU not implemented
        std::ostringstream err;
        err<<"Fatal error(code = T1). GPU version is not implemented.";
        throw std::runtime_error(exception_msg(err.str()));
    } else {
        double max;
        max=elem[0];
        
        for (size_t i=1; i<elemNum; i++)
        if (max < elem[i]) max=elem[i];
        return max;
    }
}

double  elemAbsMax(double* elem, size_t elemNum, bool ongpu){
    
    if (ongpu) {
        // GPU not implemented
        std::ostringstream err;
        err<<"Fatal error(code = T1). GPU version is not implemented.";
        throw std::runtime_error(exception_msg(err.str()));
    } else {
	
	size_t idx = 0;
        double max = fabs(elem[0]);
        
        for (size_t i=1; i<elemNum; i++)
        if (max < fabs(elem[i])){
	  max=fabs(elem[i]);
	  idx = i;
	}
        return elem[idx];
    }
}

/***** Complex version *****/
std::complex<double> getElemAt(size_t idx, std::complex<double>* elem, bool ongpu){
	return elem[idx];
}

void setElemAt(size_t idx, std::complex<double> val, std::complex<double> *elem, bool ongpu){
	elem[idx] = val;
}

void elemRand(std::complex<double>* elem, size_t N, bool ongpu){
	for(size_t i = 0; i < N; i++)
		elem[i] = std::complex<double>(((double)rand()) / RAND_MAX, ((double)rand()) / RAND_MAX); //lapack_uni01_sampler();
}

void elemCast(std::complex<double>* des, double* src, size_t N, bool des_ongpu, bool src_ongpu){
  for(size_t i = 0; i < N; i++)
    des[i] = src[i];
}
void elemCast(double* des, std::complex<double>* src, size_t N, bool des_ongpu, bool src_ongpu){
  for(size_t i = 0; i < N; i++)
    des[i] = src[i].real();
}
void setDiag(std::complex<double>* elem, std::complex<double>* diag_elem, size_t m, size_t n, size_t diag_n, bool ongpu, bool diag_ongpu){
	size_t min = m < n ? m : n;
	min = min < diag_n ? min : diag_n;
	for(size_t i = 0; i < min; i++)
		elem[i * n + i] = diag_elem[i];
}
void getDiag(std::complex<double>* elem, std::complex<double>* diag_elem, size_t m, size_t n, size_t diag_n, bool ongpu, bool diag_ongpu){
	size_t min = m < n ? m : n;
	min = min < diag_n ? min : diag_n;
	for(size_t i = 0; i < min; i++)
		diag_elem[i] = elem[i * n + i];
}

void reshapeElem(std::complex<double>* oldElem, int bondNum, size_t elemNum, size_t* offset, std::complex<double>* newElem){
  std::ostringstream err;
  err<<"Fatal error(code = T1). Please contact the developer of the uni10 library.";
  throw std::runtime_error(exception_msg(err.str()));
}


};	/* namespace uni10 */
